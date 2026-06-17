#!/usr/bin/env node
/**
 * VNNBio MCP Server
 *
 * Exposes VNNBio's interpretable pathway-constrained neural networks to Claude
 * via the Model Context Protocol. Each tool maps to a VNNBio R function; all
 * heavy computation stays in a persistent R process.
 *
 * Tool chain: load_tcga → build_pathway_map → build_architecture → train_vnn → predict → explain
 */

import { Server } from "@modelcontextprotocol/sdk/server/index.js";
import { StdioServerTransport } from "@modelcontextprotocol/sdk/server/stdio.js";
import {
  CallToolRequestSchema,
  ListToolsRequestSchema,
  type Tool,
} from "@modelcontextprotocol/sdk/types.js";
import { RBridge } from "./r-bridge.js";

// ── Tool definitions ────────────────────────────────────────────────────────
// Descriptions are written to guide Claude through correct ordering and
// parameter usage. This is half the prompt engineering for the demo.

const TOOLS: Tool[] = [
  {
    name: "load_tcga",
    description:
      "Load the bundled TCGA kidney cancer demo dataset (KIRC vs KIRP, 663 samples). " +
      "This is always the FIRST step. Returns a data reference and summary statistics. " +
      "The dataset is pre-processed with Ensembl version suffixes stripped.",
    inputSchema: {
      type: "object" as const,
      properties: {
        cancer_type: {
          type: "string",
          enum: ["KIRC_KIRP"],
          description: "Demo dataset identifier. Use 'KIRC_KIRP' for kidney renal clear cell vs papillary carcinoma.",
          default: "KIRC_KIRP",
        },
      },
      required: [],
    },
  },
  {
    name: "load_custom",
    description:
      "Load a custom SummarizedExperiment from an .rds file on disk. " +
      "The SE must have: (1) an expression assay (TPM, counts, etc.) with genes as rows, " +
      "(2) Ensembl gene IDs as rownames (version-stripped, e.g. ENSG00000141510), " +
      "(3) a label column in colData with two factor levels for binary classification. " +
      "Use this instead of load_tcga when working with your own data.",
    inputSchema: {
      type: "object" as const,
      properties: {
        path: {
          type: "string",
          description: "Absolute path to the .rds file containing a SummarizedExperiment",
        },
        label_col: {
          type: "string",
          description: "Column name in colData to use as the class label (default: 'label')",
          default: "label",
        },
      },
      required: ["path"],
    },
  },
  {
    name: "build_pathway_map",
    description:
      "Build a pathway-gene mapping from MSigDB gene set collections. " +
      "This defines which genes belong to which biological pathways, creating the " +
      "structure for the visible neural network. Must be called BEFORE build_architecture. " +
      "Pass the data_ref from load_tcga to align the mask rows to expression features. " +
      "Common categories: 'H' (Hallmark, 50 gene sets), 'C2' with subcategory 'CP:KEGG' (KEGG), " +
      "'C5' with subcategory 'GO:BP' (GO Biological Process).",
    inputSchema: {
      type: "object" as const,
      properties: {
        data_ref: {
          type: "string",
          description: "Reference key from load_tcga (e.g. 'data_1'). Aligns mask rows to expression features.",
        },
        category: {
          type: "string",
          description: "MSigDB category: 'H' (Hallmark), 'C2' (curated), 'C5' (GO), etc.",
          default: "H",
        },
        subcategory: {
          type: "string",
          description: "MSigDB subcategory (e.g. 'CP:KEGG', 'GO:BP'). Only needed for C2/C5.",
        },
        species: {
          type: "string",
          description: "Species for gene mapping",
          default: "Homo sapiens",
        },
        gene_id_type: {
          type: "string",
          enum: ["ensembl_gene", "entrez_gene", "gene_symbol"],
          description: "Gene ID type matching the expression data. Default 'ensembl_gene'.",
          default: "ensembl_gene",
        },
      },
      required: ["data_ref"],
    },
  },
  {
    name: "build_architecture",
    description:
      "Build the VNN (Visible Neural Network) architecture from a pathway map. " +
      "The architecture mirrors the biological pathway hierarchy — each hidden node " +
      "corresponds to a real pathway, making the model inherently interpretable. " +
      "Requires a map_ref from build_pathway_map. Must be called BEFORE train_vnn.",
    inputSchema: {
      type: "object" as const,
      properties: {
        map_ref: {
          type: "string",
          description: "Reference key from build_pathway_map (e.g. 'map_1')",
        },
        activation: {
          type: "string",
          enum: ["tanh", "relu", "sigmoid", "gelu", "swish"],
          description: "Activation function for hidden pathway nodes (default: 'tanh')",
          default: "tanh",
        },
        n_output: {
          type: "number",
          description: "Number of output nodes. 1 for binary classification (default).",
          default: 1,
        },
      },
      required: ["map_ref"],
    },
  },
  {
    name: "train_vnn",
    description:
      "Train the visible neural network on expression data. The model learns " +
      "pathway-level representations that can be decomposed via Shapley values. " +
      "Requires data_ref from load_tcga AND arch_ref from build_architecture. " +
      "Training uses a Julia backend for performance. Returns training metrics " +
      "including loss and top pathway importance scores. May take 1-5 minutes.",
    inputSchema: {
      type: "object" as const,
      properties: {
        data_ref: {
          type: "string",
          description: "Reference key from load_tcga (e.g. 'data_1')",
        },
        arch_ref: {
          type: "string",
          description: "Reference key from build_architecture (e.g. 'arch_1')",
        },
        label_col: {
          type: "string",
          description: "Column in colData(se) to use as the response variable (default: 'label')",
          default: "label",
        },
        epochs: {
          type: "number",
          description: "Number of training epochs (default: 50)",
          default: 50,
        },
        lr: {
          type: "number",
          description: "Learning rate (default: 0.001)",
          default: 0.001,
        },
        patience: {
          type: "number",
          description: "Early stopping patience — stop if val_loss doesn't improve for this many epochs. 0 to disable. (default: 10)",
          default: 10,
        },
      },
      required: ["data_ref", "arch_ref"],
    },
  },
  {
    name: "predict",
    description:
      "Run predictions on data using a trained VNN model. Returns predicted classes, " +
      "probabilities, and accuracy against true labels. Requires model_ref from train_vnn " +
      "and data_ref from load_tcga.",
    inputSchema: {
      type: "object" as const,
      properties: {
        model_ref: {
          type: "string",
          description: "Reference key from train_vnn (e.g. 'model_1')",
        },
        data_ref: {
          type: "string",
          description: "Reference key from load_tcga (e.g. 'data_1')",
        },
      },
      required: ["model_ref", "data_ref"],
    },
  },
  {
    name: "explain",
    description:
      "Compute per-pathway Shapley values for a specific sample using " +
      "shapleyPathwayAttribution(), revealing WHY the model made its prediction. " +
      "Each pathway gets a signed Shapley value: positive pushes toward class 1, " +
      "negative toward class 0. The decomposition exploits VNN structure for speed — " +
      "pathway activations are precomputed once, then coalitions evaluated via output " +
      "layer dot products. The efficiency axiom holds: values sum to prediction minus " +
      "baseline. Returns the top pathways ranked by |contribution|. " +
      "Requires model_ref from train_vnn and data_ref from load_tcga.",
    inputSchema: {
      type: "object" as const,
      properties: {
        model_ref: {
          type: "string",
          description: "Reference key from train_vnn (e.g. 'model_1')",
        },
        data_ref: {
          type: "string",
          description: "Reference key from load_tcga (e.g. 'data_1')",
        },
        sample_index: {
          type: "number",
          description: "1-based index of the sample to explain (default: 1)",
          default: 1,
        },
        n_perm: {
          type: "number",
          description: "Number of permutations for Shapley approximation. More = more precise but slower. (default: 200)",
          default: 200,
        },
      },
      required: ["model_ref", "data_ref"],
    },
  },
  {
    name: "visualize",
    description:
      "Generate publication-quality figures from trained VNNBio models. " +
      "Three plot types available: " +
      "(1) 'cohort_importance' — bar chart of pathway coefficients showing which pathways " +
      "the model learned as most important for classification across all patients. " +
      "(2) 'patient_shapley' — waterfall chart showing per-pathway Shapley contributions " +
      "for a specific patient, with red bars pushing toward one class and blue toward the other. " +
      "(3) 'probability_distribution' — histogram of predicted probabilities for all patients, " +
      "colored by true outcome, with a specific patient highlighted. " +
      "Requires model_ref from train_vnn and data_ref from load_tcga/load_custom. " +
      "Returns the file path of the saved PNG figure.",
    inputSchema: {
      type: "object" as const,
      properties: {
        plot_type: {
          type: "string",
          enum: ["cohort_importance", "patient_shapley", "probability_distribution"],
          description: "Type of figure to generate",
        },
        model_ref: {
          type: "string",
          description: "Reference key from train_vnn (e.g. 'model_1')",
        },
        data_ref: {
          type: "string",
          description: "Reference key from load_tcga or load_custom (e.g. 'data_1')",
        },
        sample_index: {
          type: "number",
          description: "Patient index for patient_shapley and probability_distribution plots (default: 1)",
          default: 1,
        },
        top_n: {
          type: "number",
          description: "Number of top pathways to show (default: 17)",
          default: 17,
        },
      },
      required: ["plot_type", "model_ref", "data_ref"],
    },
  },
];

// ── Server setup ────────────────────────────────────────────────────────────

const rBridge = new RBridge();

const server = new Server(
  { name: "vnnbio", version: "0.1.0" },
  { capabilities: { tools: {} } }
);

// List available tools
server.setRequestHandler(ListToolsRequestSchema, async () => ({
  tools: TOOLS,
}));

// ── Lazy R bridge — starts on first tool call, not at server boot ────────────

let bridgeReady: Promise<Record<string, unknown>> | null = null;

function ensureBridge(): Promise<Record<string, unknown>> {
  if (!bridgeReady) {
    process.stderr.write("[vnnbio-mcp] Starting R bridge (first tool call)...\n");
    bridgeReady = rBridge.start().then((info) => {
      process.stderr.write(
        `[vnnbio-mcp] R bridge ready: ${JSON.stringify(info)}\n`
      );
      return info;
    });
  }
  return bridgeReady;
}

// Dispatch tool calls to R
server.setRequestHandler(CallToolRequestSchema, async (request) => {
  const { name, arguments: args } = request.params;

  // Validate tool exists
  const tool = TOOLS.find((t) => t.name === name);
  if (!tool) {
    return {
      content: [{ type: "text", text: `Unknown tool: ${name}` }],
      isError: true,
    };
  }

  try {
    // Ensure R bridge is running before dispatching
    await ensureBridge();
    const result = await rBridge.call(name, (args ?? {}) as Record<string, unknown>);

    return {
      content: [
        {
          type: "text",
          text: JSON.stringify(result, null, 2),
        },
      ],
    };
  } catch (err) {
    const message = err instanceof Error ? err.message : String(err);
    return {
      content: [
        {
          type: "text",
          text: `VNNBio error in '${name}': ${message}`,
        },
      ],
      isError: true,
    };
  }
});

// ── Main ────────────────────────────────────────────────────────────────────

async function main(): Promise<void> {
  // Connect MCP transport FIRST so the Inspector/Claude handshake works
  // immediately. R bridge starts lazily on first tool call.
  const transport = new StdioServerTransport();
  await server.connect(transport);
  process.stderr.write("[vnnbio-mcp] MCP server connected via stdio\n");
  process.stderr.write("[vnnbio-mcp] R bridge will start on first tool call\n");

  // Clean shutdown
  const shutdown = async () => {
    process.stderr.write("[vnnbio-mcp] Shutting down...\n");
    await rBridge.stop();
    process.exit(0);
  };

  process.on("SIGINT", shutdown);
  process.on("SIGTERM", shutdown);
}

main().catch((err) => {
  process.stderr.write(`[vnnbio-mcp] Fatal error: ${err}\n`);
  process.exit(1);
});
