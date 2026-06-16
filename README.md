# VNNBio-MCP

**Give Claude eyes into biological pathways.**

An MCP server that lets Claude drive [VNNBio](https://github.com/YOUR_HANDLE/VNNBio) — interpretable, pathway-constrained neural networks for genomics. Claude chains the tools to classify cancer subtypes, then explains *which biological pathways* drove each prediction using exact Shapley decomposition.

```
You: "Classify these kidney tumors and explain the biology"

Claude: [load_tcga → build_pathway_map → build_architecture → train_vnn → predict → explain]

"The model classified sample TCGA_042 as KIRP with 98.3% confidence.
 The top contributing pathway was KEGG_CELL_CYCLE (Shapley = +0.47),
 followed by KEGG_P53_SIGNALING_PATHWAY (−0.31)..."
```

## Architecture

```
Claude ←→ MCP Server (TypeScript) ←→ Persistent R Process (VNNBio + Julia)
            stdio JSON-RPC              stdin/stdout JSON lines
```

R objects (models, data, architectures) stay in R memory. The MCP server only passes lightweight string handles (`"model_1"`, `"data_1"`).

## Requirements

- **Node.js** ≥ 18
- **R** ≥ 4.3 with packages: `VNNBio`, `jsonlite`, `SummarizedExperiment`
- **Julia** ≥ 1.9 (used by VNNBio's training backend via JuliaConnectoR)

## Quick Start

```bash
git clone https://github.com/YOUR_HANDLE/vnnbio-mcp.git
cd vnnbio-mcp
npm install
npm run build
```

### Add to Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "vnnbio": {
      "command": "node",
      "args": ["/absolute/path/to/vnnbio-mcp/dist/index.js"]
    }
  }
}
```

Restart Claude Desktop. You should see "vnnbio" in the MCP tools menu.

### Bundle Demo Data (optional)

To use real TCGA KIRC/KIRP data instead of synthetic:

```r
# prepare-demo-data.R
library(VNNBio)
# ... load your TCGA SummarizedExperiment, strip Ensembl versions ...
# rownames(se) <- sub("\\.\\d+$", "", rownames(se))
saveRDS(se, "data/tcga_kirc_kirp.rds")
```

Without bundled data the server generates a synthetic dataset so the full pipeline is testable.

## Tools

| Tool | What it does | Required inputs |
|---|---|---|
| `load_tcga` | Load KIRC vs KIRP demo data | — |
| `build_pathway_map` | Map genes → pathways via MSigDB | collection, species |
| `build_architecture` | Create pathway-constrained VNN layers | map_ref |
| `train_vnn` | Train on expression data (Julia backend) | data_ref, arch_ref |
| `predict` | Classify samples | model_ref, data_ref |
| `explain` | Per-pathway Shapley values for one sample | model_ref, data_ref, sample_index |

**Chain order:** `load_tcga` → `build_pathway_map` → `build_architecture` → `train_vnn` → `predict` / `explain`

## Demo Prompt

> Classify TCGA kidney tumor subtypes using KEGG pathways. Use the bundled KIRC vs KIRP dataset. After training, explain the prediction for sample 1 — which biological pathways mattered most?

## Testing

```bash
# Test with MCP Inspector (no Claude Desktop needed)
npm run inspect
```

## How It Works

VNNBio is a **Visible Neural Network** — the architecture mirrors the biological pathway hierarchy from MSigDB. Each hidden node *is* a pathway. This means:

1. **No black box.** Every prediction decomposes into per-pathway contributions.
2. **Exact Shapley values.** The DAG structure enables exact (not approximate) Shapley computation. The efficiency axiom holds: values sum exactly to the prediction logit.
3. **Biological interpretability.** "KEGG_CELL_CYCLE contributed +0.47 toward KIRP" is a statement a biologist can reason about.

The MCP layer lets Claude orchestrate this pipeline conversationally — load data, choose pathway databases, train, predict, and explain, all through natural language.

## License

MIT
