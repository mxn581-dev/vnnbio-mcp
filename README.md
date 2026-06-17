# VNNBio-MCP

**Give Claude eyes into biological pathways.**

An MCP server that lets Claude drive [VNNBio](https://github.com/mxn581-dev/VNNBio) — interpretable, pathway-constrained neural networks for genomics. Claude chains the tools to classify disease subtypes, then explains *which biological pathways* drove each prediction using per-pathway attribution.

```
You: "Classify these kidney tumors and explain the biology"

Claude: [load_tcga → build_pathway_map → build_architecture → train → predict → explain]

"The model classified sample 71 as KIRC with high confidence.
 The top contributing pathway was HALLMARK_HYPOXIA — consistent with
 VHL-loss driven pseudohypoxic signaling that defines clear cell RCC."
```

## Demo Results

Tested on three real-world clinical datasets from a single English prompt each:

| Disease | Dataset | Pathways | CV AUC | Key Finding |
|---|---|---|---|---|
| Kidney cancer | TCGA KIRC vs KIRP | 50 Hallmark | 0.954 | HYPOXIA drives KIRC (VHL loss) |
| Kidney cancer | TCGA KIRC vs KIRP | 186 KEGG | 0.94 | PPAR signaling drives KIRC |
| Sepsis mortality | GSE65682 ICU | 50 Hallmark | 0.69 | HYPOXIA drives death (tissue hypoperfusion) |

Every finding is biologically correct — discovered from expression data alone, then contextualized by Claude in real time.

## Architecture

```
Claude ←→ MCP Server (TypeScript) ←→ Persistent R Process (VNNBio)
            stdio JSON-RPC              stdin/stdout JSON lines
```

R objects (models, data, architectures) stay in R memory. The MCP server only passes lightweight string handles (`"model_1"`, `"data_1"`).

## Requirements

- **Node.js** ≥ 18
- **R** ≥ 4.4 with packages: `VNNBio`, `jsonlite`, `SummarizedExperiment`, `msigdbr`, `glmnet`
- **Claude Code** or **Claude Desktop** (for conversational usage)

## Quick Start

```bash
git clone https://github.com/mxn581-dev/vnnbio-mcp.git
cd vnnbio-mcp
npm install
npm run build
```

### Add to Claude Code

```bash
claude mcp add vnnbio node dist/index.js
```

Then start Claude Code and ask:

> Classify kidney tumor subtypes using Hallmark pathways. Load the demo data, build a pathway map, train a model, predict, and explain sample 1.

### Add to Claude Desktop

Edit `~/.config/Claude/claude_desktop_config.json` (Linux) or `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS):

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

### Test with MCP Inspector

```bash
npm run inspect
```

## Tools (7 total)

| Tool | What it does | Inputs |
|---|---|---|
| `load_tcga` | Load bundled KIRC vs KIRP demo data | — |
| `load_custom` | Load any SummarizedExperiment .rds file | `path`, `label_col` |
| `build_pathway_map` | Map genes → pathways via MSigDB | `data_ref`, `category`, `gene_id_type` |
| `build_architecture` | Create pathway-constrained VNN layers | `map_ref` |
| `train_vnn` | Train classifier on pathway-projected expression | `data_ref`, `arch_ref` |
| `predict` | Classify samples | `model_ref`, `data_ref` |
| `explain` | Per-pathway attribution for one sample | `model_ref`, `data_ref`, `sample_index` |

**Chain order:** `load_tcga` or `load_custom` → `build_pathway_map` → `build_architecture` → `train_vnn` → `predict` / `explain`

## Bring Your Own Data

The `load_custom` tool accepts any `.rds` file containing a `SummarizedExperiment` with:

- **Assay**: expression matrix (TPM, counts, log-counts) — genes as rows, samples as columns
- **rownames**: Gene IDs matching your chosen `gene_id_type` (Ensembl IDs or gene symbols)
- **colData**: must have a label column with two factor levels for binary classification

```r
library(SummarizedExperiment)

se <- SummarizedExperiment(
    assays  = list(tpm = your_expression_matrix),
    colData = DataFrame(label = factor(your_labels))
)
rownames(se) <- your_gene_ids
saveRDS(se, "my_dataset.rds")
```

Then in Claude: *"Load /path/to/my_dataset.rds and classify using KEGG pathways."*

## How It Works

The model projects gene expression through a **biological pathway mask** from MSigDB:

```
H = X × mask    (genes → pathway activations)
ŷ = glmnet(H)   (pathway activations → classification)
```

Each hidden node *is* a pathway. After training, per-pathway importance comes directly from the model coefficients — no post-hoc interpretation needed. For individual samples, the explanation decomposes the prediction into `activation × coefficient` per pathway, showing exactly which biological processes drove that patient's classification.

## Example: Sepsis Mortality Prediction

```
> Load sepsis_mortality.rds and explain why patient 71 died.

Claude: "The prediction is driven overwhelmingly by an elevated HYPOXIA
program (activation 734) — consistent with tissue hypoperfusion as a
hallmark of fatal sepsis — with the P53 stress response reinforcing it.
The pro-survival metabolic pathways were present but collectively
insufficient to overcome the hypoxic signal."
```

## License

MIT

## Citation

```
@software{VNNBio-MCP,
  title  = {VNNBio-MCP: Interpretable Pathway Neural Networks via Model Context Protocol},
  author = {Minh Hieu Nguyen},
  year   = {2026},
  url    = {https://github.com/mxn581-dev/vnnbio-mcp}
}
```
