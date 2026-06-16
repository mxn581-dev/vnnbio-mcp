#!/usr/bin/env Rscript
# prepare-demo-data.R
# Run once to bundle the TCGA KIRC/KIRP SummarizedExperiment for the MCP demo.
# Expects you already have the SE from your JOSS benchmark pipeline.
#
# Usage:
#   Rscript data/prepare-demo-data.R /path/to/your/tcga_kirc_kirp_se.rds
#
# Or edit EXISTING_SE_PATH below and source() this script.

suppressMessages({
  library(SummarizedExperiment)
})

EXISTING_SE_PATH <- commandArgs(trailingOnly = TRUE)[1]

if (is.na(EXISTING_SE_PATH) || !file.exists(EXISTING_SE_PATH)) {
  cat("Usage: Rscript data/prepare-demo-data.R <path-to-existing-SE.rds>\n\n")
  cat("The SE should have:\n")
  cat("  - assay 'counts' with gene expression values\n")
  cat("  - colData column 'label' with factors 'KIRC' and 'KIRP'\n")
  cat("  - rownames as Ensembl IDs (version suffixes will be stripped)\n\n")
  cat("If you don't have one, this script can build from TCGAbiolinks.\n")
  cat("Set BUILD_FROM_SCRATCH=TRUE below.\n")
  stop("No input file provided")
}

se <- readRDS(EXISTING_SE_PATH)

# Strip Ensembl version suffixes (ENSG00000141510.18 → ENSG00000141510)
rownames(se) <- sub("\\.\\d+$", "", rownames(se))

# Sanity checks
stopifnot("label" %in% colnames(colData(se)))
stopifnot(all(levels(colData(se)$label) %in% c("KIRC", "KIRP")))

cat(sprintf(
  "Bundling: %d genes × %d samples (%s)\n",
  nrow(se), ncol(se),
  paste(names(table(colData(se)$label)), table(colData(se)$label),
        sep = "=", collapse = ", ")
))

out_path <- file.path(dirname(sys.frame(1)$ofile %||% "."), "tcga_kirc_kirp.rds")
if (!interactive()) out_path <- "data/tcga_kirc_kirp.rds"

saveRDS(se, out_path)
cat("Saved to:", out_path, "\n")
