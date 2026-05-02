#!/usr/bin/env Rscript
# =============================================================================
# Step 2: Adaptive Ortholog Gene Alignment
# =============================================================================
# Core innovation: aligns HVGs across species using adaptive ranking for
# multi-copy genes (fish). Produces a gene list ready for Seurat CCA integration.
#
# Usage:
#   Rscript 02_adaptive_ortholog_alignment.R \
#     --input ../data/all_species_preprocessed.rds \
#     --n_hvg 2000 \
#     --outdir ../data
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

# Source the core algorithm
source("../R/adaptive_ortholog.R")

# ---- Parse arguments ----
option_list <- list(
  make_option("--input", type = "character", default = "../data/all_species_preprocessed.rds",
              help = "Input RDS file with list of Seurat objects"),
  make_option("--n_hvg", type = "integer", default = 2000,
              help = "Number of HVGs per species [default: %default]"),
  make_option("--outdir", type = "character", default = "../data",
              help = "Output directory [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
cat(sprintf("Loading data from: %s\n", opt$input))
seurat_list <- readRDS(opt$input)

# ---- Define species info ----
# IMPORTANT: Edit this to match your species composition
species_info <- list(
  "human" = list(copy_type = "single"),
  "mouse" = list(copy_type = "single"),
  "lizard" = list(copy_type = "single"),
  "zebrafish" = list(
    copy_type = "multi",
    copy_pattern = "^(\\w+?)([a-z]\\d*|[a-z]?)$"
  ),
  "large_yellow_croaker" = list(
    copy_type = "multi",
    copy_pattern = "^(\\w+?)([a-z]\\d*|[a-z]?)$"
  )
)

# ---- Run adaptive ortholog alignment ----
cat("\n============================================================\n")
cat("  Adaptive Ortholog Gene Alignment\n")
cat("  HVG per species:", opt$n_hvg, "\n")
cat("============================================================\n")

align_result <- adaptive_ortholog_align(
  seurat_objects = seurat_list,
  species_info   = species_info,
  n_hvg          = opt$n_hvg,
  verbose        = TRUE
)

# ---- Display summary ----
print(align_result)

# ---- Plot diagnostics ----
png(file.path(opt$outdir, "alignment_diagnostics.png"), width = 10, height = 6, units = "in", res = 150)
plots <- plot_alignment_diagnostics(align_result)
gridExtra::grid.arrange(grobs = plots, ncol = min(2, length(plots)))
dev.off()
cat(sprintf("\nDiagnostics plot saved to: %s/alignment_diagnostics.png\n", opt$outdir))

# ---- Save results ----
saveRDS(align_result, file.path(opt$outdir, "align_result.rds"))

# Write gene map to CSV
write.csv(align_result$gene_map, file.path(opt$outdir, "aligned_gene_map.csv"), row.names = FALSE)

# Write aligned gene list (for Seurat)
writeLines(align_result$aligned_genes, file.path(opt$outdir, "aligned_genes.txt"))

cat(sprintf("\nResults saved to: %s/\n", opt$outdir))
cat("  - align_result.rds        (full R object)\n")
cat("  - aligned_gene_map.csv    (gene mapping table)\n")
cat("  - aligned_genes.txt       (gene list for Seurat CCA)\n")
cat("\nReady for CCA integration!\n")
