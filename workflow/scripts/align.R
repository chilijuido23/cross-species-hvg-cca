#!/usr/bin/env Rscript
# =============================================================================
# Snakemake-compatible wrapper: Adaptive Ortholog Alignment
#
# Used by: workflow/rules/align.smk
# Accepts: --input, --species_config, --n_hvg, --outdir
# Outputs: align_result.rds, aligned_gene_map.csv, aligned_genes.txt
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

# Source core algorithm (relative to project root)
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))
][1]))
script_dir <- sub("^--file=", "", script_dir)
project_root <- normalizePath(file.path(script_dir, "..", ".."))
source(file.path(project_root, "R", "adaptive_ortholog.R"))

# ---- Parse args ----
args <- commandArgs(trailingOnly = TRUE)
arg_list <- list()
for (i in seq(1, length(args), by = 2)) {
  key <- sub("^--", "", args[i])
  val <- args[i + 1]
  if (!is.na(suppressWarnings(as.numeric(val)))) {
    val <- as.numeric(val)
  }
  arg_list[[key]] <- val
}

input_file         <- arg_list$input
species_config_file <- arg_list$species_config
n_hvg              <- arg_list$n_hvg %||% 2000
outdir             <- arg_list$outdir %||% "data"

if (is.null(input_file)) stop("Required: --input")
if (is.null(species_config_file)) stop("Required: --species_config")

# ---- Load ----
cat(sprintf("Loading: %s\n", input_file))
seurat_list <- readRDS(input_file)

species_config <- read_yaml(species_config_file)

# Convert species config to species_info format
species_info <- list()
for (sp in names(species_config)) {
  cfg <- species_config[[sp]]
  info <- list(copy_type = cfg$copy_type)
  if (cfg$copy_type == "multi" && !is.null(cfg$copy_pattern)) {
    info$copy_pattern <- cfg$copy_pattern
  }
  species_info[[sp]] <- info
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Run alignment ----
align_result <- adaptive_ortholog_align(
  seurat_objects = seurat_list,
  species_info   = species_info,
  n_hvg          = n_hvg,
  verbose        = TRUE
)

print(align_result)

# ---- Diagnostics plot ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

png(file.path(outdir, "alignment_diagnostics.png"),
    width = 10, height = 6, units = "in", res = 150)
plots <- plot_alignment_diagnostics(align_result)
gridExtra::grid.arrange(grobs = plots, ncol = min(2, length(plots)))
dev.off()

# ---- Save ----
saveRDS(align_result, file.path(outdir, "align_result.rds"))
write.csv(align_result$gene_map, file.path(outdir, "aligned_gene_map.csv"),
          row.names = FALSE)
writeLines(align_result$aligned_genes, file.path(outdir, "aligned_genes.txt"))

cat(sprintf("\nResults saved to: %s/\n", outdir))
