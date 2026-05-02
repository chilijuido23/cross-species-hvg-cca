#!/usr/bin/env Rscript
# =============================================================================
# Snakemake-compatible wrapper: Preprocessing
# 
# Used by: workflow/rules/preprocess.smk
# Accepts: --species_config, --outdir, --min_features, --max_mt_percent
# Outputs: all_species_preprocessed.rds
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

# ---- Parse Snakemake-style args ----
args <- commandArgs(trailingOnly = TRUE)
arg_list <- list()
for (i in seq(1, length(args), by = 2)) {
  key <- sub("^--", "", args[i])
  val <- args[i + 1]
  # Convert numeric args
  if (!is.na(suppressWarnings(as.numeric(val)))) {
    val <- as.numeric(val)
  }
  arg_list[[key]] <- val
}

species_config_file <- arg_list$species_config
outdir              <- arg_list$outdir %||% "data"
min_features        <- arg_list$min_features %||% 200
max_mt_percent      <- arg_list$max_mt_percent %||% 25

if (is.null(species_config_file)) {
  stop("Required argument: --species_config")
}

# ---- Load species config ----
species_config <- read_yaml(species_config_file)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Process each species ----
seurat_list <- list()

for (sp in names(species_config)) {
  cat(sprintf("\n========== Processing: %s ==========\n", sp))
  cfg <- species_config[[sp]]

  # Check data exists
  if (!dir.exists(cfg$data_dir)) {
    stop(sprintf("Data directory not found for %s: %s", sp, cfg$data_dir))
  }

  # Load 10x data
  counts <- Read10X(cfg$data_dir)
  obj <- CreateSeuratObject(
    counts = counts,
    project = sp,
    min.cells = 3,
    min.features = min_features
  )

  # QC metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = cfg$mt_pattern)

  # Filter
  obj <- subset(obj, subset = nFeature_RNA > min_features &
                         percent.mt < max_mt_percent)

  # Normalize
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)

  # Metadata
  obj$species <- sp
  obj$copy_type <- cfg$copy_type

  cat(sprintf("  Cells: %d | Genes: %d\n", ncol(obj), nrow(obj)))
  seurat_list[[sp]] <- obj

  saveRDS(obj, file.path(outdir, paste0(sp, "_preprocessed.rds")))
}

# Save combined
saveRDS(seurat_list, file.path(outdir, "all_species_preprocessed.rds"))
cat(sprintf("\nPreprocessing complete: %d species saved to %s/\n",
            length(seurat_list), outdir))
