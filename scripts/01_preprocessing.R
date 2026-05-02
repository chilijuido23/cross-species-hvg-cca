#!/usr/bin/env Rscript
# =============================================================================
# Step 1: Preprocessing — Load and QC single-cell data for each species
# =============================================================================
# This script loads 10x Genomics data for each species, performs QC filtering,
# normalization, and saves individual Seurat objects ready for alignment.
#
# Usage:
#   Rscript 01_preprocessing.R --species_config config.yaml --outdir ../data
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

# ---- Parse command-line arguments ----
option_list <- list(
  make_option("--species_config", type = "character", default = NULL,
              help = "YAML/CSV config file with species info"),
  make_option("--outdir", type = "character", default = "../data",
              help = "Output directory for processed Seurat objects [default: %default]"),
  make_option("--min_cells", type = "integer", default = 3,
              help = "Min cells expressing a gene [default: %default]"),
  make_option("--min_features", type = "integer", default = 200,
              help = "Min features per cell [default: %default]"),
  make_option("--max_mt_percent", type = "integer", default = 25,
              help = "Max mitochondrial percentage [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Species configuration ----
# If no config file provided, define inline. Edit this list to match your data.
species_config <- list(
  "human" = list(
    data_dir   = "path/to/human/filtered_feature_bc_matrix",
    copy_type  = "single",
    mt_pattern = "^MT-"
  ),
  "mouse" = list(
    data_dir   = "path/to/mouse/filtered_feature_bc_matrix",
    copy_type  = "single",
    mt_pattern = "^mt-"
  ),
  "lizard" = list(
    data_dir   = "path/to/lizard/filtered_feature_bc_matrix",
    copy_type  = "single",
    mt_pattern = "^MT-"
  ),
  "zebrafish" = list(
    data_dir   = "path/to/zebrafish/filtered_feature_bc_matrix",
    copy_type  = "multi",
    copy_pattern = "^(\\w+?)([a-z]\\d*|[a-z]?)$",
    mt_pattern = "^mt-"
  ),
  "large_yellow_croaker" = list(
    data_dir   = "path/to/croaker/filtered_feature_bc_matrix",
    copy_type  = "multi",
    copy_pattern = "^(\\w+?)([a-z]\\d*|[a-z]?)$",
    mt_pattern = "^mt-"
  )
)

# Override with config file if provided
if (!is.null(opt$species_config) && file.exists(opt$species_config)) {
  species_config <- yaml::read_yaml(opt$species_config)
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Process each species ----
seurat_list <- list()

for (sp in names(species_config)) {
  cat(sprintf("\n========== Processing: %s ==========\n", sp))
  cfg <- species_config[[sp]]

  # Load 10x data
  counts <- Read10X(cfg$data_dir)
  obj <- CreateSeuratObject(counts = counts, project = sp,
                            min.cells = opt$min_cells,
                            min.features = opt$min_features)

  # QC metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = cfg$mt_pattern)

  # Filter
  obj <- subset(obj, subset = nFeature_RNA > opt$min_features &
                         percent.mt < opt$max_mt_percent)

  # Normalize
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)

  # Add species metadata
  obj$species <- sp
  obj$copy_type <- cfg$copy_type

  cat(sprintf("  Cells retained: %d\n", ncol(obj)))
  cat(sprintf("  Genes: %d\n", nrow(obj)))

  seurat_list[[sp]] <- obj

  # Save individual object
  saveRDS(obj, file.path(opt$outdir, paste0(sp, "_preprocessed.rds")))
}

# Save combined list
saveRDS(seurat_list, file.path(opt$outdir, "all_species_preprocessed.rds"))

cat("\n========== Preprocessing complete ==========\n")
cat(sprintf("Saved %d species objects to: %s\n", length(seurat_list), opt$outdir))
