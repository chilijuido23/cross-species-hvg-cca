#!/usr/bin/env Rscript
# =============================================================================
# Step 3: Cross-Species CCA Integration with Seurat
# =============================================================================
# Uses the aligned gene list from Step 2 to perform Seurat CCA-based
# integration across all species.
#
# Usage:
#   Rscript 03_cca_integration.R \
#     --input ../data/all_species_preprocessed.rds \
#     --align_result ../data/align_result.rds \
#     --outdir ../data \
#     --n_dim 30
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

option_list <- list(
  make_option("--input", type = "character",
              default = "../data/all_species_preprocessed.rds",
              help = "Input RDS with list of Seurat objects"),
  make_option("--align_result", type = "character",
              default = "../data/align_result.rds",
              help = "Alignment result from Step 2"),
  make_option("--outdir", type = "character", default = "../data",
              help = "Output directory"),
  make_option("--n_dim", type = "integer", default = 30,
              help = "Number of CCA dimensions [default: %default]"),
  make_option("--n_hvg", type = "integer", default = 2000,
              help = "HVGs per species for FindVariableFeatures [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
cat("Loading data...\n")
seurat_list <- readRDS(opt$input)
align_result <- readRDS(opt$align_result)
aligned_genes <- align_result$aligned_genes

cat(sprintf("Aligned genes for CCA: %d\n", length(aligned_genes)))

# ---- For each species, subset to aligned genes and rename using aligned names ----
cat("\nPreparing objects for integration...\n")

seurat_prepped <- list()

for (sp in names(seurat_list)) {
  cat(sprintf("  %s:\n", sp))

  obj <- seurat_list[[sp]]
  gene_map <- align_result$gene_map

  # Map original genes -> aligned genes for this species
  sp_genes <- gene_map[[sp]]
  names(sp_genes) <- gene_map$aligned

  # Remove NAs (genes not found in this species)
  sp_genes <- sp_genes[!is.na(sp_genes)]

  cat(sprintf("    %d genes mapped to aligned set\n", length(sp_genes)))

  # Subset to mapped genes and rename
  obj <- subset(obj, features = sp_genes)

  # Rename features to aligned (uppercase) names
  # Create reverse map: original -> aligned
  rev_map <- stats::setNames(names(sp_genes), sp_genes)
  new_names <- rev_map[rownames(obj)]

  # Handle any remaining NAs
  keep <- !is.na(new_names)
  obj <- obj[keep, ]
  new_names <- new_names[keep]

  # Rename (requires rebuilding counts matrix)
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  rownames(counts) <- new_names
  obj <- CreateSeuratObject(counts = counts, meta.data = obj@meta.data)
  obj <- NormalizeData(obj, verbose = FALSE)

  obj$species <- sp
  seurat_prepped[[sp]] <- obj
  cat(sprintf("    Final: %d genes x %d cells\n", nrow(obj), ncol(obj)))
}

# ---- Merge into single object ----
cat("\nMerging species...\n")
merged <- merge(seurat_prepped[[1]],
                y = seurat_prepped[-1],
                add.cell.ids = names(seurat_prepped))

# ---- Integration pipeline ----
cat("\nRunning integration pipeline...\n")

# Split by species
obj_list <- SplitObject(merged, split.by = "species")

# Normalize and find variable features per species
obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = opt$n_hvg, verbose = FALSE)
  x
})

# Select integration features
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = opt$n_hvg)
cat(sprintf("Integration features selected: %d\n", length(features)))

# Find integration anchors (CCA)
cat("Finding integration anchors (CCA)...\n")
anchors <- FindIntegrationAnchors(
  object.list   = obj_list,
  anchor.features = features,
  reduction     = "cca",
  dims          = 1:opt$n_dim,
  verbose       = TRUE
)

# Integrate data
cat("Integrating data...\n")
integrated <- IntegrateData(
  anchorset = anchors,
  dims      = 1:opt$n_dim,
  verbose   = TRUE
)

# ---- Post-integration processing ----
DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = min(50, opt$n_dim), verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:opt$n_dim, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:opt$n_dim, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 0.5, verbose = FALSE)

# ---- Save ----
saveRDS(integrated, file.path(opt$outdir, "integrated_seurat.rds"))

cat(sprintf("\nIntegration complete!\n"))
cat(sprintf("  Integrated object: %s/integrated_seurat.rds\n", opt$outdir))
cat(sprintf("  Cells: %d | Genes: %d\n", ncol(integrated), nrow(integrated)))
cat(sprintf("  UMAP computed: %d dimensions\n", opt$n_dim))
