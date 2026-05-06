#!/usr/bin/env Rscript
# =============================================================================
# Snakemake-compatible wrapper: CCA Integration
#
# Used by: workflow/rules/integrate.smk
# Accepts: --preprocessed, --align_result, --outdir, --n_dim, --resolution
# Outputs: integrated_seurat.rds
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
})

# Source core algorithm for validation
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

preprocessed_file <- arg_list$preprocessed
align_file        <- arg_list$align_result
outdir            <- arg_list$outdir %||% "data"
n_dim             <- arg_list$n_dim %||% 30
resolution        <- arg_list$resolution %||% 0.5
integration_method <- arg_list$integration_method %||% "cca"

if (is.null(preprocessed_file)) stop("Required: --preprocessed")
if (is.null(align_file)) stop("Required: --align_result")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load ----
cat("Loading data...\n")
seurat_list <- readRDS(preprocessed_file)
align_result <- readRDS(align_file)

# ---- Validate gene consensus (PRE-CCA CHECK) ----
cat("\n========================================\n")
cat("  Pre-CCA Gene Name Validation\n")
cat("========================================\n")

# Store seurat objects in align_result for repair strategies
align_result$seurat_objects <- seurat_list

validation <- validate_gene_consensus(align_result, auto_repair = TRUE, verbose = TRUE)

if (!validation$validated && length(validation$common_genes) < 50) {
  stop(sprintf(
    "Gene consensus validation FAILED. Only %d common genes found (minimum 50 required).\n
    Check the issues above or re-run with lower n_hvg.",
    length(validation$common_genes)
  ))
}

# Use validated gene map
gene_map <- validation$gene_map
aligned_genes <- validation$common_genes
cat(sprintf("\nUsing validated gene set: %d genes for CCA\n", length(aligned_genes)))

# ---- Prepare objects ----
cat("Preparing objects for integration...\n")
seurat_prepped <- list()

for (sp in names(seurat_list)) {
  cat(sprintf("  %s: ", sp))
  obj <- seurat_list[[sp]]
  sp_genes <- gene_map[[sp]]
  names(sp_genes) <- gene_map$aligned
  sp_genes <- sp_genes[!is.na(sp_genes)]

  obj <- subset(obj, features = sp_genes)
  rev_map <- stats::setNames(names(sp_genes), sp_genes)
  new_names <- rev_map[rownames(obj)]
  keep <- !is.na(new_names)
  obj <- obj[keep, ]
  new_names <- new_names[keep]

  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  rownames(counts) <- new_names
  obj <- CreateSeuratObject(counts = counts, meta.data = obj@meta.data)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj$species <- sp

  seurat_prepped[[sp]] <- obj
  cat(sprintf("%d genes x %d cells\n", nrow(obj), ncol(obj)))
}

# ---- Merge ----
cat("Merging species...\n")
merged <- merge(seurat_prepped[[1]], y = seurat_prepped[-1],
                add.cell.ids = names(seurat_prepped))

# ---- Integration ----
cat(sprintf("Running %s integration (%d dims)...\n", toupper(integration_method), n_dim))

obj_list <- SplitObject(merged, split.by = "species")
obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  x
})

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 2000)

if (integration_method == "rpca") {
  obj_list <- lapply(obj_list, function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    x
  })
  anchors <- FindIntegrationAnchors(
    object.list = obj_list, anchor.features = features,
    reduction = "rpca", dims = 1:n_dim, verbose = TRUE
  )
} else {
  anchors <- FindIntegrationAnchors(
    object.list = obj_list, anchor.features = features,
    reduction = "cca", dims = 1:n_dim, verbose = TRUE
  )
}

integrated <- IntegrateData(anchorset = anchors, dims = 1:n_dim, verbose = TRUE)

# ---- Post-integration ----
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = min(50, n_dim), verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:n_dim, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:n_dim, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = resolution, verbose = FALSE)

saveRDS(integrated, file.path(outdir, "integrated_seurat.rds"))
cat(sprintf("\nIntegration complete: %d cells, %d clusters\n",
            ncol(integrated), length(unique(integrated$seurat_clusters))))
