#!/usr/bin/env Rscript
# =============================================================================
# Step 4: Visualization — Cross-species integration results
# =============================================================================
# Generates publication-quality figures for the integrated cross-species
# single-cell atlas.
#
# Usage:
#   Rscript 04_visualization.R \
#     --input ../data/integrated_seurat.rds \
#     --outdir ../figures
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(optparse)
})

option_list <- list(
  make_option("--input", type = "character",
              default = "../data/integrated_seurat.rds",
              help = "Integrated Seurat object from Step 3"),
  make_option("--outdir", type = "character", default = "../figures",
              help = "Output directory for figures"),
  make_option("--width", type = "integer", default = 10,
              help = "Figure width in inches"),
  make_option("--height", type = "integer", default = 8,
              help = "Figure height in inches")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
cat("Loading integrated object...\n")
integrated <- readRDS(opt$input)

# ---- Color palettes ----
species_colors <- c(
  "human"    = "#E41A1C",
  "mouse"    = "#377EB8",
  "lizard"   = "#4DAF4A",
  "zebrafish" = "#984EA3",
  "large_yellow_croaker" = "#FF7F00"
)

# ---- Figure 1: UMAP colored by species ----
cat("Generating Figure 1: UMAP by species...\n")
p1 <- DimPlot(integrated, reduction = "umap", group.by = "species",
              cols = species_colors, pt.size = 0.5) +
  ggtitle("Cross-Species Hypothalamus Atlas") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(opt$outdir, "fig1_umap_by_species.png"),
       p1, width = opt$width, height = opt$height, dpi = 150)

# ---- Figure 2: UMAP split by species (facet) ----
cat("Generating Figure 2: UMAP split view...\n")
p2 <- DimPlot(integrated, reduction = "umap", split.by = "species",
              cols = species_colors, pt.size = 0.3, ncol = 3) +
  ggtitle("Per-Species UMAP") +
  theme_minimal()

ggsave(file.path(opt$outdir, "fig2_umap_split.png"),
       p2, width = opt$width * 1.5, height = opt$height, dpi = 150)

# ---- Figure 3: UMAP colored by cluster ----
cat("Generating Figure 3: UMAP by cluster...\n")
p3 <- DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE, pt.size = 0.5) +
  ggtitle("Integrated Clusters") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(opt$outdir, "fig3_umap_by_cluster.png"),
       p3, width = opt$width, height = opt$height, dpi = 150)

# ---- Figure 4: Species composition per cluster ----
cat("Generating Figure 4: Species composition per cluster...\n")
composition <- table(integrated$seurat_clusters, integrated$species)
composition_prop <- prop.table(composition, margin = 1)

comp_df <- as.data.frame(composition_prop)
colnames(comp_df) <- c("Cluster", "Species", "Proportion")

p4 <- ggplot(comp_df, aes(x = Cluster, y = Proportion, fill = Species)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = species_colors) +
  labs(title = "Species Composition per Cluster",
       x = "Cluster", y = "Proportion") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(opt$outdir, "fig4_species_composition.png"),
       p4, width = opt$width, height = opt$height * 0.8, dpi = 150)

# ---- Figure 5: CCA colored by species ----
cat("Generating Figure 5: CCA plot...\n")
# Run CCA on the integrated object if not already present
if (!"cca" %in% names(integrated@reductions)) {
  cat("  Running CCA...\n")
  obj_list <- SplitObject(integrated, split.by = "species")
  obj_list <- lapply(obj_list, function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, nfeatures = 2000, verbose = FALSE)
    x
  })
  features <- SelectIntegrationFeatures(obj_list, nfeatures = 2000)
  obj_list <- lapply(obj_list, function(x) {
    ScaleData(x, features = features, verbose = FALSE)
  })
  integrated <- RunCCA(integrated, features = features, group.by = "species",
                       num.cc = 30, verbose = FALSE)
}

p5 <- DimPlot(integrated, reduction = "cca", group.by = "species",
              cols = species_colors, pt.size = 0.5, dims = c(1, 2)) +
  ggtitle("CCA: Cross-Species Integration") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(opt$outdir, "fig5_cca_by_species.png"),
       p5, width = opt$width, height = opt$height, dpi = 150)

# ---- Figure 6: Cell count barplot ----
cat("Generating Figure 6: Cell count...\n")
cell_counts <- as.data.frame(table(integrated$species))
colnames(cell_counts) <- c("Species", "Count")

p6 <- ggplot(cell_counts, aes(x = reorder(Species, -Count), y = Count, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = species_colors) +
  geom_text(aes(label = scales::comma(Count)), vjust = -0.5) +
  labs(title = "Cell Count per Species", x = "Species", y = "Number of Cells") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(opt$outdir, "fig6_cell_counts.png"),
       p6, width = opt$width * 0.8, height = opt$height * 0.7, dpi = 150)

# ---- Figure 7: Combined publication figure ----
cat("Generating Figure 7: Combined panel...\n")
combined <- (p1 | p3) / (p4 | p6) +
  plot_annotation(
    title = "Cross-Species Hypothalamus Single-Cell Atlas",
    subtitle = paste0("Adaptive HVG-based ortholog alignment | ",
                      ncol(integrated), " cells | ", length(unique(integrated$species)), " species"),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                  plot.subtitle = element_text(hjust = 0.5, size = 11))
  )

ggsave(file.path(opt$outdir, "fig7_combined_panel.png"),
       combined, width = opt$width * 1.8, height = opt$height * 1.3, dpi = 150)

# ---- Save integrated object with updated reductions ----
saveRDS(integrated, opt$input)

cat(sprintf("\nVisualization complete! %d figures saved to: %s/\n",
            length(list.files(opt$outdir, pattern = "\\.png$")), opt$outdir))
