#!/usr/bin/env Rscript
# =============================================================================
# Snakemake-compatible wrapper: Visualization
#
# Used by: workflow/rules/visualize.smk
# Accepts: --input, --outdir, --width, --height, --dpi
# Outputs: fig1-7 PNG files + combined panel
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

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

input_file <- arg_list$input
outdir     <- arg_list$outdir %||% "figures"
fig_width  <- arg_list$width %||% 10
fig_height <- arg_list$height %||% 8
fig_dpi    <- arg_list$dpi %||% 150

if (is.null(input_file)) stop("Required: --input")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load ----
cat("Loading integrated object...\n")
integrated <- readRDS(input_file)
species <- unique(integrated$species)

# ---- Color palette ----
species_colors <- setNames(
  RColorBrewer::brewer.pal(max(3, length(species)), "Set1"),
  species
)

# ---- Figure 1: UMAP by species ----
cat("Fig 1: UMAP by species\n")
p1 <- DimPlot(integrated, reduction = "umap", group.by = "species",
              cols = species_colors, pt.size = 0.5) +
  ggtitle("Cross-Species Hypothalamus Atlas") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(outdir, "fig1_umap_by_species.png"),
       p1, width = fig_width, height = fig_height, dpi = fig_dpi)

# ---- Figure 2: UMAP split by species ----
cat("Fig 2: UMAP split view\n")
ncol_split <- min(3, length(species))
p2 <- DimPlot(integrated, reduction = "umap", split.by = "species",
              cols = species_colors, pt.size = 0.3, ncol = ncol_split) +
  ggtitle("Per-Species UMAP Embedding") +
  theme_minimal()
ggsave(file.path(outdir, "fig2_umap_split.png"),
       p2, width = fig_width * 1.5, height = fig_height, dpi = fig_dpi)

# ---- Figure 3: UMAP by cluster ----
cat("Fig 3: UMAP by cluster\n")
p3 <- DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE, pt.size = 0.5) +
  ggtitle("Integrated Clusters") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(outdir, "fig3_umap_by_cluster.png"),
       p3, width = fig_width, height = fig_height, dpi = fig_dpi)

# ---- Figure 4: Species composition per cluster ----
cat("Fig 4: Species composition\n")
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "fig4_species_composition.png"),
       p4, width = fig_width, height = fig_height * 0.8, dpi = fig_dpi)

# ---- Figure 5: Cell count barplot ----
cat("Fig 5: Cell counts\n")
cell_counts <- as.data.frame(table(integrated$species))
colnames(cell_counts) <- c("Species", "Count")

p5 <- ggplot(cell_counts, aes(x = reorder(Species, -Count), y = Count, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = species_colors) +
  geom_text(aes(label = scales::comma(Count)), vjust = -0.5) +
  labs(title = "Cells per Species", x = "Species", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(file.path(outdir, "fig5_cell_counts.png"),
       p5, width = fig_width * 0.8, height = fig_height * 0.7, dpi = fig_dpi)

# ---- Figure 6: Combined panel ----
cat("Fig 6: Combined publication figure\n")
combined <- (p1 | p3) / (p4 | p5) +
  plot_annotation(
    title = "Cross-Species scRNA-seq Integration",
    subtitle = paste0("Adaptive HVG Ortholog Alignment | ",
                      ncol(integrated), " cells | ",
                      length(species), " species"),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                  plot.subtitle = element_text(hjust = 0.5, size = 11))
  )
ggsave(file.path(outdir, "fig6_combined_panel.png"),
       combined, width = fig_width * 1.8, height = fig_height * 1.3, dpi = fig_dpi)

cat(sprintf("\nDone! %d figures saved to %s/\n",
            length(list.files(outdir, pattern = "\\.png$")), outdir))
