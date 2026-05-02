# Data Directory

This directory contains processed data files:

- `*_preprocessed.rds` — Individual Seurat objects per species (from Step 1)
- `all_species_preprocessed.rds` — Combined list of all Seurat objects
- `align_result.rds` — Adaptive ortholog alignment result (from Step 2)
- `aligned_gene_map.csv` — Gene-level mapping across species
- `aligned_genes.txt` — Gene list for CCA integration
- `integrated_seurat.rds` — Final integrated object (from Step 3)

**Note**: Raw sequencing data should be placed in a separate location
(e.g., `~/data/raw/` or an external drive). Only processed intermediate
files go here.

This directory is gitignored. Run the pipeline scripts to populate it.
