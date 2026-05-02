#' Adaptive Ortholog Alignment for Cross-Species scRNA-seq Integration
#'
#' @description
#' Aligns genes across species for Seurat CCA integration by leveraging
#' Highly Variable Genes (HVG) and adaptive ranking for multi-copy genes.
#'
#' The core challenge: fish species have genome duplications resulting in
#' multi-copy genes (e.g., sox2a, sox2b). Standard ortholog-based alignment
#' fails because you cannot determine which paralog corresponds to the
#' single-copy mammalian ortholog.
#'
#' Our solution:
#'   1. Compute top N HVGs independently for each species
#'   2. Single-copy species: convert gene symbols to UPPERCASE
#'   3. Multi-copy species (fish): for each gene family, rank paralogs by
#'      mean expression and select the top-ranked copy
#'   4. Intersect gene sets across all species
#'   5. Output aligned gene list ready for Seurat CCA
#'
#' @param seurat_objects  Named list of Seurat objects
#' @param species_info    Named list with entries: "copy_type" ("single"|"multi"),
#'                        and for multi-copy species, "copy_pattern" (regex to
#'                        identify paralogs, e.g., "^(\\w+?)[ab]\\d*$")
#' @param n_hvg           Number of HVGs per species (default: 2000)
#' @param assay           Assay to use (default: "RNA")
#' @param verbose         Print progress messages (default: TRUE)
#'
#' @return A list with:
#'   \item{aligned_genes}{Character vector of aligned gene symbols}
#'   \item{gene_map}{Data frame mapping aligned genes to original per-species genes}
#'   \item{hvg_lists}{Per-species list of original HVGs}
#'   \item{multi_copy_choices}{For multi-copy species, which paralog was selected}
#'
#' @import Seurat
#' @importFrom Matrix rowMeans
#' @export
adaptive_ortholog_align <- function(
    seurat_objects,
    species_info,
    n_hvg = 2000,
    assay = "RNA",
    verbose = TRUE
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required. Install with: install.packages('Seurat')")
  }

  species_names <- names(seurat_objects)
  if (is.null(species_names)) {
    species_names <- paste0("species_", seq_along(seurat_objects))
    names(seurat_objects) <- species_names
  }

  # Validate species_info matches seurat_objects
  missing_info <- setdiff(species_names, names(species_info))
  if (length(missing_info) > 0) {
    stop("species_info missing entries for: ", paste(missing_info, collapse = ", "))
  }

  # ---- Step 1: Compute HVGs per species ----
  if (verbose) cat("\n=== Step 1: Computing HVGs per species ===\n")
  hvg_lists <- list()

  for (sp in species_names) {
    if (verbose) cat(sprintf("  %s: finding %d HVGs...\n", sp, n_hvg))
    obj <- seurat_objects[[sp]]

    # DefaultAssay may need setting
    if (DefaultAssay(obj) != assay) {
      DefaultAssay(obj) <- assay
    }

    # Normalize if not already done
    if (!"data" %in% names(obj[[assay]])) {
      if (verbose) cat(sprintf("    Normalizing %s...\n", sp))
      obj <- NormalizeData(obj, verbose = FALSE)
    }

    obj <- FindVariableFeatures(obj, nfeatures = n_hvg, verbose = FALSE)
    hvg_lists[[sp]] <- VariableFeatures(obj)
    if (verbose) cat(sprintf("    -> %d HVGs found\n", length(hvg_lists[[sp]])))
    seurat_objects[[sp]] <- obj
  }

  # ---- Step 2: Process gene names ----
  if (verbose) cat("\n=== Step 2: Processing gene names ===\n")

  processed_genes <- list()
  multi_copy_choices <- list()

  for (sp in species_names) {
    info <- species_info[[sp]]
    hvg <- hvg_lists[[sp]]

    if (info$copy_type == "single") {
      # Single-copy species: just uppercase
      if (verbose) cat(sprintf("  %s (single-copy): converting %d genes to uppercase\n", sp, length(hvg)))
      processed_genes[[sp]] <- toupper(hvg)
      multi_copy_choices[[sp]] <- data.frame(
        original = hvg,
        aligned = toupper(hvg),
        is_multi = FALSE,
        stringsAsFactors = FALSE
      )

    } else if (info$copy_type == "multi") {
      # Multi-copy species: adaptive ranking within gene families
      if (verbose) cat(sprintf("  %s (multi-copy): adaptive ranking of %d HVGs...\n", sp, length(hvg)))

      pattern <- info$copy_pattern
      if (is.null(pattern)) {
        # Default pattern for fish: gene root + optional letter + optional digit
        # e.g., sox2a, sox2b, sox21 -> family "sox2"
        pattern <- "^(\\w+?)([a-z]\\d*|[a-z]?)$"
      }

      # Get expression matrix
      expr_mat <- GetAssayData(seurat_objects[[sp]], assay = assay, layer = "data")
      mean_expr <- Matrix::rowMeans(expr_mat)
      mean_expr <- mean_expr[hvg]

      # Identify gene families using regex
      # Strategy: strip trailing letters/digits to find the root gene name
      gene_names <- names(mean_expr)

      # Extract root names by removing trailing paralog suffixes
      root_names <- sub("([a-z]\\d+|\\d+[a-z]?|[a-z])$", "", gene_names, perl = TRUE)

      # Group by root name
      family_map <- split(gene_names, root_names)
      multi_families <- family_map[lengths(family_map) > 1]

      if (verbose) cat(sprintf("    %d gene families with multiple copies detected\n", length(multi_families)))

      choices <- list()
      for (fam_name in names(multi_families)) {
        copies <- multi_families[[fam_name]]
        # Rank by mean expression, select top
        expr_vals <- mean_expr[copies]
        best_copy <- names(which.max(expr_vals))
        choices[[fam_name]] <- data.frame(
          family = fam_name,
          selected = best_copy,
          all_copies = paste(copies, collapse = ";"),
          max_expr = max(expr_vals),
          stringsAsFactors = FALSE
        )
        if (verbose && length(multi_families) <= 20) {
          cat(sprintf("    %s: %s selected (copies: %s)\n",
                       fam_name, best_copy, paste(copies, collapse = ", ")))
        }
      }

      # Build aligned gene list: use root name (uppercase) for alignment
      aligned_names <- character(length(gene_names))
      names(aligned_names) <- gene_names

      for (g in gene_names) {
        root <- sub("([a-z]\\d+|\\d+[a-z]?|[a-z])$", "", g, perl = TRUE)
        # If this gene is a multi-copy family, check if it's the selected one
        if (root %in% names(multi_families)) {
          best <- choices[[root]]$selected
          if (g == best) {
            aligned_names[g] <- toupper(root)  # Use root name for alignment
          } else {
            aligned_names[g] <- NA_character_  # Non-selected copies dropped
          }
        } else {
          aligned_names[g] <- toupper(g)
        }
      }

      # Filter out NAs
      aligned_names <- aligned_names[!is.na(aligned_names)]
      processed_genes[[sp]] <- unname(aligned_names)

      multi_copy_choices[[sp]] <- data.frame(
        original = names(aligned_names),
        aligned = unname(aligned_names),
        is_multi = names(aligned_names) != unname(aligned_names),
        stringsAsFactors = FALSE
      )

      if (verbose) {
        n_dropped <- length(hvg) - length(processed_genes[[sp]])
        cat(sprintf("    -> %d genes retained, %d non-selected copies dropped\n",
                     length(processed_genes[[sp]]), n_dropped))
      }
    } else {
      stop("Unknown copy_type for ", sp, ": ", info$copy_type,
           ". Must be 'single' or 'multi'.")
    }
  }

  # ---- Step 3: Intersect across species ----
  if (verbose) cat("\n=== Step 3: Intersecting gene sets ===\n")

  aligned_genes <- Reduce(intersect, processed_genes)
  if (verbose) {
    for (sp in species_names) {
      cat(sprintf("  %s: %d genes\n", sp, length(processed_genes[[sp]])))
    }
    cat(sprintf("  Intersection: %d common genes\n", length(aligned_genes)))
  }

  if (length(aligned_genes) < 50) {
    warning("Very few aligned genes (", length(aligned_genes),
            "). Consider reducing n_hvg or checking copy_pattern regex.")
  }

  # ---- Step 4: Build gene map ----
  if (verbose) cat("\n=== Step 4: Building gene map ===\n")

  gene_map <- data.frame(aligned = aligned_genes, stringsAsFactors = FALSE)
  for (sp in species_names) {
    choices_df <- multi_copy_choices[[sp]]
    # Map aligned genes back to original genes
    mapping <- stats::setNames(choices_df$original, choices_df$aligned)
    gene_map[[sp]] <- mapping[aligned_genes]
    if (verbose) {
      n_mapped <- sum(!is.na(mapping[aligned_genes]))
      cat(sprintf("  %s: %d/%d genes mapped\n", sp, n_mapped, length(aligned_genes)))
    }
  }

  # ---- Return ----
  result <- list(
    aligned_genes  = aligned_genes,
    gene_map       = gene_map,
    hvg_lists      = hvg_lists,
    multi_copy_choices = multi_copy_choices,
    species_info   = species_info,
    n_hvg          = n_hvg
  )
  class(result) <- "adaptive_ortholog"

  if (verbose) cat("\n=== Alignment complete ===\n")
  return(result)
}


#' Print summary of adaptive ortholog alignment
#' @export
print.adaptive_ortholog <- function(x, ...) {
  cat("\nAdaptive Ortholog Alignment Results\n")
  cat("====================================\n\n")
  cat(sprintf("Species integrated: %d\n", length(x$species_info)))
  cat(sprintf("HVGs per species:  %d\n", x$n_hvg))
  cat(sprintf("Final aligned genes: %d\n\n", length(x$aligned_genes)))

  cat("Per-species summary:\n")
  for (sp in names(x$species_info)) {
    n_multi <- sum(x$multi_copy_choices[[sp]]$is_multi)
    n_total <- nrow(x$multi_copy_choices[[sp]])
    cat(sprintf("  %-20s  type=%-7s  genes=%d  multi-copy=%d\n",
                 sp, x$species_info[[sp]]$copy_type, n_total, n_multi))
  }
  cat("\n")
}


#' Plot alignment diagnostics
#' @export
plot_alignment_diagnostics <- function(align_result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting")
  }

  plots <- list()

  # 1. Gene overlap UpSet-style summary
  n_species <- length(align_result$hvg_lists)
  sp_names <- names(align_result$hvg_lists)
  n_genes_vec <- sapply(align_result$hvg_lists, length)

  overlap_df <- data.frame(
    species = factor(sp_names, levels = sp_names),
    n_hvg   = n_genes_vec,
    n_aligned = c(n_genes_vec[1], rep(NA, n_species - 1)),
    stringsAsFactors = FALSE
  )
  overlap_df$n_aligned[1] <- length(align_result$aligned_genes)

  p1 <- ggplot2::ggplot(overlap_df, ggplot2::aes(x = species, y = n_hvg)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = length(align_result$aligned_genes),
                        linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::annotate("text", x = 1, y = length(align_result$aligned_genes),
                      label = paste0("Intersection: ", length(align_result$aligned_genes)),
                      vjust = -1, color = "red", size = 3.5) +
    ggplot2::labs(
      title = "Gene Count per Species and Intersection",
      x = "Species", y = "Number of Genes"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  plots$overlap <- p1

  # 2. Multi-copy gene selection summary (if any multi-copy species)
  multi_sp <- names(which(sapply(align_result$species_info, function(x) x$copy_type == "multi")))
  if (length(multi_sp) > 0) {
    all_choices <- do.call(rbind, lapply(multi_sp, function(sp) {
      df <- align_result$multi_copy_choices[[sp]]
      df$species <- sp
      df
    }))
    multi_only <- all_choices[all_choices$is_multi, ]

    if (nrow(multi_only) > 0) {
      p2 <- ggplot2::ggplot(multi_only, ggplot2::aes(x = species, fill = species)) +
        ggplot2::geom_bar() +
        ggplot2::labs(
          title = "Multi-copy Gene Families Resolved",
          x = "Species", y = "Number of Gene Families"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
      plots$multi_copy <- p2
    }
  }

  return(plots)
}
