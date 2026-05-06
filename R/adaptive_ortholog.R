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


# =============================================================================
# Gene Symbol Harmonization Module
# =============================================================================
# Ensures identical gene name casing across all species before Seurat CCA.
# This prevents the "gene not found" errors that occur when different species
# use different casing conventions (human=SOX2, mouse=Sox2, zebrafish=sox2).


#' Validate and harmonize gene names across species
#'
#' After adaptive ortholog alignment, this function:
#'   1. Checks that all species' gene maps produce identical aligned gene sets
#'   2. Reports any case mismatches with suggested fixes
#'   3. Optionally auto-repairs trivially fixable mismatches
#'   4. Returns a clean, validated gene map ready for Seurat merge
#'
#' @param align_result   Output from adaptive_ortholog_align()
#' @param auto_repair    If TRUE, attempt to fix simple casing mismatches
#' @param verbose        Print diagnostic messages
#'
#' @return List with:
#'   \item{validated}{Logical — TRUE if all species agree on gene names}
#'   \item{issues}{Data frame of any mismatches found}
#'   \item{gene_map}{Cleaned gene map (if validated)}
#'   \item{common_genes}{Final set of genes present in ALL species}
#'
#' @export
validate_gene_consensus <- function(align_result, auto_repair = TRUE, verbose = TRUE) {
  gene_map <- align_result$gene_map
  species <- names(align_result$species_info)
  aligned <- gene_map$aligned

  # ---- Check 1: All species have the same aligned gene set ----
  if (verbose) cat("\n=== Gene Consensus Validation ===\n")

  species_gene_sets <- list()
  for (sp in species) {
    sp_map <- gene_map[[sp]]
    names(sp_map) <- aligned
    sp_genes <- aligned[!is.na(sp_map)]
    species_gene_sets[[sp]] <- sp_genes
  }

  # Find genes present in all species
  common <- Reduce(intersect, species_gene_sets)

  # Find genes missing from each species
  issues <- data.frame(
    gene = character(),
    missing_in = character(),
    present_in = character(),
    suggested_fix = character(),
    stringsAsFactors = FALSE
  )

  for (sp in species) {
    missing <- setdiff(aligned, species_gene_sets[[sp]])
    if (length(missing) > 0) {
      for (g in missing) {
        # Which species HAVE this gene?
        have_it <- names(which(sapply(species_gene_sets, function(s) g %in% s)))
        issues <- rbind(issues, data.frame(
          gene = g,
          missing_in = sp,
          present_in = paste(have_it, collapse = ", "),
          suggested_fix = suggest_gene_fix(g, sp, gene_map, align_result),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # ---- Report ----
  if (verbose) {
    cat(sprintf("  Total aligned genes: %d\n", length(aligned)))
    cat(sprintf("  Common to ALL species: %d\n", length(common)))
    for (sp in species) {
      n_present <- length(species_gene_sets[[sp]])
      n_missing <- length(aligned) - n_present
      status <- if (n_missing == 0) "✓" else sprintf("✗ (%d missing)", n_missing)
      cat(sprintf("  %-25s %d genes %s\n", paste0(sp, ":"), n_present, status))
    }

    if (nrow(issues) > 0) {
      cat(sprintf("\n  ⚠ %d gene mismatches detected:\n", nrow(issues)))
      for (i in seq_len(min(10, nrow(issues)))) {
        iss <- issues[i, ]
        cat(sprintf("    %s: missing from %s (present in %s)\n",
                     iss$gene, iss$missing_in, iss$present_in))
        if (nchar(iss$suggested_fix) > 0) {
          cat(sprintf("      → fix: %s\n", iss$suggested_fix))
        }
      }
      if (nrow(issues) > 10) {
        cat(sprintf("    ... and %d more\n", nrow(issues) - 10))
      }
    } else {
      cat("\n  ✓ All genes validated — ready for CCA integration\n")
    }
  }

  # ---- Auto-repair ----
  repaired_map <- gene_map
  if (auto_repair && nrow(issues) > 0) {
    if (verbose) cat("\n  Attempting auto-repair...\n")
    for (i in seq_len(nrow(issues))) {
      iss <- issues[i, ]
      g <- iss$gene
      sp <- iss$missing_in

      # Strategy 1: Case-insensitive lookup in original HVGs
      sp_hvgs <- align_result$hvg_lists[[sp]]
      sp_hvgs_upper <- toupper(sp_hvgs)
      match_idx <- which(sp_hvgs_upper == g)
      if (length(match_idx) == 1) {
        repaired_map[[sp]][which(aligned == g)] <- sp_hvgs[match_idx]
        if (verbose) cat(sprintf("    Fixed %s in %s: %s → %s (case match)\n",
                                  g, sp, NA, sp_hvgs[match_idx]))
        next
      }

      # Strategy 2: Check if it's a multi-copy species missing a single-copy gene
      info <- align_result$species_info[[sp]]
      if (info$copy_type == "multi") {
        # Try regex-based matching against original features
        pattern <- info$copy_pattern %||% "^(\\w+?)([a-z]\\d*|[a-z]?)$"
        # The gene g is uppercase (e.g., SOX2). Try to find matching lower-case paralogs
        # by converting to lowercase and checking against HVGs
        sp_hvgs_lower <- tolower(sp_hvgs)
        g_lower <- tolower(g)
        matching_hvgs <- sp_hvgs[grepl(paste0("^", g_lower), sp_hvgs_lower)]
        if (length(matching_hvgs) > 0) {
          # Pick the one with highest mean expression
          obj <- align_result$seurat_objects %||% NULL
          if (!is.null(obj) && sp %in% names(obj)) {
            expr <- tryCatch({
              Matrix::rowMeans(GetAssayData(obj[[sp]], assay = "RNA", layer = "data"))
            }, error = function(e) NULL)
            if (!is.null(expr) && all(matching_hvgs %in% names(expr))) {
              best <- matching_hvgs[which.max(expr[matching_hvgs])]
            } else {
              best <- matching_hvgs[1]
            }
          } else {
            best <- matching_hvgs[1]
          }
          repaired_map[[sp]][which(aligned == g)] <- best
          if (verbose) cat(sprintf("    Fixed %s in %s: → %s (paralog match)\n",
                                    g, sp, best))
        }
      }
    }
  }

  # ---- Recompute common genes after repair ----
  repaired_sets <- list()
  for (sp in species) {
    sp_map <- repaired_map[[sp]]
    names(sp_map) <- aligned
    repaired_sets[[sp]] <- aligned[!is.na(sp_map)]
  }
  common_repaired <- Reduce(intersect, repaired_sets)

  result <- list(
    validated = (nrow(issues) == 0 || length(common_repaired) == length(common)),
    issues = issues,
    gene_map = repaired_map,
    common_genes = common_repaired,
    n_original = length(common),
    n_repaired = length(common_repaired)
  )

  if (verbose) {
    cat(sprintf("\n  Final consensus: %d genes\n", length(common_repaired)))
    if (length(common_repaired) < 50) {
      warning("Very few consensus genes (", length(common_repaired),
              "). CCA integration may be unstable.")
    }
  }

  return(result)
}


#' Suggest a fix for a gene mismatch
#' @keywords internal
suggest_gene_fix <- function(gene, species, gene_map, align_result) {
  # Check if this gene exists with different case in the species' HVGs
  hvg_list <- align_result$hvg_lists[[species]]
  hvg_upper <- toupper(hvg_list)

  if (gene %in% hvg_upper) {
    idx <- which(hvg_upper == gene)
    return(paste0("Use original gene name: '", hvg_list[idx],
                  "' (case mismatch — add toupper() step)"))
  }

  # Check if there's a similar gene (fuzzy match)
  distances <- utils::adist(gene, hvg_upper)
  if (min(distances) <= 2) {
    closest <- hvg_list[which.min(distances)]
    return(paste0("Closest match: '", closest,
                  "' (distance=", min(distances), " — verify manually)"))
  }

  return("")
}


#' One-shot gene harmonization: align + validate + repair
#'
#' Convenience wrapper that runs adaptive_ortholog_align() followed by
#' validate_gene_consensus() with auto-repair.
#'
#' @param ... Arguments passed to adaptive_ortholog_align()
#' @return List with align_result + validation_result
#' @export
harmonize_cross_species_genes <- function(...) {
  align_result <- adaptive_ortholog_align(...)
  validation <- validate_gene_consensus(align_result, auto_repair = TRUE, verbose = TRUE)

  list(
    align_result = align_result,
    validation = validation,
    gene_map = validation$gene_map,
    common_genes = validation$common_genes
  )
}
