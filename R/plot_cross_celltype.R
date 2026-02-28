#' Plot Cross-Cell-Type Interaction Heatmap
#'
#' @description
#' Visualises the cross-cell-type interaction structure for a gene pair as a
#' heatmap.  Each tile represents a directed cell-type pair
#' (source type \eqn{\to} neighbour type), coloured by the Pearson correlation
#' between gene A expression in the source cells and gene B expression in the
#' neighbouring cells.
#'
#' This is the primary visualisation for the trans-cellular synergy metric
#' introduced in scPairs 0.1.3.  It reveals *which* cell-type interfaces carry
#' the cross-type signal (e.g. Adora2a in T-cells correlated with Ido1 in
#' dendritic cells).
#'
#' @param object A Seurat object.
#' @param gene1 Character; first gene.
#' @param gene2 Character; second gene.
#' @param result Optional \code{scPairs_pair_result} from
#'     \code{\link{AssessGenePair}}.  If NULL, the pair is assessed internally.
#' @param assay Character; assay name.  Default: \code{DefaultAssay(object)}.
#' @param slot Character; data slot.  Default \code{"data"}.
#' @param cluster_col Character; meta.data column with cell-type labels.
#'     NULL = use \code{Idents(object)}.
#' @param neighbourhood_k Integer; k for KNN graph.  Default 20.
#' @param neighbourhood_reduction Character; reduction for KNN graph.
#'     Default \code{"pca"}.
#' @param min_cross_pairs Integer; minimum cross-type pairs per tile.
#'     Tiles with fewer pairs are greyed out.  Default 30.
#' @param min_pct_expressed Numeric; minimum percentage of cells (0-1) in a
#'     cell type that must express a gene. Default 0.01 (1%). Prevents
#'     spurious correlations with very sparse genes.
#' @param show_n Logical; annotate each tile with the number of cross-type
#'     neighbour pairs.  Default TRUE.
#' @param show_reverse Logical; if TRUE (default), show a second panel for
#'     the reverse direction (gene2 in source \eqn{\to} gene1 in neighbour).
#' @param diverging Logical; use a diverging red--white--blue colour scale
#'     centred at 0.  Default TRUE.
#' @param title Character; overall title.  NULL = auto-generated.
#'
#' @return A \code{ggplot} heatmap (or two-panel patchwork when
#'   \code{show_reverse = TRUE}) with cell-type pairs on axes and synergy
#'   enrichment encoded by colour.
#'
#' @family Section_2_Visualization
#'
#' @seealso \code{\link{AssessGenePair}} for cross-cell-type metrics,
#'   \code{\link{PlotPairDimplot}} for individual-cell UMAP display.
#'
#' @export
#'
#' @examples
#' # scpairs_testdata has clusters (seurat_clusters) and PCA already built in.
#' PlotPairCrossType(scpairs_testdata,
#'                  gene1 = "GENE3",
#'                  gene2 = "GENE4")
#'
PlotPairCrossType <- function(object,
                              gene1,
                              gene2,
                              result                  = NULL,
                              assay                   = NULL,
                              slot                    = "data",
                              cluster_col             = NULL,
                              neighbourhood_k         = 20,
                              neighbourhood_reduction = "pca",
                              min_cross_pairs         = 30,
                              min_pct_expressed       = 0.01,
                              show_n                  = TRUE,
                              show_reverse            = TRUE,
                              diverging               = TRUE,
                              title                   = NULL) {

  validated <- .validate_plot_inputs(object, gene1, gene2,
                                      assay = assay, slot = slot)
  assay <- validated$assay

  .validate_percentage(min_pct_expressed, "min_pct_expressed")

  # ------------------------------------------------------------------
  # Obtain cross-cell-type detail
  # ------------------------------------------------------------------
  if (!is.null(result) && inherits(result, "scPairs_pair_result")) {
    detail <- result$cross_celltype_detail
    global_r_ab <- result$metrics$cross_celltype_r_ab
    global_r_ba <- result$metrics$cross_celltype_r_ba
    score <- result$metrics$cross_celltype_score
  } else {
    # Compute from scratch
    mat <- .get_expression_matrix(object, features = c(gene1, gene2),
                                   assay = assay, slot = slot)
    if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
      mat <- as.matrix(mat)
    }

    embed <- tryCatch(
      Seurat::Embeddings(object, reduction = neighbourhood_reduction),
      error = function(e) {
        stop(sprintf("Reduction '%s' not found.", neighbourhood_reduction),
             call. = FALSE)
      }
    )
    dims_use <- seq_len(min(30, ncol(embed)))

    if (!is.null(cluster_col)) {
      if (!(cluster_col %in% colnames(object@meta.data))) {
        stop(sprintf("cluster_col '%s' not found in meta.data.", cluster_col),
             call. = FALSE)
      }
      cluster_ids <- as.factor(object@meta.data[[cluster_col]])
    } else {
      cluster_ids <- as.factor(Seurat::Idents(object))
    }

    cross_res <- .cross_celltype(
      x = mat[gene1, ], y = mat[gene2, ],
      cluster_ids = cluster_ids,
      embed = embed[, dims_use, drop = FALSE],
      min_bins = min_cross_pairs %/% 3,
      min_pct_expressed = min_pct_expressed
    )
    detail <- cross_res$per_celltype_pair
    global_r_ab <- cross_res$r_ab
    global_r_ba <- cross_res$r_ba
    score <- cross_res$score
  }

  # Validation
  if (is.null(detail) || nrow(detail) == 0) {
    # Provide helpful diagnostic message
    mat_check <- .get_expression_matrix(object, features = c(gene1, gene2),
                                         assay = assay, slot = slot)
    if (inherits(mat_check, "dgCMatrix") || inherits(mat_check, "dgRMatrix")) {
      mat_check <- as.matrix(mat_check)
    }
    
    pct1 <- mean(mat_check[gene1, ] > 0) * 100
    pct2 <- mean(mat_check[gene2, ] > 0) * 100
    
    msg <- sprintf(
      "No cross-cell-type pairs with sufficient data found.\n",
      "Gene expression: %s = %.1f%%, %s = %.1f%%\n",
      "Possible reasons:\n",
      "  - One or both genes are too sparse (< %.1f%% in cell types)\n",
      "  - Insufficient bins with both cell types (min_cross_pairs = %d)\n",
      "  - Only one cell type present\n",
      "Try: reducing min_cross_pairs, min_pct_expressed, or check gene expression.",
      gene1, pct1, gene2, pct2, min_pct_expressed * 100, min_cross_pairs
    )
    message(msg)
    
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::ggtitle("No cross-cell-type pairs",
                              subtitle = sprintf("%s: %.1f%% cells, %s: %.1f%% cells",
                                                 gene1, pct1, gene2, pct2)))
  }

  # ------------------------------------------------------------------
  # Build full grid (all source x neighbour type combinations)
  # ------------------------------------------------------------------
  ct_levels <- sort(unique(c(detail$source_type, detail$neighbour_type)))

  # The new micro-environment approach stores both directions in the same

  # detail data.frame: r_g1_to_g2 and r_g2_to_g1.
  # Forward panel uses r_g1_to_g2, reverse panel uses r_g2_to_g1.

  # Normalise column name for n: accept both old (n_pairs) and new (n_bins_valid)
  n_col <- if ("n_bins_valid" %in% colnames(detail)) "n_bins_valid" else "n_pairs"

  # Build reverse detail by renaming r column
  rev_detail <- detail
  if ("r_g2_to_g1" %in% colnames(detail)) {
    rev_detail$r_g1_to_g2 <- detail$r_g2_to_g1
  }

  # Compute shared colour scale limits across both panels
  all_r <- c(detail$r_g1_to_g2, rev_detail$r_g1_to_g2)
  shared_r_max <- max(abs(all_r), na.rm = TRUE)
  shared_r_max <- max(shared_r_max, 0.1)

  fwd_panel <- .build_heatmap_panel(
    detail      = detail,
    ct_levels   = ct_levels,
    gene_src    = gene1,
    gene_nbr    = gene2,
    direction   = "forward",
    global_r    = global_r_ab,
    show_n      = show_n,
    diverging   = diverging,
    n_col       = n_col,
    r_max       = shared_r_max
  )

  if (!show_reverse) {
    if (is.null(title)) {
      title <- paste0("Cross-cell-type: ", gene1, " \u2192 ", gene2)
    }
    return(
      fwd_panel +
        ggplot2::ggtitle(title,
                          subtitle = .cross_subtitle(score, global_r_ab, NULL))
    )
  }

  rev_panel <- .build_heatmap_panel(
    detail      = rev_detail,
    ct_levels   = ct_levels,
    gene_src    = gene2,
    gene_nbr    = gene1,
    direction   = "reverse",
    global_r    = global_r_ba,
    show_n      = show_n,
    diverging   = diverging,
    n_col       = n_col,
    r_max       = shared_r_max
  )

  if (is.null(title)) {
    title <- paste0("Cross-cell-type interactions: ", gene1, " & ", gene2)
  }

  # Stack vertically when many cell types to avoid squeezing
  n_types <- length(ct_levels)
  if (n_types > 6) {
    layout_ncol <- 1
  } else {
    layout_ncol <- 2
  }

  patchwork::wrap_plots(fwd_panel, rev_panel, ncol = layout_ncol) +
    patchwork::plot_annotation(
      title    = title,
      subtitle = .cross_subtitle(score, global_r_ab, global_r_ba)
    )
}


# ===========================================================================
#  Internal helpers
# ===========================================================================

#' Build a single heatmap panel for cross-cell-type correlation
#' @keywords internal
.build_heatmap_panel <- function(detail, ct_levels,
                                 gene_src, gene_nbr, direction,
                                 global_r, show_n, diverging,
                                 n_col = "n_bins_valid",
                                 r_max = NULL) {

  n_types <- length(ct_levels)

  # Full grid of all cell-type pairs (excluding self-pairs)
  grid <- expand.grid(
    source_type    = ct_levels,
    neighbour_type = ct_levels,
    stringsAsFactors = FALSE
  )
  grid <- grid[grid$source_type != grid$neighbour_type, , drop = FALSE]

  # Merge detail onto grid
  grid <- merge(grid, detail[, c("source_type", "neighbour_type",
                                  n_col, "r_g1_to_g2"), drop = FALSE],
                by = c("source_type", "neighbour_type"),
                all.x = TRUE)

  # Rename n column for consistency
  grid$n_count <- grid[[n_col]]

  # Adaptive text size for tile annotations
  text_size <- if (n_types <= 4) 2.8 else if (n_types <= 7) 2.2 else 1.8
  axis_size <- if (n_types <= 4) 9 else if (n_types <= 7) 8 else 7

  # Determine label for count column
  n_label_prefix <- if (n_col == "n_bins_valid") "bins=" else "n="

  # Tile label: only show on tiles with data
  if (show_n) {
    has_data <- !is.na(grid$n_count) & !is.na(grid$r_g1_to_g2)
    grid$label <- ifelse(
      has_data,
      paste0("r=", format(round(grid$r_g1_to_g2, 2), nsmall = 2),
             "\n", n_label_prefix, grid$n_count),
      ""
    )
  }

  # Mark tiles with data
  grid$sufficient <- !is.na(grid$n_count) & !is.na(grid$r_g1_to_g2)

  # For tiles without data, set r to NA so they appear grey
  grid$r_plot <- ifelse(grid$sufficient, grid$r_g1_to_g2, NA_real_)

  # Factor levels
  grid$source_type    <- factor(grid$source_type,    levels = ct_levels)
  grid$neighbour_type <- factor(grid$neighbour_type, levels = rev(ct_levels))

  # Build plot
  p <- ggplot2::ggplot(grid,
    ggplot2::aes(x = source_type, y = neighbour_type, fill = r_plot)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8)

  # Colour scale
  if (diverging) {
    if (is.null(r_max)) {
      r_max <- max(abs(grid$r_plot), na.rm = TRUE)
      r_max <- max(r_max, 0.1)
    }
    p <- p + ggplot2::scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-r_max, r_max),
      na.value = "grey90",
      name     = "Pearson r"
    )
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(
      option   = "D",
      na.value = "grey90",
      name     = "Pearson r"
    )
  }

  if (show_n) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label),
      size       = text_size,
      colour     = "grey20",
      lineheight = 0.85
    )
  }

  panel_title <- paste0(gene_src, " (source) \u2192 ", gene_nbr, " (neighbour)")

  p +
    ggplot2::labs(
      x = paste0("Source (", gene_src, ")"),
      y = paste0("Neighbour (", gene_nbr, ")")
    ) +
    ggplot2::ggtitle(panel_title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x    = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1,
                                              size = axis_size),
      axis.text.y    = ggplot2::element_text(size = axis_size),
      axis.title     = ggplot2::element_text(size = 9),
      plot.title     = ggplot2::element_text(size = 10, face = "bold"),
      panel.grid     = ggplot2::element_blank(),
      legend.key.width  = ggplot2::unit(0.4, "cm"),
      legend.key.height = ggplot2::unit(0.8, "cm")
    )
}


#' Format subtitle with global cross-cell-type stats
#' @keywords internal
.cross_subtitle <- function(score, r_ab, r_ba) {
  parts <- c()
  if (!is.null(score) && !is.na(score)) {
    parts <- c(parts, paste0("Score = ", round(score, 3)))
  }
  if (!is.null(r_ab) && !is.na(r_ab)) {
    parts <- c(parts, paste0("r(A\u2192B) = ", round(r_ab, 3)))
  }
  if (!is.null(r_ba) && !is.na(r_ba)) {
    parts <- c(parts, paste0("r(B\u2192A) = ", round(r_ba, 3)))
  }
  if (length(parts) == 0) return(NULL)
  paste(parts, collapse = "  |  ")
}



