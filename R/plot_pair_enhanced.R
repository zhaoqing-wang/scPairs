#' Enhanced Co-Expression Visualization with Neighbourhood Smoothing
#'
#' @description
#' A six-panel visualization that shows both raw and KNN-smoothed expression
#' for a gene pair on the UMAP (or other reduction) embedding.  The top row
#' shows raw expression (gene1, gene2, product); the bottom row shows
#' KNN-smoothed expression.  This is particularly informative for gene pairs
#' that are not co-expressed in the same cell but share neighbourhood-level
#' co-expression patterns.
#'
#' @param object A Seurat object with a dimensionality reduction.
#' @param gene1 Character; first gene.
#' @param gene2 Character; second gene.
#' @param reduction Character; reduction for plotting.  Default `"umap"`.
#' @param smooth_reduction Character; reduction for KNN graph.  Default `"pca"`.
#' @param k Integer; neighbourhood size for smoothing.  Default 20.
#' @param alpha Numeric in \[0,1\]; self-weight for smoothing.  Default 0.3.
#' @param assay Character; assay.
#' @param slot Character; data slot.
#' @param pt_size Numeric; point size.
#' @param pt_alpha Numeric; point alpha.
#' @param title Character; overall title.
#'
#' @return A combined \code{ggplot} (patchwork) with 6 panels: three showing
#'   raw expression and three showing KNN-smoothed expression, for \code{gene1},
#'   \code{gene2}, and their co-expression product.
#'
#' @family Section_2_Visualization
#'
#' @seealso \code{\link{PlotPairDimplot}} for the raw-only 3-panel variant,
#'   \code{\link{PlotPairSummary}} for a comprehensive multi-evidence summary.
#'
#' @export
#'
#' @examples
#' # scpairs_testdata has PCA (smooth_reduction) and UMAP (reduction) ready.
#' PlotPairSmoothed(scpairs_testdata, gene1 = "GENE3", gene2 = "GENE4")
#'
PlotPairSmoothed <- function(object,
                             gene1,
                             gene2,
                             reduction        = "umap",
                             smooth_reduction = "pca",
                             k                = 20,
                             alpha            = 0.3,
                             assay            = NULL,
                             slot             = "data",
                             pt_size          = 0.3,
                             pt_alpha         = 0.8,
                             title            = NULL) {

  validated <- .validate_plot_inputs(object, gene1, gene2,
                                      assay = assay, slot = slot,
                                      reduction = reduction)
  assay <- validated$assay

  embed <- Seurat::Embeddings(object, reduction = reduction)
  mat <- .prepare_expression_data(object, gene1, gene2,
                                   assay = assay, slot = slot)

  # Build KNN graph and smooth
  W <- .build_knn_graph(object, reduction = smooth_reduction, k = k)
  smoothed <- .smooth_expression(mat, W, alpha = alpha)
  if (inherits(smoothed, "sparseMatrix")) smoothed <- as.matrix(smoothed)

  common <- Reduce(intersect, list(rownames(embed), colnames(mat), colnames(smoothed)))
  embed <- embed[common, ]
  dim_names <- colnames(embed)[1:2]

  df <- data.frame(
    dim1       = embed[, 1],
    dim2       = embed[, 2],
    raw1       = mat[gene1, common],
    raw2       = mat[gene2, common],
    smooth1    = smoothed[gene1, common],
    smooth2    = smoothed[gene2, common],
    stringsAsFactors = FALSE
  )
  df$raw_product    <- df$raw1 * df$raw2
  df$smooth_product <- df$smooth1 * df$smooth2

  # Compute correlations for annotation
  raw_r <- round(stats::cor(df$raw1, df$raw2), 3)
  smooth_r <- round(stats::cor(df$smooth1, df$smooth2), 3)

  # Common theme
  base_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank(),
      axis.title  = ggplot2::element_blank(),
      plot.title  = ggplot2::element_text(size = 10, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8),
      legend.key.width = ggplot2::unit(0.3, "cm"),
      legend.key.height = ggplot2::unit(0.4, "cm"),
      legend.title = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 6)
    )

  .make_panel <- function(data, col, scale_opt, panel_title, subtitle = NULL) {
    data <- data[order(data[[col]]), ]
    mapping <- ggplot2::aes(x = dim1, y = dim2)
    mapping[["colour"]] <- as.symbol(col)
    p <- ggplot2::ggplot(data, mapping) +
      ggplot2::geom_point(size = pt_size, alpha = pt_alpha) +
      ggplot2::scale_color_viridis_c(option = scale_opt, name = "Expr") +
      ggplot2::ggtitle(panel_title, subtitle = subtitle) +
      base_theme
    p
  }

  # Row 1: Raw expression
  p1 <- .make_panel(df, "raw1", "D", gene1, "Raw")
  p2 <- .make_panel(df, "raw2", "C", gene2, "Raw")
  p3 <- .make_panel(df, "raw_product", "B",
                    paste0(gene1, " \u00D7 ", gene2),
                    paste0("Raw (r = ", raw_r, ")"))

  # Row 2: Smoothed expression
  p4 <- .make_panel(df, "smooth1", "D", gene1, "KNN-smoothed")
  p5 <- .make_panel(df, "smooth2", "C", gene2, "KNN-smoothed")
  p6 <- .make_panel(df, "smooth_product", "B",
                    paste0(gene1, " \u00D7 ", gene2),
                    paste0("Smoothed (r = ", smooth_r, ")"))

  if (is.null(title)) {
    title <- paste0("Synergy assessment: ", gene1, " & ", gene2)
  }

  patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3) +
    patchwork::plot_annotation(
      title    = title,
      subtitle = "Top: raw expression | Bottom: KNN-smoothed expression"
    )
}


#' Comprehensive Synergy Summary Plot
#'
#' @description
#' A multi-panel publication-ready figure combining:
#' 1. Raw UMAP co-expression (3 panels)
#' 2. KNN-smoothed UMAP (3 panels)
#' 3. Per-cluster expression comparison
#' 4. Metric radar/bar chart
#'
#' @param object A Seurat object.
#' @param gene1 Character; first gene.
#' @param gene2 Character; second gene.
#' @param result Optional `scPairs_pair_result` from `AssessGenePair()`.
#'     If NULL, assessment is run internally.
#' @param reduction Character; reduction for plotting.
#' @param smooth_reduction Character; reduction for KNN graph.
#' @param k Integer; KNN k for smoothing.
#' @param alpha Numeric; smoothing alpha.
#' @param assay Character; assay.
#' @param slot Character; data slot.
#' @param pt_size Numeric; point size.
#'
#' @return A combined \code{ggplot} (patchwork) with up to 10 panels:
#'   raw UMAP co-expression (3 panels), KNN-smoothed UMAP (3 panels),
#'   per-cluster expression bar chart, and metric evidence bar chart.
#'
#' @family Section_2_Visualization
#'
#' @seealso \code{\link{PlotPairSmoothed}}, \code{\link{AssessGenePair}}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' PlotPairSummary(scpairs_testdata, gene1 = "GENE3", gene2 = "GENE4")
#' }
PlotPairSummary <- function(object,
                            gene1,
                            gene2,
                            result           = NULL,
                            reduction        = "umap",
                            smooth_reduction = "pca",
                            k                = 20,
                            alpha            = 0.3,
                            assay            = NULL,
                            slot             = "data",
                            pt_size          = 0.3) {

  validated <- .validate_plot_inputs(object, gene1, gene2,
                                      assay = assay, slot = slot,
                                      reduction = reduction)
  assay <- validated$assay

  # Run assessment if not provided
  if (is.null(result)) {
    result <- AssessGenePair(object, gene1 = gene1, gene2 = gene2,
                              assay = assay, slot = slot,
                              use_neighbourhood = TRUE,
                              neighbourhood_k = k,
                              neighbourhood_reduction = smooth_reduction,
                              smooth_alpha = alpha,
                              n_perm = 0, verbose = FALSE)
  }

  # --- Panel A: Smoothed UMAP (top row) ---
  p_umap <- PlotPairSmoothed(object, gene1, gene2,
                              reduction = reduction,
                              smooth_reduction = smooth_reduction,
                              k = k, alpha = alpha,
                              assay = assay, slot = slot,
                              pt_size = pt_size, title = "")

  # --- Panel B: Per-cluster bar chart ---
  mat <- .get_expression_matrix(object, features = c(gene1, gene2),
                                 assay = assay, slot = slot)
  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }

  clusters <- Seurat::Idents(object)
  cluster_df <- data.frame(
    cluster = clusters,
    expr1 = mat[gene1, names(clusters)],
    expr2 = mat[gene2, names(clusters)],
    stringsAsFactors = FALSE
  )

  cl_means <- stats::aggregate(cbind(expr1, expr2) ~ cluster, data = cluster_df, FUN = mean)
  colnames(cl_means) <- c("cluster", gene1, gene2)
  cl_long <- tidyr::pivot_longer(cl_means, cols = -cluster,
                                  names_to = "gene", values_to = "expression")

  p_bar <- ggplot2::ggplot(cl_long,
    ggplot2::aes(x = cluster, y = expression, fill = gene)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.7, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c("#4DAF4A", "#984EA3"),
                                name = "Gene") +
    ggplot2::labs(x = NULL, y = "Mean expression") +
    ggplot2::ggtitle("Cluster-level expression") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      plot.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "bottom"
    )

  # --- Panel C: Metric comparison ---
  m <- result$metrics
  metric_df <- data.frame(
    metric = c("Pearson", "Spearman", "Biweight",
               "MI", "Ratio\nconsist.",
               if (!is.null(m$smoothed_cor)) "Smoothed\ncor" else NULL,
               if (!is.null(m$neighbourhood_score)) "Neigh.\nscore" else NULL,
               if (!is.null(m$cluster_cor)) "Cluster\ncor" else NULL),
    metric_value = c(
      abs(m$cor_pearson), abs(m$cor_spearman), abs(m$cor_biweight),
      min(m$mi_score / log(5), 1), m$ratio_consistency,
      if (!is.null(m$smoothed_cor)) abs(m$smoothed_cor) else NULL,
      if (!is.null(m$neighbourhood_score)) min(pmax(m$neighbourhood_score, 0) / 2, 1) else NULL,
      if (!is.null(m$cluster_cor) && !is.na(m$cluster_cor)) abs(m$cluster_cor) else NULL),
    stringsAsFactors = FALSE
  )

  # Colour cell-level vs neighbourhood metrics differently
  n_cell <- 5L
  n_neigh <- nrow(metric_df) - n_cell
  metric_df$type <- c(rep("Cell-level", n_cell),
                      rep("Neighbourhood", n_neigh))

  metric_df$metric <- factor(metric_df$metric, levels = metric_df$metric)

  p_metrics <- ggplot2::ggplot(metric_df,
    ggplot2::aes(x = metric, y = metric_value, fill = type)) +
    ggplot2::geom_col(width = 0.7, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c("Cell-level" = "#377EB8",
                                           "Neighbourhood" = "#E41A1C"),
                                name = "Metric type") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(x = NULL, y = "Normalised score") +
    ggplot2::ggtitle("Metric comparison",
                      subtitle = paste0("Synergy = ", round(result$synergy_score, 3),
                                        " | Confidence: ", result$confidence)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "bottom"
    )

  # --- Compose final figure ---
  bottom <- patchwork::wrap_plots(p_bar, p_metrics, ncol = 2)

  patchwork::wrap_plots(p_umap, bottom, ncol = 1, heights = c(2, 1)) +
    patchwork::plot_annotation(
      title = paste0("Synergistic gene pair: ", gene1, " & ", gene2),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold")
      )
    )
}
