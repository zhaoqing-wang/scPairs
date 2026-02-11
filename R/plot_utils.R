#' Shared Plotting Utilities for scPairs
#'
#' @description
#' Internal functions providing common validation, data preparation, and
#' theming logic used across all plot_*.R files.
#'
#' @name plot-utils
#' @keywords internal
NULL


#' Extract pairs data.frame from any scPairs result type
#'
#' All three result classes (`scPairs_result`, `scPairs_gene_result`,
#' `scPairs_pair_result`) and plain data.frames are supported.
#'
#' @param result An scPairs result object or data.frame.
#' @param require_cols Character vector of required columns.
#' @return A data.frame with the requested columns.
#' @keywords internal
.extract_pairs_df <- function(result,
                              require_cols = c("gene1", "gene2", "synergy_score")) {
  if (inherits(result, "scPairs_pair_result")) {
    # AssessGenePair now stores a single-row pair_dt in $pairs
    edges <- as.data.frame(result$pairs)
  } else if (inherits(result, "scPairs_result") ||
             inherits(result, "scPairs_gene_result")) {
    edges <- as.data.frame(result$pairs)
  } else if (is.data.frame(result)) {
    edges <- result
  } else {
    stop("Input must be an scPairs result object or a data.frame.", call. = FALSE)
  }

  missing <- setdiff(require_cols, colnames(edges))
  if (length(missing) > 0) {
    stop("Input must contain columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  edges
}


#' Validate common plot inputs
#'
#' Checks that gene names exist and reduction is available.
#'
#' @param object Seurat object.
#' @param gene1,gene2 Gene names.
#' @param assay Assay name (NULL = default).
#' @param slot Data slot.
#' @param reduction Reduction name (NULL to skip check).
#' @return A list with validated \code{assay}.
#' @keywords internal
.validate_plot_inputs <- function(object, gene1, gene2,
                                  assay = NULL, slot = "data",
                                  reduction = NULL) {
  .validate_seurat(object)
  assay <- assay %||% Seurat::DefaultAssay(object)
  .validate_features(c(gene1, gene2), object, assay)

  if (!is.null(reduction)) {
    embed <- tryCatch(
      Seurat::Embeddings(object, reduction = reduction),
      error = function(e) NULL
    )
    if (is.null(embed)) {
      stop(sprintf("Reduction '%s' not found in the Seurat object.", reduction),
           call. = FALSE)
    }
  }

  invisible(list(assay = assay))
}


#' Prepare expression data for plotting a gene pair
#'
#' Extracts expression, converts to dense, aligns with embeddings.
#'
#' @param object Seurat object.
#' @param gene1,gene2 Gene names.
#' @param assay Assay name.
#' @param slot Data slot.
#' @return Dense matrix (2 x n_cells).
#' @keywords internal
.prepare_expression_data <- function(object, gene1, gene2,
                                     assay = NULL, slot = "data") {
  assay <- assay %||% Seurat::DefaultAssay(object)
  mat <- .get_expression_matrix(object, features = c(gene1, gene2),
                                assay = assay, slot = slot)
  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }
  mat
}


#' Build a common minimal theme for panel plots
#'
#' @param axis_text_size Numeric; axis text size. NULL = hide axis text.
#' @param title_size Numeric; title text size.
#' @return A ggplot2 theme object.
#' @keywords internal
.build_panel_theme <- function(axis_text_size = NULL, title_size = 10) {
  th <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      plot.title  = ggplot2::element_text(size = title_size, face = "bold"),
      legend.key.width  = ggplot2::unit(0.3, "cm"),
      legend.key.height = ggplot2::unit(0.4, "cm"),
      legend.title = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 6)
    )

  if (is.null(axis_text_size)) {
    th <- th + ggplot2::theme(
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )
  } else {
    th <- th + ggplot2::theme(
      axis.text = ggplot2::element_text(size = axis_text_size)
    )
  }

  th
}
