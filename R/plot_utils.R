#' Shared Plotting Utilities for scPairs
#'
#' @description
#' Internal functions providing common validation, data preparation, and
#' theming logic used across all plot_*.R files.
#'
#' @name plot-utils
#' @keywords internal
NULL


#' Validate common plot inputs
#'
#' Checks that gene names exist and reduction is available.
#'
#' @param object Seurat object.
#' @param gene1,gene2 Gene names.
#' @param assay Assay name (NULL = default).
#' @param slot Data slot.
#' @param reduction Reduction name (NULL to skip check).
#' @return A list with validated \code{assay} and \code{mat} (dense matrix).
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
