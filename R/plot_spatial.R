#' Plot Spatial Co-Expression Map
#'
#' @description
#' For spatial transcriptomics data, visualises the spatial distribution of
#' two genes and their co-expression product on the tissue.  Three panels
#' are shown side by side:
#'
#' 1. Expression of gene 1.
#' 2. Expression of gene 2.
#' 3. Co-expression product (gene1 * gene2), highlighting spots where both
#'    genes are simultaneously active.
#'
#' @param object A Seurat object with spatial coordinates.
#' @param gene1 Character; first gene.
#' @param gene2 Character; second gene.
#' @param assay Character; assay.
#' @param slot Character; data slot.
#' @param pt_size Numeric; point size.
#' @param alpha Numeric; point alpha.
#' @param title Character; overall title.
#'
#' @return A combined `ggplot` (via patchwork).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' PlotPairSpatial(spatial_obj, gene1 = "CD8A", gene2 = "CD8B")
#' }
#'
PlotPairSpatial <- function(object,
                            gene1,
                            gene2,
                            assay    = NULL,
                            slot     = "data",
                            pt_size  = 1.2,
                            alpha    = 0.8,
                            title    = NULL) {

  .validate_seurat(object, require_spatial = TRUE)
  validated <- .validate_plot_inputs(object, gene1, gene2,
                                      assay = assay, slot = slot)
  assay <- validated$assay
  mat <- .prepare_expression_data(object, gene1, gene2,
                                   assay = assay, slot = slot)

  coords <- .get_spatial_coords(object)
  common <- intersect(rownames(coords), colnames(mat))
  coords <- coords[common, ]

  df <- data.frame(
    x = coords$x,
    y = coords$y,
    expr1 = mat[gene1, common],
    expr2 = mat[gene2, common],
    stringsAsFactors = FALSE
  )
  df$coexpr_product <- df$expr1 * df$expr2

  # Panel 1: gene1
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = expr1)) +
    ggplot2::geom_point(size = pt_size, alpha = alpha) +
    ggplot2::scale_color_viridis_c(option = "D", name = "Expr") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(gene1) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  # Panel 2: gene2
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = expr2)) +
    ggplot2::geom_point(size = pt_size, alpha = alpha) +
    ggplot2::scale_color_viridis_c(option = "C", name = "Expr") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(gene2) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  # Panel 3: co-expression
  p3 <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = coexpr_product)) +
    ggplot2::geom_point(size = pt_size, alpha = alpha) +
    ggplot2::scale_color_viridis_c(option = "B", name = "Co-expr") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste0(gene1, " * ", gene2)) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  if (is.null(title)) title <- paste0("Spatial co-expression: ", gene1, " & ", gene2)

  patchwork::wrap_plots(p1, p2, p3, ncol = 3) +
    patchwork::plot_annotation(title = title)
}
