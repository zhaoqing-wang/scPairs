#' Plot Gene Pair Co-Expression on UMAP / Dimensionality Reduction
#'
#' @description
#' Displays the co-expression of two genes on the UMAP (or other reduction)
#' embedding.  Three panels show: gene 1 expression, gene 2 expression, and
#' their element-wise product (co-expression intensity).  This allows visual
#' assessment of whether co-expressing cells cluster together.
#'
#' @param object A Seurat object with a dimensionality reduction.
#' @param gene1 Character; first gene.
#' @param gene2 Character; second gene.
#' @param reduction Character; reduction to use.  Default `"umap"`.
#' @param assay Character; assay.
#' @param slot Character; data slot.
#' @param pt_size Numeric; point size.
#' @param alpha Numeric; point alpha.
#' @param title Character; overall title.
#'
#' @return A combined `ggplot` (patchwork).
#'
#' @export
#'
#' @examples
#' \donttest{
#' counts <- matrix(rpois(600, 5), nrow = 20, ncol = 30,
#'   dimnames = list(paste0("Gene", 1:20), paste0("Cell", 1:30)))
#' obj <- Seurat::CreateSeuratObject(counts = counts)
#' obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#' obj[["umap"]] <- Seurat::CreateDimReducObject(
#'   embeddings = matrix(rnorm(60), ncol = 2,
#'     dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2"))),
#'   key = "UMAP_")
#'
#' PlotPairDimplot(obj, gene1 = "Gene1", gene2 = "Gene2")
#' }
#'
PlotPairDimplot <- function(object,
                            gene1,
                            gene2,
                            reduction = "umap",
                            assay     = NULL,
                            slot      = "data",
                            pt_size   = 0.5,
                            alpha     = 0.8,
                            title     = NULL) {

  validated <- .validate_plot_inputs(object, gene1, gene2,
                                      assay = assay, slot = slot,
                                      reduction = reduction)
  assay <- validated$assay

  embed <- Seurat::Embeddings(object, reduction = reduction)
  mat <- .prepare_expression_data(object, gene1, gene2,
                                   assay = assay, slot = slot)

  common <- intersect(rownames(embed), colnames(mat))
  embed <- embed[common, ]

  df <- data.frame(
    dim1 = embed[, 1],
    dim2 = embed[, 2],
    expr1 = mat[gene1, common],
    expr2 = mat[gene2, common],
    stringsAsFactors = FALSE
  )
  df$coexpr_product <- df$expr1 * df$expr2

  dim_names <- colnames(embed)[1:2]

  # Sort so high-expression cells are plotted on top
  df <- df[order(df$expr1), ]
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, colour = expr1)) +
    ggplot2::geom_point(size = pt_size, alpha = alpha) +
    ggplot2::scale_color_viridis_c(option = "D", name = "Expr") +
    ggplot2::labs(x = dim_names[1], y = dim_names[2]) +
    ggplot2::ggtitle(gene1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  df <- df[order(df$expr2), ]
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, colour = expr2)) +
    ggplot2::geom_point(size = pt_size, alpha = alpha) +
    ggplot2::scale_color_viridis_c(option = "C", name = "Expr") +
    ggplot2::labs(x = dim_names[1], y = dim_names[2]) +
    ggplot2::ggtitle(gene2) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  df <- df[order(df$coexpr_product), ]
  p3 <- ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, colour = coexpr_product)) +
    ggplot2::geom_point(size = pt_size, alpha = alpha) +
    ggplot2::scale_color_viridis_c(option = "B", name = "Co-expr") +
    ggplot2::labs(x = dim_names[1], y = dim_names[2]) +
    ggplot2::ggtitle(paste0(gene1, " * ", gene2)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  if (is.null(title)) title <- paste0("Co-expression: ", gene1, " & ", gene2)

  patchwork::wrap_plots(p1, p2, p3, ncol = 3) +
    patchwork::plot_annotation(title = title)
}


#' Violin Plot of Pair Expression Across Clusters
#'
#' @description
#' Displays side-by-side violin plots of two genes across cell clusters or
#' groups, enabling visual assessment of whether their expression patterns
#' are coordinated across populations.
#'
#' @param object Seurat object.
#' @param gene1 Character; first gene.
#' @param gene2 Character; second gene.
#' @param group_by Character; column in meta.data for grouping.  NULL = Idents.
#' @param assay Character; assay.
#' @param slot Character; data slot.
#' @param pt_size Point size for jitter (0 = no points).
#' @param title Character.
#'
#' @return A `ggplot`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' counts <- matrix(rpois(600, 5), nrow = 20, ncol = 30,
#'   dimnames = list(paste0("Gene", 1:20), paste0("Cell", 1:30)))
#' obj <- Seurat::CreateSeuratObject(counts = counts)
#' obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#' obj$group <- factor(sample(c("A", "B"), 30, replace = TRUE))
#'
#' PlotPairViolin(obj, "Gene1", "Gene2", group_by = "group")
#' }
#'
PlotPairViolin <- function(object,
                           gene1,
                           gene2,
                           group_by = NULL,
                           assay    = NULL,
                           slot     = "data",
                           pt_size  = 0,
                           title    = NULL) {

  validated <- .validate_plot_inputs(object, gene1, gene2,
                                      assay = assay, slot = slot)
  assay <- validated$assay
  mat <- .prepare_expression_data(object, gene1, gene2,
                                   assay = assay, slot = slot)

  if (!is.null(group_by)) {
    groups <- object@meta.data[[group_by]]
  } else {
    groups <- Seurat::Idents(object)
  }

  df <- data.frame(
    cluster = as.factor(groups),
    expr1   = mat[gene1, ],
    expr2   = mat[gene2, ],
    stringsAsFactors = FALSE
  )

  df_long <- tidyr::pivot_longer(df, cols = c("expr1", "expr2"),
                                  names_to = "gene", values_to = "expression")
  df_long$gene <- ifelse(df_long$gene == "expr1", gene1, gene2)

  if (is.null(title)) title <- paste0(gene1, " & ", gene2, " across clusters")

  ggplot2::ggplot(df_long, ggplot2::aes(x = cluster, y = expression,
                                         fill = gene)) +
    ggplot2::geom_violin(scale = "width", alpha = 0.7, position = ggplot2::position_dodge(0.8)) +
    {if (pt_size > 0) ggplot2::geom_jitter(size = pt_size, alpha = 0.3,
                                             position = ggplot2::position_jitterdodge(
                                               dodge.width = 0.8, jitter.width = 0.15))} +
    ggplot2::labs(x = "Cluster", y = "Expression", fill = "Gene") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}


#' Scatter Plot of Two Genes (Cell-Level)
#'
#' @description
#' Plots cell-level expression of gene1 vs. gene2 as a scatter plot, coloured
#' by cluster identity.  Marginal density curves (optional) help reveal
#' cluster-specific co-expression patterns.
#'
#' @param object Seurat object.
#' @param gene1 Character; x-axis gene.
#' @param gene2 Character; y-axis gene.
#' @param group_by Character; colour cells by this meta.data column.
#' @param assay Character; assay.
#' @param slot Character; data slot.
#' @param pt_size Numeric.
#' @param alpha Numeric.
#' @param add_density Logical; add marginal density.  Requires `ggExtra`
#'     package.
#' @param title Character.
#'
#' @return A `ggplot`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' counts <- matrix(rpois(600, 5), nrow = 20, ncol = 30,
#'   dimnames = list(paste0("Gene", 1:20), paste0("Cell", 1:30)))
#' obj <- Seurat::CreateSeuratObject(counts = counts)
#' obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#' obj$group <- factor(sample(c("A", "B"), 30, replace = TRUE))
#'
#' PlotPairScatter(obj, "Gene1", "Gene2", group_by = "group")
#' }
#'
PlotPairScatter <- function(object,
                            gene1,
                            gene2,
                            group_by    = NULL,
                            assay       = NULL,
                            slot        = "data",
                            pt_size     = 0.5,
                            alpha       = 0.6,
                            add_density = FALSE,
                            title       = NULL) {

  validated <- .validate_plot_inputs(object, gene1, gene2,
                                      assay = assay, slot = slot)
  assay <- validated$assay
  mat <- .prepare_expression_data(object, gene1, gene2,
                                   assay = assay, slot = slot)

  if (!is.null(group_by)) {
    groups <- as.factor(object@meta.data[[group_by]])
  } else {
    groups <- as.factor(Seurat::Idents(object))
  }

  df <- data.frame(
    expr1 = mat[gene1, ],
    expr2 = mat[gene2, ],
    cluster = groups,
    stringsAsFactors = FALSE
  )

  # Correlation annotation
  r <- round(stats::cor(df$expr1, df$expr2, method = "pearson"), 3)

  if (is.null(title)) title <- paste0(gene1, " vs ", gene2, "  (r = ", r, ")")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = expr1, y = expr2, colour = cluster)) +
    ggplot2::geom_point(size = pt_size, alpha = alpha) +
    ggplot2::labs(x = gene1, y = gene2, colour = "Cluster") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw()

  if (add_density && requireNamespace("ggExtra", quietly = TRUE)) {
    p <- ggExtra::ggMarginal(p, type = "density", groupFill = TRUE)
  }

  p
}
