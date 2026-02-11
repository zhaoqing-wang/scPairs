#' Plot Synergy Score Heatmap
#'
#' @description
#' Displays a symmetric heatmap of synergy scores among a set of genes.
#' Useful for visualising the overall co-expression landscape of top
#' synergistic genes or genes of interest.
#'
#' @param result An `scPairs_result` or `scPairs_gene_result` object, or a
#'     `data.frame` with `gene1`, `gene2`, `synergy_score`.
#' @param top_n Integer; include the top N genes by number of significant
#'     partnerships.  Default 30.
#' @param genes Character vector; specific genes to include.  NULL = auto.
#' @param cluster_genes Logical; cluster rows/columns by score similarity.
#'     Default TRUE.
#' @param low_color Character; colour for low scores.
#' @param high_color Character; colour for high scores.
#' @param title Character; plot title.
#'
#' @return A `ggplot` object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- FindAllPairs(seurat_obj, top_n = 500)
#' PlotPairHeatmap(result, top_n = 25)
#' }
#'
PlotPairHeatmap <- function(result,
                            top_n      = 30,
                            genes      = NULL,
                            cluster_genes = TRUE,
                            low_color  = "#F7FBFF",
                            high_color = "#08306B",
                            title      = "Gene pair synergy heatmap") {

  # Extract pairs (supports all three result classes + data.frame)
  edges <- .extract_pairs_df(result)

  # Select genes
  if (is.null(genes)) {
    # Top genes by degree (number of partnerships)
    all_genes <- c(edges$gene1, edges$gene2)
    gene_freq <- sort(table(all_genes), decreasing = TRUE)
    genes <- names(gene_freq)[seq_len(min(top_n, length(gene_freq)))]
  }

  # Filter edges to selected genes
  edges <- edges[edges$gene1 %in% genes & edges$gene2 %in% genes, ]
  if (nrow(edges) == 0) {
    message("No edges among selected genes.")
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::ggtitle("No data"))
  }

  # Build symmetric score matrix
  all_g <- sort(unique(c(edges$gene1, edges$gene2)))
  mat <- matrix(0, nrow = length(all_g), ncol = length(all_g),
                dimnames = list(all_g, all_g))

  for (i in seq_len(nrow(edges))) {
    g1 <- edges$gene1[i]
    g2 <- edges$gene2[i]
    s  <- edges$synergy_score[i]
    mat[g1, g2] <- s
    mat[g2, g1] <- s
  }
  diag(mat) <- 1

  # Optional clustering
  if (cluster_genes && nrow(mat) > 2) {
    hc <- stats::hclust(stats::dist(mat))
    ord <- hc$order
    mat <- mat[ord, ord]
  }

  # Convert to long format for ggplot
  df <- expand.grid(gene1 = rownames(mat), gene2 = colnames(mat),
                    stringsAsFactors = FALSE)
  df$synergy_score <- vapply(seq_len(nrow(df)), function(i) {
    mat[df$gene1[i], df$gene2[i]]
  }, numeric(1))

  df$gene1 <- factor(df$gene1, levels = rownames(mat))
  df$gene2 <- factor(df$gene2, levels = colnames(mat))

  ggplot2::ggplot(df, ggplot2::aes(x = gene1, y = gene2, fill = synergy_score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient(low = low_color, high = high_color,
                                  name = "Synergy\nscore") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed()
}
