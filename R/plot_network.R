#' Plot Gene Interaction Network
#'
#' @description
#' Draws a publication-ready gene interaction network from scPairs results.
#' Nodes represent genes; edges represent synergistic relationships.  Edge
#' width encodes synergy score; edge colour encodes confidence.  Node size
#' optionally reflects the number of significant partners (degree centrality).
#'
#' @param result An object of class `"scPairs_result"`, `"scPairs_gene_result"`,
#'     or a `data.frame` / `data.table` with columns `gene1`, `gene2`,
#'     `synergy_score`.
#' @param top_n Integer; show only the top N edges.  Default 50.
#' @param min_score Numeric; minimum synergy score to display an edge.
#' @param confidence Character vector; filter to these confidence levels
#'     (e.g., `c("High", "Medium")`).  NULL = no filter.
#' @param layout Character; ggraph layout algorithm.  Default `"fr"`
#'     (Fruchterman-Reingold).
#' @param node_color Character; colour for nodes.  Default `"#2C3E50"`.
#' @param edge_palette Character vector of 3 colours for confidence
#'     (High, Medium, Low).  Default blue-orange-grey scheme.
#' @param label_size Numeric; node label font size.
#' @param title Character; plot title.
#' @param show_legend Logical.
#'
#' @return A `ggplot` object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- FindAllPairs(seurat_obj, top_n = 200)
#' PlotPairNetwork(result, top_n = 50)
#'
#' # For a gene query
#' gene_result <- FindGenePairs(seurat_obj, gene = "CD3D")
#' PlotPairNetwork(gene_result, top_n = 20)
#' }
#'
PlotPairNetwork <- function(result,
                            top_n        = 50,
                            min_score    = 0,
                            confidence   = NULL,
                            layout       = "fr",
                            node_color   = "#2C3E50",
                            edge_palette = c("High" = "#E74C3C",
                                             "Medium" = "#F39C12",
                                             "Low" = "#95A5A6",
                                             "NS"  = "#D5D8DC"),
                            label_size   = 3.5,
                            title        = NULL,
                            show_legend  = TRUE) {

  # Extract pairs table
  if (inherits(result, "scPairs_result") || inherits(result, "scPairs_gene_result")) {
    edges <- as.data.frame(result$pairs)
  } else if (is.data.frame(result)) {
    edges <- result
  } else {
    stop("Input must be an scPairs result object or a data.frame.", call. = FALSE)
  }

  required <- c("gene1", "gene2", "synergy_score")
  if (!all(required %in% colnames(edges))) {
    stop("Input must contain columns: ", paste(required, collapse = ", "),
         call. = FALSE)
  }

  # Filter
  edges <- edges[edges$synergy_score >= min_score, , drop = FALSE]
  if (!is.null(confidence) && "confidence" %in% colnames(edges)) {
    edges <- edges[edges$confidence %in% confidence, , drop = FALSE]
  }

  if (nrow(edges) == 0) {
    message("No edges to plot after filtering.")
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::ggtitle("No gene pairs to display"))
  }

  # Top N
  edges <- edges[order(-edges$synergy_score), ]
  if (!is.null(top_n) && top_n < nrow(edges)) {
    edges <- edges[seq_len(top_n), ]
  }

  # Add confidence if missing

  if (!("confidence" %in% colnames(edges))) {
    edges$confidence <- "Medium"
  }

  # Build graph
  graph <- tidygraph::as_tbl_graph(
    igraph::graph_from_data_frame(
      edges[, c("gene1", "gene2", "synergy_score", "confidence")],
      directed = FALSE
    )
  )

  # Default title
  if (is.null(title)) {
    if (inherits(result, "scPairs_gene_result")) {
      title <- paste0("Gene interaction network: ", result$query_gene)
    } else {
      title <- "Synergistic gene pair network"
    }
  }

  # Plot
  p <- ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(
        width = synergy_score,
        colour = confidence
      ),
      alpha = 0.7
    ) +
    ggraph::geom_node_point(
      size = 5,
      colour = node_color,
      alpha = 0.9
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      size = label_size,
      repel = TRUE
    ) +
    ggraph::scale_edge_width_continuous(
      range = c(0.3, 3),
      name = "Synergy score"
    ) +
    ggraph::scale_edge_colour_manual(
      values = edge_palette,
      name = "Confidence"
    ) +
    ggraph::theme_graph(base_family = "") +
    ggplot2::ggtitle(title)

  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  p
}
