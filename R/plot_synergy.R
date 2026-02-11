#' Visualize Synergistic Relationship Between Gene Pairs
#'
#' @description
#' Publication-ready multi-panel visualization that integrates prior knowledge,
#' expression evidence, and neighbourhood context to show the synergistic
#' relationship between two genes.  This goes beyond co-expression to reveal
#' *why* two genes may be functionally synergistic.
#'
#' @param object A Seurat object.
#' @param gene1 Character; first gene.
#' @param gene2 Character; second gene.
#' @param prior_net Optional prior network from \code{.build_prior_network()}.
#'     If NULL, built automatically.
#' @param organism Character; "mouse" or "human".  Used if prior_net is NULL.
#' @param reduction Character; reduction for UMAP plotting.
#' @param smooth_reduction Character; reduction for KNN.
#' @param k Integer; KNN k.
#' @param alpha Numeric; smoothing alpha.
#' @param cluster_col Character; cluster column in meta.data.
#' @param assay Character; assay.
#' @param slot Character; data slot.
#' @param top_bridges Integer; maximum bridge genes to show.
#' @param pt_size Numeric; point size.
#'
#' @return A combined ggplot (patchwork) with panels showing:
#'   1. UMAP with neighbourhood synergy highlighting
#'   2. Bridge gene network
#'   3. Per-cluster expression evidence
#'   4. Metric comparison (expression + prior)
#'
#' @export
PlotPairSynergy <- function(object,
                            gene1,
                            gene2,
                            prior_net         = NULL,
                            organism          = "mouse",
                            reduction         = "umap",
                            smooth_reduction  = "pca",
                            k                 = 20,
                            alpha             = 0.3,
                            cluster_col       = NULL,
                            assay             = NULL,
                            slot              = "data",
                            top_bridges       = 10,
                            pt_size           = 0.3) {

  validated <- .validate_plot_inputs(object, gene1, gene2,
                                      assay = assay, slot = slot,
                                      reduction = reduction)
  assay <- validated$assay

  # Build prior network if not provided
  all_genes <- rownames(Seurat::GetAssayData(object, assay = assay, layer = slot))
  if (is.null(prior_net)) {
    prior_net <- tryCatch(
      .build_prior_network(organism = organism, genes = all_genes, verbose = FALSE),
      error = function(e) {
        message("[scPairs] Could not build prior network: ", conditionMessage(e))
        NULL
      }
    )
  }

  # Get expression data
  embed <- Seurat::Embeddings(object, reduction = reduction)
  mat <- .get_expression_matrix(object, features = NULL, assay = assay, slot = slot)

  # Get expression for all genes (needed for bridge checking)
  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat_dense_pair <- as.matrix(mat[c(gene1, gene2), , drop = FALSE])
  } else {
    mat_dense_pair <- mat[c(gene1, gene2), , drop = FALSE]
  }

  # Expressed genes (>5% cells expressing)
  expressed_genes <- rownames(mat)[Matrix::rowMeans(mat > 0) > 0.05]

  # Build KNN graph
  W <- .build_knn_graph(object, reduction = smooth_reduction, k = k)

  # Cluster info
  if (!is.null(cluster_col)) {
    clusters <- as.factor(object@meta.data[[cluster_col]])
  } else {
    clusters <- as.factor(Seurat::Idents(object))
  }

  common <- intersect(rownames(embed), colnames(mat))
  embed <- embed[common, ]
  dim_names <- colnames(embed)[1:2]

  x <- mat_dense_pair[gene1, common]
  y <- mat_dense_pair[gene2, common]

  # --- Panel 1: Neighbourhood synergy UMAP ---
  # Color cells by: gene A expression * neighbourhood gene B expression
  neigh_y <- as.numeric(W[common, common] %*% y)
  neigh_x <- as.numeric(W[common, common] %*% x)

  df_umap <- data.frame(
    dim1 = embed[, 1],
    dim2 = embed[, 2],
    expr1 = x,
    expr2 = y,
    neigh_synergy_ab = x * neigh_y,  # A expressed here, B in neighbours
    neigh_synergy_ba = y * neigh_x,  # B expressed here, A in neighbours
    cluster = clusters[common],
    stringsAsFactors = FALSE
  )
  df_umap$synergy_signal <- df_umap$neigh_synergy_ab + df_umap$neigh_synergy_ba

  base_theme <- .build_panel_theme()

  # UMAP panel: gene1 expression
  df_sorted <- df_umap[order(df_umap$expr1), ]
  p_g1 <- ggplot2::ggplot(df_sorted,
    ggplot2::aes(x = dim1, y = dim2, colour = expr1)) +
    ggplot2::geom_point(size = pt_size, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(option = "D", name = "Expr") +
    ggplot2::ggtitle(gene1) +
    base_theme

  # UMAP panel: gene2 expression
  df_sorted <- df_umap[order(df_umap$expr2), ]
  p_g2 <- ggplot2::ggplot(df_sorted,
    ggplot2::aes(x = dim1, y = dim2, colour = expr2)) +
    ggplot2::geom_point(size = pt_size, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(option = "C", name = "Expr") +
    ggplot2::ggtitle(gene2) +
    base_theme

  # UMAP panel: neighbourhood synergy signal
  df_sorted <- df_umap[order(df_umap$synergy_signal), ]
  p_syn <- ggplot2::ggplot(df_sorted,
    ggplot2::aes(x = dim1, y = dim2, colour = synergy_signal)) +
    ggplot2::geom_point(size = pt_size, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(option = "inferno", name = "Signal") +
    ggplot2::ggtitle("Neighbourhood synergy",
                      subtitle = paste0(gene1, " \u00D7 neigh(", gene2, ") + ",
                                        gene2, " \u00D7 neigh(", gene1, ")")) +
    base_theme

  row1 <- patchwork::wrap_plots(p_g1, p_g2, p_syn, ncol = 3)

  # --- Panel 2: Bridge gene network ---
  bridge_result <- NULL
  p_network <- NULL

  if (!is.null(prior_net) && prior_net$n_terms > 0) {
    bridge_result <- .bridge_score(gene1, gene2, prior_net, expressed_genes,
                                    top_n = top_bridges)

    p_network <- .plot_bridge_network(
      gene1, gene2,
      bridge_genes = bridge_result$bridges,
      shared_terms = bridge_result$shared_terms,
      mat = mat[, common],
      prior_net = prior_net,
      top_n = top_bridges
    )
  }

  if (is.null(p_network)) {
    # Fallback: simple scatter of co-expression
    p_network <- ggplot2::ggplot(df_umap,
      ggplot2::aes(x = expr1, y = expr2, colour = cluster)) +
      ggplot2::geom_point(size = 0.8, alpha = 0.5) +
      ggplot2::labs(x = gene1, y = gene2, colour = "Cluster") +
      ggplot2::ggtitle("Gene pair scatter") +
      ggplot2::theme_bw()
  }

  # --- Panel 3: Per-cluster expression comparison ---
  cl_df <- data.frame(
    cl = as.character(clusters[common]),
    gA = x,
    gB = y,
    stringsAsFactors = FALSE
  )
  cl_means <- stats::aggregate(cbind(gA, gB) ~ cl, data = cl_df, FUN = mean)
  colnames(cl_means) <- c("cluster", gene1, gene2)
  cl_long <- tidyr::pivot_longer(cl_means, cols = -cluster,
                                  names_to = "gene", values_to = "expression")

  p_bar <- ggplot2::ggplot(cl_long,
    ggplot2::aes(x = cluster, y = expression, fill = gene)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.8),
                       width = 0.7, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c("#4DAF4A", "#984EA3"),
                                name = "Gene") +
    ggplot2::labs(x = NULL, y = "Mean expression") +
    ggplot2::ggtitle("Cluster expression") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      plot.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "bottom"
    )

  # --- Panel 4: Metric summary with prior knowledge ---
  # Compute metrics for annotation
  met <- list()
  met$cor_pearson <- stats::cor(x, y, method = "pearson")
  met$cor_spearman <- stats::cor(x, y, method = "spearman")

  # Neighbourhood synergy
  met$neigh_synergy <- .neighbourhood_synergy(x, y, W[common, common])

  # Prior and bridge scores
  if (!is.null(prior_net) && prior_net$n_terms > 0) {
    met$prior_score <- .prior_score(gene1, gene2, prior_net)
    met$bridge_score <- bridge_result$score
    n_bridges <- length(bridge_result$bridges)
  } else {
    met$prior_score <- NA
    met$bridge_score <- NA
    n_bridges <- 0
  }

  # Smoothed correlation
  smoothed <- .smooth_expression(mat[c(gene1, gene2), common], W[common, common],
                                  alpha = alpha)
  if (inherits(smoothed, "sparseMatrix")) smoothed <- as.matrix(smoothed)
  met$smoothed_cor <- stats::cor(smoothed[gene1, ], smoothed[gene2, ])

  # Build metric bar chart
  metric_df <- data.frame(
    metric = c("Pearson", "Spearman", "Smoothed\ncor",
               "Neigh.\nsynergy",
               if (!is.na(met$prior_score)) "Prior\nknowledge" else NULL,
               if (!is.na(met$bridge_score)) "Bridge\nscore" else NULL),
    metric_value = c(
      abs(met$cor_pearson), abs(met$cor_spearman),
      abs(met$smoothed_cor),
      min(pmax(met$neigh_synergy, 0) / 2, 1),
      if (!is.na(met$prior_score)) met$prior_score else NULL,
      if (!is.na(met$bridge_score)) met$bridge_score else NULL
    ),
    type = c(rep("Expression", 3), "Neighbourhood",
             if (!is.na(met$prior_score)) "Prior" else NULL,
             if (!is.na(met$bridge_score)) "Prior" else NULL),
    stringsAsFactors = FALSE
  )
  metric_df$metric <- factor(metric_df$metric, levels = metric_df$metric)

  p_metrics <- ggplot2::ggplot(metric_df,
    ggplot2::aes(x = metric, y = metric_value, fill = type)) +
    ggplot2::geom_col(width = 0.7, alpha = 0.85) +
    ggplot2::scale_fill_manual(
      values = c("Expression" = "#377EB8",
                 "Neighbourhood" = "#E41A1C",
                 "Prior" = "#FF7F00"),
      name = "Evidence") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(x = NULL, y = "Score") +
    ggplot2::ggtitle("Multi-evidence scores",
                      subtitle = paste0("Bridges: ", n_bridges,
                                        if (!is.na(met$prior_score))
                                          paste0(" | Prior: ", round(met$prior_score, 3))
                                        else "")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "bottom"
    )

  # --- Compose using design layout for reliable rendering ---
  design <- c(
    patchwork::area(t = 1, l = 1, b = 5, r = 4),   # p_g1
    patchwork::area(t = 1, l = 5, b = 5, r = 8),   # p_g2
    patchwork::area(t = 1, l = 9, b = 5, r = 12),  # p_syn
    patchwork::area(t = 6, l = 1, b = 12, r = 5),  # p_network (wider)
    patchwork::area(t = 6, l = 6, b = 12, r = 9),  # p_bar
    patchwork::area(t = 6, l = 10, b = 12, r = 12)  # p_metrics
  )

  p_g1 + p_g2 + p_syn + p_network + p_bar + p_metrics +
    patchwork::plot_layout(design = design) +
    patchwork::plot_annotation(
      title = paste0("Synergistic gene pair: ", gene1, " & ", gene2),
      subtitle = "Top: expression & neighbourhood synergy | Bottom: prior knowledge & cluster evidence",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10, colour = "grey40")
      )
    )
}


#' Plot bridge gene network connecting a gene pair
#'
#' @param gene1,gene2 The focal gene pair.
#' @param bridge_genes Character vector of bridge genes.
#' @param shared_terms Character vector of shared GO/KEGG terms.
#' @param mat Expression matrix for computing bridge expression.
#' @param prior_net Prior network object.
#' @param top_n Maximum bridges to display.
#' @return A ggplot object.
#' @keywords internal
.plot_bridge_network <- function(gene1, gene2,
                                 bridge_genes, shared_terms,
                                 mat, prior_net, top_n = 10) {

  if (length(bridge_genes) == 0 && length(shared_terms) == 0) {
    # No prior knowledge connections found
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                         label = paste0("No prior knowledge connections\n",
                                        "found between ", gene1, " & ", gene2),
                         size = 3.5, colour = "grey50") +
      ggplot2::theme_void() +
      ggplot2::ggtitle("Bridge network")
    return(p)
  }

  # Limit bridges
  if (length(bridge_genes) > top_n) {
    bridge_genes <- bridge_genes[seq_len(top_n)]
  }

  # Build edge list
  nodes <- unique(c(gene1, gene2, bridge_genes))
  edges <- data.frame(from = character(), to = character(),
                       edge_type = character(), stringsAsFactors = FALSE)

  g2t <- prior_net$gene_to_terms

  # Direct edge between focal pair if shared terms exist
  if (length(shared_terms) > 0) {
    edges <- rbind(edges, data.frame(
      from = gene1, to = gene2, edge_type = "direct",
      stringsAsFactors = FALSE))
  }

  # Bridge edges
  for (bg in bridge_genes) {
    bg_terms <- g2t[[bg]]
    terms1 <- g2t[[gene1]]
    terms2 <- g2t[[gene2]]

    if (length(intersect(bg_terms, terms1)) > 0) {
      edges <- rbind(edges, data.frame(
        from = gene1, to = bg, edge_type = "bridge",
        stringsAsFactors = FALSE))
    }
    if (length(intersect(bg_terms, terms2)) > 0) {
      edges <- rbind(edges, data.frame(
        from = bg, to = gene2, edge_type = "bridge",
        stringsAsFactors = FALSE))
    }
  }

  if (nrow(edges) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                         label = "No prior connections found",
                         size = 3.5, colour = "grey50") +
      ggplot2::theme_void() +
      ggplot2::ggtitle("Bridge network")
    return(p)
  }

  # Compute node expression (mean across cells)
  node_expr <- vapply(nodes, function(g) {
    if (g %in% rownames(mat)) {
      mean(mat[g, ])
    } else {
      0
    }
  }, numeric(1))

  # Node type
  node_type <- ifelse(nodes %in% c(gene1, gene2), "focal", "bridge")

  # Build igraph
  g <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                      vertices = data.frame(
                                        name = nodes,
                                        expression = node_expr,
                                        node_type = node_type
                                      ))

  tg <- tidygraph::as_tbl_graph(g)

  n_shared <- length(shared_terms)
  subtitle <- if (n_shared > 0) {
    paste0(n_shared, " shared terms, ", length(bridge_genes), " bridges")
  } else {
    paste0(length(bridge_genes), " bridge genes")
  }

  ggraph::ggraph(tg, layout = "stress") +
    ggraph::geom_edge_link(
      ggplot2::aes(edge_linetype = edge_type),
      alpha = 0.6, edge_width = 0.8, edge_colour = "grey50") +
    ggraph::geom_node_point(
      ggplot2::aes(size = expression, colour = node_type)) +
    ggraph::geom_node_text(ggplot2::aes(label = name),
                            repel = TRUE, size = 3) +
    ggplot2::scale_colour_manual(
      values = c("focal" = "#E41A1C", "bridge" = "#377EB8"),
      name = "Role") +
    ggplot2::scale_size_continuous(range = c(3, 8), name = "Mean expr") +
    ggraph::scale_edge_linetype_manual(
      values = c("direct" = "solid", "bridge" = "dashed"),
      name = "Connection") +
    ggplot2::ggtitle("Prior knowledge network", subtitle = subtitle) +
    ggraph::theme_graph() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 6),
      legend.title = ggplot2::element_text(size = 7),
      legend.key.size = ggplot2::unit(0.4, "cm")
    )
}
