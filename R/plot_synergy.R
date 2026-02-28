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
#' @return A combined \code{ggplot} (patchwork) with up to 4 panels:
#' \enumerate{
#'   \item UMAP coloured by per-cell neighbourhood synergy score.
#'   \item Bridge gene network showing shared GO/KEGG pathway intermediaries.
#'   \item Per-cluster expression bar chart for both genes.
#'   \item Multi-evidence metric comparison bar chart (expression + prior).
#' }
#' Falls back gracefully when prior knowledge is unavailable (panels 2 and 4
#' are omitted).
#'
#' @family Section_2_Visualization
#'
#' @seealso \code{\link{PlotBridgeNetwork}} for a standalone bridge gene
#'   network, \code{\link{AssessGenePair}} for the underlying metrics.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Requires Bioconductor annotation packages: org.Hs.eg.db or org.Mm.eg.db
#' # and AnnotationDbi.
#' PlotPairSynergy(scpairs_testdata, gene1 = "GENE3", gene2 = "GENE4",
#'                 organism = "human")
#' }
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
  neigh_y <- as.numeric(W[common, common] %*% y)
  neigh_x <- as.numeric(W[common, common] %*% x)

  df_umap <- data.frame(
    dim1 = embed[, 1],
    dim2 = embed[, 2],
    expr1 = x,
    expr2 = y,
    neigh_synergy_ab = x * neigh_y,
    neigh_synergy_ba = y * neigh_x,
    cluster = clusters[common],
    stringsAsFactors = FALSE
  )
  df_umap$synergy_signal <- df_umap$neigh_synergy_ab + df_umap$neigh_synergy_ba

  base_theme <- .build_panel_theme()

  df_sorted <- df_umap[order(df_umap$expr1), ]
  p_g1 <- ggplot2::ggplot(df_sorted,
    ggplot2::aes(x = dim1, y = dim2, colour = expr1)) +
    ggplot2::geom_point(size = pt_size, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(option = "D", name = "Expr") +
    ggplot2::ggtitle(gene1) +
    base_theme

  df_sorted <- df_umap[order(df_umap$expr2), ]
  p_g2 <- ggplot2::ggplot(df_sorted,
    ggplot2::aes(x = dim1, y = dim2, colour = expr2)) +
    ggplot2::geom_point(size = pt_size, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(option = "C", name = "Expr") +
    ggplot2::ggtitle(gene2) +
    base_theme

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

  # --- Panel 2: Enhanced bridge gene network ---
  bridge_result <- NULL
  p_network <- NULL

  if (!is.null(prior_net) && prior_net$n_terms > 0) {
    bridge_result <- .bridge_score(gene1, gene2, prior_net, expressed_genes,
                                    top_n = top_bridges)

    p_network <- .plot_bridge_network_enhanced(
      gene1, gene2,
      bridge_genes  = bridge_result$bridges,
      shared_terms  = bridge_result$shared_terms,
      mat           = mat[, common],
      prior_net     = prior_net,
      top_n         = top_bridges,
      layout        = "auto"
    )
  }

  if (is.null(p_network)) {
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
  met <- list()
  met$cor_pearson <- stats::cor(x, y, method = "pearson")
  met$cor_spearman <- stats::cor(x, y, method = "spearman")
  met$neigh_synergy <- .neighbourhood_synergy(x, y, W[common, common])

  if (!is.null(prior_net) && prior_net$n_terms > 0) {
    met$prior_score <- .prior_score(gene1, gene2, prior_net)
    met$bridge_score <- bridge_result$score
    n_bridges <- length(bridge_result$bridges)
  } else {
    met$prior_score <- NA
    met$bridge_score <- NA
    n_bridges <- 0
  }

  smoothed <- .smooth_expression(mat[c(gene1, gene2), common], W[common, common],
                                  alpha = alpha)
  if (inherits(smoothed, "sparseMatrix")) smoothed <- as.matrix(smoothed)
  met$smoothed_cor <- stats::cor(smoothed[gene1, ], smoothed[gene2, ])

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

  # --- Compose ---
  design <- c(
    patchwork::area(t = 1, l = 1, b = 5, r = 4),
    patchwork::area(t = 1, l = 5, b = 5, r = 8),
    patchwork::area(t = 1, l = 9, b = 5, r = 12),
    patchwork::area(t = 6, l = 1, b = 12, r = 5),
    patchwork::area(t = 6, l = 6, b = 12, r = 9),
    patchwork::area(t = 6, l = 10, b = 12, r = 12)
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


#' Plot Bridge Gene Network
#'
#' @description
#' Draws a publication-ready radial bridge gene network showing the
#' prior-knowledge connections between a focal gene pair via shared GO/KEGG
#' pathway intermediaries (bridge genes).  Focal genes are placed at the
#' centre; bridge genes are arranged on a ring whose radius is inversely
#' proportional to their shared pathway count with the focal pair (more
#' shared terms => closer to centre, reflecting stronger biological relevance).
#' Solid edges connect focal genes to bridge genes: red edges originate from
#' \code{gene1}, blue edges from \code{gene2}.  Edge width encodes shared
#' term count.  Pairwise Jaccard similarity between bridge genes is overlaid
#' as thin dotted lines, whose opacity reflects similarity strength,
#' revealing functional clusters among the intermediaries.
#'
#' @param object A Seurat object.
#' @param gene1 Character; first focal gene.
#' @param gene2 Character; second focal gene.
#' @param organism Character; \code{"mouse"} or \code{"human"}.
#'     Used when \code{prior_net} is NULL.
#' @param prior_net Optional prior network object from
#'     \code{.build_prior_network()}.  If NULL, built automatically using
#'     \code{organism}.
#' @param top_bridges Integer; maximum number of bridge genes to display.
#'     Default 15.
#' @param layout Ignored; layout is always radial (kept for API compatibility).
#' @param assay Character; assay name.  NULL uses the default assay.
#' @param slot Character; data layer/slot.  Default \code{"data"}.
#' @param label_size Numeric; gene label font size.  Default 3.
#' @param pt_size_range Numeric vector of length 2; minimum and maximum node
#'     sizes for bridge genes.  Default \code{c(3, 9)}.
#' @param edge_width_range Numeric vector of length 2; minimum and maximum
#'     edge widths for focal-to-bridge connections, scaled by shared term
#'     count.  Default \code{c(0.4, 2)}.
#' @param sim_threshold Numeric (0 to 1); minimum Jaccard similarity between
#'     two bridge genes required to draw a dotted similarity edge.
#'     Default 0.05.
#' @param title Character; plot title.  NULL generates a default title.
#'
#' @return A \code{ggplot} object.  Focal genes appear as large red nodes at
#'   the centre.  Bridge genes are arranged radially, sized by node degree and
#'   coloured by mean expression.  Solid coloured edges (red = gene1,
#'   blue = gene2) connect focal genes to bridge genes, with width proportional
#'   to shared term count.  Thin dotted grey lines between bridge genes encode
#'   Jaccard pathway similarity.
#'
#' @family Section_2_Visualization
#'
#' @seealso \code{\link{PlotPairSynergy}} for the 4-panel synergy dashboard
#'   that embeds this network as one panel, \code{\link{AssessGenePair}} for
#'   extracting bridge genes programmatically.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Requires Bioconductor annotation packages (org.Hs.eg.db or org.Mm.eg.db)
#' PlotBridgeNetwork(seurat_obj, gene1 = "Adora2a", gene2 = "Ido1",
#'                   organism = "mouse")
#' }
PlotBridgeNetwork <- function(object,
                               gene1,
                               gene2,
                               organism         = "mouse",
                               prior_net        = NULL,
                               top_bridges      = 15,
                               layout           = "auto",
                               assay            = NULL,
                               slot             = "data",
                               label_size       = 3,
                               pt_size_range    = c(3, 9),
                               edge_width_range = c(0.4, 2),
                               sim_threshold    = 0.05,
                               title            = NULL) {

  validated <- .validate_plot_inputs(object, gene1, gene2,
                                      assay = assay, slot = slot,
                                      reduction = NULL)
  assay <- validated$assay

  all_genes <- rownames(Seurat::GetAssayData(object, assay = assay, layer = slot))

  if (is.null(prior_net)) {
    prior_net <- tryCatch(
      .build_prior_network(organism = organism, genes = all_genes, verbose = FALSE),
      error = function(e) {
        stop("Could not build prior network: ", conditionMessage(e), call. = FALSE)
      }
    )
  }

  mat <- .get_expression_matrix(object, features = NULL, assay = assay, slot = slot)
  expressed_genes <- rownames(mat)[Matrix::rowMeans(mat > 0) > 0.05]

  bridge_result <- .bridge_score(gene1, gene2, prior_net, expressed_genes,
                                  top_n = top_bridges)

  if (is.null(title)) {
    title <- paste0("Bridge network: ", gene1, " & ", gene2)
  }

  .plot_bridge_network_enhanced(
    gene1, gene2,
    bridge_genes     = bridge_result$bridges,
    shared_terms     = bridge_result$shared_terms,
    mat              = mat,
    prior_net        = prior_net,
    top_n            = top_bridges,
    layout           = layout,
    label_size       = label_size,
    pt_size_range    = pt_size_range,
    edge_width_range = edge_width_range,
    title            = title,
    sim_threshold    = sim_threshold
  )
}


#' Enhanced bridge gene network (internal)
#'
#' Uses MDS on the full Jaccard-distance matrix (focal pair + all bridge genes)
#' to derive angular positions that reflect inter-gene pathway similarity.
#' Radial distance is then overridden by each bridge gene's total shared term
#' count with the focal pair (more shared => closer to centre).  This gives a
#' layout where functionally similar bridge genes cluster angularly, and the
#' most important bridging intermediaries sit in the inner ring.
#'
#' @param gene1,gene2 Focal gene pair.
#' @param bridge_genes Character vector of bridge genes.
#' @param shared_terms Character vector of shared GO/KEGG terms.
#' @param mat Expression matrix (genes x cells).
#' @param prior_net Prior network object.
#' @param top_n Maximum bridges to display.
#' @param layout Ignored (kept for API compatibility); layout is MDS-radial.
#' @param label_size Numeric; gene label font size.
#' @param pt_size_range Numeric vector of length 2; node size range for bridge genes.
#' @param edge_width_range Numeric vector of length 2; edge width range for focal edges.
#' @param title Character; plot title.
#' @param sim_threshold Numeric (0 to 1); minimum Jaccard similarity between
#'   bridge genes to draw a dotted similarity edge.  Default 0.05.
#' @return A ggplot object.
#' @keywords internal
.plot_bridge_network_enhanced <- function(gene1, gene2,
                                           bridge_genes, shared_terms,
                                           mat, prior_net, top_n = 15,
                                           layout = "auto",
                                           label_size = 3,
                                           pt_size_range = c(3, 9),
                                           edge_width_range = c(0.4, 2),
                                           title = NULL,
                                           sim_threshold = 0.05) {

  if (is.null(title)) title <- "Prior knowledge network"

  if (length(bridge_genes) == 0 && length(shared_terms) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                           label = paste0("No prior knowledge connections\n",
                                          "found between ", gene1, " & ", gene2),
                           size = 3.5, colour = "grey50") +
        ggplot2::theme_void() +
        ggplot2::ggtitle(title)
    )
  }

  if (length(bridge_genes) > top_n) bridge_genes <- bridge_genes[seq_len(top_n)]

  g2t  <- prior_net$gene_to_terms
  n_bg <- length(bridge_genes)

  # ------------------------------------------------------------------ #
  #  Per-bridge: shared term counts with each focal gene                #
  # ------------------------------------------------------------------ #
  n_with_g1 <- vapply(bridge_genes, function(bg) {
    length(intersect(g2t[[bg]], g2t[[gene1]]))
  }, integer(1))
  n_with_g2 <- vapply(bridge_genes, function(bg) {
    length(intersect(g2t[[bg]], g2t[[gene2]]))
  }, integer(1))
  n_total <- n_with_g1 + n_with_g2

  # ------------------------------------------------------------------ #
  #  Build full Jaccard similarity matrix: focal pair + all bridge genes#
  #  Then MDS on (1 - Jaccard) distance to get 2-D embedding           #
  # ------------------------------------------------------------------ #
  all_g   <- c(gene1, gene2, bridge_genes)
  n_all   <- length(all_g)
  jac_mat <- matrix(0, n_all, n_all, dimnames = list(all_g, all_g))
  for (i_g in seq_len(n_all)) {
    for (j_g in seq_len(n_all)) {
      if (i_g == j_g) { jac_mat[i_g, j_g] <- 1; next }
      ti <- g2t[[all_g[i_g]]]; tj <- g2t[[all_g[j_g]]]
      u  <- length(union(ti, tj))
      jac_mat[i_g, j_g] <- if (u > 0) length(intersect(ti, tj)) / u else 0
    }
  }

  dist_mat <- stats::as.dist(1 - jac_mat)
  mds_xy   <- tryCatch(
    stats::cmdscale(dist_mat, k = 2),
    error = function(e) {
      # Fallback: random 2D positions
      matrix(stats::rnorm(n_all * 2), n_all, 2,
             dimnames = list(all_g, c("V1", "V2")))
    }
  )

  # Angular position of each bridge gene derived from MDS embedding
  bg_theta <- atan2(mds_xy[bridge_genes, 2], mds_xy[bridge_genes, 1])

  # ------------------------------------------------------------------ #
  #  Radial positions: radius overridden by n_total strength            #
  #  (MDS provides angular clustering; strength provides depth)         #
  # ------------------------------------------------------------------ #
  r_inner <- 0.60   # ring for strongest bridge genes
  r_outer <- 1.35   # ring for weakest bridge genes

  s_rng <- range(n_total)
  s_norm <- if (diff(s_rng) > 0) {
    (n_total - s_rng[1]) / diff(s_rng)
  } else {
    rep(0.5, n_bg)
  }
  radii_bg <- r_outer - s_norm * (r_outer - r_inner)

  bx <- radii_bg * cos(bg_theta)
  by <- radii_bg * sin(bg_theta)

  # Focal gene positions: derived from MDS but scaled to be near origin
  # Normalise MDS so focal genes sit within a small central radius (≤ 0.3)
  focal_mds <- mds_xy[c(gene1, gene2), , drop = FALSE]
  focal_scale <- max(sqrt(rowSums(focal_mds^2))) + 1e-9
  focal_xy    <- focal_mds / focal_scale * 0.22   # keep focals near centre

  focal_x <- focal_xy[, 1]
  focal_y <- focal_xy[, 2]

  # ------------------------------------------------------------------ #
  #  Node data frame                                                    #
  # ------------------------------------------------------------------ #
  all_nodes <- c(gene1, gene2, bridge_genes)
  node_x    <- c(focal_x, bx)
  node_y    <- c(focal_y, by)
  node_type <- c("focal", "focal", rep("bridge", n_bg))

  node_expr <- vapply(all_nodes, function(g) {
    if (g %in% rownames(mat)) mean(mat[g, ]) else 0
  }, numeric(1))

  bridge_sizes <- pt_size_range[1] + s_norm * diff(pt_size_range)
  node_size    <- c(9, 9, bridge_sizes)

  node_df <- data.frame(
    x         = node_x,
    y         = node_y,
    name      = all_nodes,
    node_type = node_type,
    mean_expr = node_expr,
    node_size = node_size,
    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------------ #
  #  Focal → bridge edges                                               #
  # ------------------------------------------------------------------ #
  max_w <- max(c(n_with_g1, n_with_g2), na.rm = TRUE)
  norm_w <- function(w) {
    edge_width_range[1] + (w / max(max_w, 1L)) * diff(edge_width_range)
  }

  focal_edges <- data.frame(
    x = numeric(), y = numeric(), xend = numeric(), yend = numeric(),
    colour = character(), lwd = numeric(), stringsAsFactors = FALSE
  )
  for (i in seq_len(n_bg)) {
    w1 <- n_with_g1[i]; w2 <- n_with_g2[i]
    if (w1 > 0)
      focal_edges <- rbind(focal_edges, data.frame(
        x = focal_x[1], y = focal_y[1], xend = bx[i], yend = by[i],
        colour = "g1", lwd = norm_w(w1), stringsAsFactors = FALSE))
    if (w2 > 0)
      focal_edges <- rbind(focal_edges, data.frame(
        x = focal_x[2], y = focal_y[2], xend = bx[i], yend = by[i],
        colour = "g2", lwd = norm_w(w2), stringsAsFactors = FALSE))
  }

  # Direct focal-pair edge
  direct_edge <- if (length(shared_terms) > 0)
    data.frame(x = focal_x[1], y = focal_y[1],
               xend = focal_x[2], yend = focal_y[2],
               stringsAsFactors = FALSE)
  else NULL

  # ------------------------------------------------------------------ #
  #  Bridge-bridge similarity edges (thin dotted)                      #
  # ------------------------------------------------------------------ #
  sim_edges <- data.frame(
    x = numeric(), y = numeric(), xend = numeric(), yend = numeric(),
    sim = numeric(), stringsAsFactors = FALSE
  )
  if (n_bg >= 2) {
    for (i in seq_len(n_bg - 1)) {
      for (j in seq(i + 1L, n_bg)) {
        ti  <- g2t[[bridge_genes[i]]]; tj <- g2t[[bridge_genes[j]]]
        u   <- length(union(ti, tj))
        if (u == 0L) next
        jac <- length(intersect(ti, tj)) / u
        if (jac >= sim_threshold)
          sim_edges <- rbind(sim_edges, data.frame(
            x = bx[i], y = by[i], xend = bx[j], yend = by[j],
            sim = jac, stringsAsFactors = FALSE))
      }
    }
  }

  # ------------------------------------------------------------------ #
  #  Plot limits                                                        #
  # ------------------------------------------------------------------ #
  pad      <- 0.55
  xy_lim   <- r_outer + pad

  expr_bridge <- node_expr[node_type == "bridge"]
  expr_lim    <- range(expr_bridge, na.rm = TRUE)

  n_shared_terms <- length(shared_terms)
  subtitle_txt <- paste0(
    n_shared_terms, if (n_shared_terms == 1) " shared term" else " shared terms",
    " \u00B7 ", n_bg, if (n_bg == 1) " bridge gene" else " bridge genes",
    if (n_bg > 0) paste0(" \u00B7 top: ", bridge_genes[which.max(n_total)]) else ""
  )

  # ------------------------------------------------------------------ #
  #  Build ggplot                                                       #
  # ------------------------------------------------------------------ #
  p <- ggplot2::ggplot() +
    # 1. Bridge-bridge similarity dotted lines (bottom layer)
    {
      if (nrow(sim_edges) > 0)
        ggplot2::geom_segment(
          data = sim_edges,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, alpha = sim),
          colour = "grey50", linewidth = 0.28, linetype = "dotted"
        )
    } +
    {
      if (nrow(sim_edges) > 0)
        ggplot2::scale_alpha_continuous(
          range = c(0.12, 0.65), name = "Bridge\nsimilarity",
          guide = ggplot2::guide_legend(
            override.aes = list(
              linewidth = 0.4, colour = "grey50", linetype = "dotted"
            ),
            keywidth  = ggplot2::unit(0.8, "cm"),
            keyheight = ggplot2::unit(0.4, "cm")
          )
        )
    } +
    # 2. Direct focal-pair edge
    {
      if (!is.null(direct_edge))
        ggplot2::geom_segment(
          data = direct_edge,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          colour = "#8B0000", linewidth = 1.1, linetype = "solid", alpha = 0.85
        )
    } +
    # 3. Focal → bridge solid edges
    {
      if (nrow(focal_edges) > 0)
        ggplot2::geom_segment(
          data = focal_edges,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                        colour = colour, linewidth = lwd),
          linetype = "solid", alpha = 0.72, lineend = "round"
        )
    } +
    ggplot2::scale_colour_manual(
      values = c("g1" = "#C0392B", "g2" = "#2471A3"),
      labels = c("g1" = gene1, "g2" = gene2),
      name   = "Connected\nto"
    ) +
    ggplot2::scale_linewidth_identity(guide = "none") +
    # 4. Bridge gene nodes (fill = mean expression)
    ggplot2::geom_point(
      data = node_df[node_df$node_type == "bridge", ],
      ggplot2::aes(x = x, y = y, fill = mean_expr, size = node_size),
      shape = 21, colour = "white", stroke = 0.55, alpha = 0.93
    ) +
    ggplot2::scale_fill_gradientn(
      colours  = c("#2166AC", "#4393C3", "#92C5DE",
                   "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
      na.value = "grey82",
      name     = "Mean expr\n(bridge)",
      limits   = if (diff(expr_lim) > 0) expr_lim else NULL,
      guide    = ggplot2::guide_colourbar(
        title.position = "top",
        barwidth  = ggplot2::unit(0.35, "cm"),
        barheight = ggplot2::unit(2.2,  "cm")
      )
    ) +
    # 5. Focal gene nodes (red, fixed, on top)
    ggplot2::geom_point(
      data = node_df[node_df$node_type == "focal", ],
      ggplot2::aes(x = x, y = y, size = node_size),
      shape = 21, fill = "#C0392B", colour = "white", stroke = 1.3, alpha = 1
    ) +
    ggplot2::scale_size_identity(guide = "none") +
    # 6. Focal gene labels (bold, red, nudged away from each other)
    ggplot2::geom_text(
      data = node_df[node_df$node_type == "focal", ],
      ggplot2::aes(
        x     = x + c(-0.13, 0.13),
        y     = y + 0.14,
        label = name
      ),
      size     = label_size + 1,
      fontface = "bold",
      colour   = "#C0392B",
      vjust    = 0
    ) +
    # 7. Bridge gene labels (repel)
    ggrepel::geom_text_repel(
      data = node_df[node_df$node_type == "bridge", ],
      ggplot2::aes(x = x, y = y, label = name),
      size           = label_size,
      colour         = "grey15",
      fontface       = "plain",
      box.padding    = ggplot2::unit(0.3,  "lines"),
      point.padding  = ggplot2::unit(0.25, "lines"),
      segment.colour = "grey65",
      segment.size   = 0.28,
      max.overlaps   = 25,
      seed           = 7L
    ) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::xlim(-xy_lim, xy_lim) +
    ggplot2::ylim(-xy_lim, xy_lim) +
    ggplot2::ggtitle(title, subtitle = subtitle_txt) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 10, face = "bold",
                                             margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40",
                                             margin = ggplot2::margin(b = 4)),
      legend.text   = ggplot2::element_text(size = 6),
      legend.title  = ggplot2::element_text(size = 7),
      legend.key.size  = ggplot2::unit(0.4, "cm"),
      legend.position  = "right",
      legend.box       = "vertical",
      legend.spacing.y = ggplot2::unit(0.25, "cm"),
      plot.margin   = ggplot2::margin(8, 8, 8, 8)
    )

  p
}
