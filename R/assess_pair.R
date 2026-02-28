#' Assess the Synergy of a Specific Gene Pair
#'
#' @description
#' Given two genes, `AssessGenePair` performs an in-depth evaluation of their
#' co-regulatory relationship.  In addition to the standard multi-metric
#' scoring, it computes:
#'
#' * **Per-cluster co-expression** -- correlation within each cell cluster.
#' * **Expression distribution overlap** -- Jaccard index of expressing cells.
#' * **Permutation-based significance** -- 999 permutations by default.
#'
#' @details
#' The `mode` parameter controls which layers are scored:
#' * `"all"` (default) -- full multi-evidence assessment.
#' * `"expression"` -- expression and neighbourhood metrics only.
#' * `"prior_only"` -- prior knowledge scores only.
#'
#' @param object A Seurat object.
#' @param gene1 Character; first gene.
#' @param gene2 Character; second gene.
#' @param assay Character; assay name.
#' @param slot Character; data slot.
#' @param cluster_col Character; cluster column.
#' @param mode Character; `"all"`, `"expression"`, or `"prior_only"`.
#' @param use_neighbourhood Logical; neighbourhood metrics.
#' @param neighbourhood_k Integer; KNN k.
#' @param neighbourhood_reduction Character; reduction for KNN.
#' @param smooth_alpha Numeric; self-weight for smoothing.
#' @param use_spatial Logical.
#' @param spatial_k Integer; spatial KNN k.
#' @param use_prior Logical; prior knowledge scores.
#' @param organism Character; `"mouse"` or `"human"`.
#' @param custom_pairs Optional data.frame of custom interactions.
#' @param n_perm Integer; permutations (default 999).
#' @param verbose Logical.
#'
#' @return A list with class `"scPairs_pair_result"`:
#' \describe{
#'   \item{`gene1`, `gene2`}{The query genes.}
#'   \item{`pairs`}{Single-row `data.table` with all metric columns and
#'     `synergy_score`, `rank`, `confidence` (same format as `FindAllPairs`
#'     output for unified downstream processing).}
#'   \item{`metrics`}{Named list of all computed metrics.}
#'   \item{`per_cluster`}{data.frame of per-cluster correlations.}
#'   \item{`synergy_score`}{Composite score.}
#'   \item{`p_value`}{Permutation p-value.}
#'   \item{`confidence`}{Categorical confidence label.}
#'   \item{`jaccard_index`}{Expression overlap Jaccard index.}
#'   \item{`has_spatial`}{Logical.}
#'   \item{`n_cells`}{Integer.}
#'   \item{`mode`}{Character.}
#' }
#'
#' @family Section_1_Discovery
#'
#' @seealso \code{\link{FindAllPairs}} for genome-wide screening,
#'   \code{\link{FindGenePairs}} for query-centric partner search,
#'   \code{\link{PlotPairSynergy}} and \code{\link{PlotPairSummary}} for
#'   visualising pair-level evidence.
#'
#' @export
#'
#' @examples
#' # Assess the injected co-expressed pair GENE3 & GENE4.
#' result <- AssessGenePair(scpairs_testdata,
#'                          gene1   = "GENE3",
#'                          gene2   = "GENE4",
#'                          mode    = "expression",
#'                          verbose = FALSE)
#' print(result)
AssessGenePair <- function(object,
                           gene1,
                           gene2,
                           assay                   = NULL,
                           slot                    = "data",
                           cluster_col             = NULL,
                           mode                    = c("all", "expression", "prior_only"),
                           use_prior               = TRUE,
                           organism                = "mouse",
                           custom_pairs            = NULL,
                           use_neighbourhood       = TRUE,
                           neighbourhood_k         = 20,
                           neighbourhood_reduction = "pca",
                           smooth_alpha            = 0.3,
                           use_spatial             = TRUE,
                           spatial_k               = 6,
                           n_perm                  = 999,
                           verbose                 = TRUE) {

  mode <- match.arg(mode)

  .validate_seurat(object)
  assay <- assay %||% Seurat::DefaultAssay(object)
  n_cells_total <- ncol(object)

  .validate_percentage(smooth_alpha, "smooth_alpha")
  if (use_neighbourhood && mode != "prior_only") {
    .validate_neighbourhood_params(neighbourhood_k, n_cells_total)
  }
  if (n_perm > 0) .validate_n_perm(n_perm)

  if (mode == "expression")  use_prior <- FALSE
  if (mode == "prior_only") {
    use_neighbourhood <- FALSE
    use_spatial       <- FALSE
  }

  .validate_features(c(gene1, gene2), object, assay)
  if (gene1 == gene2) {
    stop("gene1 and gene2 must be different genes.", call. = FALSE)
  }

  # --- Extract data ---
  mat <- .get_expression_matrix(object, features = c(gene1, gene2),
                                 assay = assay, slot = slot)
  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat_dense <- as.matrix(mat)
  } else {
    mat_dense <- mat
  }

  x <- mat_dense[gene1, ]
  y <- mat_dense[gene2, ]
  n_cells <- length(x)
  cluster_ids <- .resolve_cluster_ids(object, cluster_col)

  .msg("Assessing pair: ", gene1, " -- ", gene2, " (", n_cells, " cells)",
       verbose = verbose)

  # --- Build single-row pair table ---
  pair_dt <- data.table::data.table(gene1 = gene1, gene2 = gene2)

  # === Expression metrics (single-pair, fast) ==============================
  metrics <- list()
  if (mode != "prior_only") {
    metrics$cor_pearson  <- stats::cor(x, y, method = "pearson")
    metrics$cor_spearman <- stats::cor(x, y, method = "spearman")
    metrics$cor_biweight <- .bicor(x, y)
    metrics$mi_score     <- .mutual_info(x, y, n_bins = 5)
    metrics$ratio_consistency <- .ratio_consistency(x, y, cluster_ids)

    # Write into pair_dt for unified integration
    pair_dt[, cor_pearson := metrics$cor_pearson]
    pair_dt[, cor_spearman := metrics$cor_spearman]
    pair_dt[, cor_biweight := metrics$cor_biweight]
    pair_dt[, mi_score := metrics$mi_score]
    pair_dt[, ratio_consistency := metrics$ratio_consistency]

    # Correlation test p-values
    ct_pearson  <- stats::cor.test(x, y, method = "pearson")
    ct_spearman <- stats::cor.test(x, y, method = "spearman", exact = FALSE)
    metrics$p_pearson  <- ct_pearson$p.value
    metrics$p_spearman <- ct_spearman$p.value
  }

  # Expression overlap (Jaccard)
  expr_a <- x > 0
  expr_b <- y > 0
  intersection <- sum(expr_a & expr_b)
  union_ab     <- sum(expr_a | expr_b)
  jaccard <- if (union_ab > 0) intersection / union_ab else 0
  metrics$jaccard_index     <- jaccard
  metrics$n_coexpressing    <- intersection
  metrics$pct_coexpressing  <- intersection / n_cells

  # --- Per-cluster correlations ---
  clusters <- levels(cluster_ids)
  per_cluster <- data.frame(
    cluster      = clusters,
    n_cells      = integer(length(clusters)),
    cor_pearson  = numeric(length(clusters)),
    cor_spearman = numeric(length(clusters)),
    pct_coexpr   = numeric(length(clusters)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(clusters)) {
    idx <- which(cluster_ids == clusters[i])
    n_cl <- length(idx)
    per_cluster$n_cells[i] <- n_cl
    if (n_cl < 5 || mode == "prior_only") {
      per_cluster$cor_pearson[i]  <- NA
      per_cluster$cor_spearman[i] <- NA
      per_cluster$pct_coexpr[i]   <- NA
    } else {
      per_cluster$cor_pearson[i]  <- stats::cor(x[idx], y[idx], method = "pearson")
      per_cluster$cor_spearman[i] <- stats::cor(x[idx], y[idx], method = "spearman")
      per_cluster$pct_coexpr[i]   <- sum(x[idx] > 0 & y[idx] > 0) / n_cl
    }
  }

  # === Neighbourhood metrics (single-pair) =================================
  has_neighbourhood <- FALSE
  embed <- NULL
  W <- NULL
  if (use_neighbourhood && mode != "prior_only") {
    W <- tryCatch(
      .build_knn_graph(object, reduction = neighbourhood_reduction,
                       k = neighbourhood_k),
      error = function(e) {
        .msg("Could not build KNN graph: ", conditionMessage(e),
             ". Skipping neighbourhood metrics.", verbose = verbose)
        NULL
      }
    )

    if (!is.null(W)) {
      has_neighbourhood <- TRUE
      .msg("Computing neighbourhood-aware metrics ...", verbose = verbose)

      smoothed <- .smooth_expression(mat, W, alpha = smooth_alpha)
      if (inherits(smoothed, "sparseMatrix")) smoothed <- as.matrix(smoothed)
      sx <- smoothed[gene1, ]
      sy <- smoothed[gene2, ]
      metrics$smoothed_cor <- stats::cor(sx, sy, method = "pearson")
      pair_dt[, smoothed_cor := metrics$smoothed_cor]

      metrics$neighbourhood_score <- .neighbourhood_coexpr(x, y, W)
      pair_dt[, neighbourhood_score := metrics$neighbourhood_score]

      metrics$cluster_cor <- .cluster_cor(x, y, cluster_ids)
      pair_dt[, cluster_cor := metrics$cluster_cor]

      embed <- tryCatch(
        Seurat::Embeddings(object, reduction = neighbourhood_reduction),
        error = function(e) NULL
      )
      if (!is.null(embed)) {
        dims_use <- seq_len(min(30, ncol(embed)))
        cross_ct <- .cross_celltype(x, y, cluster_ids,
                                    embed[, dims_use, drop = FALSE])
        metrics$cross_celltype_score <- cross_ct$score
        metrics$cross_celltype_r_ab  <- cross_ct$r_ab
        metrics$cross_celltype_r_ba  <- cross_ct$r_ba
        pair_dt[, cross_celltype_score := metrics$cross_celltype_score]
      }

      metrics$neighbourhood_synergy <- .neighbourhood_synergy(x, y, W)
      pair_dt[, neighbourhood_synergy := metrics$neighbourhood_synergy]
    }
  }

  # === Prior knowledge metrics =============================================
  prior_net <- NULL
  bridge_result <- NULL
  if (use_prior) {
    all_genes <- rownames(.get_expression_matrix(object, assay = assay,
                                                  slot = slot))
    prior_net <- tryCatch(
      .build_prior_network(organism = organism, genes = all_genes,
                           custom_pairs = custom_pairs, verbose = verbose),
      error = function(e) {
        .msg("Prior knowledge not available: ", conditionMessage(e),
             verbose = verbose)
        NULL
      }
    )

    if (!is.null(prior_net) && prior_net$n_terms > 0) {
      .msg("Computing prior knowledge scores ...", verbose = verbose)
      metrics$prior_score <- .prior_score(gene1, gene2, prior_net)
      pair_dt[, prior_score := metrics$prior_score]

      expressed_genes <- all_genes[
        Matrix::rowMeans(.get_expression_matrix(object, assay = assay,
                                                slot = slot) > 0) > 0.01]
      bridge_result <- .bridge_score(gene1, gene2, prior_net, expressed_genes)
      metrics$bridge_score    <- bridge_result$score
      metrics$n_bridge_genes  <- length(bridge_result$bridges)
      metrics$bridge_genes    <- paste(utils::head(bridge_result$bridges, 20),
                                       collapse = ", ")
      metrics$shared_terms    <- paste(utils::head(bridge_result$shared_terms, 10),
                                       collapse = ", ")
      pair_dt[, bridge_score := metrics$bridge_score]
    }
  }

  # === Spatial metrics =====================================================
  has_spatial <- FALSE
  if (use_spatial && mode != "prior_only" && .has_spatial(object)) {
    has_spatial <- TRUE
    coords <- .get_spatial_coords(object)
    pair_dt <- .compute_spatial_lee(coords, mat, pair_dt, k = spatial_k,
                                    n_perm = min(n_perm, 199), verbose = verbose)
    pair_dt <- .compute_spatial_clq(coords, mat, pair_dt, k = spatial_k,
                                     verbose = verbose)
    metrics$spatial_lee_L <- pair_dt$spatial_lee_L[1]
    if ("spatial_lee_p" %in% colnames(pair_dt)) {
      metrics$spatial_lee_p <- pair_dt$spatial_lee_p[1]
    }
    metrics$spatial_clq <- pair_dt$spatial_clq[1]
  }

  # === Composite score (using mean of available normalised metrics) ========
  metric_vals <- c(
    if (!is.null(metrics$cor_pearson))  abs(metrics$cor_pearson),
    if (!is.null(metrics$cor_spearman)) abs(metrics$cor_spearman),
    if (!is.null(metrics$cor_biweight)) abs(metrics$cor_biweight),
    if (!is.null(metrics$mi_score))     min(metrics$mi_score / log(5), 1),
    if (!is.null(metrics$ratio_consistency)) metrics$ratio_consistency,
    if (has_neighbourhood && !is.null(metrics$smoothed_cor)) abs(metrics$smoothed_cor),
    if (has_neighbourhood && !is.null(metrics$neighbourhood_score))
      min(pmax(metrics$neighbourhood_score, 0) / 2, 1),
    if (has_neighbourhood && !is.na(metrics$cluster_cor %||% NA))
      abs(metrics$cluster_cor),
    if (has_neighbourhood && !is.na(metrics$cross_celltype_score %||% NA))
      metrics$cross_celltype_score,
    if (has_neighbourhood && !is.null(metrics$neighbourhood_synergy) &&
        !is.na(metrics$neighbourhood_synergy))
      min(pmax(metrics$neighbourhood_synergy, 0) / 2, 1),
    if (!is.null(metrics$prior_score) && !is.na(metrics$prior_score))
      metrics$prior_score,
    if (!is.null(metrics$bridge_score) && !is.na(metrics$bridge_score))
      metrics$bridge_score,
    if (has_spatial && !is.null(metrics$spatial_lee_L))
      abs(metrics$spatial_lee_L),
    if (has_spatial && !is.null(metrics$spatial_clq))
      min(metrics$spatial_clq / 2, 1)
  )
  synergy <- mean(metric_vals, na.rm = TRUE)

  # === Permutation p-value (vectorised) ====================================
  p_value <- NA_real_
  if (n_perm > 0 && mode != "prior_only") {
    .msg("Permutation test (", n_perm, " permutations) ...", verbose = verbose)

    perm_y_mat <- vapply(seq_len(n_perm), function(p) sample(y), numeric(n_cells))
    perm_pearson <- abs(as.numeric(stats::cor(x, perm_y_mat)))
    x_rank <- rank(x)
    perm_y_rank <- apply(perm_y_mat, 2, rank)
    perm_spearman <- abs(as.numeric(stats::cor(x_rank, perm_y_rank)))

    perm_bicor <- numeric(n_perm)
    perm_mi <- numeric(n_perm)
    perm_rc <- numeric(n_perm)
    perm_smoothed <- numeric(n_perm)
    perm_ncs <- numeric(n_perm)
    perm_clcor <- numeric(n_perm)
    perm_cross_ct <- numeric(n_perm)

    for (p in seq_len(n_perm)) {
      yp <- perm_y_mat[, p]
      perm_bicor[p] <- abs(.bicor(x, yp))
      perm_mi[p] <- min(.mutual_info(x, yp, 5) / log(5), 1)
      perm_rc[p] <- .ratio_consistency(x, yp, cluster_ids)
      if (has_neighbourhood) {
        perm_smoothed[p] <- abs(.smoothed_cor(
          rbind(x), rbind(yp), W, alpha = smooth_alpha))
        perm_ncs[p] <- min(pmax(.neighbourhood_coexpr(x, yp, W), 0) / 2, 1)
        perm_clcor[p] <- abs(.cluster_cor(x, yp, cluster_ids))
        if (!is.null(embed)) {
          ct_res <- .cross_celltype(x, yp, cluster_ids,
                                    embed[, dims_use, drop = FALSE])
          perm_cross_ct[p] <- if (!is.na(ct_res$score)) ct_res$score else 0
        }
      }
    }

    n_metrics <- 5L
    null_scores <- perm_pearson + perm_spearman + perm_bicor +
                      perm_mi + perm_rc
    if (has_neighbourhood) {
      null_scores <- null_scores + perm_smoothed + perm_ncs + perm_clcor +
                        perm_cross_ct
      n_metrics <- n_metrics + 4L
    }
    null_scores <- null_scores / n_metrics
    p_value <- (sum(null_scores >= synergy) + 1) / (n_perm + 1)
  }

  confidence <- if (!is.na(p_value)) {
    if (p_value < 0.01) "High"
    else if (p_value < 0.05) "Medium"
    else if (p_value < 0.1) "Low"
    else "NS"
  } else {
    if (synergy >= 0.6) "High"
    else if (synergy >= 0.4) "Medium"
    else if (synergy >= 0.2) "Low"
    else "NS"
  }

  # Write unified columns into pair_dt
  pair_dt[, synergy_score := synergy]
  pair_dt[, rank := 1L]
  pair_dt[, confidence := confidence]
  if (!is.na(p_value)) {
    pair_dt[, p_value := p_value]
    pair_dt[, p_adj := p_value]
  }

  # --- Cross-cell-type detail ---
  cross_celltype_detail <- if (has_neighbourhood && exists("cross_ct")) {
    cross_ct$per_celltype_pair
  } else {
    data.frame()
  }

  # --- Return ---
  structure(
    list(
      gene1                = gene1,
      gene2                = gene2,
      pairs                = pair_dt,
      metrics              = metrics,
      per_cluster          = per_cluster,
      cross_celltype_detail = cross_celltype_detail,
      bridge_genes         = if (!is.null(bridge_result)) bridge_result$bridges else character(0),
      shared_terms         = if (!is.null(bridge_result)) bridge_result$shared_terms else character(0),
      prior_net            = prior_net,
      synergy_score        = synergy,
      p_value              = p_value,
      confidence           = confidence,
      jaccard_index        = jaccard,
      has_spatial           = has_spatial,
      n_cells              = n_cells,
      mode                 = mode
    ),
    class = "scPairs_pair_result"
  )
}


#' Print method for scPairs_pair_result
#' @param x An scPairs_pair_result object.
#' @param ... Ignored.
#' @return The input object \code{x}, returned invisibly.
#' @export
#' @method print scPairs_pair_result
print.scPairs_pair_result <- function(x, ...) {
  cat("scPairs gene pair assessment\n")
  cat(sprintf("  Pair           : %s -- %s\n", x$gene1, x$gene2))
  cat(sprintf("  Cells          : %d\n", x$n_cells))
  cat(sprintf("  Synergy score  : %.4f\n", x$synergy_score))
  cat(sprintf("  Confidence     : %s\n", x$confidence))
  if (!is.na(x$p_value)) {
    cat(sprintf("  p-value        : %.4g\n", x$p_value))
  }
  cat(sprintf("  Jaccard overlap: %.3f\n", x$jaccard_index))
  if (!is.null(x$mode)) cat(sprintf("  Mode           : %s\n", x$mode))
  cat("\n  Metrics:\n")
  for (nm in names(x$metrics)) {
    val <- x$metrics[[nm]]
    if (is.numeric(val) && length(val) == 1) {
      cat(sprintf("    %-20s: %.4f\n", nm, val))
    }
  }
  if (x$has_spatial) cat("  (Spatial metrics included)\n")
  if (length(x$bridge_genes) > 0) {
    cat(sprintf("\n  Bridge genes (%d): %s\n",
                length(x$bridge_genes),
                paste(utils::head(x$bridge_genes, 10), collapse = ", ")))
  }
  if (length(x$shared_terms) > 0) {
    cat(sprintf("  Shared GO/KEGG terms: %d\n", length(x$shared_terms)))
  }
  cat("\n  Per-cluster correlations:\n")
  print(x$per_cluster, row.names = FALSE)
  if (!is.null(x$cross_celltype_detail) && nrow(x$cross_celltype_detail) > 0) {
    cat("\n  Cross-cell-type interactions:\n")
    print(x$cross_celltype_detail, row.names = FALSE)
  }
  invisible(x)
}
