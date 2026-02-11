#' Unified Core Metric Engine for scPairs
#'
#' @description
#' Shared computation backbone invoked by all three discovery functions
#' (`FindAllPairs`, `FindGenePairs`, `AssessGenePair`).
#' Accepts a pre-built pair table and computes the requested metric layers
#' (expression, neighbourhood, prior knowledge, spatial) on demand.
#'
#' Users control which layers to compute via dedicated flags.  When
#' `mode = "prior_only"`, only prior knowledge metrics are returned; when
#' `mode = "expression"`, only expression-based metrics are computed; the
#' default `mode = "all"` computes every available layer.
#'
#' @param mat Expression matrix (genes x cells, dense or sparse).
#' @param pair_dt `data.table` with columns `gene1`, `gene2`.
#' @param cluster_ids Factor of cluster assignments.
#' @param object Seurat object (needed for KNN graph / embedding).
#' @param mode Character; one of `"all"` (default), `"expression"`,
#'     `"prior_only"`.
#' @param cor_method Character vector of correlation methods.
#' @param n_mi_bins Integer; bins for mutual information.
#' @param min_cells_expressed Integer; minimum co-expressing cells.
#' @param use_neighbourhood Logical; compute neighbourhood metrics.
#' @param neighbourhood_k Integer; KNN k.
#' @param neighbourhood_reduction Character; reduction for KNN.
#' @param smooth_alpha Numeric; self-weight for KNN smoothing.
#' @param use_prior Logical; compute prior knowledge metrics.
#' @param organism Character; `"mouse"` or `"human"`.
#' @param custom_pairs Optional data.frame of custom interactions.
#' @param use_spatial Logical; compute spatial metrics.
#' @param spatial_k Integer; spatial neighbourhood k.
#' @param n_perm Integer; permutations for score integration.
#' @param weights Named numeric; metric weights.
#' @param verbose Logical.
#'
#' @return A list with:
#' \describe{
#'   \item{`pair_dt`}{Updated `data.table` with all metric columns.}
#'   \item{`prior_net`}{Prior network object (if computed).}
#'   \item{`W`}{KNN weight matrix (if computed).}
#'   \item{`has_spatial`}{Logical.}
#'   \item{`has_neighbourhood`}{Logical.}
#' }
#'
#' @keywords internal
.compute_pair_metrics <- function(mat,
                                  pair_dt,
                                  cluster_ids,
                                  object,
                                  mode                    = c("all", "expression", "prior_only"),
                                  cor_method              = c("pearson", "spearman", "biweight"),
                                  n_mi_bins               = 5,
                                  min_cells_expressed     = 10,
                                  use_neighbourhood       = TRUE,
                                  neighbourhood_k         = 20,
                                  neighbourhood_reduction = "pca",
                                  smooth_alpha            = 0.3,
                                  use_prior               = TRUE,
                                  organism                = "mouse",
                                  custom_pairs            = NULL,
                                  use_spatial             = TRUE,
                                  spatial_k               = 6,
                                  n_perm                  = 0,
                                  weights                 = NULL,
                                  verbose                 = TRUE) {

  mode <- match.arg(mode)

  has_spatial       <- FALSE
  has_neighbourhood <- FALSE
  prior_net         <- NULL
  W                 <- NULL

  features <- unique(c(pair_dt$gene1, pair_dt$gene2))

  # === Prior-only mode: skip expression metrics ============================
  if (mode == "prior_only") {
    use_neighbourhood <- FALSE
    use_spatial       <- FALSE
  }

  # === Expression metrics ==================================================
  # Skip when caller has already computed expression metrics (empty cor_method)
  has_expr_work <- length(cor_method) > 0 || n_mi_bins > 0
  if (mode != "prior_only" && nrow(pair_dt) > 0 && has_expr_work) {

    # Dense conversion
    if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
      mat_dense <- as.matrix(mat)
    } else {
      mat_dense <- mat
    }

    n_genes <- length(features)
    gene_names <- rownames(mat_dense)

    # Pearson
    if ("pearson" %in% cor_method && !("cor_pearson" %in% colnames(pair_dt))) {
      .msg("  Pearson correlation ...", verbose = verbose)
      cor_mat <- stats::cor(t(mat_dense[intersect(features, gene_names), , drop = FALSE]),
                            method = "pearson")
      pair_dt[, cor_pearson := .extract_pair_vals(cor_mat, pair_dt)]
    }

    # Spearman
    if ("spearman" %in% cor_method && !("cor_spearman" %in% colnames(pair_dt))) {
      .msg("  Spearman correlation ...", verbose = verbose)
      rank_mat <- t(apply(mat_dense[intersect(features, gene_names), , drop = FALSE], 1, rank))
      cor_sp <- stats::cor(t(rank_mat), method = "pearson")
      pair_dt[, cor_spearman := .extract_pair_vals(cor_sp, pair_dt)]
    }

    # Biweight
    if ("biweight" %in% cor_method && !("cor_biweight" %in% colnames(pair_dt))) {
      .msg("  Biweight midcorrelation ...", verbose = verbose)
      bicor_mat <- .bicor_matrix(mat_dense[intersect(features, gene_names), , drop = FALSE])
      pair_dt[, cor_biweight := .extract_pair_vals(bicor_mat, pair_dt)]
    }

    # MI
    if (n_mi_bins > 0 && !("mi_score" %in% colnames(pair_dt))) {
      .msg("  Mutual information ...", verbose = verbose)
      sub_mat <- mat_dense[intersect(features, gene_names), , drop = FALSE]
      sub_genes <- rownames(sub_mat)
      gi <- match(pair_dt$gene1, sub_genes)
      gj <- match(pair_dt$gene2, sub_genes)
      valid <- !is.na(gi) & !is.na(gj)
      pidx <- rbind(gi[valid], gj[valid])
      mi_v <- rep(NA_real_, nrow(pair_dt))
      if (ncol(pidx) > 0) {
        mi_v[valid] <- .mutual_info_batch(sub_mat, pidx, n_bins = n_mi_bins)
      }
      pair_dt[, mi_score := mi_v]
    }

    # Ratio consistency (only if not already present)
    if (!is.null(cluster_ids) && !("ratio_consistency" %in% colnames(pair_dt))) {
      .msg("  Ratio consistency ...", verbose = verbose)
      sub_mat <- mat_dense[intersect(features, gene_names), , drop = FALSE]
      sub_genes <- rownames(sub_mat)
      gi <- match(pair_dt$gene1, sub_genes)
      gj <- match(pair_dt$gene2, sub_genes)
      valid <- !is.na(gi) & !is.na(gj)
      pidx <- rbind(gi[valid], gj[valid])
      rc_v <- rep(NA_real_, nrow(pair_dt))
      if (ncol(pidx) > 0) {
        rc_v[valid] <- .ratio_consistency_batch(sub_mat, pidx, cluster_ids)
      }
      pair_dt[, ratio_consistency := rc_v]
    }
  }

  # === Neighbourhood metrics ===============================================
  if (use_neighbourhood && mode != "prior_only" && nrow(pair_dt) > 0) {
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

      .msg("  KNN-smoothed correlation ...", verbose = verbose)
      pair_dt[, smoothed_cor := .smoothed_cor_batch(
        mat, pair_dt, W, alpha = smooth_alpha, method = "pearson")]

      .msg("  Neighbourhood co-expression score ...", verbose = verbose)
      pair_dt[, neighbourhood_score := .neighbourhood_coexpr_batch(
        mat, pair_dt, W)]

      .msg("  Cluster-level correlation ...", verbose = verbose)
      pair_dt[, cluster_cor := .cluster_cor_batch(mat, cluster_ids, pair_dt)]

      .msg("  Cross-cell-type interaction score ...", verbose = verbose)
      embed <- tryCatch(
        Seurat::Embeddings(object, reduction = neighbourhood_reduction),
        error = function(e) NULL
      )
      if (!is.null(embed)) {
        dims_use <- seq_len(min(30, ncol(embed)))
        pair_dt[, cross_celltype_score := .cross_celltype_batch(
          mat, pair_dt, cluster_ids, embed[, dims_use, drop = FALSE])]
      }

      .msg("  Neighbourhood synergy score ...", verbose = verbose)
      pair_dt[, neighbourhood_synergy := .neighbourhood_synergy_batch(
        mat, pair_dt, W)]
    }
  }

  # === Prior knowledge metrics =============================================
  if (use_prior && nrow(pair_dt) > 0) {
    prior_net <- tryCatch(
      .build_prior_network(organism = organism, genes = features,
                           custom_pairs = custom_pairs, verbose = verbose),
      error = function(e) {
        .msg("Prior knowledge not available: ", conditionMessage(e),
             ". Continuing without prior scores.", verbose = verbose)
        NULL
      }
    )

    if (!is.null(prior_net) && prior_net$n_terms > 0) {
      .msg("Computing prior knowledge scores ...", verbose = verbose)
      pair_dt[, prior_score := .prior_score_batch(pair_dt, prior_net)]

      expressed_genes <- features[features %in% rownames(mat)]
      bridge_res <- .bridge_score_batch(pair_dt, prior_net, expressed_genes)
      pair_dt[, bridge_score := bridge_res$scores]
    }
  }

  # === Spatial metrics =====================================================
  if (use_spatial && mode != "prior_only" && .has_spatial(object) &&
      nrow(pair_dt) > 0) {
    has_spatial <- TRUE
    .msg("Spatial data detected. Computing spatial metrics ...", verbose = verbose)
    coords <- .get_spatial_coords(object)

    pair_dt <- .compute_spatial_lee(
      coords  = coords, mat = mat, pair_dt = pair_dt,
      k = spatial_k, n_perm = min(n_perm, 199), verbose = verbose
    )
    pair_dt <- .compute_spatial_clq(
      coords = coords, mat = mat, pair_dt = pair_dt,
      k = spatial_k, verbose = verbose
    )
  }

  # === Score integration ===================================================
  if (nrow(pair_dt) > 0) {
    pair_dt <- .integrate_scores(
      pair_dt     = pair_dt,
      weights     = weights,
      n_perm      = n_perm,
      mat         = mat,
      cluster_ids = cluster_ids,
      coords      = if (has_spatial) .get_spatial_coords(object) else NULL,
      spatial_k   = spatial_k,
      verbose     = verbose
    )
    data.table::setorder(pair_dt, rank)
  }

  list(
    pair_dt           = pair_dt,
    prior_net         = prior_net,
    W                 = W,
    has_spatial        = has_spatial,
    has_neighbourhood  = has_neighbourhood
  )
}


#' Extract pairwise values from a symmetric matrix given a pair table
#' @keywords internal
.extract_pair_vals <- function(sym_mat, pair_dt) {
  gene_names <- rownames(sym_mat)
  gi <- match(pair_dt$gene1, gene_names)
  gj <- match(pair_dt$gene2, gene_names)
  valid <- !is.na(gi) & !is.na(gj)
  vals <- rep(NA_real_, nrow(pair_dt))
  vals[valid] <- sym_mat[cbind(gi[valid], gj[valid])]
  vals
}
