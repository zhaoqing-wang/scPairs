#' Integrate multi-evidence scores into a composite synergy score
#'
#' Combines co-expression metrics (correlation, mutual information, ratio
#' consistency) and spatial metrics (Lee's L, CLQ) into a single synergy score.
#' Scores are rank-normalised to \[0, 1\] before weighted combination so that
#' different metric scales do not dominate.
#'
#' An empirical p-value is obtained by permuting cell labels.
#'
#' @param pair_dt `data.table` from `.compute_coexpression` and optionally
#'     spatial functions, containing any subset of metric columns.
#' @param weights Named numeric vector of metric weights.  Names must match
#'     column names in `pair_dt`.  Any missing metric is silently skipped.
#'     Default weights:
#'     * `cor_pearson = 1`, `cor_spearman = 1`, `cor_biweight = 1.5`
#'     * `mi_score = 1`, `ratio_consistency = 1.2`
#'     * `spatial_lee_L = 1.5`, `spatial_clq = 1.2`
#' @param n_perm Integer; permutations for composite score p-value (0 to skip).
#' @param mat Expression matrix (needed for permutation null).
#' @param cluster_ids Cluster factor (for permutation null).
#' @param coords Spatial coordinates (for spatial permutation null).
#' @param spatial_k KNN k for spatial metrics in permutation null.
#' @param verbose Logical.
#'
#' @return `pair_dt` with added columns:
#'   * `synergy_score` – composite score in \[0, 1\].
#'   * `p_value` – permutation-based p-value (if `n_perm > 0`).
#'   * `p_adj` – BH-adjusted p-value.
#'   * `rank` – rank by synergy score (1 = strongest).
#'   * `confidence` – categorical label: "High" (p_adj < 0.01),
#'     "Medium" (< 0.05), "Low" (< 0.1), "NS" otherwise.
#'
#' @keywords internal
.integrate_scores <- function(pair_dt,
                              weights     = NULL,
                              n_perm      = 0,
                              mat         = NULL,
                              cluster_ids = NULL,
                              coords      = NULL,
                              spatial_k   = 6,
                              verbose     = TRUE) {

  if (nrow(pair_dt) == 0) return(pair_dt)

  # Default weights
  default_w <- c(
    cor_pearson       = 1.0,
    cor_spearman      = 1.0,
    cor_biweight      = 1.5,
    mi_score          = 1.0,
    ratio_consistency = 1.2,
    spatial_lee_L     = 1.5,
    spatial_clq       = 1.2
  )

  if (is.null(weights)) {
    weights <- default_w
  } else {
    # Fill in any missing weights from defaults
    for (nm in names(default_w)) {
      if (!(nm %in% names(weights))) weights[nm] <- default_w[nm]
    }
  }

  # Identify available metric columns
  metric_cols <- intersect(names(weights), colnames(pair_dt))
  metric_cols <- metric_cols[vapply(metric_cols, function(col) {
    is.numeric(pair_dt[[col]]) && !all(is.na(pair_dt[[col]]))
  }, logical(1))]

  if (length(metric_cols) == 0) {
    warning("No valid metric columns found; synergy_score set to NA.",
            call. = FALSE)
    pair_dt[, synergy_score := NA_real_]
    pair_dt[, rank := NA_integer_]
    pair_dt[, confidence := "NS"]
    return(pair_dt)
  }

  .msg("Integrating ", length(metric_cols), " metrics: ",
       paste(metric_cols, collapse = ", "), verbose = verbose)

  # --- Rank-normalise each metric to [0, 1] --------------------------------
  # For correlations / Lee's L, higher absolute value = stronger; use absolute
  # value for rank normalisation.  For CLQ, higher = stronger co-location.
  # For MI, higher = more dependent.

  abs_metrics <- c("cor_pearson", "cor_spearman", "cor_biweight", "spatial_lee_L")
  raw_metrics <- c("mi_score", "ratio_consistency", "spatial_clq")

  rank_norm <- function(x) {
    r <- rank(x, na.last = "keep", ties.method = "average")
    (r - 1) / (sum(!is.na(r)) - 1)
  }

  normed <- data.table::copy(pair_dt)
  for (col in metric_cols) {
    vals <- pair_dt[[col]]
    if (col %in% abs_metrics) vals <- abs(vals)
    data.table::set(normed, j = col, value = rank_norm(vals))
  }

  # --- Weighted combination ------------------------------------------------
  w_vec <- weights[metric_cols]
  w_vec <- w_vec / sum(w_vec)  # normalise to sum = 1

  score <- rep(0, nrow(normed))
  for (col in metric_cols) {
    vals <- normed[[col]]
    vals[is.na(vals)] <- 0
    score <- score + w_vec[col] * vals
  }

  pair_dt[, synergy_score := score]

  # --- Permutation p-values (optional) -------------------------------------
  if (n_perm > 0 && !is.null(mat)) {
    .msg("Permutation test for composite score (", n_perm, " permutations) ...",
         verbose = verbose)

    features <- unique(c(pair_dt$gene1, pair_dt$gene2))
    features <- intersect(features, rownames(mat))

    null_scores <- matrix(0, nrow = nrow(pair_dt), ncol = n_perm)

    for (p in seq_len(n_perm)) {
      # Permute cell labels (columns)
      perm_idx <- sample(ncol(mat))
      mat_perm <- mat[, perm_idx, drop = FALSE]
      colnames(mat_perm) <- colnames(mat)

      perm_dt <- .compute_coexpression(
        mat         = mat_perm,
        features    = features,
        cluster_ids = if (!is.null(cluster_ids)) cluster_ids[perm_idx] else NULL,
        cor_method  = intersect(c("pearson", "spearman", "biweight"),
                                gsub("^cor_", "", metric_cols)),
        n_mi_bins   = if ("mi_score" %in% metric_cols) 5 else 0,
        min_cells_expressed = 0,
        verbose     = FALSE
      )

      # Merge on gene1+gene2 to align
      perm_dt <- merge(pair_dt[, .(gene1, gene2)], perm_dt,
                       by = c("gene1", "gene2"), all.x = TRUE)

      # Rank-normalise and combine
      perm_score <- rep(0, nrow(perm_dt))
      for (col in intersect(metric_cols, colnames(perm_dt))) {
        vals <- perm_dt[[col]]
        if (col %in% abs_metrics) vals <- abs(vals)
        vals_normed <- rank_norm(vals)
        vals_normed[is.na(vals_normed)] <- 0
        perm_score <- perm_score + w_vec[col] * vals_normed
      }
      null_scores[, p] <- perm_score
    }

    # Empirical p-value
    obs <- pair_dt$synergy_score
    pvals <- vapply(seq_len(nrow(pair_dt)), function(i) {
      (sum(null_scores[i, ] >= obs[i]) + 1) / (n_perm + 1)
    }, numeric(1))

    pair_dt[, p_value := pvals]
    pair_dt[, p_adj := stats::p.adjust(pvals, method = "BH")]
  }

  # --- Rank and confidence -------------------------------------------------
  pair_dt[, rank := frank(-synergy_score, ties.method = "min")]

  if ("p_adj" %in% colnames(pair_dt)) {
    pair_dt[, confidence := ifelse(
      p_adj < 0.01, "High",
      ifelse(p_adj < 0.05, "Medium",
             ifelse(p_adj < 0.1, "Low", "NS"))
    )]
  } else {
    # Without permutation test, use top-percentile heuristic
    pair_dt[, confidence := ifelse(
      synergy_score >= stats::quantile(synergy_score, 0.95, na.rm = TRUE), "High",
      ifelse(synergy_score >= stats::quantile(synergy_score, 0.80, na.rm = TRUE), "Medium",
             ifelse(synergy_score >= stats::quantile(synergy_score, 0.50, na.rm = TRUE), "Low", "NS"))
    )]
  }

  pair_dt
}
