#' Compute pairwise co-expression metrics for a set of genes
#'
#' Calculates Pearson correlation, Spearman correlation, biweight
#' midcorrelation, mutual information (discretised), and expression ratio
#' consistency across cell clusters.  Designed for speed: uses vectorised
#' matrix algebra where possible and avoids per-pair loops for correlation.
#'
#' @param mat Numeric (or sparse) matrix, genes in rows, cells in columns.
#'     Should be log-normalised expression values.
#' @param features Character vector of gene names to include (must be rownames
#'     of `mat`).
#' @param cluster_ids Factor or character vector of cluster assignments,
#'     length = ncol(mat).  Used for ratio-consistency calculation.  NULL to
#'     skip that metric.
#' @param cor_method Character vector of correlation types to compute.
#'     Any subset of `c("pearson", "spearman", "biweight")`.
#' @param n_mi_bins Integer; number of bins for mutual information
#'     discretisation. Set to 0 to skip MI.
#' @param min_cells_expressed Integer; minimum number of cells where *both*
#'     genes of a pair must be expressed (> 0) to retain the pair.
#' @param verbose Logical; print progress messages.
#'
#' @return A `data.table` with columns: `gene1`, `gene2`, plus any computed
#'     metric columns (`cor_pearson`, `cor_spearman`, `cor_biweight`,
#'     `mi_score`, `ratio_consistency`).
#'
#' @details
#' **Biweight midcorrelation** (Langfelder & Horvath, 2012) is a robust
#' alternative to Pearson correlation that down-weights outlier observations,
#' which is particularly valuable for noisy single-cell data.
#'
#' **Ratio consistency** measures whether the expression ratio of two genes is
#' stable across clusters: for each cluster the log-ratio median is computed,
#' and consistency = 1 - CoV of those medians (bounded to \[0, 1\]).  High
#' values indicate the two genes maintain a fixed stoichiometric relationship
#' across cell populations -- a hallmark of genuine co-regulation.
#'
#' **Mutual information** captures non-linear dependencies missed by
#' correlation.  Expression values are discretised into equal-frequency bins
#' and MI is estimated via the plug-in estimator.
#'
#' @keywords internal
.compute_coexpression <- function(mat,
                                  features,
                                  cluster_ids = NULL,
                                  cor_method  = c("pearson", "spearman", "biweight"),
                                  n_mi_bins   = 5,
                                  min_cells_expressed = 10,
                                  verbose     = TRUE) {

  # --- Subset matrix to requested features ----------------------------------
  features <- intersect(features, rownames(mat))
  if (length(features) < 2) {
    stop("At least 2 valid features required for co-expression analysis.",
         call. = FALSE)
  }
  mat <- mat[features, , drop = FALSE]

  # Convert sparse to dense for correlation (faster for moderate gene counts)
  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat_dense <- as.matrix(mat)
  } else {
    mat_dense <- mat
  }

  n_genes <- nrow(mat_dense)
  n_cells <- ncol(mat_dense)
  gene_names <- rownames(mat_dense)

  .msg("Computing co-expression for ", n_genes, " genes across ",
       n_cells, " cells ...", verbose = verbose)


  # --- Generate all unique pairs -------------------------------------------
  pair_idx <- utils::combn(n_genes, 2)
  n_pairs <- ncol(pair_idx)

  result <- data.table::data.table(
    gene1 = gene_names[pair_idx[1, ]],
    gene2 = gene_names[pair_idx[2, ]]
  )

  # --- Filter by minimum co-expression ------------------------------------
  if (min_cells_expressed > 0) {
    co_expressed <- vapply(seq_len(n_pairs), function(k) {
      sum(mat_dense[pair_idx[1, k], ] > 0 & mat_dense[pair_idx[2, k], ] > 0)
    }, integer(1))
    keep <- co_expressed >= min_cells_expressed
    result <- result[keep, ]
    pair_idx <- pair_idx[, keep, drop = FALSE]
    n_pairs <- ncol(pair_idx)
    .msg("  Retained ", n_pairs, " pairs with >= ", min_cells_expressed,
         " co-expressing cells.", verbose = verbose)
    if (n_pairs == 0) return(result)
  }

  # --- Pearson correlation (vectorised) ------------------------------------
  if ("pearson" %in% cor_method) {
    .msg("  Pearson correlation ...", verbose = verbose)
    cor_mat <- stats::cor(t(mat_dense), method = "pearson")
    result[, cor_pearson := cor_mat[cbind(pair_idx[1, ], pair_idx[2, ])]]
  }

  # --- Spearman correlation (vectorised) -----------------------------------
  if ("spearman" %in% cor_method) {
    .msg("  Spearman correlation ...", verbose = verbose)
    rank_mat <- t(apply(mat_dense, 1, rank))
    cor_sp <- stats::cor(t(rank_mat), method = "pearson")  # Spearman = Pearson of ranks
    result[, cor_spearman := cor_sp[cbind(pair_idx[1, ], pair_idx[2, ])]]
  }

  # --- Biweight midcorrelation ---------------------------------------------
  if ("biweight" %in% cor_method) {
    .msg("  Biweight midcorrelation ...", verbose = verbose)
    result[, cor_biweight := vapply(seq_len(n_pairs), function(k) {
      .bicor(mat_dense[pair_idx[1, k], ], mat_dense[pair_idx[2, k], ])
    }, numeric(1))]
  }

  # --- Mutual information --------------------------------------------------
  if (n_mi_bins > 0) {
    .msg("  Mutual information ...", verbose = verbose)
    result[, mi_score := vapply(seq_len(n_pairs), function(k) {
      .mutual_info(mat_dense[pair_idx[1, k], ],
                   mat_dense[pair_idx[2, k], ],
                   n_bins = n_mi_bins)
    }, numeric(1))]
  }

  # --- Ratio consistency across clusters -----------------------------------
  if (!is.null(cluster_ids)) {
    .msg("  Ratio consistency ...", verbose = verbose)
    cluster_ids <- as.factor(cluster_ids)
    result[, ratio_consistency := vapply(seq_len(n_pairs), function(k) {
      .ratio_consistency(mat_dense[pair_idx[1, k], ],
                         mat_dense[pair_idx[2, k], ],
                         cluster_ids)
    }, numeric(1))]
  }

  result
}


# ---- Low-level helper functions (not exported) ----------------------------

#' Biweight midcorrelation between two numeric vectors
#'
#' Implements the Tukey biweight robust correlation (Langfelder & Horvath 2012).
#' @param x,y Numeric vectors of equal length.
#' @return Scalar biweight midcorrelation.
#' @keywords internal
.bicor <- function(x, y) {
  # Median and MAD
  mx <- stats::median(x)
  my <- stats::median(y)
  mad_x <- stats::median(abs(x - mx)) * 1.4826  # consistent MAD
  mad_y <- stats::median(abs(y - my)) * 1.4826

  # For sparse scRNA-seq data where median is 0 and MAD is 0,
  # fall back to IQR-based scale or Pearson on non-zero cells
  if (mad_x < .Machine$double.eps || mad_y < .Machine$double.eps) {
    # Fallback: compute Pearson on cells where both genes are expressed
    nonzero <- x > 0 & y > 0
    if (sum(nonzero) < 5) return(0)
    
    # Check for zero variance (constant vectors)
    if (stats::sd(x[nonzero]) < .Machine$double.eps || 
        stats::sd(y[nonzero]) < .Machine$double.eps) {
      return(0)
    }
    
    cor_val <- stats::cor(x[nonzero], y[nonzero], method = "pearson")
    return(if (is.na(cor_val)) 0 else cor_val)
  }

  # Standardised residuals
  u <- (x - mx) / (9 * mad_x)
  v <- (y - my) / (9 * mad_y)

  # Indicator: |u| < 1, |v| < 1
  Iu <- abs(u) < 1
  Iv <- abs(v) < 1

  # Biweight kernel
  a <- (x - mx) * (1 - u^2)^2 * Iu
  b <- (y - my) * (1 - v^2)^2 * Iv

  denom <- sqrt(sum(a^2) * sum(b^2))
  if (denom < .Machine$double.eps) return(0)

  sum(a * b) / denom
}


#' Plug-in mutual information estimator with equal-frequency binning
#'
#' @param x,y Numeric vectors.
#' @param n_bins Number of bins.
#' @return Scalar MI in nats.
#' @keywords internal
.mutual_info <- function(x, y, n_bins = 5) {
  n <- length(x)
  if (n < n_bins * 2) return(0)

  # For sparse data, first separate zero/non-zero, then bin the non-zeros
  .safe_bin <- function(v, nb) {
    # If most values are 0, create a "zero" bin + bins for non-zeros
    pct_zero <- sum(v == 0) / length(v)
    if (pct_zero > 0.5 && sum(v > 0) >= nb) {
      bins <- rep(1L, length(v))
      nz <- v > 0
      nz_vals <- v[nz]
      nz_brks <- unique(stats::quantile(nz_vals,
                                         probs = seq(0, 1, length.out = nb)))
      if (length(nz_brks) >= 2) {
        bins[nz] <- as.integer(cut(nz_vals, breaks = nz_brks,
                                    include.lowest = TRUE)) + 1L
      }
      return(bins)
    }
    brks <- unique(stats::quantile(v, probs = seq(0, 1, length.out = nb + 1)))
    if (length(brks) < 2) return(rep(1L, length(v)))
    as.integer(cut(v, breaks = brks, include.lowest = TRUE))
  }

  bx <- .safe_bin(x, n_bins)
  by <- .safe_bin(y, n_bins)

  # Joint and marginal counts
  joint <- table(bx, by)
  pxy <- joint / n
  px <- rowSums(joint) / n
  py <- colSums(joint) / n

  # MI = sum p(x,y) * log(p(x,y) / (p(x)*p(y)))
  mi <- 0
  for (i in seq_along(px)) {
    for (j in seq_along(py)) {
      if (pxy[i, j] > 0 && px[i] > 0 && py[j] > 0) {
        mi <- mi + pxy[i, j] * log(pxy[i, j] / (px[i] * py[j]))
      }
    }
  }
  max(mi, 0)
}


#' Expression ratio consistency across clusters
#'
#' For each cluster, computes the median log2(expr_g1 + 1) - log2(expr_g2 + 1).
#' Returns 1 - CoV of cluster medians, bounded to \[0, 1\].
#'
#' @param x,y Numeric vectors (expression of gene 1, gene 2).
#' @param clusters Factor of cluster assignments.
#' @return Scalar in \[0, 1\].
#' @keywords internal
.ratio_consistency <- function(x, y, clusters) {
  # Only consider cells where at least one gene is expressed
  expressed <- (x > 0) | (y > 0)
  if (sum(expressed) < 10) return(0)

  log_ratio <- log2(x + 1) - log2(y + 1)
  cluster_medians <- tapply(log_ratio[expressed], clusters[expressed], stats::median)
  cluster_medians <- cluster_medians[!is.na(cluster_medians)]

  if (length(cluster_medians) < 2) return(0)

  mn <- mean(cluster_medians)
  if (abs(mn) < .Machine$double.eps) {
    # If mean ratio is ~0, check if sd is also small (= consistent at 0)
    if (stats::sd(cluster_medians) < 0.1) return(0.8)
    return(0)
  }

  cov <- stats::sd(cluster_medians) / abs(mn)
  # Bounded consistency score
  max(1 - cov, 0)
}
