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
#' **Performance (v0.1.1):** Co-expression filtering, biweight midcorrelation,
#' mutual information, and ratio consistency are all vectorised where possible.
#' The co-expression filter uses a single matrix crossproduct instead of
#' per-pair loops.  Biweight midcorrelation uses a fast vectorised
#' implementation that processes all pairs at once via matrix operations.
#' Mutual information uses pre-computed bin matrices.
#' Ratio consistency is vectorised over clusters using matrix operations.
#'
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

  # --- Filter by minimum co-expression (vectorised via crossproduct) --------
  if (min_cells_expressed > 0) {
    # Binary matrix: gene expressed (> 0) in each cell
    bin_mat <- mat_dense > 0
    storage.mode(bin_mat) <- "double"
    # tcrossprod gives the count of cells where both genes are expressed
    co_count_mat <- tcrossprod(bin_mat)
    co_expressed <- co_count_mat[cbind(pair_idx[1, ], pair_idx[2, ])]
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
    cor_sp <- stats::cor(t(rank_mat), method = "pearson")
    result[, cor_spearman := cor_sp[cbind(pair_idx[1, ], pair_idx[2, ])]]
  }

  # --- Biweight midcorrelation (batch-vectorised) ---------------------------
  if ("biweight" %in% cor_method) {
    .msg("  Biweight midcorrelation ...", verbose = verbose)
    bicor_mat <- .bicor_matrix(mat_dense)
    result[, cor_biweight := bicor_mat[cbind(pair_idx[1, ], pair_idx[2, ])]]
  }

  # --- Mutual information (batch-vectorised) --------------------------------
  if (n_mi_bins > 0) {
    .msg("  Mutual information ...", verbose = verbose)
    mi_vals <- .mutual_info_batch(mat_dense, pair_idx, n_bins = n_mi_bins)
    result[, mi_score := mi_vals]
  }

  # --- Ratio consistency across clusters (batch-vectorised) -----------------
  if (!is.null(cluster_ids)) {
    .msg("  Ratio consistency ...", verbose = verbose)
    rc_vals <- .ratio_consistency_batch(mat_dense, pair_idx, cluster_ids)
    result[, ratio_consistency := rc_vals]
  }

  result
}


# ---- Low-level helper functions (not exported) ----------------------------

#' Compute biweight midcorrelation matrix for all gene pairs
#'
#' Vectorised implementation that processes all genes at once using matrix
#' operations.  For genes where the MAD is zero (common in sparse scRNA-seq),
#' falls back to Pearson on non-zero cells.
#'
#' @param mat Dense numeric matrix, genes in rows, cells in columns.
#' @return Symmetric numeric matrix of biweight midcorrelations.
#' @keywords internal
.bicor_matrix <- function(mat) {
  n_genes <- nrow(mat)
  n_cells <- ncol(mat)

  # Compute medians and MADs for all genes at once
  medians <- apply(mat, 1, stats::median)
  residuals <- mat - medians                    # gene × cell
  abs_resid <- abs(residuals)
  mads <- apply(abs_resid, 1, stats::median) * 1.4826  # consistent MAD

  # Identify genes with zero MAD (need fallback)
  zero_mad <- mads < .Machine$double.eps

  # --- For genes with non-zero MAD: batch biweight computation ---
  # Standardised residuals u_ij = (x_ij - median_i) / (9 * MAD_i)
  # Biweight kernel: a_ij = (x_ij - median_i) * (1 - u_ij^2)^2 * I(|u_ij| < 1)
  if (any(!zero_mad)) {
    safe_mads <- mads
    safe_mads[zero_mad] <- 1  # placeholder to avoid division by zero
    u <- residuals / (9 * safe_mads)             # gene × cell
    within_range <- abs(u) < 1                   # logical gene × cell
    a <- residuals * (1 - u^2)^2 * within_range  # gene × cell
    # Zero out rows with zero MAD (will be handled by fallback)
    a[zero_mad, ] <- 0
    # Denominator: sqrt(sum(a_i^2)) for each gene
    row_norms <- sqrt(rowSums(a^2))
    row_norms[row_norms < .Machine$double.eps] <- 1  # avoid div-by-zero
    # Normalise
    a_normed <- a / row_norms  # gene × cell
    # Biweight correlation matrix = a_normed %*% t(a_normed)
    bicor_mat <- tcrossprod(a_normed)
  } else {
    bicor_mat <- matrix(0, n_genes, n_genes)
  }

  # --- Fallback for zero-MAD genes: Pearson on non-zero cells ---
  zero_idx <- which(zero_mad)
  if (length(zero_idx) > 0) {
    for (gi in zero_idx) {
      for (gj in seq_len(n_genes)) {
        if (gi == gj) {
          bicor_mat[gi, gj] <- 1
          next
        }
        nonzero <- mat[gi, ] > 0 & mat[gj, ] > 0
        if (sum(nonzero) < 5) {
          bicor_mat[gi, gj] <- 0
          bicor_mat[gj, gi] <- 0
          next
        }
        x_nz <- mat[gi, nonzero]
        y_nz <- mat[gj, nonzero]
        if (stats::sd(x_nz) < .Machine$double.eps ||
            stats::sd(y_nz) < .Machine$double.eps) {
          bicor_mat[gi, gj] <- 0
          bicor_mat[gj, gi] <- 0
          next
        }
        val <- stats::cor(x_nz, y_nz, method = "pearson")
        if (is.na(val)) val <- 0
        bicor_mat[gi, gj] <- val
        bicor_mat[gj, gi] <- val
      }
    }
  }

  diag(bicor_mat) <- 1
  rownames(bicor_mat) <- rownames(mat)
  colnames(bicor_mat) <- rownames(mat)
  bicor_mat
}


#' Biweight midcorrelation between two numeric vectors
#'
#' Implements the Tukey biweight robust correlation (Langfelder & Horvath 2012).
#' Retained for single-pair computation in \code{AssessGenePair()}.
#'
#' @param x,y Numeric vectors of equal length.
#' @return Scalar biweight midcorrelation.
#' @keywords internal
.bicor <- function(x, y) {
  # Median and MAD
  mx <- stats::median(x)
  my <- stats::median(y)
  mad_x <- stats::median(abs(x - mx)) * 1.4826
  mad_y <- stats::median(abs(y - my)) * 1.4826

  if (mad_x < .Machine$double.eps || mad_y < .Machine$double.eps) {
    nonzero <- x > 0 & y > 0
    if (sum(nonzero) < 5) return(0)
    if (stats::sd(x[nonzero]) < .Machine$double.eps ||
        stats::sd(y[nonzero]) < .Machine$double.eps) {
      return(0)
    }
    cor_val <- stats::cor(x[nonzero], y[nonzero], method = "pearson")
    return(if (is.na(cor_val)) 0 else cor_val)
  }

  u <- (x - mx) / (9 * mad_x)
  v <- (y - my) / (9 * mad_y)
  Iu <- abs(u) < 1
  Iv <- abs(v) < 1
  a <- (x - mx) * (1 - u^2)^2 * Iu
  b <- (y - my) * (1 - v^2)^2 * Iv
  denom <- sqrt(sum(a^2) * sum(b^2))
  if (denom < .Machine$double.eps) return(0)
  sum(a * b) / denom
}


#' Batch mutual information computation for all gene pairs
#'
#' Vectorised implementation that pre-computes bin assignments for all genes
#' once, then uses fast table lookups for each pair.
#'
#' @param mat Dense numeric matrix, genes in rows, cells in columns.
#' @param pair_idx 2-row integer matrix of pair indices (from combn).
#' @param n_bins Number of bins for discretisation.
#' @return Numeric vector of MI values for each pair.
#' @keywords internal
.mutual_info_batch <- function(mat, pair_idx, n_bins = 5) {
  n_genes <- nrow(mat)
  n_cells <- ncol(mat)
  n_pairs <- ncol(pair_idx)

  if (n_cells < n_bins * 2) return(rep(0, n_pairs))

  # Pre-compute bin assignments for all genes (genes × cells integer matrix)
  bin_mat <- matrix(1L, nrow = n_genes, ncol = n_cells)
  for (g in seq_len(n_genes)) {
    v <- mat[g, ]
    pct_zero <- sum(v == 0) / n_cells
    if (pct_zero > 0.5 && sum(v > 0) >= n_bins) {
      nz <- v > 0
      nz_vals <- v[nz]
      nz_brks <- unique(stats::quantile(nz_vals,
                                         probs = seq(0, 1, length.out = n_bins)))
      if (length(nz_brks) >= 2) {
        bin_mat[g, nz] <- as.integer(cut(nz_vals, breaks = nz_brks,
                                          include.lowest = TRUE)) + 1L
      }
    } else {
      brks <- unique(stats::quantile(v, probs = seq(0, 1, length.out = n_bins + 1)))
      if (length(brks) >= 2) {
        bin_mat[g, ] <- as.integer(cut(v, breaks = brks, include.lowest = TRUE))
      }
    }
  }

  # Compute MI for each pair using pre-binned values
  mi_vals <- numeric(n_pairs)
  for (k in seq_len(n_pairs)) {
    bx <- bin_mat[pair_idx[1, k], ]
    by <- bin_mat[pair_idx[2, k], ]
    joint <- table(bx, by)
    pxy <- joint / n_cells
    px <- rowSums(joint) / n_cells
    py <- colSums(joint) / n_cells

    # Vectorised MI computation
    outer_prod <- outer(px, py)
    valid <- pxy > 0 & outer_prod > 0
    mi <- sum(pxy[valid] * log(pxy[valid] / outer_prod[valid]))
    mi_vals[k] <- max(mi, 0)
  }

  mi_vals
}


#' Plug-in mutual information estimator with equal-frequency binning
#'
#' Single-pair version retained for \code{AssessGenePair()}.
#'
#' @param x,y Numeric vectors.
#' @param n_bins Number of bins.
#' @return Scalar MI in nats.
#' @keywords internal
.mutual_info <- function(x, y, n_bins = 5) {
  n <- length(x)
  if (n < n_bins * 2) return(0)

  .safe_bin <- function(v, nb) {
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

  joint <- table(bx, by)
  pxy <- joint / n
  px <- rowSums(joint) / n
  py <- colSums(joint) / n

  # Vectorised MI (avoid double loop)
  outer_prod <- outer(px, py)
  valid <- pxy > 0 & outer_prod > 0
  mi <- sum(pxy[valid] * log(pxy[valid] / outer_prod[valid]))
  max(mi, 0)
}


#' Batch ratio consistency computation for all gene pairs
#'
#' Vectorised implementation that computes log-ratios and cluster medians
#' using matrix operations instead of per-pair tapply calls.
#'
#' @param mat Dense numeric matrix, genes in rows, cells in columns.
#' @param pair_idx 2-row integer matrix of pair indices.
#' @param cluster_ids Factor of cluster assignments.
#' @return Numeric vector of ratio consistency values.
#' @keywords internal
.ratio_consistency_batch <- function(mat, pair_idx, cluster_ids) {
  cluster_ids <- as.factor(cluster_ids)
  clusters <- levels(cluster_ids)
  n_clusters <- length(clusters)
  n_pairs <- ncol(pair_idx)

  if (n_clusters < 2) return(rep(0, n_pairs))

  # Pre-compute log2(x + 1) for all genes
  log_mat <- log2(mat + 1)

  # Build cluster membership matrix (cells × clusters) for fast group means
  # Using split indices for efficient cluster-level operations
  cl_indices <- split(seq_along(cluster_ids), cluster_ids)

  rc_vals <- numeric(n_pairs)
  for (k in seq_len(n_pairs)) {
    g1_log <- log_mat[pair_idx[1, k], ]
    g2_log <- log_mat[pair_idx[2, k], ]
    log_ratio <- g1_log - g2_log

    # Filter to cells where at least one gene is expressed
    expressed <- mat[pair_idx[1, k], ] > 0 | mat[pair_idx[2, k], ] > 0
    if (sum(expressed) < 10) {
      rc_vals[k] <- 0
      next
    }

    # Compute cluster medians
    cluster_medians <- numeric(n_clusters)
    valid_clusters <- 0L
    for (ci in seq_len(n_clusters)) {
      idx <- cl_indices[[ci]]
      idx_expr <- idx[expressed[idx]]
      if (length(idx_expr) > 0) {
        valid_clusters <- valid_clusters + 1L
        cluster_medians[valid_clusters] <- stats::median(log_ratio[idx_expr])
      }
    }

    if (valid_clusters < 2) {
      rc_vals[k] <- 0
      next
    }

    cluster_medians <- cluster_medians[seq_len(valid_clusters)]
    mn <- mean(cluster_medians)

    if (abs(mn) < .Machine$double.eps) {
      rc_vals[k] <- if (stats::sd(cluster_medians) < 0.1) 0.8 else 0
      next
    }

    cov_val <- stats::sd(cluster_medians) / abs(mn)
    rc_vals[k] <- max(1 - cov_val, 0)
  }

  rc_vals
}


#' Expression ratio consistency across clusters
#'
#' Single-pair version retained for \code{AssessGenePair()}.
#' For each cluster, computes the median log2(expr_g1 + 1) - log2(expr_g2 + 1).
#' Returns 1 - CoV of cluster medians, bounded to \[0, 1\].
#'
#' @param x,y Numeric vectors (expression of gene 1, gene 2).
#' @param clusters Factor of cluster assignments.
#' @return Scalar in \[0, 1\].
#' @keywords internal
.ratio_consistency <- function(x, y, clusters) {
  expressed <- (x > 0) | (y > 0)
  if (sum(expressed) < 10) return(0)

  log_ratio <- log2(x + 1) - log2(y + 1)
  cluster_medians <- tapply(log_ratio[expressed], clusters[expressed], stats::median)
  cluster_medians <- cluster_medians[!is.na(cluster_medians)]

  if (length(cluster_medians) < 2) return(0)

  mn <- mean(cluster_medians)
  if (abs(mn) < .Machine$double.eps) {
    if (stats::sd(cluster_medians) < 0.1) return(0.8)
    return(0)
  }

  cov <- stats::sd(cluster_medians) / abs(mn)
  max(1 - cov, 0)
}
