#' Compute bivariate spatial autocorrelation (Lee's L) for gene pairs
#'
#' Lee's L statistic (Lee, 2001) generalises Moran's I to the bivariate case,
#' measuring spatial co-variation of two variables simultaneously.  A positive L
#' means that nearby locations tend to have similar *joint* expression patterns;
#' a negative L indicates spatial segregation.
#'
#' For efficiency, a k-nearest-neighbour spatial weights matrix is used rather
#' than a full distance matrix.
#'
#' @param coords Data.frame with columns `x`, `y` and rownames = cell barcodes.
#' @param mat Expression matrix (genes x cells).
#' @param pair_dt `data.table` with columns `gene1`, `gene2`.
#' @param k Integer; number of spatial nearest neighbours (default 6 for
#'     Visium hexagonal grids; 15 for general ST).
#' @param n_perm Integer; number of permutations for empirical p-values.
#'     Set to 0 to skip.
#' @param verbose Logical.
#'
#' @return The input `pair_dt` with added columns `spatial_lee_L` and
#'     optionally `spatial_lee_p`.
#'
#' @details
#' **Lee's L statistic** is defined as:
#' \deqn{L(x,y) = \frac{n}{S_0} \cdot \frac{(Wx)^\top (Wy)}{(\sum x_i^2)^{1/2} (\sum y_i^2)^{1/2}}}
#' where W is a row-standardised spatial weight matrix and S_0 = n (under
#' row-standardisation).  The statistic ranges from -1 to 1.
#'
#' @references
#' Lee, S.-I. (2001). Developing a bivariate spatial association measure: An
#' integration of Pearson's r and Moran's I. *Journal of Geographical Systems*,
#' 3(4), 369-385.
#'
#' @keywords internal
.compute_spatial_lee <- function(coords,
                                 mat,
                                 pair_dt,
                                 k       = 6,
                                 n_perm  = 199,
                                 verbose = TRUE) {

  if (nrow(pair_dt) == 0) return(pair_dt)

  # Align cells
  common_cells <- intersect(rownames(coords), colnames(mat))
  if (length(common_cells) < 20) {
    warning("Fewer than 20 cells with both spatial coords and expression. ",
            "Skipping spatial analysis.", call. = FALSE)
    pair_dt[, spatial_lee_L := NA_real_]
    if (n_perm > 0) pair_dt[, spatial_lee_p := NA_real_]
    return(pair_dt)
  }

  coords <- coords[common_cells, ]
  mat <- mat[, common_cells, drop = FALSE]
  n <- length(common_cells)

  # Convert sparse to dense if needed
  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }

  .msg("Computing spatial weight matrix (k=", k, ", n=", n, " spots) ...",
       verbose = verbose)

  # KNN spatial weight matrix
  k_use <- min(k, n - 1)
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(as.matrix(coords[, c("x", "y")]), k = k_use + 1)
    nn_idx <- nn$nn.idx[, -1, drop = FALSE]  # exclude self
  } else {
    # Fallback: brute-force KNN
    dist_mat <- as.matrix(stats::dist(coords[, c("x", "y")]))
    nn_idx <- t(apply(dist_mat, 1, function(row) {
      order(row)[2:(k_use + 1)]
    }))
  }

  # Build row-standardised weight matrix (sparse-like via list)
  # W[i,] = 1/k for each neighbour j, 0 otherwise
  # Wx = for each cell i, mean expression of its neighbours
  .spatial_lag <- function(vec) {
    lag_vec <- numeric(n)
    for (i in seq_len(n)) {
      lag_vec[i] <- mean(vec[nn_idx[i, ]])
    }
    lag_vec
  }

  .lee_L <- function(x_vec, y_vec) {
    # Centre
    x_c <- x_vec - mean(x_vec)
    y_c <- y_vec - mean(y_vec)
    Wx <- .spatial_lag(x_c)
    Wy <- .spatial_lag(y_c)
    denom <- sqrt(sum(x_c^2)) * sqrt(sum(y_c^2))
    if (denom < .Machine$double.eps) return(0)
    n * sum(Wx * Wy) / (n * denom)
  }

  n_pairs <- nrow(pair_dt)
  .msg("Computing Lee's L for ", n_pairs, " gene pairs ...", verbose = verbose)

  # Observed Lee's L
  lee_vals <- vapply(seq_len(n_pairs), function(k_i) {
    g1 <- pair_dt$gene1[k_i]
    g2 <- pair_dt$gene2[k_i]
    if (!(g1 %in% rownames(mat)) || !(g2 %in% rownames(mat))) return(NA_real_)
    .lee_L(mat[g1, ], mat[g2, ])
  }, numeric(1))

  pair_dt[, spatial_lee_L := lee_vals]

  # Permutation p-values
  if (n_perm > 0) {
    .msg("Permutation test (", n_perm, " permutations) ...", verbose = verbose)
    pvals <- vapply(seq_len(n_pairs), function(k_i) {
      if (is.na(lee_vals[k_i])) return(NA_real_)
      g1 <- pair_dt$gene1[k_i]
      g2 <- pair_dt$gene2[k_i]
      x_vec <- mat[g1, ]
      y_vec <- mat[g2, ]
      obs <- lee_vals[k_i]

      null_dist <- vapply(seq_len(n_perm), function(p) {
        .lee_L(x_vec, sample(y_vec))
      }, numeric(1))

      # Two-sided pseudo p-value
      (sum(abs(null_dist) >= abs(obs)) + 1) / (n_perm + 1)
    }, numeric(1))

    pair_dt[, spatial_lee_p := pvals]
  }

  pair_dt
}
