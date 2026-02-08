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
#' **Performance (v0.1.1):** The spatial lag is computed via sparse matrix
#' multiplication (\code{W \%*\% t(mat)}) for all genes at once, replacing the
#' previous per-pair R-level loop.  Lee's L for all pairs is then computed via
#' vectorised column inner products.  Permutation testing is similarly
#' vectorised: the entire gene expression matrix is permuted and all pairs
#' evaluated per permutation.
#'
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
    nn_idx <- nn$nn.idx[, -1, drop = FALSE]
  } else {
    dist_mat <- as.matrix(stats::dist(coords[, c("x", "y")]))
    nn_idx <- t(apply(dist_mat, 1, function(row) {
      order(row)[2:(k_use + 1)]
    }))
  }

  # --- Build sparse row-standardised weight matrix ---------------------------
  # W[i,j] = 1/k for each neighbour j of cell i
  row_i <- rep(seq_len(n), each = k_use)
  col_j <- as.integer(t(nn_idx))
  vals  <- rep(1 / k_use, length(row_i))
  W <- Matrix::sparseMatrix(i = row_i, j = col_j, x = vals,
                            dims = c(n, n))

  # --- Vectorised Lee's L computation for all pairs at once -----------------
  # Centre each gene across cells
  gene_means <- rowMeans(mat)
  mat_c <- mat - gene_means  # centred (genes × cells)

  # Spatial lag for all genes: Wx = W %*% t(mat_c) → cells × genes
  # Then transpose to get genes × cells
  Wmat <- t(as.matrix(W %*% t(mat_c)))  # genes × cells

  # Denominator: sqrt(sum(x_c^2)) for each gene
  gene_norms <- sqrt(rowSums(mat_c^2))
  gene_norms[gene_norms < .Machine$double.eps] <- 1

  n_pairs <- nrow(pair_dt)
  .msg("Computing Lee's L for ", n_pairs, " gene pairs ...", verbose = verbose)

  # Map gene names to row indices
  gene_idx <- stats::setNames(seq_len(nrow(mat)), rownames(mat))

  # Compute Lee's L for all pairs using vectorised inner products
  g1_idx <- gene_idx[pair_dt$gene1]
  g2_idx <- gene_idx[pair_dt$gene2]

  # Check for missing genes
  valid <- !is.na(g1_idx) & !is.na(g2_idx)

  lee_vals <- rep(NA_real_, n_pairs)
  if (any(valid)) {
    # Batch compute: L = n * sum(Wx * Wy) / (n * ||x|| * ||y||)
    #               = sum(Wx * Wy) / (||x|| * ||y||)
    # Use column-wise operations on the index vectors
    vi <- g1_idx[valid]
    vj <- g2_idx[valid]
    # For each valid pair, compute sum of element-wise product of spatial lags
    # Vectorised via rowSums on indexed matrices
    Wx_mat <- Wmat[vi, , drop = FALSE]  # n_valid × cells
    Wy_mat <- Wmat[vj, , drop = FALSE]  # n_valid × cells
    numerator <- rowSums(Wx_mat * Wy_mat)
    denominator <- gene_norms[vi] * gene_norms[vj]
    lee_vals[valid] <- numerator / denominator
  }

  pair_dt[, spatial_lee_L := lee_vals]

  # --- Permutation p-values (vectorised) ------------------------------------
  if (n_perm > 0) {
    .msg("Permutation test (", n_perm, " permutations) ...", verbose = verbose)

    valid_idx <- which(valid)
    n_valid <- length(valid_idx)
    if (n_valid > 0) {
      null_counts <- integer(n_valid)
      obs_abs <- abs(lee_vals[valid_idx])

      for (p in seq_len(n_perm)) {
        # Permute cell order for one variable (gene2 only, more efficient)
        perm <- sample(n)
        mat_c_perm <- mat_c[, perm, drop = FALSE]
        Wmat_perm <- t(as.matrix(W %*% t(mat_c_perm)))

        # Recompute Lee's L for all valid pairs
        Wy_perm <- Wmat_perm[vj, , drop = FALSE]
        numer_perm <- rowSums(Wx_mat * Wy_perm)
        # Denominator uses original norms (permuting doesn't change norms)
        lee_perm <- numer_perm / denominator
        null_counts <- null_counts + as.integer(abs(lee_perm) >= obs_abs)
      }

      pvals <- rep(NA_real_, n_pairs)
      pvals[valid_idx] <- (null_counts + 1) / (n_perm + 1)
      pair_dt[, spatial_lee_p := pvals]
    } else {
      pair_dt[, spatial_lee_p := NA_real_]
    }
  }

  pair_dt
}
