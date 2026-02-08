#' Spatial Transcriptomics Metrics for scPairs
#'
#' @description
#' This module combines bivariate spatial autocorrelation (Lee's L) and
#' co-location quotient (CLQ) metrics for spatial transcriptomics data.
#' Both metrics measure different aspects of spatial co-expression:
#'
#' * **Lee's L** captures spatial co-variation of continuous expression values.
#' * **CLQ** measures whether expressing cells are spatially attracted.
#'
#' @name spatial-metrics
#' @keywords internal
NULL


# ===========================================================================
#  Shared: Build spatial KNN weight matrix
# ===========================================================================

#' Build a spatial KNN weight matrix
#'
#' Shared helper used by both Lee's L and CLQ computations.
#'
#' @param coords Data.frame with \code{x}, \code{y} columns.
#' @param k Integer; number of spatial nearest neighbours.
#' @param row_standardise Logical; if TRUE, return row-standardised weights.
#'     If FALSE, return binary indicator matrix.
#' @return A sparse matrix (n x n).
#' @keywords internal
.build_spatial_knn <- function(coords, k, row_standardise = TRUE) {
  n <- nrow(coords)
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

  row_i <- rep(seq_len(n), each = k_use)
  col_j <- as.integer(t(nn_idx))

  if (row_standardise) {
    vals <- rep(1 / k_use, length(row_i))
  } else {
    vals <- rep(1, length(row_i))
  }

  Matrix::sparseMatrix(i = row_i, j = col_j, x = vals, dims = c(n, n))
}


#' Align spatial coordinates and expression matrix
#'
#' @param coords Data.frame with x, y columns, rownames = cell barcodes.
#' @param mat Expression matrix (genes x cells).
#' @param min_cells Integer; minimum cells required.
#' @return A list with aligned \code{coords}, \code{mat}, and \code{n}.
#'     Returns NULL if insufficient cells.
#' @keywords internal
.align_spatial_data <- function(coords, mat, min_cells = 20) {
  common_cells <- intersect(rownames(coords), colnames(mat))
  if (length(common_cells) < min_cells) {
    return(NULL)
  }

  coords <- coords[common_cells, ]
  mat <- mat[, common_cells, drop = FALSE]

  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }

  list(coords = coords, mat = mat, n = length(common_cells))
}


# ===========================================================================
#  1. Lee's L statistic (bivariate spatial autocorrelation)
# ===========================================================================

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
#' @param coords Data.frame with columns \code{x}, \code{y} and
#'     rownames = cell barcodes.
#' @param mat Expression matrix (genes x cells).
#' @param pair_dt \code{data.table} with columns \code{gene1}, \code{gene2}.
#' @param k Integer; number of spatial nearest neighbours (default 6 for
#'     Visium hexagonal grids; 15 for general ST).
#' @param n_perm Integer; number of permutations for empirical p-values.
#'     Set to 0 to skip.
#' @param verbose Logical.
#'
#' @return The input \code{pair_dt} with added columns \code{spatial_lee_L}
#'     and optionally \code{spatial_lee_p}.
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

  aligned <- .align_spatial_data(coords, mat, min_cells = 20)
  if (is.null(aligned)) {
    warning("Fewer than 20 cells with both spatial coords and expression. ",
            "Skipping spatial analysis.", call. = FALSE)
    pair_dt[, spatial_lee_L := NA_real_]
    if (n_perm > 0) pair_dt[, spatial_lee_p := NA_real_]
    return(pair_dt)
  }

  coords <- aligned$coords
  mat <- aligned$mat
  n <- aligned$n

  .msg("Computing spatial weight matrix (k=", k, ", n=", n, " spots) ...",
       verbose = verbose)

  # KNN spatial weight matrix (row-standardised)
  W <- .build_spatial_knn(coords, k = k, row_standardise = TRUE)

  # --- Vectorised Lee's L computation for all pairs at once -----------------
  # Centre each gene across cells
  gene_means <- rowMeans(mat)
  mat_c <- mat - gene_means  # centred (genes x cells)

  # Spatial lag for all genes: Wx = W %*% t(mat_c) -> cells x genes
  # Then transpose to get genes x cells
  Wmat <- t(as.matrix(W %*% t(mat_c)))  # genes x cells

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
    vi <- g1_idx[valid]
    vj <- g2_idx[valid]
    Wx_mat <- Wmat[vi, , drop = FALSE]
    Wy_mat <- Wmat[vj, , drop = FALSE]
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
      vi <- g1_idx[valid_idx]
      vj <- g2_idx[valid_idx]
      Wx_mat <- Wmat[vi, , drop = FALSE]
      denominator <- gene_norms[vi] * gene_norms[vj]

      for (p in seq_len(n_perm)) {
        perm <- sample(n)
        mat_c_perm <- mat_c[, perm, drop = FALSE]
        Wmat_perm <- t(as.matrix(W %*% t(mat_c_perm)))

        Wy_perm <- Wmat_perm[vj, , drop = FALSE]
        numer_perm <- rowSums(Wx_mat * Wy_perm)
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


# ===========================================================================
#  2. Co-location quotient (CLQ)
# ===========================================================================

#' Compute co-location quotient (CLQ) for gene pairs in spatial data
#'
#' The co-location quotient measures whether cells expressing gene A are
#' disproportionately located near cells expressing gene B, relative to a
#' random spatial arrangement.  CLQ > 1 indicates spatial co-location (the
#' two expression patterns are spatially attracted); CLQ < 1 indicates spatial
#' segregation; CLQ = 1 indicates random spatial mixing.
#'
#' @details
#' For each cell *i* expressing gene A, we count how many of its *k* nearest
#' neighbours express gene B, and compare to the global proportion of cells
#' expressing gene B.  The CLQ is the ratio of observed-to-expected proportions:
#'
#' \deqn{CLQ_{A \to B} = \frac{1}{N_A} \sum_{i \in A} \frac{n_{iB}/k}{N_B / N}}
#'
#' where \eqn{N_A} and \eqn{N_B} are numbers of cells expressing A and B,
#' \eqn{n_{iB}} is the number of B-expressing neighbours of cell *i*, *k* is
#' the neighbourhood size, and *N* is the total number of cells.
#'
#' The symmetric CLQ is computed as the geometric mean of \eqn{CLQ_{A \to B}}
#' and \eqn{CLQ_{B \to A}}.
#'
#' **Performance (v0.1.1):** The neighbour expression counts are computed via
#' matrix multiplication.
#'
#' @param coords Data.frame with \code{x}, \code{y} columns,
#'     rownames = cell barcodes.
#' @param mat Expression matrix (genes x cells).
#' @param pair_dt \code{data.table} with columns \code{gene1}, \code{gene2}.
#' @param k Integer; neighbourhood size.
#' @param expr_threshold Numeric; threshold above which a gene is considered
#'     "expressed" in a cell.
#' @param verbose Logical.
#'
#' @return The input \code{pair_dt} with added column \code{spatial_clq}.
#'
#' @references
#' Leslie, T.F. & Kronenfeld, B.J. (2011). The colocation quotient: a new
#' measure of spatial association between categorical subsets of points.
#' *Geographical Analysis*, 43(3), 306-326.
#'
#' @keywords internal
.compute_spatial_clq <- function(coords,
                                  mat,
                                  pair_dt,
                                  k              = 6,
                                  expr_threshold = 0,
                                  verbose        = TRUE) {

  if (nrow(pair_dt) == 0) return(pair_dt)

  aligned <- .align_spatial_data(coords, mat, min_cells = 20)
  if (is.null(aligned)) {
    pair_dt[, spatial_clq := NA_real_]
    return(pair_dt)
  }

  coords <- aligned$coords
  mat <- aligned$mat
  n <- aligned$n
  k_use <- min(k, n - 1)

  .msg("Computing CLQ for ", nrow(pair_dt), " pairs ...", verbose = verbose)

  # Build spatial KNN (binary, not row-standardised)
  N_sparse <- .build_spatial_knn(coords, k = k, row_standardise = FALSE)

  # Pre-compute binary expression for all genes in pair_dt
  unique_genes <- unique(c(pair_dt$gene1, pair_dt$gene2))
  unique_genes <- intersect(unique_genes, rownames(mat))

  expr_binary <- mat[unique_genes, , drop = FALSE] > expr_threshold
  storage.mode(expr_binary) <- "double"

  # Neighbour expression counts for all genes at once
  neigh_counts <- as.matrix(expr_binary %*% Matrix::t(N_sparse))

  # Global expression counts per gene
  n_expr <- rowSums(expr_binary)
  expr_idx <- stats::setNames(seq_along(unique_genes), unique_genes)

  n_pairs <- nrow(pair_dt)
  clq_vals <- numeric(n_pairs)

  for (pi in seq_len(n_pairs)) {
    g1 <- pair_dt$gene1[pi]
    g2 <- pair_dt$gene2[pi]
    if (!(g1 %in% unique_genes) || !(g2 %in% unique_genes)) {
      clq_vals[pi] <- NA_real_
      next
    }

    idx_a <- expr_idx[g1]
    idx_b <- expr_idx[g2]
    n_a <- n_expr[idx_a]
    n_b <- n_expr[idx_b]

    if (n_a == 0 || n_b == 0) {
      clq_vals[pi] <- NA_real_
      next
    }

    prop_global_b <- n_b / n
    prop_global_a <- n_a / n

    cells_a <- which(expr_binary[idx_a, ] > 0)
    clq_ab <- mean(neigh_counts[idx_b, cells_a] / k_use) / prop_global_b

    cells_b <- which(expr_binary[idx_b, ] > 0)
    clq_ba <- mean(neigh_counts[idx_a, cells_b] / k_use) / prop_global_a

    clq_vals[pi] <- sqrt(clq_ab * clq_ba)
  }

  pair_dt[, spatial_clq := clq_vals]
  pair_dt
}
