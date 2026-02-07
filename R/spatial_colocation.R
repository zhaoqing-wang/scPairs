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
#' A binary expression threshold is applied: a cell "expresses" a gene if its
#' normalised expression exceeds `expr_threshold` (default: 0, i.e. any
#' non-zero count).
#'
#' @param coords Data.frame with `x`, `y` columns, rownames = cell barcodes.
#' @param mat Expression matrix (genes x cells).
#' @param pair_dt `data.table` with columns `gene1`, `gene2`.
#' @param k Integer; neighbourhood size.
#' @param expr_threshold Numeric; threshold above which a gene is considered
#'     "expressed" in a cell.
#' @param verbose Logical.
#'
#' @return The input `pair_dt` with added column `spatial_clq`.
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

  common_cells <- intersect(rownames(coords), colnames(mat))
  if (length(common_cells) < 20) {
    pair_dt[, spatial_clq := NA_real_]
    return(pair_dt)
  }

  coords <- coords[common_cells, ]
  mat <- mat[, common_cells, drop = FALSE]
  n <- length(common_cells)

  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }

  k_use <- min(k, n - 1)

  # KNN
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(as.matrix(coords[, c("x", "y")]), k = k_use + 1)
    nn_idx <- nn$nn.idx[, -1, drop = FALSE]
  } else {
    dist_mat <- as.matrix(stats::dist(coords[, c("x", "y")]))
    nn_idx <- t(apply(dist_mat, 1, function(row) order(row)[2:(k_use + 1)]))
  }

  .msg("Computing CLQ for ", nrow(pair_dt), " pairs ...", verbose = verbose)

  clq_vals <- vapply(seq_len(nrow(pair_dt)), function(pi) {
    g1 <- pair_dt$gene1[pi]
    g2 <- pair_dt$gene2[pi]
    if (!(g1 %in% rownames(mat)) || !(g2 %in% rownames(mat))) return(NA_real_)

    expr_a <- mat[g1, ] > expr_threshold
    expr_b <- mat[g2, ] > expr_threshold

    n_a <- sum(expr_a)
    n_b <- sum(expr_b)
    if (n_a == 0 || n_b == 0) return(NA_real_)

    # CLQ A->B
    cells_a <- which(expr_a)
    prop_global_b <- n_b / n
    if (prop_global_b == 0) return(NA_real_)
    clq_ab <- mean(vapply(cells_a, function(i) {
      (sum(expr_b[nn_idx[i, ]]) / k_use) / prop_global_b
    }, numeric(1)))

    # CLQ B->A
    cells_b <- which(expr_b)
    prop_global_a <- n_a / n
    if (prop_global_a == 0) return(NA_real_)
    clq_ba <- mean(vapply(cells_b, function(i) {
      (sum(expr_a[nn_idx[i, ]]) / k_use) / prop_global_a
    }, numeric(1)))

    # Symmetric CLQ = geometric mean
    sqrt(clq_ab * clq_ba)
  }, numeric(1))

  pair_dt[, spatial_clq := clq_vals]
  pair_dt
}
