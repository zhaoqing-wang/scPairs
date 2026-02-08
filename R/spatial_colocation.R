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
#' **Performance (v0.1.1):** The neighbour expression counts are computed via
#' matrix multiplication.  For each gene, the binary expression vector is
#' multiplied against the neighbour indicator matrix to get per-cell counts of
#' expressing neighbours, replacing the previous nested vapply loops.
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

  # --- Build sparse neighbour indicator matrix (n × n) ----------------------
  # N[i,j] = 1 if j is a neighbour of i
  row_i <- rep(seq_len(n), each = k_use)
  col_j <- as.integer(t(nn_idx))
  N_sparse <- Matrix::sparseMatrix(i = row_i, j = col_j, x = 1,
                                   dims = c(n, n))

  # Map gene names to row indices
  gene_idx <- stats::setNames(seq_len(nrow(mat)), rownames(mat))

  # Pre-compute binary expression for all genes in pair_dt
  unique_genes <- unique(c(pair_dt$gene1, pair_dt$gene2))
  unique_genes <- intersect(unique_genes, rownames(mat))

  # Binary expression matrix (genes × cells): TRUE if expressed
  expr_binary <- mat[unique_genes, , drop = FALSE] > expr_threshold
  storage.mode(expr_binary) <- "double"

  # Neighbour expression counts for all genes at once:
  # For each gene g and cell i, count how many neighbours of i express g
  # = expr_binary %*% t(N_sparse) → genes × cells
  # But N_sparse is cells × cells, so:
  # neigh_counts[g, i] = sum_{j in N(i)} expr_binary[g, j]
  # = (expr_binary %*% t(N_sparse))[g, i]
  neigh_counts <- as.matrix(expr_binary %*% Matrix::t(N_sparse))

  # Global expression counts per gene
  n_expr <- rowSums(expr_binary)

  # Gene name to row index in expr_binary
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

    # CLQ A→B: for cells expressing A, mean(neighbour_B_count / k) / prop_global_B
    cells_a <- which(expr_binary[idx_a, ] > 0)
    clq_ab <- mean(neigh_counts[idx_b, cells_a] / k_use) / prop_global_b

    # CLQ B→A: for cells expressing B, mean(neighbour_A_count / k) / prop_global_A
    cells_b <- which(expr_binary[idx_b, ] > 0)
    clq_ba <- mean(neigh_counts[idx_a, cells_b] / k_use) / prop_global_a

    clq_vals[pi] <- sqrt(clq_ab * clq_ba)
  }

  pair_dt[, spatial_clq := clq_vals]
  pair_dt
}
