#' Compute neighbourhood-aware co-expression metrics
#'
#' Standard single-cell co-expression methods require two genes to be
#' detected in the *same* cell.
#' Due to the inherent sparsity of scRNA-seq (dropout, low capture efficiency)
#' and the stochastic nature of transcription, many genuinely synergistic gene
#' pairs are rarely or never co-detected in the same cell -- especially when one
#' or both genes are expressed at low levels.
#'
#' This module addresses the problem by leveraging **neighbourhood information**
#' from the cell-cell similarity graph (typically derived from a UMAP / PCA
#' embedding or spatial coordinates).  Instead of asking "are these two genes
#' detected in the same cell?", we ask:
#'
#' \itemize{
#'   \item **KNN-smoothed correlation** -- After smoothing each gene's expression
#'     over its k nearest neighbours, do the smoothed profiles correlate?
#'   \item **Neighbourhood co-expression score** -- For cells expressing gene A,
#'     are their neighbours enriched for gene B expression (and vice versa)?
#'   \item **Cluster-level (pseudo-bulk) correlation** -- At the cluster
#'     level (aggregated expression), do the two genes correlate?
#' }
#'
#' These three metrics form a complementary evidence layer that is particularly
#' powerful for detecting synergistic pairs where cell-level co-expression is
#' absent due to technical or biological sparsity.
#'
#' @name neighbourhood-metrics
#' @keywords internal
NULL


# ===========================================================================
#  1. KNN graph construction (from PCA / reduction embedding)
# ===========================================================================

#' Build a KNN graph from a Seurat reduction embedding
#'
#' @param object Seurat object.
#' @param reduction Character; reduction to use (default "pca").
#' @param k Integer; number of nearest neighbours (default 20).
#' @param dims Integer vector; dimensions to use (default 1:30).
#' @return A sparse row-standardised weight matrix (n_cells x n_cells).
#' @keywords internal
.build_knn_graph <- function(object,
                             reduction = "pca",
                             k         = 20,
                             dims      = 1:30) {

  embed <- tryCatch(
    Seurat::Embeddings(object, reduction = reduction),
    error = function(e) NULL
  )

  if (is.null(embed)) {
    stop(sprintf("Reduction '%s' not found. Cannot build KNN graph.", reduction),
         call. = FALSE)
  }

  # Limit to available dims
  dims <- dims[dims <= ncol(embed)]
  embed <- embed[, dims, drop = FALSE]
  n <- nrow(embed)

  k_use <- min(k, n - 1)

  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(embed, k = k_use + 1)
    nn_idx <- nn$nn.idx[, -1, drop = FALSE]
  } else {
    dist_mat <- as.matrix(stats::dist(embed))
    nn_idx <- t(apply(dist_mat, 1, function(row) {
      order(row)[2:(k_use + 1)]
    }))
  }

  # Sparse row-standardised weight matrix
  row_i <- rep(seq_len(n), each = k_use)
  col_j <- as.integer(t(nn_idx))
  vals  <- rep(1 / k_use, length(row_i))

  W <- Matrix::sparseMatrix(i = row_i, j = col_j, x = vals, dims = c(n, n))
  rownames(W) <- rownames(embed)
  colnames(W) <- rownames(embed)
  W
}


# ===========================================================================
#  2. KNN-smoothed expression
# ===========================================================================

#' Smooth expression vectors using a KNN weight matrix
#'
#' For each cell, the smoothed expression is a weighted average of its
#' neighbours' expression values:
#' \deqn{\tilde{x}_i = \alpha \cdot x_i + (1 - \alpha) \cdot \sum_{j \in N(i)} w_{ij} x_j}
#'
#' @param mat Expression matrix (genes x cells).
#' @param W Sparse weight matrix (n_cells x n_cells), row-standardised.
#' @param alpha Numeric in \[0, 1\]; self-weight. Default 0.3.
#' @return Smoothed expression matrix (same dimensions).
#' @keywords internal
.smooth_expression <- function(mat, W, alpha = 0.3) {
  # Ensure alignment
  common <- intersect(colnames(mat), rownames(W))
  if (length(common) == 0) {
    return(matrix(numeric(0), nrow = nrow(mat), ncol = 0,
                  dimnames = list(rownames(mat), character(0))))
  }
  mat <- mat[, common, drop = FALSE]
  W   <- W[common, common, drop = FALSE]

  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }

  # Neighbour average: mat %*% t(W) gives genes x cells
  neigh_avg <- mat %*% Matrix::t(W)
  if (inherits(neigh_avg, "sparseMatrix")) {
    neigh_avg <- as.matrix(neigh_avg)
  }

  smoothed <- alpha * mat + (1 - alpha) * neigh_avg
  smoothed
}


# ===========================================================================
#  3. KNN-smoothed correlation (batch)
# ===========================================================================

#' Compute correlation on KNN-smoothed expression for all gene pairs
#'
#' @param mat Expression matrix (genes x cells).
#' @param pair_dt data.table with gene1, gene2 columns.
#' @param W KNN weight matrix.
#' @param alpha Self-weight for smoothing.
#' @param method Correlation method ("pearson" or "spearman").
#' @return Numeric vector of smoothed correlations.
#' @keywords internal
.smoothed_cor_batch <- function(mat, pair_dt, W,
                                alpha  = 0.3,
                                method = "pearson") {
  smoothed <- .smooth_expression(mat, W, alpha = alpha)

  # Ensure dense base matrix (cor() does not accept dgeMatrix)
  if (inherits(smoothed, "Matrix") || inherits(smoothed, "dgeMatrix")) {
    smoothed <- as.matrix(smoothed)
  }

  if (is.null(dim(smoothed))) {
    smoothed <- matrix(smoothed, nrow = nrow(mat), ncol = length(smoothed))
    rownames(smoothed) <- rownames(mat)
  }

  n_pairs <- nrow(pair_dt)
  if (nrow(smoothed) < 2 || ncol(smoothed) < 2) {
    return(rep(NA_real_, n_pairs))
  }

  gene_names <- rownames(smoothed)
  gene_idx <- stats::setNames(seq_along(gene_names), gene_names)

  cors <- rep(NA_real_, n_pairs)

  # Build full correlation matrix for efficiency when many pairs
  if (method == "spearman") {
    smoothed <- t(apply(smoothed, 1, rank))
  }

  cor_mat <- stats::cor(t(smoothed), method = "pearson")

  g1_idx <- gene_idx[pair_dt$gene1]
  g2_idx <- gene_idx[pair_dt$gene2]
  valid <- !is.na(g1_idx) & !is.na(g2_idx)

  cors[valid] <- cor_mat[cbind(g1_idx[valid], g2_idx[valid])]
  cors
}


#' Compute smoothed correlation for a single gene pair
#' @keywords internal
.smoothed_cor <- function(x, y, W, alpha = 0.3) {
  mat <- rbind(x, y)
  rownames(mat) <- c("g1", "g2")
  smoothed <- .smooth_expression(mat, W, alpha = alpha)
  stats::cor(smoothed["g1", ], smoothed["g2", ], method = "pearson")
}


# ===========================================================================
#  4. Neighbourhood co-expression score
# ===========================================================================

#' Neighbourhood co-expression score for gene pairs
#'
#' For each cell expressing gene A, we compute the fraction of its k
#' neighbours that express gene B (and vice versa).  The neighbourhood
#' co-expression score is the geometric mean of the two directional
#' enrichments:
#'
#' \deqn{NCS_{A \to B} = \frac{\text{mean}(\text{neigh\_B\_frac for A-cells})}
#'                              {\text{global\_frac\_B}}}
#' \deqn{NCS = \sqrt{NCS_{A \to B} \cdot NCS_{B \to A}}}
#'
#' A score > 1 indicates that expressing cells of one gene tend to have
#' neighbours expressing the other gene more than expected by chance.
#'
#' @param mat Expression matrix (genes x cells, dense).
#' @param pair_dt data.table with gene1, gene2.
#' @param W KNN weight matrix (row-standardised).
#' @param expr_threshold Numeric; expression above this is "expressed".
#' @return Numeric vector of NCS values.
#' @keywords internal
.neighbourhood_coexpr_batch <- function(mat, pair_dt, W,
                                        expr_threshold = 0) {
  common <- intersect(colnames(mat), rownames(W))
  mat <- mat[, common, drop = FALSE]
  W   <- W[common, common, drop = FALSE]
  n   <- length(common)

  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }

  # Binary expression matrix
  unique_genes <- unique(c(pair_dt$gene1, pair_dt$gene2))
  unique_genes <- intersect(unique_genes, rownames(mat))
  expr_binary <- mat[unique_genes, , drop = FALSE] > expr_threshold
  storage.mode(expr_binary) <- "double"

  # Get k from W (non-zero entries per row)
  k <- max(Matrix::rowSums(W > 0))

  # Build unweighted neighbour matrix (binary, as sparse CsparseMatrix)
  N_binary <- methods::as(W > 0, "CsparseMatrix")

  # Neighbour expression counts: for each gene g and cell i,
  # count how many neighbours of i express g
  neigh_counts <- as.matrix(expr_binary %*% Matrix::t(N_binary))

  # Global expression fractions
  n_expr <- rowSums(expr_binary)

  expr_idx <- stats::setNames(seq_along(unique_genes), unique_genes)

  n_pairs <- nrow(pair_dt)
  ncs_vals <- numeric(n_pairs)

  for (pi in seq_len(n_pairs)) {
    g1 <- pair_dt$gene1[pi]
    g2 <- pair_dt$gene2[pi]

    if (!(g1 %in% unique_genes) || !(g2 %in% unique_genes)) {
      ncs_vals[pi] <- NA_real_
      next
    }

    idx_a <- expr_idx[g1]
    idx_b <- expr_idx[g2]
    n_a <- n_expr[idx_a]
    n_b <- n_expr[idx_b]

    if (n_a == 0 || n_b == 0) {
      ncs_vals[pi] <- NA_real_
      next
    }

    frac_global_b <- n_b / n
    frac_global_a <- n_a / n

    # A->B: for cells expressing A, mean fraction of neighbours expressing B
    cells_a <- which(expr_binary[idx_a, ] > 0)
    if (length(cells_a) > 0 && frac_global_b > 0) {
      ncs_ab <- mean(neigh_counts[idx_b, cells_a] / k) / frac_global_b
    } else {
      ncs_ab <- 0
    }

    # B->A: for cells expressing B, mean fraction of neighbours expressing A
    cells_b <- which(expr_binary[idx_b, ] > 0)
    if (length(cells_b) > 0 && frac_global_a > 0) {
      ncs_ba <- mean(neigh_counts[idx_a, cells_b] / k) / frac_global_a
    } else {
      ncs_ba <- 0
    }

    # Geometric mean
    if (ncs_ab > 0 && ncs_ba > 0) {
      ncs_vals[pi] <- sqrt(ncs_ab * ncs_ba)
    } else {
      ncs_vals[pi] <- 0
    }
  }

  ncs_vals
}


#' Single-pair neighbourhood co-expression score
#' @keywords internal
.neighbourhood_coexpr <- function(x, y, W, k = NULL) {
  n <- length(x)
  if (is.null(k)) k <- max(Matrix::rowSums(W > 0))

  # Build binary neighbour matrix
  N_binary <- methods::as(W > 0, "CsparseMatrix")

  expr_a <- as.double(x > 0)
  expr_b <- as.double(y > 0)

  n_a <- sum(expr_a)
  n_b <- sum(expr_b)

  if (n_a == 0 || n_b == 0) return(NA_real_)

  # Neighbour counts
  neigh_b <- as.numeric(N_binary %*% expr_b)
  neigh_a <- as.numeric(N_binary %*% expr_a)

  frac_global_b <- n_b / n
  frac_global_a <- n_a / n

  cells_a <- which(expr_a > 0)
  cells_b <- which(expr_b > 0)

  ncs_ab <- mean(neigh_b[cells_a] / k) / frac_global_b
  ncs_ba <- mean(neigh_a[cells_b] / k) / frac_global_a

  if (ncs_ab > 0 && ncs_ba > 0) sqrt(ncs_ab * ncs_ba) else 0
}


# ===========================================================================
#  5. Cluster-level (pseudo-bulk) correlation
# ===========================================================================

#' Cluster-level (pseudo-bulk) correlation for gene pairs
#'
#' Computes the Pearson correlation of cluster-level mean expression.
#' This captures patterns where two genes tend to be expressed in the same
#' cell populations, even if they are never co-detected in the same cell.
#'
#' @param mat Expression matrix (genes x cells).
#' @param cluster_ids Factor of cluster assignments.
#' @param pair_dt data.table with gene1, gene2.
#' @return Numeric vector of cluster-level correlations.
#' @keywords internal
.cluster_cor_batch <- function(mat, cluster_ids, pair_dt) {
  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }

  clusters <- levels(as.factor(cluster_ids))
  n_clusters <- length(clusters)

  if (n_clusters < 3) return(rep(NA_real_, nrow(pair_dt)))

  # Build cluster-mean matrix (genes x clusters)
  cl_indices <- split(seq_along(cluster_ids), cluster_ids)

  unique_genes <- unique(c(pair_dt$gene1, pair_dt$gene2))
  unique_genes <- intersect(unique_genes, rownames(mat))
  sub_mat <- mat[unique_genes, , drop = FALSE]

  cl_means <- matrix(0, nrow = length(unique_genes), ncol = n_clusters,
                     dimnames = list(unique_genes, clusters))

  for (ci in seq_len(n_clusters)) {
    idx <- cl_indices[[ci]]
    if (length(idx) > 0) {
      cl_means[, ci] <- rowMeans(sub_mat[, idx, drop = FALSE])
    }
  }

  # Correlation matrix of cluster means
  cl_cor <- stats::cor(t(cl_means), method = "pearson")

  gene_idx <- stats::setNames(seq_along(unique_genes), unique_genes)
  n_pairs <- nrow(pair_dt)
  cors <- rep(NA_real_, n_pairs)

  g1_idx <- gene_idx[pair_dt$gene1]
  g2_idx <- gene_idx[pair_dt$gene2]
  valid <- !is.na(g1_idx) & !is.na(g2_idx)

  cors[valid] <- cl_cor[cbind(g1_idx[valid], g2_idx[valid])]
  cors
}


#' Single-pair cluster-level correlation
#' @keywords internal
.cluster_cor <- function(x, y, cluster_ids) {
  clusters <- levels(as.factor(cluster_ids))
  n_cl <- length(clusters)
  if (n_cl < 3) return(NA_real_)

  cl_mean_x <- tapply(x, cluster_ids, mean)
  cl_mean_y <- tapply(y, cluster_ids, mean)

  # Remove clusters with no expression in either
  valid <- !is.na(cl_mean_x) & !is.na(cl_mean_y)
  if (sum(valid) < 3) return(NA_real_)

  stats::cor(cl_mean_x[valid], cl_mean_y[valid], method = "pearson")
}
