#' Cross-cell-type interaction metrics
#'
#' Standard co-expression metrics measure whether two genes are expressed
#' together in the **same** cell.  However, many biologically important
#' interactions operate **across** cell types -- gene A expressed in one cell
#' type signals to (or synergises with) gene B expressed in a neighbouring cell
#' of a different type.  Classic examples include Adora2a--Ido1 trans-cellular
#' signalling and ligand--receptor pairs.
#'
#' This module provides a **cross-cell-type interaction score** that detects
#' such trans-cellular synergies.  For every ordered pair of cell types
#' (type_i, type_j where i != j), we measure:
#'
#' \enumerate{
#'   \item The mean expression of gene A in cells of type_i whose neighbours
#'         include cells of type_j.
#'   \item The mean expression of gene B in those type_j neighbour cells.
#'   \item The correlation of these boundary expression levels across all
#'         type_i -> type_j neighbour pairs.
#' }
#'
#' The final **cross-cell-type score** aggregates this signal over all
#' cell-type pairs, capturing whether gene A in one cell type co-varies with
#' gene B in a neighbouring cell type (and vice versa).
#'
#' @name cross-celltype-metrics
#' @keywords internal
NULL


# ===========================================================================
#  1. Cross-cell-type interaction score (batch)
# ===========================================================================

#' Compute cross-cell-type interaction scores for gene pairs
#'
#' For each gene pair (A, B), this metric measures whether expression of
#' gene A in cells of one type correlates with expression of gene B in
#' neighbouring cells of a *different* type.
#'
#' @section Algorithm:
#'
#' 1. Build a KNN graph (from PCA embedding or spatial coordinates).
#' 2. For each cell \eqn{i} of type \eqn{c_i}, identify its neighbours
#'    that belong to a *different* cell type \eqn{c_j \neq c_i}.
#' 3. For each such cross-type neighbour \eqn{j}, record the pair
#'    \eqn{(expr_A(i),\; expr_B(j))}.
#' 4. Compute the Pearson correlation of these cross-type expression
#'    pairs.  This gives the directed score \eqn{A \to B} ("gene A in
#'    the source cell, gene B in the neighbour").
#' 5. Repeat for the reverse direction \eqn{B \to A}.
#' 6. The final score is the geometric mean of the absolute values:
#'    \deqn{cross\_celltype\_score = \sqrt{|r_{A \to B}| \times |r_{B \to A}|}}
#'
#' A high score indicates that expression of gene A in one cell type
#' *systematically* co-varies with expression of gene B in a
#' neighbouring cell type, consistent with trans-cellular signalling
#' or paracrine co-regulation.
#'
#' @param mat Expression matrix (genes x cells, dense or sparse).
#' @param pair_dt \code{data.table} with columns \code{gene1}, \code{gene2}.
#' @param W Sparse KNN weight matrix (n_cells x n_cells), row-standardised.
#' @param cluster_ids Factor of cell-type / cluster assignments (one per cell).
#' @param min_cross_pairs Integer; minimum number of cross-type neighbour
#'   pairs required to compute a correlation.
#'   Pairs with fewer observations get NA.  Default 30.
#'
#' @return Numeric vector of cross-cell-type interaction scores (one per row
#'   of \code{pair_dt}).
#'
#' @keywords internal
.cross_celltype_batch <- function(mat, pair_dt, W, cluster_ids,
                                  min_cross_pairs = 30) {

  common <- intersect(colnames(mat), rownames(W))
  mat <- mat[, common, drop = FALSE]
  W   <- W[common, common, drop = FALSE]
  n   <- length(common)

  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat <- as.matrix(mat)
  }

  # Align cluster_ids with common cells
  if (!is.null(names(cluster_ids))) {
    cluster_ids <- cluster_ids[common]
  } else {
    # Assume same order as original Seurat columns
    cluster_ids <- cluster_ids[seq_len(n)]
  }
  cluster_ids <- as.integer(as.factor(cluster_ids))

  # --- Build cross-type neighbour pairs ------------------------------------
  # For each cell i, identify all neighbours j where cluster(j) != cluster(i)
  # W is row-standardised sparse; non-zero entries indicate neighbours
  W_binary <- methods::as(W > 0, "CsparseMatrix")

  # Extract (i, j) pairs from sparse matrix
  # dgCMatrix stores in compressed column format; convert to triplet for rows
  W_triplet <- methods::as(W_binary, "TsparseMatrix")
  row_idx <- W_triplet@i + 1L    # 0-indexed -> 1-indexed
  col_idx <- W_triplet@j + 1L

  # Keep only cross-type pairs: cluster(i) != cluster(j)
  cross_mask <- cluster_ids[row_idx] != cluster_ids[col_idx]
  src_cells <- row_idx[cross_mask]   # "source" cell (type A)
  nbr_cells <- col_idx[cross_mask]   # "neighbour" cell (type B)

  n_cross <- length(src_cells)

  if (n_cross < min_cross_pairs) {
    return(rep(NA_real_, nrow(pair_dt)))
  }

  # --- Pre-extract expression for cross-type pairs -------------------------
  unique_genes <- unique(c(pair_dt$gene1, pair_dt$gene2))
  unique_genes <- intersect(unique_genes, rownames(mat))
  gene_idx <- stats::setNames(seq_along(unique_genes), unique_genes)

  # Expression at source cells and neighbour cells for all genes
  # These are matrices: n_unique_genes x n_cross_pairs
  expr_src <- mat[unique_genes, src_cells, drop = FALSE]
  expr_nbr <- mat[unique_genes, nbr_cells, drop = FALSE]

  # --- Compute per-pair cross-type correlations ----------------------------
  n_pairs <- nrow(pair_dt)
  scores <- numeric(n_pairs)

  for (pi in seq_len(n_pairs)) {
    g1 <- pair_dt$gene1[pi]
    g2 <- pair_dt$gene2[pi]

    if (!(g1 %in% unique_genes) || !(g2 %in% unique_genes)) {
      scores[pi] <- NA_real_
      next
    }

    gi1 <- gene_idx[g1]
    gi2 <- gene_idx[g2]

    # Direction A->B: gene1 in source cell, gene2 in neighbour cell
    a_src <- expr_src[gi1, ]   # gene1 expression at source cells
    b_nbr <- expr_nbr[gi2, ]  # gene2 expression at neighbour cells

    # Direction B->A: gene2 in source cell, gene1 in neighbour cell
    b_src <- expr_src[gi2, ]
    a_nbr <- expr_nbr[gi1, ]

    # Filter to pairs with some expression (avoid all-zero correlations)
    valid_ab <- (a_src > 0) | (b_nbr > 0)
    valid_ba <- (b_src > 0) | (a_nbr > 0)

    n_ab <- sum(valid_ab)
    n_ba <- sum(valid_ba)

    if (n_ab < min_cross_pairs || n_ba < min_cross_pairs) {
      scores[pi] <- NA_real_
      next
    }

    # Pearson correlations for both directions
    r_ab <- stats::cor(a_src[valid_ab], b_nbr[valid_ab], method = "pearson")
    r_ba <- stats::cor(b_src[valid_ba], a_nbr[valid_ba], method = "pearson")

    if (is.na(r_ab) || is.na(r_ba)) {
      scores[pi] <- NA_real_
      next
    }

    # Geometric mean of absolute correlations
    abs_ab <- abs(r_ab)
    abs_ba <- abs(r_ba)

    if (abs_ab > 0 && abs_ba > 0) {
      scores[pi] <- sqrt(abs_ab * abs_ba)
    } else {
      scores[pi] <- 0
    }
  }

  scores
}


# ===========================================================================
#  2. Cross-cell-type interaction score (single pair)
# ===========================================================================

#' Cross-cell-type interaction score for a single gene pair
#'
#' Measures whether gene A expressed in one cell type correlates with gene B
#' expressed in neighbouring cells of a different type.
#'
#' @param x Numeric vector; expression of gene1 (length n_cells).
#' @param y Numeric vector; expression of gene2 (length n_cells).
#' @param W Sparse KNN weight matrix (n_cells x n_cells).
#' @param cluster_ids Factor/integer of cluster assignments.
#' @param min_cross_pairs Minimum cross-type pairs for correlation.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{score}}{Geometric mean of |r(A->B)| and |r(B->A)|.}
#'     \item{\code{r_ab}}{Pearson r for gene1 in source -> gene2 in neighbour.}
#'     \item{\code{r_ba}}{Pearson r for gene2 in source -> gene1 in neighbour.}
#'     \item{\code{n_cross_pairs}}{Number of cross-type neighbour pairs used.}
#'     \item{\code{per_celltype_pair}}{data.frame of per-celltype-pair results
#'       (only for the single-pair function).}
#'   }
#'
#' @keywords internal
.cross_celltype <- function(x, y, W, cluster_ids, min_cross_pairs = 30) {

  n <- length(x)
  cluster_ids <- as.factor(cluster_ids)
  cl_int <- as.integer(cluster_ids)
  cl_levels <- levels(cluster_ids)

  # Build cross-type neighbour pairs

  W_binary <- methods::as(W > 0, "CsparseMatrix")
  W_triplet <- methods::as(W_binary, "TsparseMatrix")
  row_idx <- W_triplet@i + 1L
  col_idx <- W_triplet@j + 1L

  cross_mask <- cl_int[row_idx] != cl_int[col_idx]
  src_cells <- row_idx[cross_mask]
  nbr_cells <- col_idx[cross_mask]
  n_cross <- length(src_cells)

  if (n_cross < min_cross_pairs) {
    return(list(
      score = NA_real_, r_ab = NA_real_, r_ba = NA_real_,
      n_cross_pairs = n_cross, per_celltype_pair = data.frame()
    ))
  }

  # --- Global cross-type correlations ---
  # A->B: gene1 at source, gene2 at neighbour
  a_src <- x[src_cells]
  b_nbr <- y[nbr_cells]
  valid_ab <- (a_src > 0) | (b_nbr > 0)

  # B->A: gene2 at source, gene1 at neighbour
  b_src <- y[src_cells]
  a_nbr <- x[nbr_cells]
  valid_ba <- (b_src > 0) | (a_nbr > 0)

  r_ab <- if (sum(valid_ab) >= min_cross_pairs) {
    stats::cor(a_src[valid_ab], b_nbr[valid_ab], method = "pearson")
  } else {
    NA_real_
  }

  r_ba <- if (sum(valid_ba) >= min_cross_pairs) {
    stats::cor(b_src[valid_ba], a_nbr[valid_ba], method = "pearson")
  } else {
    NA_real_
  }

  # Composite score
  score <- if (!is.na(r_ab) && !is.na(r_ba) && abs(r_ab) > 0 && abs(r_ba) > 0) {
    sqrt(abs(r_ab) * abs(r_ba))
  } else if (!is.na(r_ab) && !is.na(r_ba)) {
    0
  } else {
    NA_real_
  }

  # --- Per cell-type pair breakdown ---
  src_type <- cl_int[src_cells]
  nbr_type <- cl_int[nbr_cells]
  pair_key <- paste0(src_type, "->", nbr_type)
  unique_keys <- unique(pair_key)

  ct_results <- data.frame(
    source_type     = character(0),
    neighbour_type  = character(0),
    n_pairs         = integer(0),
    r_g1_to_g2      = numeric(0),
    stringsAsFactors = FALSE
  )

  for (uk in unique_keys) {
    mask <- pair_key == uk
    n_k <- sum(mask)
    if (n_k < max(10, min_cross_pairs %/% 3)) next

    parts <- strsplit(uk, "->")[[1]]
    s_type <- cl_levels[as.integer(parts[1])]
    n_type <- cl_levels[as.integer(parts[2])]

    r_k <- tryCatch(
      stats::cor(x[src_cells[mask]], y[nbr_cells[mask]], method = "pearson"),
      error = function(e) NA_real_
    )

    ct_results <- rbind(ct_results, data.frame(
      source_type    = s_type,
      neighbour_type = n_type,
      n_pairs        = n_k,
      r_g1_to_g2     = r_k,
      stringsAsFactors = FALSE
    ))
  }

  list(
    score            = score,
    r_ab             = r_ab,
    r_ba             = r_ba,
    n_cross_pairs    = n_cross,
    per_celltype_pair = ct_results
  )
}
