# Helper functions and test data creation for scPairs tests

#' Create a small test Seurat object for scPairs testing
#'
#' Generates a synthetic Seurat object with known co-expression patterns:
#' - GENE1 and GENE2 are strongly positively correlated in cluster 1
#' - GENE3 and GENE4 are correlated across all clusters
#' - Other genes are noise
#'
#' @param n_cells Integer; total cells.
#' @param n_genes Integer; total genes.
#' @param n_clusters Integer; number of clusters.
#' @return A Seurat object with PCA and UMAP.
#' @keywords internal
create_test_seurat <- function(n_cells = 200, n_genes = 30, n_clusters = 3) {
  suppressWarnings(suppressMessages({
    set.seed(42)

    # Base expression
    expr_mat <- matrix(
      abs(rnorm(n_cells * n_genes, mean = 0.5, sd = 0.5)),
      nrow = n_genes, ncol = n_cells
    )
    rownames(expr_mat) <- paste0("GENE", 1:n_genes)
    colnames(expr_mat) <- paste0("CELL", 1:n_cells)

    clusters <- rep(1:n_clusters, length.out = n_cells)

    # Inject co-expression: GENE1 & GENE2 correlated in cluster 1
    cl1 <- which(clusters == 1)
    shared_signal <- abs(rnorm(length(cl1), mean = 3, sd = 0.5))
    expr_mat["GENE1", cl1] <- shared_signal + rnorm(length(cl1), 0, 0.3)
    expr_mat["GENE2", cl1] <- shared_signal * 0.8 + rnorm(length(cl1), 0, 0.3)

    # GENE3 & GENE4 globally correlated
    global_signal <- abs(rnorm(n_cells, mean = 2, sd = 1))
    expr_mat["GENE3", ] <- global_signal + rnorm(n_cells, 0, 0.2)
    expr_mat["GENE4", ] <- global_signal * 1.1 + rnorm(n_cells, 0, 0.2)

    # Ensure non-negative
    expr_mat[expr_mat < 0] <- 0

    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat package required for tests")
    }

    sce <- Seurat::CreateSeuratObject(counts = expr_mat, project = "test_scPairs")
    sce <- Seurat::NormalizeData(sce, verbose = FALSE)
    sce <- Seurat::FindVariableFeatures(sce, selection.method = "vst",
                                         nfeatures = min(2000, n_genes),
                                         verbose = FALSE)
    sce <- Seurat::ScaleData(sce,
                              features = Seurat::VariableFeatures(sce),
                              verbose = FALSE)

    sce$seurat_clusters <- factor(clusters)
    Seurat::Idents(sce) <- sce$seurat_clusters

    var_features <- Seurat::VariableFeatures(sce)
    if (length(var_features) > 0) {
      npcs_use <- min(10, n_genes - 1, length(var_features))
      sce <- Seurat::RunPCA(sce, features = var_features,
                             verbose = FALSE, npcs = npcs_use)
      sce <- Seurat::RunUMAP(sce, dims = 1:npcs_use, verbose = FALSE)
    }

    return(sce)
  }))
}


#' Create a test Seurat object with spatial coordinates
#' @keywords internal
create_test_spatial_seurat <- function(n_cells = 100, n_genes = 20) {
  sce <- create_test_seurat(n_cells = n_cells, n_genes = n_genes, n_clusters = 3)

  # Add spatial coordinates to meta.data
  set.seed(42)
  sce$x_centroid <- runif(n_cells, 0, 100)
  sce$y_centroid <- runif(n_cells, 0, 100)

  sce
}


#' Skip tests if Seurat is not available
#' @keywords internal
skip_if_no_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    testthat::skip("Seurat not available")
  }
}

#' Skip on CRAN
#' @keywords internal
skip_if_cran <- function() {
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    testthat::skip("Skipping on CRAN")
  }
}
