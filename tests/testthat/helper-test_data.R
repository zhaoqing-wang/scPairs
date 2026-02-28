# Helper functions for scPairs tests
#
# Most tests use the built-in `scpairs_testdata` object (100 cells × 20 genes,
# 3 clusters, PCA + UMAP, with GENE3/GENE4 globally co-expressed and
# GENE1/GENE2 cluster-1-specific).  Use `scpairs_testdata` directly in tests
# that require a standard Seurat object.
#
# `create_test_seurat()` is retained for tests that need a custom configuration
# (e.g. a single cluster, an unusually small dataset, or a spatial object).


#' Create a small custom Seurat object for edge-case tests
#'
#' @param n_cells Integer; total cells.
#' @param n_genes Integer; total genes.
#' @param n_clusters Integer; number of clusters.
#' @return A Seurat object with PCA and UMAP reductions.
#' @keywords internal
create_test_seurat <- function(n_cells = 100, n_genes = 20, n_clusters = 3) {
  suppressWarnings(suppressMessages({
    expr_mat <- matrix(
      abs(rnorm(n_cells * n_genes, mean = 0.5, sd = 0.5)),
      nrow = n_genes, ncol = n_cells,
      dimnames = list(paste0("GENE", seq_len(n_genes)),
                      paste0("CELL", seq_len(n_cells)))
    )

    clusters <- rep(seq_len(n_clusters), length.out = n_cells)

    # Inject co-expression: GENE3 & GENE4 globally
    global_signal <- abs(rnorm(n_cells, mean = 2, sd = 1))
    expr_mat["GENE3", ] <- global_signal + rnorm(n_cells, 0, 0.2)
    if (n_genes >= 4)
      expr_mat["GENE4", ] <- global_signal * 1.1 + rnorm(n_cells, 0, 0.2)

    # GENE1 & GENE2 cluster-1-specific
    cl1 <- which(clusters == 1)
    shared <- abs(rnorm(length(cl1), mean = 3, sd = 0.5))
    expr_mat["GENE1", cl1] <- shared + rnorm(length(cl1), 0, 0.3)
    if (n_genes >= 2)
      expr_mat["GENE2", cl1] <- shared * 0.8 + rnorm(length(cl1), 0, 0.3)

    expr_mat[expr_mat < 0] <- 0

    sce <- Seurat::CreateSeuratObject(counts = expr_mat, project = "test_scPairs")
    sce <- Seurat::NormalizeData(sce, verbose = FALSE)
    sce <- Seurat::FindVariableFeatures(sce,
                                        selection.method = "vst",
                                        nfeatures        = min(2000, n_genes),
                                        verbose          = FALSE)
    sce <- Seurat::ScaleData(sce,
                             features = Seurat::VariableFeatures(sce),
                             verbose  = FALSE)

    sce$seurat_clusters <- factor(clusters)
    Seurat::Idents(sce) <- sce$seurat_clusters

    var_features <- Seurat::VariableFeatures(sce)
    if (length(var_features) > 0) {
      npcs <- min(10, n_genes - 1, length(var_features))
      sce  <- Seurat::RunPCA(sce, features = var_features,
                              npcs = npcs, verbose = FALSE)
      sce  <- Seurat::RunUMAP(sce, dims = seq_len(npcs), verbose = FALSE)
    }

    sce
  }))
}


#' Create a test Seurat object with synthetic spatial coordinates
#' @keywords internal
create_test_spatial_seurat <- function(n_cells = 100, n_genes = 20) {
  sce <- create_test_seurat(n_cells = n_cells, n_genes = n_genes, n_clusters = 3)
  sce$x_centroid <- runif(n_cells, 0, 100)
  sce$y_centroid <- runif(n_cells, 0, 100)
  sce
}


#' Skip the current test if Seurat is not installed
#' @keywords internal
skip_if_no_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE))
    testthat::skip("Seurat not available")
}


#' Skip the current test when running on CRAN
#' @keywords internal
skip_if_cran <- function() {
  if (!identical(Sys.getenv("NOT_CRAN"), "true"))
    testthat::skip("Skipping on CRAN")
}
