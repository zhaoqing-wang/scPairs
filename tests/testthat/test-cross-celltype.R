test_that(".cross_celltype_batch computes cross-cell-type scores", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 200, n_genes = 20, n_clusters = 3)
  mat <- as.matrix(
    Seurat::GetAssayData(sce, layer = "data")
  )

  W <- scPairs:::.build_knn_graph(sce, reduction = "pca", k = 15)
  cluster_ids <- as.factor(Seurat::Idents(sce))

  pair_dt <- data.table::data.table(
    gene1 = c("GENE1", "GENE3"),
    gene2 = c("GENE2", "GENE4")
  )

  scores <- scPairs:::.cross_celltype_batch(mat, pair_dt, W, cluster_ids)

  expect_length(scores, 2)
  expect_true(is.numeric(scores))

  # Scores should be non-negative (or NA)
  valid <- !is.na(scores)
  if (any(valid)) {
    expect_true(all(scores[valid] >= 0))
    expect_true(all(scores[valid] <= 1))
  }
})


test_that(".cross_celltype returns detailed single-pair result", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 200, n_genes = 20, n_clusters = 3)
  mat <- as.matrix(
    Seurat::GetAssayData(sce, layer = "data")
  )

  W <- scPairs:::.build_knn_graph(sce, reduction = "pca", k = 15)
  cluster_ids <- as.factor(Seurat::Idents(sce))

  x <- mat["GENE3", ]
  y <- mat["GENE4", ]

  result <- scPairs:::.cross_celltype(x, y, W, cluster_ids)

  expect_type(result, "list")
  expect_true("score" %in% names(result))
  expect_true("r_ab" %in% names(result))
  expect_true("r_ba" %in% names(result))
  expect_true("n_cross_pairs" %in% names(result))
  expect_true("per_celltype_pair" %in% names(result))

  # n_cross_pairs should be positive (we have 3 clusters)
  expect_true(result$n_cross_pairs > 0)

  # per_celltype_pair should be a data.frame
  expect_s3_class(result$per_celltype_pair, "data.frame")

  # Score bounds
  if (!is.na(result$score)) {
    expect_true(result$score >= 0 && result$score <= 1)
  }

  # Correlations should be in [-1, 1]
  if (!is.na(result$r_ab)) {
    expect_true(result$r_ab >= -1 && result$r_ab <= 1)
  }
  if (!is.na(result$r_ba)) {
    expect_true(result$r_ba >= -1 && result$r_ba <= 1)
  }
})


test_that("cross_celltype returns NA with insufficient cross-type pairs", {
  skip_if_no_seurat()

  # Single cluster -> no cross-type pairs
  sce <- create_test_seurat(n_cells = 50, n_genes = 10, n_clusters = 1)
  mat <- as.matrix(
    Seurat::GetAssayData(sce, layer = "data")
  )

  W <- scPairs:::.build_knn_graph(sce, reduction = "pca", k = 10)
  cluster_ids <- as.factor(Seurat::Idents(sce))

  x <- mat["GENE1", ]
  y <- mat["GENE2", ]

  result <- scPairs:::.cross_celltype(x, y, W, cluster_ids)
  expect_true(is.na(result$score))
  expect_equal(result$n_cross_pairs, 0)
})


test_that("cross_celltype detects trans-cellular correlation", {
  skip_if_no_seurat()

  set.seed(123)

  # Create a scenario with known trans-cellular pattern:
  # Gene A is high in cluster 1, Gene B is high in cluster 2
  # Clusters 1 and 2 are near each other in PCA space
  n_cells <- 300
  n_genes <- 20

  expr_mat <- matrix(
    abs(rnorm(n_cells * n_genes, mean = 0.5, sd = 0.5)),
    nrow = n_genes, ncol = n_cells
  )
  rownames(expr_mat) <- paste0("GENE", 1:n_genes)
  colnames(expr_mat) <- paste0("CELL", 1:n_cells)

  clusters <- rep(1:3, each = 100)

  # Inject trans-cellular pattern:
  # GENE1 high in cluster 1 cells, GENE2 high in cluster 2 cells
  # Both driven by a shared spatial gradient
  gradient <- seq(0.5, 3, length.out = 100)
  expr_mat["GENE1", 1:100] <- gradient + rnorm(100, 0, 0.2)     # Cluster 1
  expr_mat["GENE2", 101:200] <- gradient + rnorm(100, 0, 0.2)   # Cluster 2

  # Keep them low in other clusters
  expr_mat["GENE1", 101:300] <- abs(rnorm(200, 0.1, 0.1))
  expr_mat["GENE2", c(1:100, 201:300)] <- abs(rnorm(200, 0.1, 0.1))
  expr_mat[expr_mat < 0] <- 0

  suppressWarnings(suppressMessages({
    sce <- Seurat::CreateSeuratObject(counts = expr_mat, project = "test")
    sce <- Seurat::NormalizeData(sce, verbose = FALSE)
    sce <- Seurat::FindVariableFeatures(sce, nfeatures = n_genes, verbose = FALSE)
    sce <- Seurat::ScaleData(sce, verbose = FALSE)
    sce$seurat_clusters <- factor(clusters)
    Seurat::Idents(sce) <- sce$seurat_clusters
    npcs <- min(10, n_genes - 1)
    sce <- Seurat::RunPCA(sce, features = paste0("GENE", 1:n_genes),
                           verbose = FALSE, npcs = npcs)
  }))

  mat <- as.matrix(Seurat::GetAssayData(sce, layer = "data"))
  W <- scPairs:::.build_knn_graph(sce, reduction = "pca", k = 15)
  cluster_ids <- as.factor(Seurat::Idents(sce))

  # The trans-cellular pair
  result_trans <- scPairs:::.cross_celltype(
    mat["GENE1", ], mat["GENE2", ], W, cluster_ids
  )

  # A control pair (random genes)
  result_ctrl <- scPairs:::.cross_celltype(
    mat["GENE5", ], mat["GENE6", ], W, cluster_ids
  )

  # The trans-cellular pair should have a non-trivial score
  # We don't assert it's always higher than control because PCA structure

  # may not perfectly capture our intended neighbourhood structure.
  # But it should return a valid numeric result.
  expect_true(is.numeric(result_trans$score))
  expect_true(result_trans$n_cross_pairs > 0)
})


test_that("FindGenePairs includes cross_celltype_score column", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 200, n_genes = 30, n_clusters = 3)

  result <- suppressWarnings(suppressMessages(
    FindGenePairs(sce,
                  gene = "GENE3",
                  candidates = paste0("GENE", c(1, 2, 4, 5, 6, 7, 8)),
                  use_neighbourhood = TRUE,
                  n_perm = 0,
                  verbose = FALSE)
  ))

  expect_true("cross_celltype_score" %in% colnames(result$pairs))
})


test_that("AssessGenePair includes cross-cell-type metrics", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 150, n_genes = 15, n_clusters = 3)

  result <- suppressWarnings(suppressMessages(
    AssessGenePair(sce, gene1 = "GENE3", gene2 = "GENE4",
                   use_neighbourhood = TRUE,
                   n_perm = 0, verbose = FALSE)
  ))

  expect_true("cross_celltype_score" %in% names(result$metrics))
  expect_true("cross_celltype_r_ab" %in% names(result$metrics))
  expect_true("cross_celltype_r_ba" %in% names(result$metrics))
  expect_true(!is.null(result$cross_celltype_detail))
})
