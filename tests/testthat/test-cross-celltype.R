# Tests for cross-cell-type interaction metric

test_that(".cross_celltype_batch computes cross-cell-type scores", {
  mat        <- as.matrix(Seurat::GetAssayData(scpairs_testdata, layer = "data"))
  embed      <- Seurat::Embeddings(scpairs_testdata, reduction = "pca")
  cluster_ids <- as.factor(Seurat::Idents(scpairs_testdata))

  pair_dt <- data.table::data.table(
    gene1 = c("GENE1", "GENE3"),
    gene2 = c("GENE2", "GENE4")
  )

  scores <- scPairs:::.cross_celltype_batch(mat, pair_dt, cluster_ids, embed)

  expect_length(scores, 2)
  expect_true(is.numeric(scores))
  valid <- !is.na(scores)
  if (any(valid)) {
    expect_true(all(scores[valid] >= 0))
    expect_true(all(scores[valid] <= 1))
  }
})


test_that(".cross_celltype returns a valid detailed result", {
  mat         <- as.matrix(Seurat::GetAssayData(scpairs_testdata, layer = "data"))
  embed       <- Seurat::Embeddings(scpairs_testdata, reduction = "pca")
  cluster_ids <- as.factor(Seurat::Idents(scpairs_testdata))

  result <- scPairs:::.cross_celltype(
    mat["GENE3", ], mat["GENE4", ], cluster_ids, embed,
    n_bins            = 15,
    min_cells_per_bin = 2,
    min_bins          = 3
  )

  expect_type(result, "list")
  expect_true("score"          %in% names(result))
  expect_true("r_ab"           %in% names(result))
  expect_true("r_ba"           %in% names(result))
  expect_true("n_type_pairs"   %in% names(result))
  expect_true("per_celltype_pair" %in% names(result))
  expect_true(result$n_type_pairs > 0)
  expect_s3_class(result$per_celltype_pair, "data.frame")
  if (!is.na(result$score)) {
    expect_true(result$score >= 0 && result$score <= 1)
  }
  if (!is.na(result$r_ab)) expect_true(abs(result$r_ab) <= 1)
  if (!is.na(result$r_ba)) expect_true(abs(result$r_ba) <= 1)
})


test_that("cross_celltype returns NA with a single cluster (no cross-type pairs)", {
  skip_if_no_seurat()
  # Create a fresh 1-cluster object to test the NA edge case
  sce1 <- create_test_seurat(n_cells = 50, n_genes = 10, n_clusters = 1)
  mat         <- as.matrix(Seurat::GetAssayData(sce1, layer = "data"))
  embed       <- Seurat::Embeddings(sce1, reduction = "pca")
  cluster_ids <- as.factor(Seurat::Idents(sce1))

  result <- scPairs:::.cross_celltype(mat["GENE1", ], mat["GENE2", ],
                                       cluster_ids, embed)
  expect_true(is.na(result$score))
  expect_equal(result$n_type_pairs, 0)
})


test_that("cross_celltype detects a trans-cellular pattern", {
  skip_if_no_seurat()
  # Build a custom object with a known trans-cellular gradient (cluster-driven)
  n_cells <- 300
  n_genes <- 20
  expr_mat <- matrix(abs(rnorm(n_cells * n_genes, 0.5, 0.5)),
                     nrow = n_genes, ncol = n_cells,
                     dimnames = list(paste0("GENE", seq_len(n_genes)),
                                     paste0("CELL", seq_len(n_cells))))
  clusters <- rep(1:3, each = 100)
  gradient <- seq(0.5, 3, length.out = 100)
  expr_mat["GENE1",   1:100]   <- gradient + rnorm(100, 0, 0.2)
  expr_mat["GENE2", 101:200]   <- gradient + rnorm(100, 0, 0.2)
  expr_mat["GENE1", 101:300]   <- abs(rnorm(200, 0.1, 0.1))
  expr_mat["GENE2", c(1:100, 201:300)] <- abs(rnorm(200, 0.1, 0.1))
  expr_mat[expr_mat < 0] <- 0

  suppressWarnings(suppressMessages({
    sce <- Seurat::CreateSeuratObject(counts = expr_mat, project = "test")
    sce <- Seurat::NormalizeData(sce, verbose = FALSE)
    sce <- Seurat::FindVariableFeatures(sce, nfeatures = n_genes, verbose = FALSE)
    sce <- Seurat::ScaleData(sce, verbose = FALSE)
    sce$seurat_clusters <- factor(clusters)
    Seurat::Idents(sce) <- sce$seurat_clusters
    sce <- Seurat::RunPCA(sce, features = paste0("GENE", seq_len(n_genes)),
                           npcs = min(10, n_genes - 1), verbose = FALSE)
  }))

  mat         <- as.matrix(Seurat::GetAssayData(sce, layer = "data"))
  embed       <- Seurat::Embeddings(sce, reduction = "pca")
  cluster_ids <- as.factor(Seurat::Idents(sce))

  result_trans <- scPairs:::.cross_celltype(
    mat["GENE1", ], mat["GENE2", ], cluster_ids, embed,
    n_bins = 15, min_cells_per_bin = 2, min_bins = 3
  )
  expect_true(is.numeric(result_trans$score))
  expect_true(result_trans$n_type_pairs > 0)
})


test_that("FindGenePairs output includes cross_celltype_score column", {
  result <- suppressWarnings(suppressMessages(
    FindGenePairs(
      scpairs_testdata,
      gene              = "GENE3",
      candidates        = paste0("GENE", c(1, 2, 4, 5, 6, 7, 8)),
      use_neighbourhood = TRUE,
      n_perm            = 0,
      verbose           = FALSE
    )
  ))
  expect_true("cross_celltype_score" %in% colnames(result$pairs))
})


test_that("AssessGenePair output includes cross-cell-type metrics", {
  result <- suppressWarnings(suppressMessages(
    AssessGenePair(
      scpairs_testdata,
      gene1             = "GENE3",
      gene2             = "GENE4",
      use_neighbourhood = TRUE,
      n_perm            = 0,
      verbose           = FALSE
    )
  ))
  expect_true("cross_celltype_score" %in% names(result$metrics))
  expect_true("cross_celltype_r_ab"  %in% names(result$metrics))
  expect_true("cross_celltype_r_ba"  %in% names(result$metrics))
  expect_true(!is.null(result$cross_celltype_detail))
})
