test_that(".bicor computes correct biweight midcorrelation", {
  set.seed(42)
  x <- rnorm(100)
  y <- x + rnorm(100, sd = 0.1)

  bc <- scPairs:::.bicor(x, y)
  expect_true(is.numeric(bc))
  expect_true(bc > 0.8)  # should be strongly positive
  expect_true(bc <= 1)

  # Constant vector
  bc_const <- scPairs:::.bicor(rep(1, 100), rnorm(100))
  expect_equal(bc_const, 0)
})


test_that(".mutual_info is non-negative", {
  set.seed(42)
  x <- rnorm(200)
  y <- x^2 + rnorm(200, sd = 0.1)

  mi <- scPairs:::.mutual_info(x, y, n_bins = 5)
  expect_true(mi >= 0)

  # Independent variables should have low MI
  mi_indep <- scPairs:::.mutual_info(rnorm(200), rnorm(200), n_bins = 5)
  expect_true(mi > mi_indep)
})


test_that(".ratio_consistency returns bounded value", {
  set.seed(42)
  n <- 100
  x <- rnorm(n, 2, 0.1)
  y <- x * 1.5 + rnorm(n, 0, 0.05)  # constant ratio
  cl <- factor(rep(1:3, length.out = n))

  rc <- scPairs:::.ratio_consistency(x, y, cl)
  expect_true(rc >= 0 && rc <= 1)
  expect_true(rc > 0.5)  # should be high for constant ratio
})


test_that(".validate_seurat catches non-Seurat input", {
  expect_error(scPairs:::.validate_seurat("not a seurat"),
               "Seurat object")
  expect_error(scPairs:::.validate_seurat(list()),
               "Seurat object")
})


test_that(".has_spatial detects spatial columns", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 50, n_genes = 10)
  expect_false(scPairs:::.has_spatial(sce))

  sce$x_centroid <- runif(50)
  sce$y_centroid <- runif(50)
  expect_true(scPairs:::.has_spatial(sce))
})
