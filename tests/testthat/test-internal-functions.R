# Tests for core mathematical helper functions (R/coexpression.R, R/validation.R)

test_that(".bicor computes biweight midcorrelation correctly", {
  x <- rnorm(100)
  y <- x + rnorm(100, sd = 0.1)

  bc <- scPairs:::.bicor(x, y)
  expect_true(is.numeric(bc) && length(bc) == 1)
  expect_true(bc > 0.8)   # strongly positive for correlated pair
  expect_true(bc <= 1)

  # Constant input vector returns 0 (no variance)
  expect_equal(scPairs:::.bicor(rep(1, 100), rnorm(100)), 0)
})


test_that(".mutual_info is non-negative and higher for dependent variables", {
  x <- rnorm(200)
  y <- x^2 + rnorm(200, sd = 0.1)

  mi_dep   <- scPairs:::.mutual_info(x, y,         n_bins = 5)
  mi_indep <- scPairs:::.mutual_info(x, rnorm(200), n_bins = 5)

  expect_true(mi_dep   >= 0)
  expect_true(mi_indep >= 0)
  expect_true(mi_dep   > mi_indep)
})


test_that(".ratio_consistency returns a value in [0, 1] and is high for constant-ratio data", {
  n  <- 100
  x  <- rnorm(n, 2, 0.1)
  y  <- x * 1.5 + rnorm(n, 0, 0.05)   # near-constant ratio
  cl <- factor(rep(1:3, length.out = n))

  rc <- scPairs:::.ratio_consistency(x, y, cl)
  expect_true(rc >= 0 && rc <= 1)
  expect_true(rc > 0.5)
})


test_that(".validate_seurat rejects non-Seurat objects", {
  expect_error(scPairs:::.validate_seurat("not a seurat"), "Seurat object")
  expect_error(scPairs:::.validate_seurat(list()),          "Seurat object")
})


test_that(".has_spatial detects spatial coordinate columns", {
  expect_false(scPairs:::.has_spatial(scpairs_testdata))

  # Add spatial coords to a copy and recheck
  sce_sp <- scpairs_testdata
  sce_sp$x_centroid <- runif(ncol(sce_sp))
  sce_sp$y_centroid <- runif(ncol(sce_sp))
  expect_true(scPairs:::.has_spatial(sce_sp))
})
