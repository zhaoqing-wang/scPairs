# Tests for R/validation.R helper functions

test_that(".validate_min_cells_expressed rejects invalid input", {
  expect_error(
    scPairs:::.validate_min_cells_expressed("abc", 100),
    "single numeric"
  )
  expect_error(
    scPairs:::.validate_min_cells_expressed(-5, 100),
    "non-negative"
  )
})

test_that(".validate_min_cells_expressed warns for large values", {
  expect_warning(
    scPairs:::.validate_min_cells_expressed(50, 100),
    "very high"
  )
  # Should pass silently for reasonable values
  expect_silent(
    scPairs:::.validate_min_cells_expressed(5, 100)
  )
})


test_that(".validate_neighbourhood_params rejects bad k", {
  expect_error(
    scPairs:::.validate_neighbourhood_params("not_a_number", 100),
    "single numeric"
  )
  expect_error(
    scPairs:::.validate_neighbourhood_params(1, 100),
    ">= 2"
  )
})

test_that(".validate_neighbourhood_params warns for large k", {
  expect_warning(
    scPairs:::.validate_neighbourhood_params(60, 100),
    "more than half"
  )
  expect_silent(
    scPairs:::.validate_neighbourhood_params(20, 100)
  )
})


test_that(".validate_percentage rejects out-of-range values", {
  expect_error(
    scPairs:::.validate_percentage(-0.1, "alpha"),
    "between 0 and 1"
  )
  expect_error(
    scPairs:::.validate_percentage(1.5, "alpha"),
    "between 0 and 1"
  )
  expect_silent(
    scPairs:::.validate_percentage(0.5, "alpha")
  )
})


test_that(".validate_cor_method rejects invalid methods", {
  expect_error(
    scPairs:::.validate_cor_method(c("pearson", "invalid")),
    "Invalid"
  )
  expect_error(
    scPairs:::.validate_cor_method(character(0)),
    "specified"
  )
  expect_silent(
    scPairs:::.validate_cor_method(c("pearson", "spearman"))
  )
})


test_that(".validate_n_perm rejects negative values", {
  expect_error(
    scPairs:::.validate_n_perm(-1),
    "non-negative"
  )
})

test_that(".validate_n_perm warns for small permutation counts", {
  expect_warning(
    scPairs:::.validate_n_perm(50),
    "< 100"
  )
  expect_silent(
    scPairs:::.validate_n_perm(999)
  )
})


test_that(".validate_features detects missing genes", {
  skip_if_no_seurat()

  sce <- create_test_seurat(n_cells = 50, n_genes = 10)

  # All missing -> error

  expect_error(
    scPairs:::.validate_features(c("NONEXISTENT1", "NONEXISTENT2"), sce),
    "None of the requested features"
  )

  # Some missing -> warning, returns intersection
  expect_warning(
    result <- scPairs:::.validate_features(c("GENE1", "NONEXISTENT"), sce),
    "not found"
  )
  expect_equal(result, "GENE1")

  # All present -> silent
  expect_silent(
    result2 <- scPairs:::.validate_features(c("GENE1", "GENE2"), sce)
  )
  expect_equal(sort(result2), c("GENE1", "GENE2"))
})


test_that(".validate_binning_params detects infeasible parameters", {
  # Warn when parameters are too strict
  expect_warning(
    scPairs:::.validate_binning_params(
      n_bins = 100, min_cells_per_bin = 5,
      min_bins = 8, n_cells = 50, n_types = 3
    ),
    "too strict"
  )

  # Pass for reasonable params
  expect_silent(
    scPairs:::.validate_binning_params(
      n_bins = 10, min_cells_per_bin = 3,
      min_bins = 5, n_cells = 500, n_types = 3
    )
  )
})
