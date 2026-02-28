# scPairs Tests

This directory contains the test suite for the **scPairs** package, built on the `testthat` framework.

## Test Structure

```
tests/
├── testthat.R                           # Main test runner
└── testthat/
    ├── helper-test_data.R               # Shared helpers and custom data-creation utilities
    ├── test-find-all-pairs.R            # Tests for FindAllPairs()
    ├── test-find-gene-pairs.R           # Tests for FindGenePairs()
    ├── test-assess-pair.R               # Tests for AssessGenePair()
    ├── test-cross-celltype.R            # Tests for cross-cell-type metrics
    ├── test-consolidated-functions.R    # Tests for unified computation engine
    ├── test-internal-functions.R        # Tests for core mathematical helpers
    └── test-validation.R               # Tests for input-validation helpers
```

## Running Tests

### Run all tests
```r
devtools::test()
```

### Run a specific test file
```r
devtools::test(filter = "find-all-pairs")
devtools::test(filter = "assess-pair")
devtools::test(filter = "validation")
```

### Run with coverage report
```r
covr::package_coverage()
```

### Run as part of package check
```r
devtools::check()
```

---

## Test Categories

### 1. `FindAllPairs` Tests (`test-find-all-pairs.R`)

Validates global pair-discovery output and behavior:

- Output class is `scPairs_result` with required columns (`gene1`, `gene2`, `synergy_score`, `rank`, `confidence`)
- **Known-pair detection**: GENE3/GENE4 (injected co-expression, r ≈ 0.89) appears in the top 10 pairs
- `top_n` parameter correctly limits the number of returned pairs
- Error is raised when fewer than 2 features are provided

### 2. `FindGenePairs` Tests (`test-find-gene-pairs.R`)

Validates query-centric partner search:

- Output class is `scPairs_gene_result`; all rows have `gene1 == query_gene`
- **Ranking validation**: GENE4 is identified as the top partner for GENE3
- Error is raised for a non-existent query gene

### 3. `AssessGenePair` Tests (`test-assess-pair.R`)

Validates in-depth single-pair assessment:

- Output class is `scPairs_pair_result` with all required fields
- **Score ordering**: GENE3/GENE4 achieves a higher synergy score than a random noise pair (GENE3/GENE10)
- Permutation p-value is in (0, 1] and confidence is one of `"High"`, `"Medium"`, `"Low"`, `"NS"`
- Error is raised when `gene1 == gene2`

### 4. Cross-Cell-Type Tests (`test-cross-celltype.R`)

Validates trans-cellular interaction metrics:

- `.cross_celltype_batch()` returns numeric scores in [0, 1] (or NA)
- `.cross_celltype()` returns a complete list with score, r_ab, r_ba, n_type_pairs, and per_celltype_pair
- Single-cluster data returns NA (no cross-type pairs exist)
- A deliberately injected trans-cellular gradient produces a non-trivial score
- `FindGenePairs()` and `AssessGenePair()` both include `cross_celltype_score` in their output

### 5. Consolidated Engine Tests (`test-consolidated-functions.R`)

Validates the unified computation backend and schema:

- `.compute_bicor()`, `.compute_mi()`, `.compute_ratio_consistency()` dispatch correctly in both single-pair and batch modes, returning values consistent with the underlying primitive functions
- `.compute_cross_celltype()` dispatches correctly in both modes
- Schema constants (`ABS_METRICS`, `RAW_METRICS`, `DEFAULT_WEIGHTS`, `RESULT_COLUMNS`) are self-consistent

### 6. Internal Function Tests (`test-internal-functions.R`)

Validates core mathematical helpers:

- `.bicor()`: strongly positive for correlated pairs; returns 0 for constant inputs; bounded ≤ 1
- `.mutual_info()`: non-negative; higher for dependent than independent variables
- `.ratio_consistency()`: bounded in [0, 1]; high (> 0.5) for near-constant-ratio data
- `.validate_seurat()`: rejects non-Seurat objects with an informative error
- `.has_spatial()`: correctly detects `x_centroid`/`y_centroid` metadata columns

### 7. Input-Validation Tests (`test-validation.R`)

Validates all `.validate_*()` helper functions:

| Function | Tests |
|---|---|
| `.validate_min_cells_expressed()` | Rejects non-numeric / negative; warns for > 50 % of cells |
| `.validate_neighbourhood_params()` | Rejects non-numeric / k < 2; warns when k > 50 % of cells |
| `.validate_percentage()` | Rejects values outside [0, 1] |
| `.validate_cor_method()` | Rejects invalid or empty method vectors |
| `.validate_n_perm()` | Rejects negative; warns for < 100 permutations |
| `.validate_features()` | Errors on all-missing; warns + returns intersection on partial match |
| `.validate_binning_params()` | Warns for infeasible binning parameters |

---

## Test Data

Most tests use the **built-in dataset `scpairs_testdata`** directly:

```r
# Available as soon as the package is loaded
scpairs_testdata        # 100 cells × 20 genes, 3 clusters, PCA + UMAP
```

Two co-expression patterns are injected:

| Pair | Pattern | Pearson r |
|---|---|---|
| GENE3 & GENE4 | Global (all cells) | ≈ 0.89 |
| GENE1 & GENE2 | Cluster-1 specific | moderate |

`create_test_seurat()` in `helper-test_data.R` is retained for edge-case tests
that require a custom configuration (e.g., single-cluster, unusual gene count,
spatial coordinates).

---

## Test Helpers (`helper-test_data.R`)

| Function | Purpose |
|---|---|
| `create_test_seurat(n_cells, n_genes, n_clusters)` | Custom Seurat object with injected co-expression + PCA + UMAP |
| `create_test_spatial_seurat(n_cells, n_genes)` | Extends `create_test_seurat()` with synthetic spatial coordinates |
| `skip_if_no_seurat()` | Skips the current test if Seurat is not installed |
| `skip_if_cran()` | Skips the current test in CRAN check environments |

---

## Expected Test Behavior

| Status | Meaning |
|---|---|
| ✅ Pass | All assertions satisfied |
| ⚠️ Warning | Expected package warnings (e.g., `n_perm < 100`, `NaNs produced` in batch mode) |
| ⏭ Skip | Optional dependency absent or CRAN environment detected |
| ❌ Fail | Unexpected failure — should be investigated and fixed before submission |

---

## Test Coverage Targets

| Area | Target |
|---|---|
| Main discovery functions | 100 % |
| Visualization functions | ≥ 80 % (examples exercised via `\donttest{}`) |
| Input validation | ≥ 95 % |
| Internal mathematical helpers | ≥ 90 % |
| Cross-cell-type metrics | ≥ 85 % |
| **Overall** | **> 80 %** |

---

## Adding New Tests

1. Create `test-<feature>.R` in `tests/testthat/`
2. Use `scpairs_testdata` as the default test object
3. Add `skip_if_cran()` for computationally expensive operations
4. Follow this template:

```r
test_that("<Function> does <expected behaviour>", {
  # Use built-in data where possible
  result <- SomeFunction(scpairs_testdata, ...)

  expect_s3_class(result, "expected_class")
  expect_true("expected_column" %in% colnames(result$pairs))
  expect_true(result$some_value >= 0 && result$some_value <= 1)
})
```

---

## Continuous Integration

Tests are run on:

- **Local development**: `devtools::test()`
- **Package check**: `R CMD check --as-cran`
- **Before every CRAN submission**

---

## Contact

For test-related issues:
- File an issue: <https://github.com/zhaoqing-wang/scPairs/issues>
- Email: <zhaoqingwang@mail.sdu.edu.cn>
