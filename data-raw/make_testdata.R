## data-raw/make_testdata.R
## Run this script once to regenerate the built-in test dataset.
## Do NOT call set.seed() inside package functions; it lives only here.

set.seed(7391)

suppressPackageStartupMessages({
  library(Seurat)
})

n_cells    <- 100
n_genes    <- 20
n_clusters <- 3

# ------------------------------------------------------------------
# 1. Base count matrix (Poisson noise)
# ------------------------------------------------------------------
counts <- matrix(
  rpois(n_genes * n_cells, lambda = 3),
  nrow = n_genes, ncol = n_cells
)
rownames(counts) <- paste0("GENE", seq_len(n_genes))
colnames(counts) <- sprintf("CELL%03d", seq_len(n_cells))

clusters <- rep(seq_len(n_clusters), length.out = n_cells)

# ------------------------------------------------------------------
# 2. Inject co-expression patterns
# ------------------------------------------------------------------

# GENE3 & GENE4: strongly correlated across ALL cells (global signal)
global_signal <- abs(rnorm(n_cells, mean = 8, sd = 2))
counts["GENE3", ] <- round(global_signal + rnorm(n_cells, 0, 0.5))
counts["GENE4", ] <- round(global_signal * 1.1 + rnorm(n_cells, 0, 0.5))

# GENE1 & GENE2: moderately correlated within cluster 1 only
cl1_idx <- which(clusters == 1)
local_signal <- abs(rnorm(length(cl1_idx), mean = 6, sd = 2))
counts["GENE1", cl1_idx] <- round(local_signal + rnorm(length(cl1_idx), 0, 1))
counts["GENE2", cl1_idx] <- round(local_signal * 0.9 + rnorm(length(cl1_idx), 0, 1))

# Ensure all counts non-negative
counts[counts < 0] <- 0

# ------------------------------------------------------------------
# 3. Build Seurat object with full pipeline
# ------------------------------------------------------------------
scpairs_testdata <- Seurat::CreateSeuratObject(counts = counts,
                                               project = "scPairs_example",
                                               min.cells = 0,
                                               min.features = 0)
scpairs_testdata <- Seurat::NormalizeData(scpairs_testdata, verbose = FALSE)
scpairs_testdata <- Seurat::FindVariableFeatures(scpairs_testdata,
                                                 selection.method = "vst",
                                                 nfeatures        = n_genes,
                                                 verbose          = FALSE)
scpairs_testdata <- Seurat::ScaleData(scpairs_testdata,
                                      features = Seurat::VariableFeatures(scpairs_testdata),
                                      verbose  = FALSE)

# PCA on all variable genes (5 PCs sufficient for 20 genes)
npcs <- min(5, n_genes - 1)
scpairs_testdata <- Seurat::RunPCA(scpairs_testdata,
                                   features = Seurat::VariableFeatures(scpairs_testdata),
                                   npcs     = npcs,
                                   verbose  = FALSE)

# UMAP from PCA
scpairs_testdata <- Seurat::RunUMAP(scpairs_testdata,
                                    dims    = seq_len(npcs),
                                    verbose = FALSE)

# Cluster metadata
scpairs_testdata$seurat_clusters <- factor(clusters)
Seurat::Idents(scpairs_testdata) <- "seurat_clusters"

# ------------------------------------------------------------------
# 4. Verify before saving
# ------------------------------------------------------------------
stopifnot(
  ncol(scpairs_testdata) == n_cells,
  nrow(scpairs_testdata) == n_genes,
  "umap"  %in% names(scpairs_testdata@reductions),
  "pca"   %in% names(scpairs_testdata@reductions),
  "seurat_clusters" %in% colnames(scpairs_testdata@meta.data)
)

cat("scpairs_testdata: ",
    nrow(scpairs_testdata), "genes x",
    ncol(scpairs_testdata), "cells\n")
cat("Reductions:", paste(names(scpairs_testdata@reductions), collapse = ", "), "\n")
cat("Clusters:", levels(scpairs_testdata$seurat_clusters), "\n")

norm_mat <- Seurat::GetAssayData(scpairs_testdata, layer = "data")
r34 <- cor(as.numeric(norm_mat["GENE3", ]),
           as.numeric(norm_mat["GENE4", ]))
cat(sprintf("GENE3-GENE4 Pearson r = %.3f (expected > 0.8)\n", r34))

# ------------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------------
usethis::use_data(scpairs_testdata, overwrite = TRUE, compress = "xz")
