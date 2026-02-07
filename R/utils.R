#' Validate that input is a Seurat object
#'
#' @param object Object to validate.
#' @param require_spatial Logical; if TRUE, require spatial coordinates.
#' @return TRUE invisibly; stops with informative error otherwise.
#' @keywords internal
.validate_seurat <- function(object, require_spatial = FALSE) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.", call. = FALSE)
  }

  assay <- Seurat::DefaultAssay(object)
  if (is.null(assay)) {
    stop("Seurat object has no default assay set.", call. = FALSE)
  }

  if (require_spatial && !.has_spatial(object)) {
    stop(
      "Spatial coordinates are required but not found in the Seurat object.\n",
      "  Ensure the object contains an image / spatial coordinates slot.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Detect whether a Seurat object contains spatial information
#'
#' @param object A Seurat object.
#' @return Logical scalar.
#' @keywords internal
.has_spatial <- function(object) {
  imgs <- tryCatch(SeuratObject::Images(object), error = function(e) character(0))
  if (length(imgs) > 0) return(TRUE)


  # Fallback: check for common spatial coordinate columns in meta.data
  spatial_cols <- c("x", "y", "imagerow", "imagecol",
                    "x_centroid", "y_centroid", "array_row", "array_col")
  any(tolower(spatial_cols) %in% tolower(colnames(object@meta.data)))
}


#' Extract normalised expression matrix (genes x cells, dense or sparse)
#'
#' For speed-critical operations the matrix is kept sparse (dgCMatrix) when
#' possible.  Only the `data` slot (log-normalised) is used by default.
#'
#' @param object Seurat object.
#' @param features Character vector of gene names.  NULL = all genes.
#' @param assay Character; assay name.  NULL = default assay.
#' @param slot Character; slot to pull ("data", "counts", "scale.data").
#' @return A (possibly sparse) matrix with genes in rows, cells in columns.
#' @keywords internal
.get_expression_matrix <- function(object,
                                   features = NULL,
                                   assay    = NULL,
                                   slot     = "data") {
  assay <- assay %||% Seurat::DefaultAssay(object)
  # Seurat v5 uses 'layer' instead of 'slot'
  mat <- tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = slot),
    error = function(e) {
      Seurat::GetAssayData(object, assay = assay, slot = slot)
    }
  )

  if (!is.null(features)) {
    features <- intersect(features, rownames(mat))
    if (length(features) == 0) {
      stop("None of the requested features found in the expression matrix.",
           call. = FALSE)
    }
    mat <- mat[features, , drop = FALSE]
  }

  mat
}


#' Extract spatial coordinates from a Seurat object
#'
#' Attempts `GetTissueCoordinates()` first; falls back to common meta.data
#' column patterns.
#'
#' @param object Seurat object.
#' @return A data.frame with columns `x` and `y`, rownames = cell barcodes.
#' @keywords internal
.get_spatial_coords <- function(object) {
  # Try GetTissueCoordinates (Seurat v4/v5)
  coords <- tryCatch({
    imgs <- SeuratObject::Images(object)
    if (length(imgs) == 0) stop("no images")
    cc <- Seurat::GetTissueCoordinates(object, image = imgs[1])
    if (ncol(cc) >= 2) cc else stop("bad dims")
  }, error = function(e) NULL)

  if (!is.null(coords)) {
    colnames(coords)[1:2] <- c("x", "y")
    return(coords)
  }

  # Fallback: probe meta.data for coordinate columns
  md <- object@meta.data
  x_candidates <- c("x", "x_centroid", "imagerow", "array_row")
  y_candidates <- c("y", "y_centroid", "imagecol", "array_col")

  x_col <- x_candidates[tolower(x_candidates) %in% tolower(colnames(md))][1]
  y_col <- y_candidates[tolower(y_candidates) %in% tolower(colnames(md))][1]

  if (is.na(x_col) || is.na(y_col)) {
    stop("Cannot locate spatial coordinates in the Seurat object.", call. = FALSE)
  }

  data.frame(
    x = md[[x_col]],
    y = md[[y_col]],
    row.names = rownames(md)
  )
}


#' Select highly-variable genes or top-expressed genes for analysis
#'
#' When the number of genes is very large we pre-filter to a tractable set.
#' Priority order: user-supplied genes > Seurat VariableFeatures > top genes
#' by mean expression.
#'
#' @param object Seurat object.
#' @param features Character vector of gene names; NULL for auto-selection.
#' @param n_top Integer; maximum number of genes to select.
#' @param assay Character; assay name.
#' @return Character vector of gene names.
#' @keywords internal
.select_features <- function(object,
                             features = NULL,
                             n_top    = 2000,
                             assay    = NULL) {
  assay <- assay %||% Seurat::DefaultAssay(object)
  all_genes <- rownames(tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = "data"),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = "data")
  ))

  if (!is.null(features)) {
    found <- intersect(features, all_genes)
    if (length(found) == 0) {
      stop("None of the supplied features found in the assay.", call. = FALSE)
    }
    if (length(found) < length(features)) {
      missing <- setdiff(features, all_genes)
      warning(sprintf("%d features not found and removed: %s",
                      length(missing),
                      paste(utils::head(missing, 5), collapse = ", ")),
              call. = FALSE)
    }
    return(found)
  }

  # Use Seurat variable features if available
  vf <- tryCatch(Seurat::VariableFeatures(object, assay = assay),
                 error = function(e) character(0))
  if (length(vf) > 0) {
    return(utils::head(vf, n_top))
  }

  # Fallback: top genes by mean expression
  mat <- tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = "data"),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = "data")
  )
  gene_means <- Matrix::rowMeans(mat)
  names(sort(gene_means, decreasing = TRUE))[seq_len(min(n_top, length(gene_means)))]
}


#' Null-coalescing operator (if not already imported from Seurat)
#'
#' @param a First value.
#' @param b Second value (returned if \code{a} is NULL).
#' @return \code{a} if not NULL, otherwise \code{b}.
#' @keywords internal
#' @name grapes-or-or-grapes
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Fast row-wise variance for sparse matrices
#'
#' @param x Sparse or dense matrix (genes x cells).
#' @return Named numeric vector of variances.
#' @keywords internal
.row_vars <- function(x) {
  rm <- Matrix::rowMeans(x)
  n <- ncol(x)
  # E[X^2] - (E[X])^2, with Bessel correction
  if (inherits(x, "dgCMatrix") || inherits(x, "dgRMatrix")) {
    rsq <- Matrix::rowSums(x^2)
  } else {
    rsq <- rowSums(x^2)
  }
  vars <- (rsq - n * rm^2) / (n - 1)
  names(vars) <- rownames(x)
  vars
}


#' Print a progress message (respects verbose flag)
#'
#' @param ... Parts of the message (passed to `paste0`).
#' @param verbose Logical.
#' @keywords internal
.msg <- function(..., verbose = TRUE) {
  if (verbose) message(paste0("[scPairs] ", ...))
}
