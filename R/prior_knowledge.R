#' Prior Knowledge Integration for Synergistic Gene Pair Discovery
#'
#' @description
#' Integrates biological prior knowledge databases to shift the scoring
#' from pure co-expression towards functional synergy.  Three complementary
#' evidence layers are provided:
#'
#' 1. **Prior interaction score** -- direct gene-gene functional annotation
#'    overlap from GO, KEGG, and curated interaction databases.
#' 2. **Pathway bridge score** -- indirect synergy through shared intermediary
#'    genes that connect the pair in biological pathways AND are expressed in
#'    the current dataset.
#' 3. **Neighbourhood synergy score** -- directional enrichment of gene B in
#'    the biological neighbourhood of gene A-expressing cells, capturing
#'    paracrine/juxtacrine interactions even when the two genes are never
#'    co-detected in the same cell.
#'
#' These scores are designed to up-weight gene pairs with biological plausibility
#' for functional synergy, even when cell-level co-expression is weak or absent.
#'
#' @name prior-knowledge
#' @keywords internal
NULL


# ===========================================================================
#  1. Prior knowledge database construction
# ===========================================================================

#' Build a prior knowledge gene interaction network
#'
#' Constructs a gene-gene interaction list from available annotation sources.
#' The function queries GO (Biological Process), KEGG pathways, and optionally
#' user-supplied interaction databases.  It returns a list structure that can
#' be used for scoring gene pairs.
#'
#' @param organism Character; organism identifier for annotation lookup.
#'     Supported: "mouse" (Mus musculus) or "human" (Homo sapiens).
#' @param genes Character vector of gene symbols to include (typically the
#'     features in the Seurat object).  Restricts the network to relevant genes.
#' @param sources Character vector; knowledge sources to use.  Any subset of
#'     \code{c("GO", "KEGG", "custom")}.  Default: \code{c("GO", "KEGG")}.
#' @param custom_pairs Optional data.frame with columns \code{gene1}, \code{gene2}
#'     (and optionally \code{source}, \code{weight}) for user-supplied interactions.
#'     This can include CellChatDB, CellPhoneDB, SCENIC regulon targets, etc.
#' @param min_genes Integer; minimum number of genes in a GO/KEGG term to be
#'     included (avoids overly broad terms).  Default 5.
#' @param max_genes Integer; maximum genes per term (avoids overly broad terms
#'     like "protein binding").  Default 500.
#' @param verbose Logical.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{gene_sets}}{Named list of character vectors; each element is a
#'     functional term/pathway mapping to its member genes.}
#'   \item{\code{gene_to_terms}}{Named list; for each gene, the set of terms
#'     it belongs to.}
#'   \item{\code{interactions}}{data.table with columns \code{gene1}, \code{gene2},
#'     \code{source}, \code{n_shared_terms}, \code{jaccard}.}
#'   \item{\code{organism}}{Character.}
#'   \item{\code{n_genes}}{Integer; number of genes covered.}
#'   \item{\code{n_terms}}{Integer; number of functional terms.}
#' }
#'
#' @keywords internal
.build_prior_network <- function(organism    = "mouse",
                                 genes       = NULL,
                                 sources     = c("GO", "KEGG"),
                                 custom_pairs = NULL,
                                 min_genes   = 5,
                                 max_genes   = 500,
                                 verbose     = TRUE) {

  gene_sets <- list()

  # --- GO Biological Process ---
  if ("GO" %in% sources) {
    .msg("Building prior network: GO Biological Process ...", verbose = verbose)
    go_sets <- .get_go_gene_sets(organism, genes, min_genes, max_genes)
    gene_sets <- c(gene_sets, go_sets)
  }

  # --- KEGG pathways ---
  if ("KEGG" %in% sources) {
    .msg("Building prior network: KEGG pathways ...", verbose = verbose)
    kegg_sets <- .get_kegg_gene_sets(organism, genes, min_genes, max_genes)
    gene_sets <- c(gene_sets, kegg_sets)
  }

  if (length(gene_sets) == 0 && is.null(custom_pairs)) {
    .msg("No annotation sources available. Using expression-only mode.",
         verbose = verbose)
    return(list(
      gene_sets     = list(),
      gene_to_terms = list(),
      interactions  = data.table::data.table(
        gene1 = character(), gene2 = character(),
        source = character(), n_shared_terms = integer(),
        jaccard = numeric()),
      organism = organism,
      n_genes  = 0L,
      n_terms  = 0L
    ))
  }

  # Filter gene sets by size
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= min_genes &
                         vapply(gene_sets, length, integer(1)) <= max_genes]

  # Restrict to genes present in data

  if (!is.null(genes)) {
    gene_sets <- lapply(gene_sets, function(gs) intersect(gs, genes))
    gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= 2]
  }

  .msg("  ", length(gene_sets), " functional terms retained.", verbose = verbose)

  # Build gene-to-terms mapping
  gene_to_terms <- list()
  for (term_name in names(gene_sets)) {
    for (g in gene_sets[[term_name]]) {
      gene_to_terms[[g]] <- c(gene_to_terms[[g]], term_name)
    }
  }

  # Add custom interactions
  custom_dt <- data.table::data.table(
    gene1 = character(), gene2 = character(),
    source = character(), weight = numeric())

  if (!is.null(custom_pairs) && nrow(custom_pairs) > 0) {
    .msg("Adding ", nrow(custom_pairs), " custom interactions ...", verbose = verbose)
    cp <- data.table::as.data.table(custom_pairs)
    if (!"source" %in% colnames(cp)) cp[, source := "custom"]
    if (!"weight" %in% colnames(cp)) cp[, weight := 1.0]
    # Restrict to genes in data
    if (!is.null(genes)) {
      cp <- cp[gene1 %in% genes & gene2 %in% genes]
    }
    custom_dt <- cp[, .(gene1, gene2, source, weight)]
  }

  list(
    gene_sets     = gene_sets,
    gene_to_terms = gene_to_terms,
    custom_pairs  = custom_dt,
    organism      = organism,
    n_genes       = length(gene_to_terms),
    n_terms       = length(gene_sets)
  )
}


# ===========================================================================
#  2. GO gene sets (using org.Mm.eg.db / org.Hs.eg.db if available)
# ===========================================================================

#' Retrieve GO Biological Process gene sets
#' @keywords internal
.get_go_gene_sets <- function(organism, genes = NULL,
                              min_genes = 5, max_genes = 500) {
  org_pkg <- switch(organism,
    mouse = "org.Mm.eg.db",
    human = "org.Hs.eg.db",
    stop("Unsupported organism: ", organism, ". Use 'mouse' or 'human'.",
         call. = FALSE)
  )

  if (!requireNamespace(org_pkg, quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    message("[scPairs] ", org_pkg, " or AnnotationDbi not available. ",
            "Skipping GO annotations.\n",
            "  Install with: BiocManager::install(c('", org_pkg, "', 'AnnotationDbi'))")
    return(list())
  }

  db <- getExportedValue(org_pkg, org_pkg)

  # Map gene symbols to GO BP terms
  # Use tryCatch to handle potential database issues gracefully
  go_map <- tryCatch({
    AnnotationDbi::select(db,
      keys = AnnotationDbi::keys(db, keytype = "SYMBOL"),
      columns = c("SYMBOL", "GO", "ONTOLOGY"),
      keytype = "SYMBOL"
    )
  }, error = function(e) {
    message("[scPairs] Could not query GO annotations: ", conditionMessage(e))
    return(data.frame())
  })

  if (nrow(go_map) == 0) return(list())

  # Filter to Biological Process only
  go_map <- go_map[!is.na(go_map$ONTOLOGY) & go_map$ONTOLOGY == "BP", ]
  go_map <- go_map[!is.na(go_map$GO) & !is.na(go_map$SYMBOL), ]

  # Filter to relevant genes
  if (!is.null(genes)) {
    go_map <- go_map[go_map$SYMBOL %in% genes, ]
  }

  # Build gene sets
  go_terms <- split(go_map$SYMBOL, go_map$GO)
  go_terms <- lapply(go_terms, unique)

  # Filter by size
  sizes <- vapply(go_terms, length, integer(1))
  go_terms <- go_terms[sizes >= min_genes & sizes <= max_genes]

  # Prefix names for clarity
  names(go_terms) <- paste0("GO:", names(go_terms))

  go_terms
}


#' Retrieve KEGG pathway gene sets
#' @keywords internal
.get_kegg_gene_sets <- function(organism, genes = NULL,
                                min_genes = 5, max_genes = 500) {
  org_code <- switch(organism,
    mouse = "mmu",
    human = "hsa",
    stop("Unsupported organism: ", organism, call. = FALSE)
  )

  org_pkg <- switch(organism,
    mouse = "org.Mm.eg.db",
    human = "org.Hs.eg.db"
  )

  if (!requireNamespace(org_pkg, quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    message("[scPairs] ", org_pkg, " or AnnotationDbi not available. ",
            "Skipping KEGG annotations.")
    return(list())
  }

  db <- getExportedValue(org_pkg, org_pkg)

  # Map symbols to KEGG pathways via Entrez IDs
  kegg_map <- tryCatch({
    AnnotationDbi::select(db,
      keys = AnnotationDbi::keys(db, keytype = "SYMBOL"),
      columns = c("SYMBOL", "PATH"),
      keytype = "SYMBOL"
    )
  }, error = function(e) {
    message("[scPairs] Could not query KEGG annotations: ", conditionMessage(e))
    return(data.frame())
  })

  if (nrow(kegg_map) == 0) return(list())

  kegg_map <- kegg_map[!is.na(kegg_map$PATH) & !is.na(kegg_map$SYMBOL), ]

  if (!is.null(genes)) {
    kegg_map <- kegg_map[kegg_map$SYMBOL %in% genes, ]
  }

  kegg_sets <- split(kegg_map$SYMBOL, kegg_map$PATH)
  kegg_sets <- lapply(kegg_sets, unique)

  sizes <- vapply(kegg_sets, length, integer(1))
  kegg_sets <- kegg_sets[sizes >= min_genes & sizes <= max_genes]

  names(kegg_sets) <- paste0("KEGG:", names(kegg_sets))

  kegg_sets
}


# ===========================================================================
#  3. Prior interaction score (Jaccard similarity of functional annotations)
# ===========================================================================

#' Compute prior interaction score for gene pairs
#'
#' For each gene pair, computes the Jaccard similarity of their functional
#' annotation sets (GO terms + KEGG pathways).  A high Jaccard indicates that
#' the two genes participate in many of the same biological processes --
#' evidence for functional relatedness beyond mere co-expression.
#'
#' @param pair_dt data.table with gene1, gene2 columns.
#' @param prior_net Prior network from \code{.build_prior_network()}.
#' @return Numeric vector of prior scores in \[0, 1\].
#' @keywords internal
.prior_score_batch <- function(pair_dt, prior_net) {
  g2t <- prior_net$gene_to_terms
  custom_dt <- prior_net$custom_pairs
  n_pairs <- nrow(pair_dt)
  scores <- numeric(n_pairs)

  for (i in seq_len(n_pairs)) {
    g1 <- pair_dt$gene1[i]
    g2 <- pair_dt$gene2[i]

    # Annotation-based Jaccard
    terms1 <- g2t[[g1]]
    terms2 <- g2t[[g2]]

    jaccard <- 0
    if (length(terms1) > 0 && length(terms2) > 0) {
      inter <- length(intersect(terms1, terms2))
      union <- length(union(terms1, terms2))
      if (union > 0) jaccard <- inter / union
    }

    # Custom interaction bonus
    custom_bonus <- 0
    if (nrow(custom_dt) > 0) {
      match <- custom_dt[(gene1 == g1 & gene2 == g2) |
                         (gene1 == g2 & gene2 == g1)]
      if (nrow(match) > 0) {
        custom_bonus <- max(match$weight)
      }
    }

    # Combine: Jaccard + custom bonus, capped at 1
    scores[i] <- min(jaccard + custom_bonus * 0.5, 1)
  }

  scores
}


#' Compute prior interaction score for a single gene pair
#' @keywords internal
.prior_score <- function(gene1, gene2, prior_net) {
  g2t <- prior_net$gene_to_terms
  terms1 <- g2t[[gene1]]
  terms2 <- g2t[[gene2]]

  jaccard <- 0
  if (length(terms1) > 0 && length(terms2) > 0) {
    inter <- length(intersect(terms1, terms2))
    union <- length(union(terms1, terms2))
    if (union > 0) jaccard <- inter / union
  }

  custom_bonus <- 0
  if (nrow(prior_net$custom_pairs) > 0) {
    cp <- prior_net$custom_pairs
    match <- cp[(gene1 == !!gene1 & gene2 == !!gene2) |
                (gene1 == !!gene2 & gene2 == !!gene1)]
    if (nrow(match) > 0) custom_bonus <- max(match$weight)
  }

  min(jaccard + custom_bonus * 0.5, 1)
}


# ===========================================================================
#  4. Pathway bridge score
# ===========================================================================

#' Compute pathway bridge score for gene pairs
#'
#' For a gene pair (A, B), identifies intermediate "bridge" genes C such that
#' C shares functional annotations with BOTH A and B, AND C is expressed in
#' the current dataset.  The bridge score reflects the strength of indirect
#' connectivity:
#'
#' \deqn{bridge\_score = \frac{n\_bridges}{\sqrt{|terms_A| \cdot |terms_B|}}}
#'
#' This captures synergistic relationships where two genes are not directly
#' co-annotated but are connected through shared pathway intermediaries.
#'
#' @param pair_dt data.table with gene1, gene2 columns.
#' @param prior_net Prior network from \code{.build_prior_network()}.
#' @param expressed_genes Character vector of genes expressed in the dataset
#'     (bridge genes must be expressed to be relevant).
#' @param min_shared Integer; minimum shared terms between a bridge gene
#'     and each member of the pair.  Default 1.
#' @return A list with:
#' \describe{
#'   \item{\code{scores}}{Numeric vector of bridge scores.}
#'   \item{\code{bridges}}{List of character vectors; bridge genes for each pair.}
#' }
#' @keywords internal
.bridge_score_batch <- function(pair_dt, prior_net, expressed_genes,
                                min_shared = 1) {
  g2t <- prior_net$gene_to_terms
  gene_sets <- prior_net$gene_sets
  n_pairs <- nrow(pair_dt)
  scores <- numeric(n_pairs)
  bridge_list <- vector("list", n_pairs)

  for (i in seq_len(n_pairs)) {
    g1 <- pair_dt$gene1[i]
    g2 <- pair_dt$gene2[i]

    terms1 <- g2t[[g1]]
    terms2 <- g2t[[g2]]

    if (length(terms1) == 0 || length(terms2) == 0) {
      scores[i] <- 0
      bridge_list[[i]] <- character(0)
      next
    }

    # Find all genes in terms shared with g1 and in terms shared with g2
    genes_in_g1_terms <- unique(unlist(gene_sets[terms1], use.names = FALSE))
    genes_in_g2_terms <- unique(unlist(gene_sets[terms2], use.names = FALSE))

    # Bridge candidates: genes in both gene sets, expressed, and not g1/g2
    bridges <- intersect(genes_in_g1_terms, genes_in_g2_terms)
    bridges <- setdiff(bridges, c(g1, g2))
    bridges <- intersect(bridges, expressed_genes)

    # Optionally filter by minimum shared terms with each gene
    if (min_shared > 1 && length(bridges) > 0) {
      bridges <- bridges[vapply(bridges, function(bg) {
        bg_terms <- g2t[[bg]]
        length(intersect(bg_terms, terms1)) >= min_shared &&
          length(intersect(bg_terms, terms2)) >= min_shared
      }, logical(1))]
    }

    n_bridges <- length(bridges)
    denominator <- sqrt(length(terms1) * length(terms2))
    scores[i] <- if (denominator > 0) min(n_bridges / denominator, 1) else 0
    bridge_list[[i]] <- bridges
  }

  list(scores = scores, bridges = bridge_list)
}


#' Compute bridge score and bridge genes for a single pair
#' @keywords internal
.bridge_score <- function(gene1, gene2, prior_net, expressed_genes,
                          min_shared = 1, top_n = 20) {
  g2t <- prior_net$gene_to_terms
  gene_sets <- prior_net$gene_sets

  terms1 <- g2t[[gene1]]
  terms2 <- g2t[[gene2]]

  if (length(terms1) == 0 || length(terms2) == 0) {
    return(list(score = 0, bridges = character(0),
                shared_terms = character(0)))
  }

  # Shared terms (direct overlap)
  shared_terms <- intersect(terms1, terms2)

  # Find bridge genes
  genes_in_g1_terms <- unique(unlist(gene_sets[terms1], use.names = FALSE))
  genes_in_g2_terms <- unique(unlist(gene_sets[terms2], use.names = FALSE))

  bridges <- intersect(genes_in_g1_terms, genes_in_g2_terms)
  bridges <- setdiff(bridges, c(gene1, gene2))
  bridges <- intersect(bridges, expressed_genes)

  # Score each bridge by its connectivity strength
  if (length(bridges) > 0) {
    bridge_strength <- vapply(bridges, function(bg) {
      bg_terms <- g2t[[bg]]
      s1 <- length(intersect(bg_terms, terms1))
      s2 <- length(intersect(bg_terms, terms2))
      sqrt(s1 * s2)
    }, numeric(1))

    # Sort by strength
    ord <- order(bridge_strength, decreasing = TRUE)
    bridges <- bridges[ord]
    bridge_strength <- bridge_strength[ord]

    if (length(bridges) > top_n) {
      bridges <- bridges[seq_len(top_n)]
    }
  }

  denominator <- sqrt(length(terms1) * length(terms2))
  score <- if (denominator > 0) min(length(bridges) / denominator, 1) else 0

  list(score = score, bridges = bridges, shared_terms = shared_terms)
}


# ===========================================================================
#  5. Neighbourhood synergy score (directional paracrine enrichment)
# ===========================================================================

#' Compute neighbourhood synergy score for gene pairs
#'
#' Unlike the existing neighbourhood co-expression score (which measures whether
#' A- and B-expressing cells tend to share neighbours), this metric explicitly
#' quantifies **directional** neighbourhood enrichment:
#'
#' For cells highly expressing gene A (top quartile), is gene B's expression
#' in their neighbourhood significantly higher than expected?  The score is:
#'
#' \deqn{NS = \frac{\text{mean neigh expr B for top-A cells}}
#'                  {\text{mean neigh expr B for all cells}}}
#'
#' This captures paracrine-like interactions where one gene's product in one
#' cell influences the expression of another gene in nearby cells.
#'
#' @param mat Expression matrix (genes x cells).
#' @param pair_dt data.table with gene1, gene2.
#' @param W KNN weight matrix.
#' @return Numeric vector of neighbourhood synergy scores.
#' @keywords internal
.neighbourhood_synergy_batch <- function(mat, pair_dt, W) {
  common <- intersect(colnames(mat), rownames(W))
  mat <- mat[, common, drop = FALSE]
  W <- W[common, common, drop = FALSE]

  if (inherits(mat, "dgCMatrix") || inherits(mat, "dgRMatrix")) {
    mat_dense <- as.matrix(mat)
  } else {
    mat_dense <- mat
  }

  n_pairs <- nrow(pair_dt)
  ns_vals <- numeric(n_pairs)

  # Pre-compute neighbourhood expression for all genes:
  # neigh_expr[g, i] = weighted average of gene g in cell i's neighbourhood
  unique_genes <- unique(c(pair_dt$gene1, pair_dt$gene2))
  unique_genes <- intersect(unique_genes, rownames(mat_dense))
  sub_mat <- mat_dense[unique_genes, , drop = FALSE]

  # Neighbourhood average: mat %*% t(W) (KNN-smoothed without self)
  neigh_expr <- as.matrix(sub_mat %*% Matrix::t(W))

  gene_idx <- stats::setNames(seq_along(unique_genes), unique_genes)

  for (i in seq_len(n_pairs)) {
    g1 <- pair_dt$gene1[i]
    g2 <- pair_dt$gene2[i]

    if (!(g1 %in% unique_genes) || !(g2 %in% unique_genes)) {
      ns_vals[i] <- NA_real_
      next
    }

    idx1 <- gene_idx[g1]
    idx2 <- gene_idx[g2]

    expr_a <- sub_mat[idx1, ]
    neigh_b <- neigh_expr[idx2, ]

    # Cells expressing gene A (top quartile of non-zero)
    nonzero_a <- expr_a[expr_a > 0]
    if (length(nonzero_a) < 5) {
      ns_vals[i] <- NA_real_
      next
    }

    thresh_a <- stats::quantile(nonzero_a, 0.5)
    high_a_cells <- which(expr_a >= thresh_a)

    # Mean neighbourhood expr B for high-A cells vs all cells
    mean_neigh_b_highA <- mean(neigh_b[high_a_cells])
    mean_neigh_b_all   <- mean(neigh_b)

    if (mean_neigh_b_all > 0) {
      ratio_ab <- mean_neigh_b_highA / mean_neigh_b_all
    } else {
      ratio_ab <- 1
    }

    # Symmetric: also check B -> A direction
    expr_b <- sub_mat[idx2, ]
    neigh_a <- neigh_expr[idx1, ]

    nonzero_b <- expr_b[expr_b > 0]
    if (length(nonzero_b) < 5) {
      ns_vals[i] <- max(ratio_ab - 1, 0)  # Use unidirectional score
      next
    }

    thresh_b <- stats::quantile(nonzero_b, 0.5)
    high_b_cells <- which(expr_b >= thresh_b)

    mean_neigh_a_highB <- mean(neigh_a[high_b_cells])
    mean_neigh_a_all   <- mean(neigh_a)

    if (mean_neigh_a_all > 0) {
      ratio_ba <- mean_neigh_a_highB / mean_neigh_a_all
    } else {
      ratio_ba <- 1
    }

    # Geometric mean of directional enrichments (minus 1 baseline)
    if (ratio_ab > 0 && ratio_ba > 0) {
      ns_vals[i] <- sqrt(pmax(ratio_ab, 1) * pmax(ratio_ba, 1)) - 1
    } else {
      ns_vals[i] <- 0
    }
  }

  ns_vals
}


#' Single-pair neighbourhood synergy score
#' @keywords internal
.neighbourhood_synergy <- function(x, y, W) {
  n <- length(x)

  # Neighbourhood expression
  neigh_y <- as.numeric(W %*% y)
  neigh_x <- as.numeric(W %*% x)

  # A -> B direction
  nonzero_x <- x[x > 0]
  if (length(nonzero_x) < 5) return(NA_real_)

  thresh_x <- stats::quantile(nonzero_x, 0.5)
  high_x <- which(x >= thresh_x)

  mean_ny_highx <- mean(neigh_y[high_x])
  mean_ny_all   <- mean(neigh_y)
  ratio_ab <- if (mean_ny_all > 0) mean_ny_highx / mean_ny_all else 1

  # B -> A direction
  nonzero_y <- y[y > 0]
  if (length(nonzero_y) < 5) return(max(ratio_ab - 1, 0))

  thresh_y <- stats::quantile(nonzero_y, 0.5)
  high_y <- which(y >= thresh_y)

  mean_nx_highy <- mean(neigh_x[high_y])
  mean_nx_all   <- mean(neigh_x)
  ratio_ba <- if (mean_nx_all > 0) mean_nx_highy / mean_nx_all else 1

  if (ratio_ab > 0 && ratio_ba > 0) {
    sqrt(pmax(ratio_ab, 1) * pmax(ratio_ba, 1)) - 1
  } else {
    0
  }
}
