#' Load Curated Cell Death Gene Sets
#'
#' Load and prepare cell death pathway gene sets for ORA and GSEA analysis.
#' This function returns a gene list compatible with clusterProfiler functions.
#' The 13 recommended representative gene sets are automatically merged in
#' (see \code{representative_pathways} for their curation metadata).
#'
#' @param .test_df Internal parameter for unit testing. A data.frame with
#'   \code{term} and \code{gene} columns; when provided, the built-in dataset
#'   and the representative gene sets are NOT used.
#' @param .test_rep Internal parameter for unit testing. A named list of gene
#'   vectors used in place of the built-in representative gene sets; when
#'   provided together with \code{.test_df}, the merge step is still exercised
#' @return A named list of gene vectors, each element representing a cell death pathway.
#' @export
#' @examples
#' gene_sets <- load_cell_death_genes()
#' length(gene_sets)
#' names(gene_sets)
load_cell_death_genes <- function(.test_df = NULL, .test_rep = NULL) {

  # ===================== 1. Load built-in rda dataset of the package =====================
  if (!is.null(.test_df)) {
    df <- .test_df
  } else {
    utils::data("cd_Genesets", envir = environment())
    df <- cd_Genesets
  }

  # ===================== 2. Mandatory check of column names =====================
  if (!all(c("term", "gene") %in% colnames(df))) {
    stop("Error: Gene set data must contain two columns: term and gene!")
  }

  # ===================== 3. Clean gene symbols =====================
  df$gene <- trimws(as.character(df$gene))
  df$term <- trimws(as.character(df$term))

  # ===================== 4. Filter out NA values and empty strings =====================
  df <- df[!is.na(df$gene) & df$gene != "", ]
  df <- df[!is.na(df$term) & df$term != "", ]

  # ===================== 5. Convert to list format =====================
  gene_list <- split(df$gene, df$term)

  # ===================== 6. Filter out pathways with insufficient genes =====================
  gene_list <- gene_list[sapply(gene_list, length) >= 3]

  # ===================== 7. Merge recommended representative gene sets =====================
  # Skipped only in pure test mode (.test_df given without .test_rep).
  # Existing pathways with the same name are never overwritten.
  if (is.null(.test_df) || !is.null(.test_rep)) {
    rep_sets <- if (!is.null(.test_rep)) .test_rep else representative_genes
    missing_sets <- rep_sets[!names(rep_sets) %in% names(gene_list)]
    missing_sets <- missing_sets[sapply(missing_sets, length) >= 3]
    gene_list <- c(gene_list, missing_sets)
  }

  return(gene_list)
}
