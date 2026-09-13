#' Cell Death Pathway GSEA Analysis
#'
#' This function performs GSEA using user-defined cell death gene sets.
#' Support two modes:
#' 1. Global mode: ranked gene list vs all cell death pathways
#' 2. Targeted mode: ranked gene list vs one specific pathway
#'
#' @param geneList A named numeric vector, sorted descendingly.
#' @param term_select Character. Either an **exact gene set name**
#' (e.g. "Disulfidptosis_Up_36747082") used directly, or a **cell death
#' type name** (e.g. "Disulfidptosis"): for the 13 types with curated
#' recommendations, the recommended representative gene set is used
#' automatically (see \code{\link{get_representative_pathway}}); for other
#' types with a single gene set, that gene set is used with a notice.
#' Default NULL uses all cell death pathways (global mode).
#' @param use_recommended Logical. Only takes effect when \code{term_select}
#' is a cell death type with a curated recommendation: TRUE (default) uses
#' only the recommended representative gene set; FALSE runs against ALL
#' gene sets under that death type.
#' @param pvalueCutoff P-value cutoff.
#' @param verbose Print messages.
#' @param savefile Whether to save the results as a CSV file (default: FALSE)
#' @param filename File name (default: "celldeath_gsea_result.csv")
#' @return A gseaResult object.
#' @seealso \code{\link{get_representative_pathway}} for recommendation details.
#' @export
#' @importFrom clusterProfiler GSEA
#' @importFrom utils write.csv
#' @examples
#' # Prepare ranked gene vector (descending)
#' gene_rank <- c(GPX4=2.5, ACSL4=2.2, SLC7A11=1.8, BAX=-1.2, CASP3=-1.9)
#'
#' # 1. Global GSEA: vs all cell death pathways
#' gsea_all <- celldeath_gsea(gene_rank)
#'
#' # 2. Targeted GSEA by death type: use the recommended gene set (default)
#' gsea_rec <- celldeath_gsea(gene_rank, term_select = "Ferroptosis")
#'
#' # 3. Targeted GSEA by death type: run against ALL Ferroptosis gene sets
#' gsea_all_type <- celldeath_gsea(gene_rank, term_select = "Ferroptosis",
#'                                  use_recommended = FALSE)
#'
#' # 4. Targeted GSEA by exact gene set name
#' gsea_exact <- celldeath_gsea(gene_rank, term_select = "Disulfidptosis_Up_36747082")
celldeath_gsea <- function(
    geneList,
    term_select = NULL,
    use_recommended = TRUE,
    pvalueCutoff = 1,
    verbose = FALSE,
    savefile = FALSE,
    filename = "celldeath_gsea_result.csv"
) {

  # =========== Input validation ============
  # 1. Must be numeric vector
  if (!is.numeric(geneList)) {
    stop("'geneList' must be a numeric vector.", call. = FALSE)
  }
  # 2. Must have non-empty names
  if (is.null(names(geneList)) || any(names(geneList) == "")) {
    stop("'geneList' must be a named numeric vector.", call. = FALSE)
  }
  # 3. Must not contain NA values
  if (any(is.na(geneList))) {
    stop("'geneList' cannot contain missing values (NA).", call. = FALSE)
  }
  # 4. Must be sorted in decreasing order
  if (is.unsorted(rev(geneList))) {
    stop("'geneList' must be sorted in decreasing order.", call. = FALSE)
  }

  # ========== Load cell death gene sets ===========
  death_gene_list_all <- load_cell_death_genes()

  # ============ Targeted mode: resolve user-specified pathway ============
  death_gene_list <- .resolve_term_select(term_select, death_gene_list_all,
                                          use_recommended = use_recommended)

  # ========== Convert to TERM2GENE format required by GSEA ==========
  term2gene <- data.frame(
    term = rep(names(death_gene_list), vapply(death_gene_list, length, integer(1))),
    gene = unlist(death_gene_list, use.names = FALSE),
    stringsAsFactors = FALSE
  )

  # ============ Run GSEA =============
  gsea_res <- clusterProfiler::GSEA(
    geneList = geneList,
    TERM2GENE = term2gene,
    pvalueCutoff = pvalueCutoff,
    verbose = verbose,
    pAdjustMethod = "BH"
  )

  # ========== Friendly message for empty results ============
  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
    if (is.null(term_select)) {
      message("No significantly enriched cell death pathways detected. Try relaxing pvalueCutoff and retry!")
    } else {
      message("No significant GSEA enrichment detected for pathway [",
                     paste(names(death_gene_list), collapse = ", "), "]!")
    }
    return(NULL)
  }

  # ========== Export CSV results ==========
  if (savefile) {
    result_df <- as.data.frame(gsea_res)
    utils::write.csv(result_df, file = filename, row.names = FALSE, fileEncoding = "UTF-8")
    message("Results have been saved to: ", filename)
  }

  return(gsea_res)
}


#' Multiple Group Cell Death Pathway GSEA Analysis
#'
#' This function performs GSEA for multiple gene lists.
#' Support two modes:
#' 1. Global mode: all groups vs all cell death pathways
#' 2. Targeted mode: all groups vs one specific pathway
#'
#' @param geneList_list A named list of numeric vectors, each vector is sorted descendingly.
#' @param term_select Character. Either an **exact gene set name**
#' (e.g. "Disulfidptosis_Up_36747082") used directly, or a **cell death
#' type name** (e.g. "Disulfidptosis"): for the 13 types with curated
#' recommendations, the recommended representative gene set is used
#' automatically (see \code{\link{get_representative_pathway}}); for other
#' types with a single gene set, that gene set is used with a notice.
#' Default NULL uses all cell death pathways (global mode).
#' @param use_recommended Logical. Only takes effect when \code{term_select}
#' is a cell death type with a curated recommendation: TRUE (default) uses
#' only the recommended representative gene set; FALSE runs against ALL
#' gene sets under that death type.
#' @param pvalueCutoff P-value cutoff.
#' @param verbose Print messages.
#' @param savefile Whether to save all results as CSV files (default: FALSE)
#' @return A named list of gseaResult objects.
#' @seealso \code{\link{get_representative_pathway}} for recommendation details.
#' @export
#' @importFrom clusterProfiler GSEA
#' @importFrom utils write.csv
#' @examples
#' # Multi-group ranked list
#' group1 <- c(GPX4=2.5, ACSL4=2.2, NUBPL=1.8, BAX=-1.2)
#' group2 <- c(GSDMD=2.8, SLC3A2=2.0, GPX4=-1.5)
#' gene_list <- list(Control = group1, Treatment = group2)
#'
#' # 1. Global mode: all pathways
#' res_all <- celldeath_gsea_multiple(gene_list)
#'
#' # 2. Targeted mode by death type: use the recommended gene set (default)
#' res_rec <- celldeath_gsea_multiple(gene_list, term_select = "Ferroptosis")
#'
#' # 3. Targeted mode by death type: run against ALL Ferroptosis gene sets
#' res_all_type <- celldeath_gsea_multiple(gene_list, term_select = "Ferroptosis",
#'                                          use_recommended = FALSE)
#'
#' # 4. Targeted mode by exact gene set name
#' res_exact <- celldeath_gsea_multiple(gene_list,
#'                                       term_select = "Disulfidptosis_Up_36747082")
celldeath_gsea_multiple <- function(
    geneList_list,
    term_select = NULL,
    use_recommended = TRUE,
    pvalueCutoff = 1,
    verbose = FALSE,
    savefile = FALSE
) {
  # ========== Input validity check ==========
  # 1. Check if input is a named list
  if (!is.list(geneList_list)) {
    stop("'geneList_list' must be a list.", call. = FALSE)
  }
  if (is.null(names(geneList_list)) || any(names(geneList_list) == "")) {
    stop("'geneList_list' must be a NAMED list (each group must have a name).", call. = FALSE)
  }

  # 2. Validate each group geneList (reuse single-group validation logic)
  for (group_name in names(geneList_list)) {
    gl <- geneList_list[[group_name]]
    if (!is.numeric(gl)) {
      stop("Group ", group_name, ": 'geneList' must be a numeric vector.", call. = FALSE)
    }
    if (is.null(names(gl)) || any(names(gl) == "")) {
      stop("Group ", group_name, ": 'geneList' must be a named numeric vector.", call. = FALSE)
    }
    if (any(is.na(gl))) {
      stop("Group ", group_name, ": 'geneList' cannot contain missing values (NA).", call. = FALSE)
    }
    if (is.unsorted(rev(gl))) {
      stop("Group ", group_name, ": 'geneList' must be sorted in decreasing order.", call. = FALSE)
    }
  }

  # ========== Load cell death gene sets + targeted pathway resolution ==========
  death_gene_list_all <- load_cell_death_genes()

  # Targeted mode: resolve user-specified pathway
  death_gene_list <- .resolve_term_select(term_select, death_gene_list_all,
                                          use_recommended = use_recommended)

  # Convert to TERM2GENE format required by GSEA
  term2gene <- data.frame(
    term = rep(names(death_gene_list), vapply(death_gene_list, length, integer(1))),
    gene = unlist(death_gene_list, use.names = FALSE),
    stringsAsFactors = FALSE
  )

  # ========== Loop for multi-group GSEA ==========
  gsea_result_list <- list()
  for (group_name in names(geneList_list)) {
    message("Running GSEA for group: ", group_name, " ...")
    # Get ranked gene list for current group
    current_geneList <- geneList_list[[group_name]]

    # Skip groups with no overlapping genes
    if (!any(names(current_geneList) %in% term2gene$gene)) {
      message(group_name, ": no genes overlap with the selected pathway(s), skipped.")
      gsea_result_list[group_name] <- list(NULL)
      next
    }

    # Run GSEA
    gsea_res <- clusterProfiler::GSEA(
      geneList = current_geneList,
      TERM2GENE = term2gene,
      pvalueCutoff = pvalueCutoff,
      verbose = verbose,
      pAdjustMethod = "BH"
    )

    # Store result
    gsea_result_list[[group_name]] <- gsea_res

    # Friendly message for empty result
    if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
      if (is.null(term_select)) {
        message(group_name, ": No significantly enriched cell death pathways detected")
      } else {
        message(group_name, ": No significant enrichment detected for pathway [",
                paste(names(death_gene_list), collapse = ", "), "]")
      }
    }

    # Save CSV file (saved directly in current working directory)
    if (savefile) {
      filename <- paste0(group_name, "_celldeath_gsea.csv")
      result_df <- as.data.frame(gsea_res)
      utils::write.csv(result_df, file = filename, row.names = FALSE, fileEncoding = "UTF-8")
      message("Result saved to: ", filename)
    }
  }

  message("All groups GSEA completed!")
  return(gsea_result_list)
}
