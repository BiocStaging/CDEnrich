#' ORA Enrichment Analysis for Cell Death Pathways (Single Group)
#'
#' Perform Over-Representation Analysis (ORA) for a list of DEGs
#' based on curated cell death gene sets.
#' Support two modes:
#' 1. Global mode: one group vs all cell death pathways
#' 2. Targeted mode: one group vs one specific pathway
#'
#' @param deg A character vector of input genes (gene symbols)
#' for enrichment analysis.
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
#' @param universe Character vector of background genes. Default NULL uses
#' all genes in the curated cell death collection. For formal analysis,
#' pass all genes detected in your experiment (e.g. all expressed genes).
#' @param pvalueCutoff P-value cutoff for significance (default: 0.05).
#' @param qvalueCutoff Q-value cutoff for filtering false positives (default: 0.2).
#' @param minGSSize Minimum gene set size to retain (default: 3).
#' @param savefile Whether to save the results as a CSV file (default: FALSE)
#' @param filename File name (default: "celldeath_ORA_result.csv")
#' @return An \code{enrichResult} object for visualization and result export.
#' @seealso \code{\link{get_representative_pathway}} for recommendation details.
#' @export
#' @importFrom clusterProfiler enricher
#' @importFrom utils write.csv
#' @examples
#' deg <- c("GPX4","ACSL4","SLC7A11","GSDMD","CASP3","BAX")
#'
#' # 1. Global ORA: vs all cell death pathways
#' enrich_res <- celldeath_enrich(deg)
#'
#' # 2. Targeted ORA by death type: use the recommended gene set (default)
#' enrich_res_rec <- celldeath_enrich(deg, term_select = "Ferroptosis")
#'
#' # 3. Targeted ORA by death type: run against ALL Ferroptosis gene sets
#' enrich_res_all <- celldeath_enrich(deg, term_select = "Ferroptosis",
#'                                    use_recommended = FALSE)
#'
#' # 4. Targeted ORA by exact gene set name
#' enrich_res_exact <- celldeath_enrich(deg, term_select = "Disulfidptosis_Up_36747082")
celldeath_enrich <- function(
    deg,
    term_select = NULL,
    use_recommended = TRUE,
    universe = NULL,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 3,
    savefile = FALSE,
    filename = "celldeath_ORA_result.csv"
) {
  # ========== Input validity check ==========
  if (!is.character(deg)) {
    stop("Parameter deg must be a character vector of gene symbols, e.g. c('GPX4','ACSL4')!")
  }
  # Remove empty strings and duplicate genes
  deg <- unique(deg[deg != ""])
  # Check sufficient gene quantity
  if (length(deg) < 3) {
    stop("The number of differential genes cannot be less than 3! Please add more genes and retry.")
  }

  # ========== Load cell death gene sets ==========
  death_gene_list_all <- load_cell_death_genes()

  # ========== Background universe (computed BEFORE pathway filtering) ==========
  # Default: all curated cell death genes.
  # Best practice: users should pass all genes detected in their experiment.
  if (is.null(universe)) {
    universe <- unique(unlist(death_gene_list_all, use.names = FALSE))
  }

  # ========== Targeted mode: resolve user-specified pathway ==========
  death_gene_list <- .resolve_term_select(term_select, death_gene_list_all,
                                          use_recommended = use_recommended)

  # ========== Convert to TERM2GENE format ==========
  term <- rep(names(death_gene_list), vapply(death_gene_list, length, integer(1)))
  gene <- unlist(death_gene_list, use.names = FALSE)

  term2gene <- data.frame(
    Term = term,
    Gene = gene,
    stringsAsFactors = FALSE
  )

  # clusterProfiler otherwise intersects the universe with genes represented
  # in TERM2GENE, shrinking targeted-mode background to the selected pathway.
  old_force_universe <- getOption("enrichment_force_universe", FALSE)
  options(enrichment_force_universe = TRUE)
  on.exit(options(enrichment_force_universe = old_force_universe), add = TRUE)

  # ========== ORA enrichment (core call of clusterProfiler::enricher) ==========
  enrich_result <- clusterProfiler::enricher(
    gene = deg,
    universe = universe,
    TERM2GENE = term2gene,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize,
    pAdjustMethod = "BH"
  )

  # Check for significant enrichment results (prompt user if no hits)
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) {
    message("No significantly enriched cell death pathways detected. Try adjusting pvalueCutoff or checking your differential genes!")
    return(NULL)
  }

  # ========== Export CSV ==========
  if (savefile) {
    result_df <- as.data.frame(enrich_result)
    utils::write.csv(result_df, file = filename, row.names = FALSE, fileEncoding = "UTF-8")
    message("Results have been saved to:", filename)
  }

  return(enrich_result)
}


#' ORA Enrichment Analysis for Cell Death Pathways (Multi-Group Comparison)
#'
#' Perform multi-group over-representation analysis (ORA)
#' to compare cell death pathway enrichment across multiple gene lists.
#' Support two modes:
#' 1. Global mode: all groups vs all cell death pathways
#' 2. Targeted mode: all groups vs one specific pathway
#'
#' @param deg_list_list A named list of gene vectors, where each element
#' corresponds to a group of differentially expressed genes (gene symbols).
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
#' @param universe Character vector of background genes. Default NULL uses
#' all genes in the curated cell death collection. For formal analysis,
#' pass all genes detected in your experiment (e.g. all expressed genes).
#' @param pvalueCutoff P-value cutoff for significance (default: 0.05).
#' @param qvalueCutoff Q-value cutoff for false positive control (default: 0.2).
#' @param minGSSize Minimum gene set size to retain (default: 3).
#' @param savefile Whether to save the results as a CSV file (default: FALSE)
#' @param filename File name (default: "cell_death_ORA_result.csv")
#' @return A \code{compareClusterResult} object for visualization.
#' Returns \code{NULL} if no significant pathways are identified.
#' @export
#' @importFrom clusterProfiler compareCluster enricher
#' @examples
#' group1 <- c("GPX4","ACSL4","SLC7A11","GSDMD")
#' group2 <- c("BAX","CASP3","RIPK1","FDX1")
#' deg_list <- list(Control = group1, Treatment = group2)
#'
#' # 1. Global mode: all groups vs all cell death pathways
#' compare_res_all <- celldeath_compare_enrich(deg_list)
#'
#' # 2. Targeted mode by death type: use the recommended gene set (default)
#' compare_res_rec <- celldeath_compare_enrich(deg_list, term_select = "Ferroptosis")
#'
#' # 3. Targeted mode by death type: run against ALL Ferroptosis gene sets
#' compare_res_type <- celldeath_compare_enrich(deg_list, term_select = "Ferroptosis",
#'                                              use_recommended = FALSE)
#'
#' # 4. Targeted mode by exact gene set name
#' compare_res_exact <- celldeath_compare_enrich(deg_list,
#'                                               term_select = "Disulfidptosis_Up_36747082")
celldeath_compare_enrich <- function(
    deg_list_list,
    term_select = NULL,
    use_recommended = TRUE,
    universe = NULL,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 3,
    savefile = FALSE,
    filename = "celldeath_compareORA_result.csv"
) {
  # ========== Input validity check ==========
  # Check if input is a list
  if (!is.list(deg_list_list)) {
    stop("Parameter deg_list_list must be a list! Each element should be a character vector of differential genes for one group.")
  }
  # Check if list has group names
  if (is.null(names(deg_list_list)) || any(names(deg_list_list) == "")) {
    stop("The list must be named! Names correspond to group labels (e.g. Control/Treatment).")
  }
  # Validate gene sets per group
  for (group_name in names(deg_list_list)) {
    deg_vec <- deg_list_list[[group_name]]
    # Check if it is a character vector
    if (!is.character(deg_vec)) {
      stop("Genes in group [", group_name, "] must be a character vector of gene symbols!")
    }
    # Remove empty strings and duplicate genes
    deg_vec <- unique(deg_vec[deg_vec != ""])
    # Check gene count threshold
    if (length(deg_vec) < 3) {
      stop("The number of differential genes in group [", group_name, "] cannot be less than 3!")
    }
    # Replace cleaned gene list
    deg_list_list[[group_name]] <- deg_vec
  }

  # ========== Load built-in gene sets + support single pathway filtering ==========
  death_gene_list_all <- load_cell_death_genes()

  # ========== Background universe (computed BEFORE pathway filtering) ==========
  # Default: all curated cell death genes.
  # Best practice: users should pass all genes detected in their experiment.
  if (is.null(universe)) {
    universe <- unique(unlist(death_gene_list_all, use.names = FALSE))
  }

  # Targeted mode: resolve user-specified pathway
  death_gene_list <- .resolve_term_select(term_select, death_gene_list_all,
                                          use_recommended = use_recommended)

  # Convert to TERM2GENE format
  term <- rep(names(death_gene_list), vapply(death_gene_list, length, integer(1)))
  gene <- unlist(death_gene_list, use.names = FALSE)
  term2gene <- data.frame(
    Term = term,
    Gene = gene,
    stringsAsFactors = FALSE
  )

  # Construct long-format data frame required by compareCluster (gene + group)
  gene_group_df <- data.frame(
    gene = unlist(deg_list_list, use.names = FALSE),
    group = rep(names(deg_list_list), vapply(deg_list_list, length, integer(1))),
    stringsAsFactors = FALSE
  )

  # Preserve the full background in targeted comparisons as well.
  old_force_universe <- getOption("enrichment_force_universe", FALSE)
  options(enrichment_force_universe = TRUE)
  on.exit(options(enrichment_force_universe = old_force_universe), add = TRUE)

  # Run multi-group comparative enrichment analysis
  compare_result <- clusterProfiler::compareCluster(
    gene ~ group,                # Formula: gene column ~ group column
    data = gene_group_df,
    fun = function(gene) {
      clusterProfiler::enricher(
        gene = gene,
        universe = universe,
        TERM2GENE = term2gene,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        minGSSize = minGSSize,
        pAdjustMethod = "BH"
      )
    }
  )

  # Check analysis results
  if (is.null(compare_result) || nrow(as.data.frame(compare_result)) == 0) {
    if (is.null(term_select)) {
      msg <- "No significantly enriched cell death pathways across multiple groups detected. Try relaxing pvalueCutoff and retry!"
    } else {
      msg <- paste0("No significant enrichment detected for pathway [", names(death_gene_list), "] in multi-group comparison!")
    }
    message(msg)
    return(NULL)
  }

  # ========== Export CSV ==========
  if (savefile) {
    result_df <- as.data.frame(compare_result)
    utils::write.csv(result_df, file = filename, row.names = FALSE, fileEncoding = "UTF-8")
    message("Results have been saved to:", filename)
  }


  return(compare_result)
}
