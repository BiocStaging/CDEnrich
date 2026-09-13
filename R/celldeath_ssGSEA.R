#' ssGSEA Analysis for Cell Death Pathways
#'
#' Perform single-sample Gene Set Enrichment Analysis (ssGSEA) on an
#' expression matrix based on curated cell death gene sets.
#' Unlike ORA/GSEA, ssGSEA computes an enrichment score for each sample
#' on each pathway, enabling sample-level comparison of cell death activity.
#' Support two modes:
#' 1. Global mode: expression matrix vs all cell death pathways
#' 2. Targeted mode: expression matrix vs one specific pathway
#'
#' Multi-group comparison is supported via the \code{group} parameter.
#' All samples are scored together on the same expression background to
#' guarantee score comparability across groups.
#'
#' @param expr A numeric matrix or data.frame of expression values,
#' with genes (gene symbols) as rows and samples as columns.
#' @param group An optional character or factor vector of group labels,
#' with length equal to the number of samples (columns of \code{expr}).
#' If provided, group information is attached to the result for downstream
#' differential analysis and visualization. Default NULL.
#' @param term_select Character. Either an **exact gene set name**
#' (e.g. "Disulfidptosis_Up_36747082") used directly, or a **cell death
#' type name** (e.g. "Disulfidptosis"): for the 13 types with curated
#' recommendations, the recommended representative gene set is used
#' automatically (see \code{\link{get_representative_pathway}}); for other
#' types with a single gene set, that gene set is used with a notice.
#' Default NULL uses all cell death pathways (global mode).
#' @param use_recommended Logical. Only takes effect when \code{term_select}
#' is a cell death type with a curated recommendation: TRUE (default) uses
#' only the recommended representative gene set; FALSE scores ALL
#' gene sets under that death type.
#' @param min_size Minimum gene set size to retain (default: 3).
#' @param normalize Logical. Whether to normalize ssGSEA scores by the
#' difference between the maximum and minimum scores (default: TRUE).
#' @param savefile Whether to save the results as a CSV file (default: FALSE)
#' @param filename File name (default: "celldeath_ssgsea_score.csv")
#' @return A named list of class \code{celldeath_ssgsea} with two elements:
#' \itemize{
#'   \item \code{score}: a numeric matrix of ssGSEA scores (pathways x samples);
#'   \item \code{group}: a factor of sample group labels (NULL if not provided).
#' }
#' @seealso \code{\link{get_representative_pathway}} for recommendation details.
#' @export
#' @importFrom GSVA ssgseaParam gsva
#' @importFrom utils write.csv
#' @examples
#' # expr: gene x sample expression matrix (rows = genes, columns = samples)
#' gene_sets <- load_cell_death_genes()
#' all_genes <- unique(unlist(gene_sets))
#' set.seed(123)
#' expr <- matrix(rnorm(length(all_genes) * 10),
#'                nrow = length(all_genes), ncol = 10,
#'                dimnames = list(all_genes, paste0("Sample", 1:10)))
#' group <- rep(c("Control", "Treatment"), each = 5)
#'
#' # 1. Global mode: all cell death pathways, no grouping
#' ssgsea_res <- celldeath_ssgsea(expr)
#'
#' # 2. Multi-group mode: with group labels for downstream comparison
#' ssgsea_res_group <- celldeath_ssgsea(expr, group = group)
#'
#' # 3. Targeted mode by death type: use the recommended gene set (default)
#' ssgsea_ferro <- celldeath_ssgsea(expr, group = group, term_select = "Ferroptosis")
#'
#' # 4. Targeted mode by death type: score ALL Ferroptosis gene sets
#' ssgsea_ferro_all <- celldeath_ssgsea(expr, group = group, term_select = "Ferroptosis",
#'                                      use_recommended = FALSE)
#'
#' # 5. Targeted mode by exact gene set name
#' ssgsea_exact <- celldeath_ssgsea(expr, group = group,
#'                                  term_select = "Disulfidptosis_Up_36747082")
celldeath_ssgsea <- function(
    expr,
    group = NULL,
    term_select = NULL,
    use_recommended = TRUE,
    min_size = 3,
    normalize = TRUE,
    savefile = FALSE,
    filename = "celldeath_ssgsea_score.csv"
) {
  # ========== Input validity check ==========
  # Convert data.frame to matrix
  if (is.data.frame(expr)) {
    expr <- as.matrix(expr)
  }
  # Must be a numeric matrix
  if (!is.matrix(expr) || !is.numeric(expr)) {
    stop("Parameter expr must be a numeric matrix or data.frame (genes x samples)!")
  }
  # Must have gene symbols as rownames
  if (is.null(rownames(expr)) || any(rownames(expr) == "")) {
    stop("The expression matrix must have gene symbols as rownames!")
  }
  # Must have sample names as colnames
  if (is.null(colnames(expr)) || any(colnames(expr) == "")) {
    stop("The expression matrix must have sample names as colnames!")
  }
  # Gene names must be unique
  if (anyDuplicated(rownames(expr))) {
    stop("Duplicated gene symbols detected in rownames! Please deduplicate before running ssGSEA (e.g. keep the probe with highest mean expression).")
  }
  # Must not contain NA values
  if (any(is.na(expr))) {
    stop("The expression matrix cannot contain missing values (NA). Please impute or remove them first!")
  }

  # ========== Validate group labels ==========
  if (!is.null(group)) {
    # Check length matches sample number
    if (length(group) != ncol(expr)) {
      stop(paste0(
        "The length of group (", length(group),
        ") does not match the number of samples (", ncol(expr), ")!"
      ))
    }
    # Check no empty labels or NA
    if (any(is.na(group)) || any(group == "")) {
      stop("Group labels cannot contain NA or empty strings!")
    }
    # Check at least 2 groups
    if (length(unique(group)) < 2) {
      stop("At least 2 distinct groups are required when group is provided!")
    }
    # Convert to factor for stable downstream handling
    group <- factor(group)
    # Each group should contain at least 2 samples
    if (any(table(group) < 2)) {
      warning("Some groups contain fewer than 2 samples. Differential analysis may be unreliable!")
    }
  }

  # ========== Load cell death gene sets ==========
  death_gene_list <- load_cell_death_genes()

  # ========== Targeted mode: resolve user-specified pathway ==========
  death_gene_list <- .resolve_term_select(term_select, death_gene_list,
                                          use_recommended = use_recommended)

  # ========== Filter gene sets by size ==========
  death_gene_list <- death_gene_list[vapply(death_gene_list, length, integer(1)) >= min_size]
  if (length(death_gene_list) == 0) {
    stop("No cell death pathway left after min_size filtering! Try lowering min_size.")
  }

  # ========== Check gene overlap between matrix and gene sets ==========
  overlap_n <- length(intersect(rownames(expr), unique(unlist(death_gene_list))))
  if (overlap_n == 0) {
    stop("No overlap between expression matrix rownames and cell death genes! Please check gene symbol format.")
  }
  if (overlap_n < min_size) {
    warning(
      "Only ", overlap_n, " cell death genes found in the expression matrix. ",
      "Results may be unreliable!"
    )
  }

  # ========== Run ssGSEA (all samples scored on the same background) ==========
  ssgsea_param <- GSVA::ssgseaParam(
    exprData = expr,
    geneSets = death_gene_list,
    minSize = min_size,
    normalize = normalize
  )
  ssgsea_score <- GSVA::gsva(ssgsea_param, verbose = FALSE)

  # ========== Assemble result object ==========
  result <- list(
    score = ssgsea_score,
    group = group
  )
  class(result) <- "celldeath_ssgsea"

  # ========== Export CSV ==========
  if (savefile) {
    result_df <- data.frame(
      Pathway = rownames(ssgsea_score),
      as.data.frame(ssgsea_score, check.names = FALSE),
      check.names = FALSE
    )
    utils::write.csv(result_df, file = filename, row.names = FALSE, fileEncoding = "UTF-8")
    message("Results have been saved to: ", filename)
  }

  return(result)
}


#' Differential Analysis of ssGSEA Scores Between Groups
#'
#' Perform group-wise differential analysis on cell death pathway
#' ssGSEA scores. Wilcoxon rank-sum test is used for two groups,
#' and Kruskal-Wallis test for three or more groups.
#'
#' @param ssgsea_result A \code{celldeath_ssgsea} object returned by
#' \code{\link{celldeath_ssgsea}} (must have non-NULL \code{group}).
#' @param pAdjustMethod Multiple testing correction method passed to
#' \code{\link[stats]{p.adjust}} (default: "BH").
#' @param savefile Whether to save the results as a CSV file (default: FALSE)
#' @param filename File name (default: "celldeath_ssgsea_diff.csv")
#' @return A data.frame with per-pathway statistics, sorted by p-value:
#' pathway name, mean score per group, test statistic, p-value and
#' adjusted p-value.
#' @export
#' @importFrom stats wilcox.test kruskal.test p.adjust
#' @importFrom utils write.csv
#' @examples
#' # expr: gene x sample expression matrix (rows = genes, columns = samples)
#' gene_sets <- load_cell_death_genes()
#' all_genes <- unique(unlist(gene_sets))
#' set.seed(123)
#' expr <- matrix(rnorm(length(all_genes) * 10),
#'                nrow = length(all_genes), ncol = 10,
#'                dimnames = list(all_genes, paste0("Sample", 1:10)))
#' group <- rep(c("Control", "Treatment"), each = 5)
#'
#' ssgsea_res <- celldeath_ssgsea(expr, group = group)
#' diff_res <- celldeath_ssgsea_diff(ssgsea_res)
#' head(diff_res)
celldeath_ssgsea_diff <- function(
    ssgsea_result,
    pAdjustMethod = "BH",
    savefile = FALSE,
    filename = "celldeath_ssgsea_diff.csv"
) {
  # ========== Input validity check ==========
  if (!inherits(ssgsea_result, "celldeath_ssgsea")) {
    stop("Parameter ssgsea_result must be the output of celldeath_ssgsea()!")
  }
  if (is.null(ssgsea_result$group)) {
    stop("No group information found! Please provide the group parameter when running celldeath_ssgsea().")
  }

  score_mat <- ssgsea_result$score
  group <- ssgsea_result$group
  group_levels <- levels(group)
  n_group <- length(group_levels)

  # ========== Per-pathway differential test ==========
  diff_list <- lapply(rownames(score_mat), function(pathway) {
    score_vec <- score_mat[pathway, ]

    # Mean score per group
    group_means <- tapply(score_vec, group, mean)

    # Select test according to group number
    if (n_group == 2) {
      test_res <- stats::wilcox.test(score_vec[group == group_levels[1]],
                                     score_vec[group == group_levels[2]],
                                     exact = FALSE)
      stat_value <- unname(test_res$statistic)
      stat_name <- "W"
    } else {
      test_res <- stats::kruskal.test(score_vec, group)
      stat_value <- unname(test_res$statistic)
      stat_name <- "Kruskal-Wallis chi-squared"
    }

    data.frame(
      Pathway = pathway,
      t(data.frame(group_means, check.names = FALSE)),
      Statistic = stat_value,
      P_value = test_res$p.value,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  # Merge per-pathway results
  diff_df <- do.call(rbind, diff_list)

  # ========== Multiple testing correction ==========
  diff_df$P_adjust <- stats::p.adjust(diff_df$P_value, method = pAdjustMethod)

  # Sort by raw p-value
  diff_df <- diff_df[order(diff_df$P_value), ]
  rownames(diff_df) <- NULL

  # ========== Export CSV ==========
  if (savefile) {
    utils::write.csv(diff_df, file = filename, row.names = FALSE, fileEncoding = "UTF-8")
    message("Results have been saved to: ", filename)
  }

  return(diff_df)
}


#' Print Method for celldeath_ssgsea Objects
#'
#' @param x A \code{celldeath_ssgsea} object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.celldeath_ssgsea <- function(x, ...) {
  cat("Cell death ssGSEA result\n")
  cat("  Pathways scored:", nrow(x$score), "\n")
  cat("  Samples scored :", ncol(x$score), "\n")
  if (!is.null(x$group)) {
    cat("  Groups         :", paste(levels(x$group), collapse = ", "), "\n")
    cat("  Samples/group  :", paste(paste0(levels(x$group), "=", table(x$group)), collapse = ", "), "\n")
  } else {
    cat("  Groups         : none (run with group = ... for multi-group comparison)\n")
  }
  invisible(x)
}
