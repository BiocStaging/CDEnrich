#' Get the Recommended Representative Gene Set for a Cell Death Type
#'
#' Look up the curated recommendation for a cell death type that has
#' multiple PMID-based gene sets. 13 of the cell death types in this
#' package have a recommended representative gene set.
#'
#' The lookup table can be overridden via
#' \code{options(CDEnrich.representative_pathways = <data.frame>)}
#' (mainly for testing); the override replaces the built-in table
#' completely for the duration it is set.
#'
#' @param death_type Character(1). A cell death type name, e.g. "Disulfidptosis".
#'   Matching is exact (case-sensitive).
#' @return A one-row data.frame with the recommendation metadata
#'   (pathway name, PMID, source, etc.), or NULL if the death type has
#'   no recommendation.
#' @export
#' @examples
#' get_representative_pathway("Disulfidptosis")
#' get_representative_pathway("Mitosis")   # returns NULL
get_representative_pathway <- function(death_type) {
  # Input validation with a readable error message
  if (!is.character(death_type) || length(death_type) != 1 || is.na(death_type)) {
    stop("Error: death_type must be a single cell death type name (character of length 1)!")
  }

  # Option override (mainly for tests) replaces the built-in table entirely
  tbl <- getOption("CDEnrich.representative_pathways",
                   default = representative_pathways)

  # Exact match on the Type column
  hit <- tbl[tbl$Type == death_type, , drop = FALSE]
  if (nrow(hit) == 0) return(NULL)
  hit[1, , drop = FALSE]
}



#' Resolve User-Specified Pathway Selection (Internal)
#'
#' Resolve the \code{term_select} argument shared by analysis and
#' visualization functions. Handles four cases:
#' 1. Exact gene set name -> use it directly;
#' 2. Cell death type with a curated recommendation -> use the recommended
#'    representative gene set (default), or ALL gene sets of this type
#'    when \code{use_recommended = FALSE};
#' 3. Cell death type with a single gene set (no recommendation) ->
#'    use it and inform the user that only one gene set exists;
#' 4. Unknown name -> stop and list all available gene sets.
#'
#' @param term_select Character(1) or NULL. User-specified pathway name
#'   or cell death type. NULL returns \code{gene_list} unchanged.
#' @param gene_list Named list of gene vectors, as returned by
#'   \code{\link{load_cell_death_genes}}.
#' @param use_recommended Logical. When \code{term_select} is a death type
#'   with a recommendation: TRUE (default) uses only the recommended
#'   representative gene set; FALSE uses ALL gene sets under this type.
#' @return A subset of \code{gene_list} containing the resolved gene set(s).
#' @keywords internal
#' @noRd
.resolve_term_select <- function(term_select, gene_list, use_recommended = TRUE) {

  # NULL -> global mode, return the full list unchanged
  if (is.null(term_select)) {
    return(gene_list)
  }

  # Defensive check: only a single pathway/type is supported
  if (!is.character(term_select) || length(term_select) != 1) {
    stop("Error: term_select must be a single pathway name (character of length 1)!")
  }

  # ---- Case 1: exact gene set name -> use directly ----
  if (term_select %in% names(gene_list)) {
    return(gene_list[term_select])
  }

  # Match gene sets whose names start with "<type>_" (e.g. "Ferroptosis_")
  matched <- names(gene_list)[startsWith(names(gene_list), paste0(term_select, "_"))]

  # ---- Case 2: death type with a curated recommendation ----
  rec <- get_representative_pathway(term_select)
  if (!is.null(rec) && rec$Pathway %in% names(gene_list)) {

    if (use_recommended) {
      # 2a. Use ONLY the recommended representative gene set
      message(sprintf(
        paste0(
          "'%s' has %d curated gene sets. ",
          "Using the recommended representative gene set:\n",
          "  %s (PMID: %s)\n",
          "To run against ALL gene sets of this type, set use_recommended = FALSE; ",
          "to use a specific one, pass its exact name via term_select."
        ),
        term_select, length(matched), rec$Pathway, rec$PMID
      ))
      return(gene_list[rec$Pathway])

    } else {
      # 2b. Use ALL gene sets under this death type
      message(sprintf(
        paste0(
          "'%s' has %d curated gene sets. use_recommended = FALSE, ",
          "running against ALL of them:\n  %s\n",
          "(Recommended representative set: %s, PMID: %s)"
        ),
        term_select, length(matched), paste(matched, collapse = "\n  "),
        rec$Pathway, rec$PMID
      ))
      return(gene_list[matched])
    }
  }

  # ---- Case 3: death type without recommendation ----
  if (length(matched) == 1) {
    # Single gene set available (from GeneCards) -> use it and inform the user
    message(sprintf(
      paste0(
        "'%s' has only one curated gene set (from GeneCards):\n",
        "  %s\n",
        "Using it for the analysis."
      ),
      term_select, matched
    ))
    return(gene_list[matched])
  }

  if (length(matched) > 1) {
    # Multiple gene sets but no recommendation -> ask for an exact name
    stop(sprintf(
      paste0(
        "'%s' has multiple gene sets but no recommended one yet. ",
        "Please specify one exactly:\n  %s"
      ),
      term_select, paste(matched, collapse = "\n  ")
    ), call. = FALSE)
  }

  # ---- Case 4: unknown name -> stop and list all available gene sets ----
  stop(paste0("Error: Pathway name does not exist! Available pathways: ",
              paste(names(gene_list), collapse = ", ")))
}


#' Resolve a Pathway Name Against Result IDs (Internal)
#'
#' Thin wrapper around \code{.resolve_term_select} for plotting functions:
#' the candidate names are the pathway IDs present in the result object
#' rather than the full gene set collection. Accepts an exact pathway ID
#' or a cell death type name (resolved to the recommended set if present
#' in the result).
#'
#' @param pathway Character(1). Pathway ID or cell death type name.
#' @param result_ids Character vector of pathway IDs present in the result.
#' @param use_recommended Logical. Passed to \code{.resolve_term_select}.
#' @return Character vector of resolved pathway ID(s).
#' @keywords internal
#' @noRd
.resolve_plot_pathway <- function(pathway, result_ids, use_recommended = TRUE) {
  dummy_list <- stats::setNames(as.list(result_ids), result_ids)
  names(.resolve_term_select(pathway, dummy_list, use_recommended = use_recommended))
}


# Wrap long pathway names at underscores, keeping the underscores visible
.wrap_pathway_labels <- function(x, width = 22) {
  vapply(x, function(s) {
    if (nchar(s) <= width) return(s)
    parts <- strsplit(s, "_", fixed = TRUE)[[1]]
    lines <- parts[1]
    for (part in parts[-1]) {
      current <- lines[length(lines)]
      if (nchar(current) + nchar(part) + 1 <= width) {
        lines[length(lines)] <- paste0(current, "_", part)
      } else {
        lines <- c(lines, part)
      }
    }
    paste0(lines, c(rep("_\n", length(lines) - 1), ""), collapse = "")
  }, character(1), USE.NAMES = FALSE)
}
