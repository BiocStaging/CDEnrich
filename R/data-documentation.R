#' Curated Data Included in CDEnrich
#'
#' Data objects used by CDEnrich to provide curated cell death pathway
#' gene sets and representative pathway recommendations.
#'
#' \itemize{
#'   \item \code{\link{cd_Genesets}}: curated cell death gene sets.
#'   \item \code{\link{representative_genes}}: representative gene sets
#'   for the 13 cell death types with curated recommendations.
#'   \item \code{\link{representative_pathways}}: metadata for the
#'   representative pathway recommendations.
#' }
#'
#' @return This help page is an overview of the curated datasets shipped
#'   with CDEnrich; see the linked pages for each dataset.
#'
#' @name CDEnrich-data
#' @keywords datasets
NULL

#' Curated Cell Death Pathway Gene Sets
#'
#' Curated cell death pathway gene sets collected from published
#' literature, GeneCards and pathway databases. This is the underlying
#' table loaded by \code{\link{load_cell_death_genes}}.
#'
#' @docType data
#' @name cd_Genesets
#' @usage data(cd_Genesets)
#' @format A data frame with two columns:
#' \describe{
#'   \item{term}{Pathway (gene set) name, e.g.
#'   \code{Ferroptosis_FerrDb V2(36305834)}.}
#'   \item{gene}{Gene symbol belonging to the pathway.}
#' }
#' @source Curated from public cell death literature and databases;
#' see \code{data-raw/} for the construction scripts.
#' @keywords datasets
"cd_Genesets"

#' Representative Gene Sets for Cell Death Types
#'
#' The recommended representative gene set for each of the 13 cell death
#' types with curated recommendations. These are merged into the result
#' of \code{\link{load_cell_death_genes}} and are used as the default
#' gene set in targeted analysis modes.
#'
#' @docType data
#' @name representative_genes
#' @usage data(representative_genes)
#' @format A named list of character vectors; names are pathway names
#'   (e.g. \code{Ferroptosis_FerrDb V2(36305834)}) and each element is a
#'   vector of gene symbols.
#' @source Curated from the seminal publication of each cell death type;
#' see \code{\link{representative_pathways}} for the curation metadata.
#' @keywords datasets
"representative_genes"

#' Representative Pathway Recommendation Metadata
#'
#' Curated recommendation table recording, for each cell death type with
#' multiple PMID-based gene sets, which gene set is the recommended
#' representative. Used by \code{\link{get_representative_pathway}}.
#'
#' @docType data
#' @name representative_pathways
#' @usage data(representative_pathways)
#' @format A data frame with one row per recommended cell death type
#'   and the following columns:
#' \describe{
#'   \item{Type}{Cell death type name, e.g. \code{Ferroptosis}.}
#'   \item{PMID}{PubMed ID of the source publication.}
#'   \item{IF}{Journal impact factor of the source publication.}
#'   \item{Cite}{Citation count of the source publication.}
#'   \item{Pathway}{Name of the recommended representative gene set,
#'   e.g. \code{Ferroptosis_FerrDb V2(36305834)}.}
#'   \item{Source}{Where the gene set was obtained from in the
#'   publication, e.g. \code{Supplementary-table1}.}
#'   \item{Description}{Free-text description of how the gene set was
#'   curated.}
#'   \item{Notes}{Additional curation notes (may be empty).}
#' }
#' @source Curated manually from the primary literature;
#' see \code{data-raw/representative_pathways.csv}.
#' @keywords datasets
"representative_pathways"
