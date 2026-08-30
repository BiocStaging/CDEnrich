# data-raw/build_representative_pathways.R
# Run once to convert the csv into data/*.rda; rerun after updating the csv

raw <- read.csv(
  "data-raw/representative_pathways.csv",
  stringsAsFactors = FALSE,
  fileEncoding = "GBK",   # change to "UTF-8" if your csv is UTF-8 encoded
  check.names = FALSE
)

meta_cols <- c("Type", "PMID", "IF", "Cite", "Pathway",
               "Source", "Description", "Notes")

# 1) Metadata table: keep a record of the rationale for each recommendation
representative_pathways <- as.data.frame(raw[, meta_cols], stringsAsFactors = FALSE)
# Coerce PMID to character to accommodate non-numeric values such as "GeneCards"
representative_pathways$PMID <- as.character(representative_pathways$PMID)
# Trim whitespace from key matching columns
representative_pathways$Type <- trimws(representative_pathways$Type)
representative_pathways$Pathway <- trimws(representative_pathways$Pathway)

# 2) Gene sets: collect non-empty genes row-wise starting from column 9
# Select gene columns by POSITION (not by name) to avoid issues with
# empty trailing header cells exported from Excel
gene_cols <- setdiff(seq_along(raw), match(meta_cols, names(raw)))
representative_genes <- lapply(seq_len(nrow(raw)), function(i) {
  genes <- trimws(as.character(unlist(raw[i, gene_cols], use.names = FALSE)))
  unique(genes[!is.na(genes) & genes != ""])
})
names(representative_genes) <- representative_pathways$Pathway

usethis::use_data(representative_pathways, representative_genes, overwrite = TRUE)
