# data-raw/build_cd_genesets.R
# Run once to convert the corrected csv into data/cd_Genesets.rda;
# rerun this script whenever the source csv is updated

# Read the wide-format source csv:
# column 1 = pathway name, remaining columns = genes (one pathway per row)
raw <- read.csv(
  "data-raw/celldeathDATA.csv",
  stringsAsFactors = FALSE,
  fileEncoding = "GBK",
  check.names = FALSE
)

# Reshape from wide to long format (term + gene),
# trimming whitespace (some cells contain leading spaces)
pathways <- trimws(raw[[1]])
gene_mat <- as.matrix(raw[, -1, drop = FALSE])

cd_Genesets <- data.frame(
  term = rep(pathways, each = ncol(gene_mat)),   # <-- times 改为 each
  gene = trimws(as.vector(t(gene_mat))),
  stringsAsFactors = FALSE
)

# Remove empty cells (rows have unequal numbers of genes) and duplicates
cd_Genesets <- cd_Genesets[!is.na(cd_Genesets$gene) & cd_Genesets$gene != "", ]
cd_Genesets <- unique(cd_Genesets)

# Sanity check: the object name MUST stay "cd_Genesets",
# because load_cell_death_genes() loads it via data("cd_Genesets")
stopifnot(all(c("term", "gene") %in% colnames(cd_Genesets)))

# Save as package data, overwriting the old .rda
usethis::use_data(cd_Genesets, overwrite = TRUE)
