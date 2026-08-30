# Build reproducible demo data for CDEnrich documentation and vignettes

devtools::load_all(".")

set.seed(123)

# ===================== Shared gene universe =====================

gene_sets <- load_cell_death_genes()
ferroptosis_id <- get_representative_pathway("Ferroptosis")$Pathway
ferroptosis_genes <- unique(gene_sets[[ferroptosis_id]])

background_genes <- paste0("BackgroundGene_", seq_len(3000))
demo_all_genes <- unique(c(ferroptosis_genes, background_genes))

# ===================== Ferroptosis ORA demo data =====================

stopifnot(
  is.character(ferroptosis_id),
  length(ferroptosis_id) == 1L,
  ferroptosis_id %in% names(gene_sets)
)

ferroptosis_genes <- unique(as.character(gene_sets[[ferroptosis_id]]))
ferroptosis_genes <- ferroptosis_genes[
  !is.na(ferroptosis_genes) &
    nzchar(ferroptosis_genes)
]

if (length(ferroptosis_genes) < 3L) {
  stop("The recommended Ferroptosis gene set contains fewer than 3 valid genes.")
}

# Use up to 100 real genes from the recommended Ferroptosis gene set.
demo_deg <- head(
  ferroptosis_genes,
  min(100L, length(ferroptosis_genes))
)

# Use a separate object name for the targeted ORA example.
demo_targeted_deg <- demo_deg

# Two-group ORA demonstration data
demo_deg_list <- list(
  Control = head(
    ferroptosis_genes,
    min(15L, length(ferroptosis_genes))
  ),
  Treatment = demo_deg
)

# ===================== Apoptosis GSEA demo data =====================

# Build a ranked list with a strong Apoptosis signal using the
# representative Apoptosis gene set (PMID 37968457)

apoptosis_id <- get_representative_pathway("Apoptosis")$Pathway

stopifnot(
  is.character(apoptosis_id),
  length(apoptosis_id) == 1L,
  apoptosis_id %in% names(gene_sets)
)

apoptosis_genes <- unique(as.character(gene_sets[[apoptosis_id]]))
apoptosis_genes <- apoptosis_genes[
  !is.na(apoptosis_genes) & nzchar(apoptosis_genes)
]

if (length(apoptosis_genes) < 3L) {
  stop("The representative Apoptosis gene set contains fewer than 3 valid genes.")
}

demo_all_genes_apop <- unique(c(apoptosis_genes, background_genes))

demo_rank_apop <- setNames(
  rnorm(length(demo_all_genes_apop)),
  demo_all_genes_apop
)

# Place all Apoptosis genes at the top of the ranked list
demo_rank_apop[apoptosis_genes] <- seq(
  from = 8,
  to = 3,
  length.out = length(apoptosis_genes)
)

demo_rank_apop <- sort(demo_rank_apop, decreasing = TRUE)

# ===================== Multi-group GSEA demo data =====================

# Control: no intentionally introduced Apoptosis signal
demo_rank_control <- sort(
  setNames(rnorm(length(demo_all_genes_apop)), demo_all_genes_apop),
  decreasing = TRUE
)

# Treatment: same strong Apoptosis signal as the single-group example
demo_rank_treatment <- demo_rank_apop

demo_rank_apop_list <- list(
  Control = demo_rank_control,
  Treatment = demo_rank_treatment
)

# ===================== ssGSEA demo data =====================

demo_group <- factor(rep(c("Control", "Treatment"), each = 6))

demo_expr <- matrix(
  rnorm(length(demo_all_genes) * length(demo_group), mean = 8, sd = 1.5),
  nrow = length(demo_all_genes),
  dimnames = list(demo_all_genes, paste0("Sample", seq_along(demo_group)))
)

# Create higher Ferroptosis activity in the Treatment group
demo_expr[ferroptosis_genes, demo_group == "Treatment"] <-
  demo_expr[ferroptosis_genes, demo_group == "Treatment"] + 2

# Save all public demo objects into data/*.rda
usethis::use_data(
  demo_all_genes,
  demo_deg,
  demo_targeted_deg,
  demo_deg_list,
  demo_rank_apop,
  demo_rank_apop_list,
  demo_expr,
  demo_group,
  overwrite = TRUE
)
