library(testthat)

# Mock gene set collection shared by resolver tests
# - "Ferroptosis": two gene sets (tests recommendation logic)
# - "Pyroptosis":  a single gene set (tests the "only one set" branch)
# - "Necroptosis": two gene sets (tests the "no recommendation yet" branch)
mock_gene_list <- list(
  Ferroptosis_FerrDb_V2_36305834 = c("GPX4", "ACSL4", "SLC7A11"),
  Ferroptosis_ALI_37149072       = c("FTH1", "TFRC", "HMOX1"),
  Pyroptosis_GeneCards           = c("GSDMD", "CASP1", "IL1B"),
  Necroptosis_MLKL_12345678      = c("MLKL", "RIPK3", "RIPK1"),
  Necroptosis_PanCancer_23456789 = c("RIPK3", "CYLD", "FADD")
)

# ===================== Test get_representative_pathway =====================

test_that("get_representative_pathway returns a one-row data.frame", {
  type_with_rec <- representative_pathways$Type[1]
  res <- get_representative_pathway(type_with_rec)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_equal(res$Type, type_with_rec)
  expect_true(all(c("Pathway", "PMID") %in% colnames(res)))
})

test_that("get_representative_pathway uses exact case-sensitive matching", {
  type_with_rec <- representative_pathways$Type[1]

  # Exact spelling should match
  res <- get_representative_pathway(type_with_rec)
  expect_s3_class(res, "data.frame")
  expect_equal(res$Type, type_with_rec)

  # Different case should not match
  expect_null(get_representative_pathway(tolower(type_with_rec)))
  expect_null(get_representative_pathway(toupper(type_with_rec)))
})

test_that("get_representative_pathway does not trim whitespace", {
  type_with_rec <- representative_pathways$Type[1]

  expect_null(
    get_representative_pathway(paste0(" ", type_with_rec))
  )

  expect_null(
    get_representative_pathway(paste0(type_with_rec, " "))
  )
})

test_that("get_representative_pathway returns NULL for an unknown type", {
  expect_null(
    get_representative_pathway("DefinitelyNotACellDeathType")
  )
})

test_that("get_representative_pathway validates death_type", {
  expect_error(
    get_representative_pathway(123),
    "death_type must be a single cell death type name"
  )

  expect_error(
    get_representative_pathway(c("Ferroptosis", "Pyroptosis")),
    "death_type must be a single cell death type name"
  )

  expect_error(
    get_representative_pathway(NA_character_),
    "death_type must be a single cell death type name"
  )
})

test_that("get_representative_pathway honours a valid option override", {
  opt_tbl <- data.frame(
    Type = "Customdeath",
    Pathway = "Customdeath_Mock_99999999",
    PMID = "99999999",
    stringsAsFactors = FALSE
  )

  withr::local_options(
    CDEnrich.representative_pathways = opt_tbl
  )

  res <- get_representative_pathway("Customdeath")

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_equal(res$Type, "Customdeath")
  expect_equal(res$Pathway, "Customdeath_Mock_99999999")
  expect_equal(res$PMID, "99999999")
})

test_that("option override completely replaces the built-in table", {
  opt_tbl <- data.frame(
    Type = "Customdeath",
    Pathway = "Customdeath_Mock_99999999",
    PMID = "99999999",
    stringsAsFactors = FALSE
  )

  withr::local_options(
    CDEnrich.representative_pathways = opt_tbl
  )

  # Ferroptosis exists in the built-in table but not in the override.
  # Because the override replaces the built-in table completely,
  # the lookup should return NULL.
  expect_null(
    get_representative_pathway("Ferroptosis")
  )
})

test_that("Pathway prefix does not replace an exact Type match", {
  opt_tbl <- data.frame(
    Type = "MockType",
    Pathway = "Mockdeath_Set_00000001",
    PMID = "00000001",
    stringsAsFactors = FALSE
  )

  withr::local_options(
    CDEnrich.representative_pathways = opt_tbl
  )

  # The query matches Pathway only, not Type.
  # The current implementation uses exact Type matching,
  # so it should return NULL.
  expect_null(
    get_representative_pathway("Mockdeath")
  )

  # The exact Type value still works.
  res <- get_representative_pathway("MockType")
  expect_equal(res$Type, "MockType")
  expect_equal(res$Pathway, "Mockdeath_Set_00000001")
})

# ===================== Test .resolve_term_select =====================

test_that(".resolve_term_select NULL returns the full gene list unchanged (global mode)", {
  res <- .resolve_term_select(NULL, mock_gene_list)
  expect_identical(res, mock_gene_list)
})

test_that(".resolve_term_select validates that term_select is a single character string", {
  expect_error(
    .resolve_term_select(123, mock_gene_list),
    "term_select must be a single pathway name"
  )
  expect_error(
    .resolve_term_select(c("Ferroptosis", "Pyroptosis"), mock_gene_list),
    "term_select must be a single pathway name"
  )
})

# ---- Case 1: exact gene set name is used directly ----
test_that(".resolve_term_select uses an exact gene set name directly without any message", {
  expect_no_message(
    res <- .resolve_term_select("Ferroptosis_ALI_37149072", mock_gene_list)
  )
  expect_equal(names(res), "Ferroptosis_ALI_37149072")
  expect_equal(res[[1]], c("FTH1", "TFRC", "HMOX1"))
})

# ---- Case 2a: death type with recommendation, use_recommended = TRUE ----
test_that(".resolve_term_select death type resolves to the recommended gene set by default", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) {
      data.frame(Type = death_type,
                 Pathway = "Ferroptosis_FerrDb_V2_36305834",
                 PMID = "36305834",
                 stringsAsFactors = FALSE)
    },
    .package = "CDEnrich"
  )

  expect_message(
    res <- .resolve_term_select("Ferroptosis", mock_gene_list),
    "recommended representative gene set"
  )
  expect_equal(names(res), "Ferroptosis_FerrDb_V2_36305834")
  expect_equal(res[[1]], c("GPX4", "ACSL4", "SLC7A11"))
})

# ---- Case 2b: death type with recommendation, use_recommended = FALSE ----
test_that(".resolve_term_select use_recommended=FALSE returns ALL gene sets of the type", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) {
      data.frame(Type = death_type,
                 Pathway = "Ferroptosis_FerrDb_V2_36305834",
                 PMID = "36305834",
                 stringsAsFactors = FALSE)
    },
    .package = "CDEnrich"
  )

  expect_message(
    res <- .resolve_term_select("Ferroptosis", mock_gene_list,
                                use_recommended = FALSE),
    "running against ALL"
  )
  expect_equal(sort(names(res)),
               c("Ferroptosis_ALI_37149072", "Ferroptosis_FerrDb_V2_36305834"))
})

# ---- Case 2 edge: recommendation exists but its pathway is missing from the list ----
test_that(".resolve_term_select errors when the recommended set is absent from the gene list", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) {
      data.frame(Type = death_type,
                 Pathway = "Ferroptosis_Missing_00000000",
                 PMID = "00000000",
                 stringsAsFactors = FALSE)
    },
    .package = "CDEnrich"
  )

  # The invalid recommendation is ignored, falling through to the
  # "multiple sets but no usable recommendation" error
  expect_error(
    .resolve_term_select("Ferroptosis", mock_gene_list),
    "has multiple gene sets but no recommended one yet"
  )
})

# ---- Case 3: death type with a single gene set and no recommendation ----
test_that(".resolve_term_select single-set death type is used with an informative message", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) NULL,
    .package = "CDEnrich"
  )

  expect_message(
    res <- .resolve_term_select("Pyroptosis", mock_gene_list),
    "has only one curated gene set"
  )
  expect_equal(names(res), "Pyroptosis_GeneCards")
  expect_equal(res[[1]], c("GSDMD", "CASP1", "IL1B"))
})

# ---- Case 3: death type with multiple gene sets but no recommendation ----
test_that(".resolve_term_select multi-set death type without recommendation asks for an exact name", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) NULL,
    .package = "CDEnrich"
  )

  expect_error(
    .resolve_term_select("Necroptosis", mock_gene_list),
    "has multiple gene sets but no recommended one yet"
  )
  # Error message lists the candidate gene set names
  expect_error(
    .resolve_term_select("Necroptosis", mock_gene_list),
    "Necroptosis_MLKL_12345678"
  )
})

# ---- Case 4: unknown name ----
test_that(".resolve_term_select unknown name errors and lists available gene sets", {
  expect_error(
    .resolve_term_select("FakePathway", mock_gene_list),
    "Pathway name does not exist"
  )
  expect_error(
    .resolve_term_select("FakePathway", mock_gene_list),
    "Ferroptosis_FerrDb_V2_36305834"
  )
})

# ---- Prefix matching does not confuse similar type names ----
test_that(".resolve_term_select only matches gene sets starting with '<type>_'", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) NULL,
    .package = "CDEnrich"
  )

  # "Necroptosis" must NOT match "Ferroptosis_*" sets etc.; and a bare
  # substring like "ptosis" should not match anything (unknown -> error)
  expect_error(
    .resolve_term_select("ptosis", mock_gene_list),
    "Pathway name does not exist"
  )
})

# ===================== Test .resolve_plot_pathway =====================

test_that(".resolve_plot_pathway returns an exact result ID directly", {
  result_ids <- c("Ferroptosis_FerrDb_V2_36305834",
                  "Ferroptosis_ALI_37149072",
                  "Apoptosis_PanCancer_37968457")

  expect_no_message(
    res <- .resolve_plot_pathway("Apoptosis_PanCancer_37968457", result_ids)
  )
  expect_equal(res, "Apoptosis_PanCancer_37968457")
  expect_type(res, "character")
})

test_that(".resolve_plot_pathway death type resolves to the recommended set present in the result", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) {
      data.frame(Type = death_type,
                 Pathway = "Ferroptosis_FerrDb_V2_36305834",
                 PMID = "36305834",
                 stringsAsFactors = FALSE)
    },
    .package = "CDEnrich"
  )
  result_ids <- c("Ferroptosis_FerrDb_V2_36305834",
                  "Ferroptosis_ALI_37149072",
                  "Apoptosis_PanCancer_37968457")

  expect_message(
    res <- .resolve_plot_pathway("Ferroptosis", result_ids),
    "recommended representative gene set"
  )
  expect_equal(res, "Ferroptosis_FerrDb_V2_36305834")

  # use_recommended = FALSE resolves to ALL sets of the type in the result
  expect_message(
    res_all <- .resolve_plot_pathway("Ferroptosis", result_ids,
                                     use_recommended = FALSE),
    "running against ALL"
  )
  expect_equal(sort(res_all),
               c("Ferroptosis_ALI_37149072", "Ferroptosis_FerrDb_V2_36305834"))
})

test_that(".resolve_plot_pathway errors for a pathway absent from the result", {
  result_ids <- c("Ferroptosis_FerrDb_V2_36305834", "Apoptosis_PanCancer_37968457")
  expect_error(
    .resolve_plot_pathway("NotExistPathway", result_ids),
    "Pathway name does not exist"
  )
})

test_that(".resolve_plot_pathway delegates resolution to .resolve_term_select", {
  resolve_calls <- list()
  local_mocked_bindings(
    .resolve_term_select = function(term_select, gene_list, use_recommended) {
      resolve_calls <<- append(resolve_calls,
                               list(list(term_select = term_select,
                                         candidate_ids = names(gene_list),
                                         use_recommended = use_recommended)))
      gene_list
    },
    .package = "CDEnrich"
  )
  result_ids <- c("Ferroptosis_FerrDb_V2_36305834", "Apoptosis_PanCancer_37968457")

  res <- .resolve_plot_pathway("Ferroptosis", result_ids, use_recommended = FALSE)

  expect_length(resolve_calls, 1)
  expect_equal(resolve_calls[[1]]$term_select, "Ferroptosis")
  expect_false(resolve_calls[[1]]$use_recommended)
  # Candidate names handed to the resolver are exactly the result IDs
  expect_equal(resolve_calls[[1]]$candidate_ids, result_ids)
  # Return value is the names of whatever the resolver returned
  expect_equal(res, result_ids)
})

# ===================== Test .wrap_pathway_labels =====================

test_that(".wrap_pathway_labels leaves short names unchanged", {
  expect_equal(.wrap_pathway_labels("Pyroptosis_GeneCards", width = 22),
               "Pyroptosis_GeneCards")
  expect_equal(.wrap_pathway_labels("Short", width = 22), "Short")
  expect_equal(.wrap_pathway_labels("", width = 22), "")
})

test_that(".wrap_pathway_labels wraps long names after underscores and keeps them", {
  wrapped <- .wrap_pathway_labels("Necroptosis_Glioma_35719999", width = 22)
  expect_equal(wrapped, "Necroptosis_Glioma_\n35719999")

  # Underscores are preserved: same count before and after wrapping
  orig <- "Necroptosis_Glioma_35719999"
  expect_equal(
    lengths(regmatches(wrapped, gregexpr("_", wrapped, fixed = TRUE))),
    lengths(regmatches(orig, gregexpr("_", orig, fixed = TRUE)))
  )
})

test_that(".wrap_pathway_labels keeps each wrapped line within width (+1 trailing underscore)", {
  wrapped <- .wrap_pathway_labels(
    "Immunogenic_Cell_Death_HCC_Prognosis_40918463", width = 22
  )
  lines <- strsplit(wrapped, "\n", fixed = TRUE)[[1]]
  expect_true(length(lines) > 1)
  # Each displayed line may carry one trailing separator underscore
  expect_true(all(nchar(lines) <= 22 + 1))
})

test_that(".wrap_pathway_labels returns overlong single-token names unchanged", {
  # No underscore to break at -> cannot wrap, returned as-is
  s <- "Supercalifragilisticexpialidocious"
  expect_equal(.wrap_pathway_labels(s, width = 10), s)
})

test_that(".wrap_pathway_labels is vectorized and unnamed", {
  input <- c("Short", "Necroptosis_Glioma_35719999")
  res <- .wrap_pathway_labels(input, width = 22)

  expect_type(res, "character")
  expect_length(res, 2)
  expect_null(names(res))
  expect_equal(res[1], "Short")
  expect_equal(res[2], "Necroptosis_Glioma_\n35719999")
})

test_that(".wrap_pathway_labels handles multiple wraps in one name", {
  wrapped <- .wrap_pathway_labels(
    "Ferroptosis_EGFRWTLUAD_Prognosis_39966805", width = 20
  )
  lines <- strsplit(wrapped, "\n", fixed = TRUE)[[1]]
  expect_true(length(lines) >= 2)
  # Reassembling the lines restores the original name exactly
  expect_equal(gsub("\n", "", wrapped, fixed = TRUE),
               "Ferroptosis_EGFRWTLUAD_Prognosis_39966805")
})

# ===================== Data integrity (real package data, no mocks) =====================

test_that("every recommended pathway exists in the curated gene set collection", {
  all_sets <- load_cell_death_genes()
  expect_true(all(representative_pathways$Pathway %in% names(all_sets)),
              info = paste("Missing gene sets:",
                           paste(setdiff(representative_pathways$Pathway,
                                         names(all_sets)),
                                 collapse = ", ")))
})

test_that("the recommendation table covers the 13 documented cell death types", {
  expect_equal(nrow(representative_pathways), 13)
  # No duplicated death types: one recommendation per type
  expect_equal(length(unique(representative_pathways$Type)),
               nrow(representative_pathways))
})

test_that("every recommended gene set is large enough for analysis (>= 3 genes)", {
  all_sets <- load_cell_death_genes()
  rec_sizes <- lengths(all_sets[representative_pathways$Pathway])
  expect_true(all(rec_sizes >= 3),
              info = paste("Too-small recommended sets:",
                           paste(names(rec_sizes)[rec_sizes < 3],
                                 collapse = ", ")))
})
