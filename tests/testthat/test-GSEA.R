library(testthat)

# Mock internal gene set function
local_mocked_bindings(
  load_cell_death_genes = function() {
    list(
      Ferroptosis = c("GPX4", "ACSL4", "SLC7A11"),
      Pyroptosis  = c("GSDMD", "CASP1", "IL1B"),
      Apoptosis   = c("CASP3", "BAX", "BCL2")
    )
  },
  .package = "CDEnrich"
)

# ===================== Test celldeath_gsea (single group) =====================

test_that("celldeath_gsea passes correct arguments to GSEA", {
  # Prepare valid geneList
  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1, SLC7A11 = 1.8, GSDMD = 1.5, CASP3 = 1.2)

  captured_args <- list()
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, pvalueCutoff, verbose, pAdjustMethod, ...) {
      captured_args <<- as.list(environment())
      data.frame(ID = "Ferroptosis", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )

  result <- celldeath_gsea(geneList, pvalueCutoff = 0.05, verbose = FALSE)

  expect_equal(captured_args$geneList, geneList)
  expect_true(is.data.frame(captured_args$TERM2GENE))
  expect_equal(captured_args$pvalueCutoff, 0.05)
  expect_equal(captured_args$verbose, FALSE)
  expect_equal(captured_args$pAdjustMethod, "BH")
})

test_that("celldeath_gsea passes correct arguments to GSEA (global mode term_select=NULL)", {
  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1, SLC7A11 = 1.8, GSDMD = 1.5, CASP3 = 1.2)

  captured_args <- list()
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, pvalueCutoff, verbose, pAdjustMethod, ...) {
      captured_args <<- as.list(environment())
      data.frame(ID = "Ferroptosis", pvalue = 0.01, NES = 1.6)
    },
    .package = "clusterProfiler"
  )

  result <- celldeath_gsea(geneList, pvalueCutoff = 0.05, verbose = FALSE)

  expect_equal(captured_args$geneList, geneList)
  expect_true(is.data.frame(captured_args$TERM2GENE))
  expect_equal(captured_args$pvalueCutoff, 0.05)
  expect_equal(captured_args$verbose, FALSE)
  expect_equal(captured_args$pAdjustMethod, "BH")
})

test_that("celldeath_gsea term_select specifies single pathway and filters TERM2GENE", {
  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1)
  captured_t2g <- NULL
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "Ferroptosis", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  celldeath_gsea(geneList, term_select = "Ferroptosis")
  expect_equal(unique(captured_t2g$term), "Ferroptosis")
})

# ---- New: term_select as an EXACT gene set name (used directly) ----
test_that("celldeath_gsea term_select as exact gene set name is used directly", {
  local_mocked_bindings(
    load_cell_death_genes = function() {
      list(
        Ferroptosis_Up_36747082   = c("GPX4", "ACSL4", "SLC7A11"),
        Ferroptosis_Down_36747082 = c("FTH1", "TFRC"),
        Pyroptosis                = c("GSDMD", "CASP1", "IL1B")
      )
    },
    .package = "CDEnrich"
  )
  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1, SLC7A11 = 1.8)
  captured_t2g <- NULL
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  celldeath_gsea(geneList, term_select = "Ferroptosis_Up_36747082")
  expect_equal(unique(captured_t2g$term), "Ferroptosis_Up_36747082")
})

# ---- New: use_recommended controls gene sets for a multi-set death type
#           (integration through .resolve_term_select) ----
test_that("celldeath_gsea use_recommended controls gene sets for a multi-set death type", {
  local_mocked_bindings(
    load_cell_death_genes = function() {
      list(
        Ferroptosis_Up_36747082   = c("GPX4", "ACSL4", "SLC7A11"),
        Ferroptosis_Down_36747082 = c("FTH1", "TFRC"),
        Pyroptosis                = c("GSDMD", "CASP1", "IL1B")
      )
    },

    get_representative_pathway = function(death_type) {
      list(Pathway = "Ferroptosis_Up_36747082", PMID = "36747082")
    },
    .package = "CDEnrich"
  )
  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1, SLC7A11 = 1.8, FTH1 = 1.5, TFRC = 1.2, GSDMD = 1.0)
  captured_t2g <- NULL
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )

  # use_recommended = TRUE (default): only the recommended representative gene set
  celldeath_gsea(geneList, term_select = "Ferroptosis", use_recommended = TRUE)
  expect_equal(sort(unique(captured_t2g$term)), "Ferroptosis_Up_36747082")

  # use_recommended = FALSE: ALL gene sets under Ferroptosis
  celldeath_gsea(geneList, term_select = "Ferroptosis", use_recommended = FALSE)
  expect_equal(sort(unique(captured_t2g$term)),
               c("Ferroptosis_Down_36747082", "Ferroptosis_Up_36747082"))
})

# ---- New: .resolve_term_select is delegated the right arguments ----
test_that("celldeath_gsea delegates term_select/use_recommended to .resolve_term_select", {
  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1)
  resolve_calls <- list()
  local_mocked_bindings(
    .resolve_term_select = function(term_select, death_gene_list_all, use_recommended) {
      resolve_calls <<- append(resolve_calls,
                               list(list(term_select = term_select,
                                         use_recommended = use_recommended,
                                         n_sets = length(death_gene_list_all))))
      list(Resolved = c("GPX4", "ACSL4"))
    },
    .package = "CDEnrich"
  )
  local_mocked_bindings(
    GSEA = function(...) data.frame(ID = "x", pvalue = 0.01),
    .package = "clusterProfiler"
  )
  celldeath_gsea(geneList, term_select = "Ferroptosis", use_recommended = FALSE)
  expect_length(resolve_calls, 1)
  expect_equal(resolve_calls[[1]]$term_select, "Ferroptosis")
  expect_false(resolve_calls[[1]]$use_recommended)
  expect_equal(resolve_calls[[1]]$n_sets, 3)
})

# ---- New: default pvalueCutoff (1) and verbose (FALSE) forwarded to GSEA ----
test_that("celldeath_gsea passes default pvalueCutoff=1/verbose=FALSE to GSEA when not specified", {
  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1, SLC7A11 = 1.8)
  captured <- list()
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, pvalueCutoff, verbose, pAdjustMethod, ...) {
      captured$args <<- as.list(environment())
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  celldeath_gsea(geneList)
  expect_equal(captured$args$pvalueCutoff, 1)
  expect_equal(captured$args$verbose, FALSE)
})

# ---- New: verbose=TRUE is forwarded to GSEA ----
test_that("celldeath_gsea forwards verbose=TRUE to GSEA", {
  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1)
  captured <- list()
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, pvalueCutoff, verbose, pAdjustMethod, ...) {
      captured$verbose <<- verbose
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  celldeath_gsea(geneList, verbose = TRUE)
  expect_true(captured$verbose)
})

test_that("celldeath_gsea throws error when specified pathway does not exist", {
  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1)
  expect_error(
    celldeath_gsea(geneList, term_select = "FakePath"),
    "Pathway name does not exist"
  )
})

test_that("celldeath_gsea input validation: illegal inputs trigger errors", {
  # Not numeric vector
  expect_error(celldeath_gsea(c("a", "b")), "numeric vector")
  # No names
  expect_error(celldeath_gsea(c(1, 2, 3)), "named numeric vector")
  # Incomplete names (empty name exists)
  expect_error(celldeath_gsea(c(a = 1, 2, b = 3)), "named numeric vector")
  # Contains NA
  expect_error(celldeath_gsea(c(A = 1, B = NA, C = 2)), "missing values")
  # Not sorted decreasing
  expect_error(celldeath_gsea(c(A = 1, B = 3, C = 2)), "decreasing order")
  # Valid input should not error (mock GSEA to avoid real computation)
  local_mocked_bindings(GSEA = function(...) data.frame(), .package = "clusterProfiler")
  geneList_ok <- c(GeneA = 3, GeneB = 2, GeneC = 1)
  expect_no_error(celldeath_gsea(geneList_ok))
})

test_that("celldeath_gsea returns NULL and prints distinct messages when no significant results", {
  local_mocked_bindings(
    GSEA = function(...) data.frame(),
    .package = "clusterProfiler"
  )
  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1)
  # Global mode
  expect_message(
    res_global <- celldeath_gsea(geneList),
    "No significantly enriched cell death pathways detected"
  )
  expect_null(res_global)

  # Targeted single pathway mode has different message text
  expect_message(
    res_target <- celldeath_gsea(geneList, term_select = "Ferroptosis"),
    "No significant GSEA enrichment detected for pathway"
  )
  expect_null(res_target)
})

test_that("celldeath_gsea correctly writes CSV when savefile=TRUE", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  fake_df <- data.frame(ID = "FakePathway", pvalue = 0.001)
  local_mocked_bindings(
    GSEA = function(...) fake_df,
    .package = "clusterProfiler"
  )

  geneList <- c(GPX4 = 2.5, ACSL4 = 2.1)
  result <- celldeath_gsea(geneList, savefile = TRUE, filename = tmp)

  expect_true(file.exists(tmp))
  csv_data <- read.csv(tmp)
  expect_equal(csv_data$ID, "FakePathway")
})

# ===================== Test celldeath_gsea_multiple (multi-group) =====================

test_that("celldeath_gsea_multiple properly calls GSEA for each group (global mode)", {
  geneList_list <- list(
    Group1 = c(GPX4 = 2.5, ACSL4 = 2.1, SLC7A11 = 1.8),
    Group2 = c(GSDMD = 1.5, CASP1 = 1.3, IL1B = 1.1)
  )

  call_count <- 0
  captured_geneLists <- list()

  local_mocked_bindings(
    GSEA = function(geneList, ...) {
      call_count <<- call_count + 1
      captured_geneLists[[call_count]] <<- geneList
      data.frame(ID = paste0("Pathway_", call_count))
    },
    .package = "clusterProfiler"
  )

  result_list <- celldeath_gsea_multiple(geneList_list, pvalueCutoff = 0.05)

  expect_equal(call_count, 2)
  expect_equal(captured_geneLists[[1]], geneList_list$Group1)
  expect_equal(captured_geneLists[[2]], geneList_list$Group2)
  expect_length(result_list, 2)
})

test_that("celldeath_gsea_multiple term_select filters for single pathway", {
  geneList_list <- list(G1 = c(GPX4=2,ACSL4=1), G2 = c(GPX4=1.5))
  captured_t2g <- NULL
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE,...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "Ferroptosis")
    },
    .package = "clusterProfiler"
  )
  celldeath_gsea_multiple(geneList_list, term_select = "Ferroptosis")
  expect_equal(unique(captured_t2g$term), "Ferroptosis")
})

# ---- New: term_select as exact gene set name for multi-group ----
test_that("celldeath_gsea_multiple term_select as exact gene set name is used directly", {
  local_mocked_bindings(
    load_cell_death_genes = function() {
      list(
        Ferroptosis_Up_36747082   = c("GPX4","ACSL4","SLC7A11"),
        Ferroptosis_Down_36747082 = c("FTH1","TFRC"),
        Pyroptosis                = c("GSDMD","CASP1","IL1B")
      )
    },
    .package = "CDEnrich"
  )
  geneList_list <- list(G1 = c(GPX4=2,ACSL4=1.5,SLC7A11=1), G2 = c(GSDMD=2,FTH1=1))
  captured_t2g <- NULL
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "x")
    },
    .package = "clusterProfiler"
  )
  celldeath_gsea_multiple(geneList_list, term_select = "Ferroptosis_Up_36747082")
  expect_equal(unique(captured_t2g$term), "Ferroptosis_Up_36747082")
})

# ---- New: use_recommended controls gene sets for a multi-set death type ----
test_that("celldeath_gsea_multiple use_recommended controls gene sets for a multi-set death type", {
  geneList_list <- list(G1 = c(GPX4=2,ACSL4=1,SLC7A11=0.5), G2 = c(FTH1=2,TFRC=1))
  local_mocked_bindings(
    .resolve_term_select = function(term_select, death_gene_list_all, use_recommended) {
      if (use_recommended) {
        list(Ferroptosis_Reco = c("GPX4", "ACSL4"))
      } else {
        list(Ferroptosis_Reco = c("GPX4", "ACSL4"),
             Ferroptosis_Other = c("SLC7A11", "GSDMD"))
      }
    },
    .package = "CDEnrich"
  )
  captured_t2g <- NULL
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "x")
    },
    .package = "clusterProfiler"
  )

  celldeath_gsea_multiple(geneList_list, term_select = "Ferroptosis")
  expect_equal(unique(captured_t2g$term), "Ferroptosis_Reco")

  celldeath_gsea_multiple(geneList_list, term_select = "Ferroptosis", use_recommended = FALSE)
  expect_true(all(c("Ferroptosis_Reco", "Ferroptosis_Other") %in% captured_t2g$term))
})

# ---- New: .resolve_term_select is delegated the right arguments (multi-group) ----
test_that("celldeath_gsea_multiple delegates term_select/use_recommended to .resolve_term_select", {
  geneList_list <- list(G1 = c(GPX4=2,ACSL4=1))
  resolve_calls <- list()
  local_mocked_bindings(
    .resolve_term_select = function(term_select, death_gene_list_all, use_recommended) {
      resolve_calls <<- append(resolve_calls,
                               list(list(term_select = term_select,
                                         use_recommended = use_recommended,
                                         n_sets = length(death_gene_list_all))))
      list(Resolved = c("GPX4", "ACSL4"))
    },
    .package = "CDEnrich"
  )
  local_mocked_bindings(
    GSEA = function(...) data.frame(ID = "x"),
    .package = "clusterProfiler"
  )
  celldeath_gsea_multiple(geneList_list, term_select = "Ferroptosis", use_recommended = FALSE)
  expect_length(resolve_calls, 1)
  expect_equal(resolve_calls[[1]]$term_select, "Ferroptosis")
  expect_false(resolve_calls[[1]]$use_recommended)
  expect_equal(resolve_calls[[1]]$n_sets, 3)
})

# ---- New: default pvalueCutoff=1 forwarded to GSEA (multi-group) ----
test_that("celldeath_gsea_multiple passes default pvalueCutoff=1 to GSEA", {
  geneList_list <- list(G1 = c(GPX4=2,ACSL4=1), G2 = c(GSDMD=1.5,CASP1=1))
  captured <- list()
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, pvalueCutoff, verbose, pAdjustMethod, ...) {
      captured$pvalueCutoff <<- pvalueCutoff
      data.frame(ID = "x")
    },
    .package = "clusterProfiler"
  )
  celldeath_gsea_multiple(geneList_list)
  expect_equal(captured$pvalueCutoff, 1)
})

# ---- New: verbose=TRUE forwarded to GSEA (multi-group) ----
test_that("celldeath_gsea_multiple forwards verbose=TRUE to GSEA", {
  geneList_list <- list(G1 = c(GPX4=2,ACSL4=1))
  captured <- list()
  local_mocked_bindings(
    GSEA = function(geneList, TERM2GENE, pvalueCutoff, verbose, pAdjustMethod, ...) {
      captured$verbose <<- verbose
      data.frame(ID = "x")
    },
    .package = "clusterProfiler"
  )
  celldeath_gsea_multiple(geneList_list, verbose = TRUE)
  expect_true(captured$verbose)
})

test_that("celldeath_gsea_multiple throws error for non-existing pathway", {
  geneList_list <- list(G1 = c(GPX4=2,ACSL4=1))
  expect_error(
    celldeath_gsea_multiple(geneList_list, term_select = "InvalidPath"),
    "Pathway name does not exist"
  )
})

test_that("celldeath_gsea_multiple input validation", {
  # Not a list
  expect_error(celldeath_gsea_multiple("not a list"), "must be a list")
  # Unnamed list
  expect_error(celldeath_gsea_multiple(list(c(1,2), c(3,4))), "NAMED list")
  # Invalid geneList inside group
  expect_error(
    celldeath_gsea_multiple(list(A = c(1,2,3), B = c(a=1, b=2))),
    "Group A.*named numeric vector"
  )
  expect_error(
    celldeath_gsea_multiple(list(A = c(A=1, B=3, C=2), B = c(D=1))),
    "Group A.*decreasing order"
  )
})

test_that("celldeath_gsea_multiple saves CSV for each group when savefile=TRUE", {
  geneList_list <- list(
    GroupA = c(GPX4 = 2.5, ACSL4 = 2.1),
    GroupB = c(GSDMD = 1.5, CASP1 = 1.3)
  )

  local_mocked_bindings(
    GSEA = function(geneList, ...) {
      data.frame(ID = if (geneList[1] > 2) "PathA" else "PathB")
    },
    .package = "clusterProfiler"
  )

  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  old_wd <- getwd()
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  result <- celldeath_gsea_multiple(geneList_list, savefile = TRUE)

  expect_true(file.exists("GroupA_celldeath_gsea.csv"))
  expect_true(file.exists("GroupB_celldeath_gsea.csv"))
  csv_a <- read.csv("GroupA_celldeath_gsea.csv")
  csv_b <- read.csv("GroupB_celldeath_gsea.csv")
  expect_equal(csv_a$ID, "PathA")
  expect_equal(csv_b$ID, "PathB")
})

test_that("celldeath_gsea_multiple stores empty results and prints message for groups with no enrichment", {
  local_mocked_bindings(
    GSEA = function(...) data.frame(),
    .package = "clusterProfiler"
  )
  # Real genes: overlap exists, but mocked GSEA returns nothing
  geneList_list <- list(G1 = c(GPX4 = 3, ACSL4 = 2), G2 = c(GSDMD = 2, CASP1 = 1))
  result <- celldeath_gsea_multiple(geneList_list)
  expect_length(result, 2)
  expect_equal(nrow(result$G1), 0)
  expect_equal(nrow(result$G2), 0)
})

test_that("celldeath_gsea_multiple skips groups with zero overlapping genes", {
  local_mocked_bindings(
    GSEA = function(...) data.frame(),
    .package = "clusterProfiler"
  )
  geneList_list <- list(G1 = c(Gene1 = 3, Gene2 = 2), G2 = c(GPX4 = 2, ACSL4 = 1))
  expect_message(
    result <- celldeath_gsea_multiple(geneList_list),
    "no genes overlap"
  )
  expect_length(result, 2)      # G1 kept as a NULL placeholder
  expect_null(result$G1)
  expect_equal(nrow(result$G2), 0)
})
