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

# ===================== Test celldeath_enrich =====================

test_that("celldeath_enrich passes correct parameters to enricher", {
  deg <- c("GPX4", "ACSL4", "CASP3")
  captured_args <- list()

  local_mocked_bindings(
    enricher = function(gene, universe, TERM2GENE, pvalueCutoff, qvalueCutoff, minGSSize, pAdjustMethod) {
      captured_args <<- as.list(environment())
      # Return a data frame directly, real environment requires enrichResult object
      data.frame(ID = "test", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )

  result <- celldeath_enrich(deg, pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 3)

  expect_equal(captured_args$gene, deg)
  expect_true(is.data.frame(captured_args$TERM2GENE))
  expect_equal(captured_args$pvalueCutoff, 0.05)
  expect_equal(captured_args$qvalueCutoff, 0.2)
  expect_equal(captured_args$minGSSize, 3)
  expect_equal(captured_args$pAdjustMethod, "BH")
  # Default universe is the full curated collection (all 9 genes)
  expect_true(all(c("GPX4", "ACSL4", "SLC7A11", "GSDMD", "CASP1",
                    "IL1B", "CASP3", "BAX", "BCL2") %in% captured_args$universe))
})

test_that("celldeath_enrich term_select specifies single valid pathway", {
  deg <- c("GPX4", "ACSL4", "CASP3")
  captured_t2g <- NULL
  local_mocked_bindings(
    enricher = function(gene, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "Ferroptosis", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  res <- celldeath_enrich(deg, term_select = "Ferroptosis")
  # TERM2GENE should only contain Ferroptosis
  expect_equal(unique(captured_t2g$Term), "Ferroptosis")
})

# ---- term_select as an EXACT gene set name (used directly) ----
test_that("celldeath_enrich term_select as exact gene set name is used directly", {
  # Use a collection where a death type has multiple gene sets, so that
  # the exact-name path can be distinguished from the death-type path.
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
  deg <- c("GPX4", "ACSL4", "SLC7A11")
  captured_t2g <- NULL
  local_mocked_bindings(
    enricher = function(gene, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  celldeath_enrich(deg, term_select = "Ferroptosis_Up_36747082")
  expect_equal(unique(captured_t2g$Term), "Ferroptosis_Up_36747082")
})

# ---- use_recommended controls which gene sets are used for a
#      multi-gene-set death type (integration through .resolve_term_select) ----
test_that("celldeath_enrich use_recommended controls gene sets for a multi-set death type", {
  local_mocked_bindings(
    load_cell_death_genes = function() {
      list(
        Ferroptosis_Up_36747082   = c("GPX4", "ACSL4", "SLC7A11"),
        Ferroptosis_Down_36747082 = c("FTH1", "TFRC"),
        Pyroptosis                = c("GSDMD", "CASP1", "IL1B")
      )
    },
    get_representative_pathway = function(death_type) {
      list(Pathway = "Ferroptosis_Up_36747082")
    },
    .package = "CDEnrich"
  )
  deg <- c("GPX4", "ACSL4", "SLC7A11", "FTH1", "TFRC", "GSDMD", "CASP1", "IL1B")
  captured_t2g <- NULL
  local_mocked_bindings(
    enricher = function(gene, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )

  # use_recommended = TRUE (default): only the recommended representative gene set
  celldeath_enrich(deg, term_select = "Ferroptosis", use_recommended = TRUE)
  expect_equal(sort(unique(captured_t2g$Term)), "Ferroptosis_Up_36747082")

  # use_recommended = FALSE: ALL gene sets under Ferroptosis
  celldeath_enrich(deg, term_select = "Ferroptosis", use_recommended = FALSE)
  expect_equal(sort(unique(captured_t2g$Term)),
               c("Ferroptosis_Down_36747082", "Ferroptosis_Up_36747082"))
})

# ---- .resolve_term_select is delegated the right arguments ----
test_that("celldeath_enrich delegates term_select/use_recommended to .resolve_term_select", {
  deg <- c("GPX4", "ACSL4", "SLC7A11")
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
    enricher = function(...) data.frame(ID = "x", pvalue = 0.01),
    .package = "clusterProfiler"
  )
  celldeath_enrich(deg, term_select = "Ferroptosis", use_recommended = FALSE)
  expect_length(resolve_calls, 1)
  expect_equal(resolve_calls[[1]]$term_select, "Ferroptosis")
  expect_false(resolve_calls[[1]]$use_recommended)
  # The full curated collection is handed to the resolver (3 sets in default mock)
  expect_equal(resolve_calls[[1]]$n_sets, 3)
})

# ---- universe defaults to the full collection even in targeted mode,
#      and can be overridden by the user ----
test_that("celldeath_enrich keeps full curated collection as universe in targeted mode", {
  deg <- c("GPX4", "ACSL4", "SLC7A11")
  captured <- list()
  local_mocked_bindings(
    enricher = function(gene, universe, TERM2GENE, ...) {
      captured$universe <<- universe
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  celldeath_enrich(deg, term_select = "Ferroptosis")   # targeted mode
  # universe must include genes from OTHER pathways too (computed before filtering)
  expect_true(all(c("GSDMD", "CASP1", "IL1B", "CASP3", "BAX", "BCL2") %in% captured$universe))
})

test_that("celldeath_enrich honours a user-supplied universe", {
  deg <- c("GPX4", "ACSL4", "SLC7A11")
  captured <- list()
  local_mocked_bindings(
    enricher = function(gene, universe, TERM2GENE, ...) {
      captured$universe <<- universe
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  celldeath_enrich(deg, universe = c("GPX4", "ACSL4", "SLC7A11", "TP53"))
  expect_equal(sort(captured$universe), c("ACSL4", "GPX4", "SLC7A11", "TP53"))
})

# ---- NEW: enrichment_force_universe is TRUE while enricher runs and restored after ----
test_that("celldeath_enrich forces full universe inside enricher and restores the option", {
  deg <- c("GPX4", "ACSL4", "SLC7A11")
  captured <- list()
  local_mocked_bindings(
    enricher = function(...) {
      captured$force <<- getOption("enrichment_force_universe", FALSE)
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )

  old_value <- getOption("enrichment_force_universe", FALSE)
  celldeath_enrich(deg, term_select = "Ferroptosis")

  # Option is TRUE while enricher runs (prevents universe shrinking in targeted mode)
  expect_true(captured$force)
  # Option is restored to its previous value after the call
  expect_identical(getOption("enrichment_force_universe", FALSE), old_value)
})

test_that("celldeath_enrich throws error for non-existing pathway name via term_select", {
  deg <- c("GPX4", "ACSL4", "CASP3")
  expect_error(
    celldeath_enrich(deg, term_select = "FakePathway"),
    "Pathway name does not exist"
  )
})

test_that("celldeath_enrich performs input validation", {
  expect_error(celldeath_enrich(123), "must be a character vector")
  expect_error(celldeath_enrich(c("A", "B")), "cannot be less than 3")
  expect_error(celldeath_enrich(c("", "", "GPX4")), "cannot be less than 3")
})

test_that("celldeath_enrich returns NULL and shows message when no significant results", {
  local_mocked_bindings(
    enricher = function(...) NULL,
    .package = "clusterProfiler"
  )
  expect_message(
    result <- celldeath_enrich(c("GPX4", "ACSL4", "SLC7A11")),
    "No significantly enriched cell death pathways detected"
  )
  expect_null(result)
})

test_that("celldeath_enrich correctly writes CSV when savefile=TRUE", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  fake_df <- data.frame(ID = "FakePathway", pvalue = 0.001)
  local_mocked_bindings(
    enricher = function(...) fake_df,
    .package = "clusterProfiler"
  )

  result <- celldeath_enrich(c("GPX4", "ACSL4", "SLC7A11"),
                             savefile = TRUE,
                             filename = tmp)

  expect_true(file.exists(tmp))
  csv_data <- read.csv(tmp)
  expect_equal(csv_data$ID, "FakePathway")
})

# ===================== Test celldeath_compare_enrich =====================

test_that("celldeath_compare_enrich properly calls compareCluster", {
  deg_list <- list(GroupA = c("GPX4","ACSL4","CASP3"),
                   GroupB = c("GSDMD","IL1B","BAX"))
  captured_fun <- NULL

  local_mocked_bindings(
    compareCluster = function(formula, data, fun, ...) {
      captured_fun <<- fun
      # Return data frame directly
      data.frame(Cluster = "A", ID = "B")
    },
    enricher = function(...) data.frame(ID = "test", pvalue = 0.01),
    .package = "clusterProfiler"
  )

  result <- celldeath_compare_enrich(deg_list)

  expect_type(captured_fun, "closure")
  test_res <- captured_fun(c("GPX4", "ACSL4"))
  expect_true(is.data.frame(test_res))
  expect_true(is.data.frame(result))
})

# ---- compare passes the documented default cutoffs (1, 1) into enricher ----
test_that("celldeath_compare_enrich passes default pvalueCutoff=1/qvalueCutoff=1 to enricher", {
  deg_list <- list(G1 = c("GPX4","ACSL4","SLC7A11"),
                   G2 = c("GSDMD","CASP1","IL1B"))
  captured <- list()
  local_mocked_bindings(
    compareCluster = function(formula, data, fun, ...) {
      fun(c("GPX4", "ACSL4"))   # exercise the inner enricher
      data.frame(Cluster = "G1", ID = "x")
    },
    enricher = function(gene, universe, TERM2GENE, pvalueCutoff, qvalueCutoff, minGSSize, pAdjustMethod) {
      captured$args <<- as.list(environment())
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  celldeath_compare_enrich(deg_list)
  expect_equal(captured$args$pvalueCutoff, 1)
  expect_equal(captured$args$qvalueCutoff, 1)
  expect_equal(captured$args$pAdjustMethod, "BH")
})

test_that("celldeath_compare_enrich term_select specifies single pathway", {
  deg_list <- list(G1 = c("GPX4","ACSL4","SLC7A11"), G2 = c("GPX4","ACSL4","CASP3"))
  captured_t2g <- NULL

  local_mocked_bindings(
    compareCluster = function(formula, data, fun, ...) {
      # Capture TERM2GENE inside mocked enricher
      captured_fun_inner <- fun
      captured_fun_inner(c("GPX4","ACSL4","SLC7A11"))
      data.frame(Cluster = "G1", ID = "Ferroptosis")
    },
    enricher = function(gene, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "Ferroptosis", pvalue = 0.02)
    },
    .package = "clusterProfiler"
  )

  res <- celldeath_compare_enrich(deg_list, term_select = "Ferroptosis")
  expect_equal(unique(captured_t2g$Term), "Ferroptosis")
})

# ---- compare term_select as exact gene set name ----
test_that("celldeath_compare_enrich term_select as exact gene set name is used directly", {
  deg_list <- list(G1 = c("GPX4","ACSL4","SLC7A11"), G2 = c("GSDMD","CASP1","IL1B"))
  captured_t2g <- NULL
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
  local_mocked_bindings(
    compareCluster = function(formula, data, fun, ...) {
      fun(c("GPX4","ACSL4","SLC7A11"))
      data.frame(Cluster = "G1", ID = "x")
    },
    enricher = function(gene, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  celldeath_compare_enrich(deg_list, term_select = "Ferroptosis_Up_36747082")
  expect_equal(unique(captured_t2g$Term), "Ferroptosis_Up_36747082")
})

# ---- compare use_recommended controls gene sets for a multi-set death type ----
test_that("celldeath_compare_enrich use_recommended controls gene sets for a multi-set death type", {
  deg_list <- list(G1 = c("GPX4","ACSL4","SLC7A11"), G2 = c("GSDMD","CASP1","IL1B"))
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
    compareCluster = function(formula, data, fun, ...) {
      fun(c("GPX4", "ACSL4", "SLC7A11"))
      fun(c("GSDMD","CASP1","IL1B"))
      data.frame(Cluster = c("G1","G2"), ID = c("x","y"))
    },
    enricher = function(gene, TERM2GENE, ...) {
      captured_t2g <<- TERM2GENE
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )

  # use_recommended = TRUE
  celldeath_compare_enrich(deg_list, term_select = "Ferroptosis")
  expect_equal(sort(unique(captured_t2g$Term)), "Ferroptosis_Reco")

  # use_recommended = FALSE
  celldeath_compare_enrich(deg_list, term_select = "Ferroptosis", use_recommended = FALSE)
  expect_equal(sort(unique(captured_t2g$Term)),
               sort(c("Ferroptosis_Reco", "Ferroptosis_Other")))
})

# ---- compare keeps full curated collection as universe in targeted mode ----
test_that("celldeath_compare_enrich keeps full curated collection as universe in targeted mode", {
  deg_list <- list(G1 = c("GPX4","ACSL4","SLC7A11"), G2 = c("GSDMD","CASP1","IL1B"))
  captured <- list()
  local_mocked_bindings(
    compareCluster = function(formula, data, fun, ...) {
      fun(c("GPX4", "ACSL4", "SLC7A11"))
      data.frame(Cluster = "G1", ID = "x")
    },
    enricher = function(gene, universe, TERM2GENE, ...) {
      captured$universe <<- universe
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )
  celldeath_compare_enrich(deg_list, term_select = "Ferroptosis")
  expect_true(all(c("GSDMD", "CASP1", "IL1B") %in% captured$universe))
})

# ---- NEW: enrichment_force_universe is TRUE while the inner enricher runs and restored after ----
test_that("celldeath_compare_enrich forces full universe inside enricher and restores the option", {
  deg_list <- list(G1 = c("GPX4","ACSL4","SLC7A11"), G2 = c("GSDMD","CASP1","IL1B"))
  captured <- list()
  local_mocked_bindings(
    compareCluster = function(formula, data, fun, ...) {
      fun(c("GPX4", "ACSL4", "SLC7A11"))   # exercise the inner enricher
      data.frame(Cluster = "G1", ID = "x")
    },
    enricher = function(...) {
      captured$force <<- getOption("enrichment_force_universe", FALSE)
      data.frame(ID = "x", pvalue = 0.01)
    },
    .package = "clusterProfiler"
  )

  old_value <- getOption("enrichment_force_universe", FALSE)
  celldeath_compare_enrich(deg_list, term_select = "Ferroptosis")

  # Option is TRUE while enricher runs (prevents universe shrinking in targeted mode)
  expect_true(captured$force)
  # Option is restored to its previous value after the call
  expect_identical(getOption("enrichment_force_universe", FALSE), old_value)
})

test_that("celldeath_compare_enrich throws error for invalid pathway name", {
  deg_list <- list(G1 = c("GPX4","ACSL4","SLC7A11"))
  expect_error(
    celldeath_compare_enrich(deg_list, term_select = "NotExistPathway"),
    "Pathway name does not exist"
  )
})

test_that("celldeath_compare_enrich shows message for targeted mode with no significant results", {
  deg_list <- list(G1 = c("GPX4","ACSL4","SLC7A11"))
  local_mocked_bindings(
    compareCluster = function(...) NULL,
    .package = "clusterProfiler"
  )
  expect_message(
    res <- celldeath_compare_enrich(deg_list, term_select = "Ferroptosis"),
    "No significant enrichment detected for pathway"
  )
  expect_null(res)
})

test_that("celldeath_compare_enrich performs input validation", {
  expect_error(celldeath_compare_enrich("not list"), "must be a list")
  expect_error(celldeath_compare_enrich(list(c("A","B"))), "must be named")
  expect_error(
    celldeath_compare_enrich(list(G1 = c("A","B"), G2 = c("C","D","E"))),
    "cannot be less than 3"
  )
})

test_that("celldeath_compare_enrich returns NULL for global mode with no significant results", {
  local_mocked_bindings(
    compareCluster = function(...) NULL,
    .package = "clusterProfiler"
  )
  deg_list <- list(G1 = c("GPX4","ACSL4","SLC7A11"), G2 = c("BAX","CASP3","GSDMD"))
  expect_message(
    result <- celldeath_compare_enrich(deg_list),
    "No significantly enriched cell death pathways across multiple groups detected"
  )
  expect_null(result)
})

# ---- compare writes CSV when savefile=TRUE ----
test_that("celldeath_compare_enrich correctly writes CSV when savefile=TRUE", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  local_mocked_bindings(
    compareCluster = function(...) data.frame(Cluster = "G1", ID = "Ferroptosis", pvalue = 0.01),
    .package = "clusterProfiler"
  )

  result <- celldeath_compare_enrich(
    list(G1 = c("GPX4","ACSL4","SLC7A11"), G2 = c("GSDMD","CASP1","IL1B")),
    savefile = TRUE,
    filename = tmp
  )

  expect_true(file.exists(tmp))
  csv_data <- read.csv(tmp)
  expect_true("Ferroptosis" %in% csv_data$ID)
})
