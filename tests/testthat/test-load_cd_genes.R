library(testthat)

# ===================== Standard behavior =====================

test_that("load_cell_death_genes loads built-in gene sets and returns a standard named list", {
  result <- load_cell_death_genes()

  # Return type checks
  expect_type(result, "list")
  expect_true(!is.null(names(result)))

  # Each pathway is a character vector of genes
  expect_true(all(sapply(result, is.character)))

  # Filtering rule: only pathways with >= 3 genes are kept
  expect_true(all(sapply(result, length) >= 3))
})

test_that("load_cell_death_genes automatically trims leading and trailing whitespace in pathway names and gene symbols", {
  test_df <- data.frame(
    term = c(" Ferroptosis ", "Ferroptosis", "Pyroptosis "),
    gene = c(" GPX4 ", "ACSL4", " GSDMD "),
    stringsAsFactors = FALSE
  )
  res <- load_cell_death_genes(.test_df = test_df)

  expect_false(any(grepl("^\\s|\\s$", names(res))))
  expect_false(any(grepl("^\\s|\\s$", unlist(res))))
})

test_that("load_cell_death_genes automatically filters rows containing NA or empty strings", {
  test_df <- data.frame(
    term = c("Ferroptosis", NA, "", "Pyroptosis"),
    gene = c("", "GPX4", "ACSL4", NA),
    stringsAsFactors = FALSE
  )
  res <- load_cell_death_genes(.test_df = test_df)

  expect_length(res, 0)
})

test_that("load_cell_death_genes automatically removes pathways with fewer than 3 genes", {
  test_df <- data.frame(
    term = c("PathShort", "PathShort", "PathOK", "PathOK", "PathOK"),
    gene = c("G1", "G2", "G3", "G4", "G5"),
    stringsAsFactors = FALSE
  )
  res <- load_cell_death_genes(.test_df = test_df)

  expect_named(res, "PathOK")
  expect_equal(res$PathOK, c("G3", "G4", "G5"))
})

# ===================== Invalid input =====================

test_that("load_cell_death_genes throws an error when required columns are missing from the dataset", {
  # Missing term column
  bad_df1 <- data.frame(gene = c("A", "B"), other = 1:2, stringsAsFactors = FALSE)
  expect_error(load_cell_death_genes(.test_df = bad_df1), "must contain two columns")

  # Missing gene column
  bad_df2 <- data.frame(term = c("A", "B"), other = 1:2, stringsAsFactors = FALSE)
  expect_error(load_cell_death_genes(.test_df = bad_df2), "must contain two columns")
})

# ===================== Merging of recommended representative gene sets =====================

test_that("load_cell_death_genes merges recommended representative gene sets", {
  test_df <- data.frame(term = rep("BuiltInPath", 3),
                        gene = c("B1", "B2", "B3"),
                        stringsAsFactors = FALSE)
  rep_sets <- list(RepNewPath = c("R1", "R2", "R3"))

  res <- load_cell_death_genes(.test_df = test_df, .test_rep = rep_sets)
  expect_named(res, c("BuiltInPath", "RepNewPath"))
  expect_equal(res$RepNewPath, c("R1", "R2", "R3"))
})

test_that("load_cell_death_genes never overwrites existing pathways when merging", {
  test_df <- data.frame(term = rep("SharedPath", 3),
                        gene = c("ORIG1", "ORIG2", "ORIG3"),
                        stringsAsFactors = FALSE)
  rep_sets <- list(SharedPath = c("FAKE1", "FAKE2", "FAKE3"),
                   NewRepPath = c("N1", "N2", "N3"))

  res <- load_cell_death_genes(.test_df = test_df, .test_rep = rep_sets)
  # Existing pathway keeps its original genes (not overwritten)
  expect_equal(res$SharedPath, c("ORIG1", "ORIG2", "ORIG3"))
  # Brand-new representative set is merged in
  expect_equal(res$NewRepPath, c("N1", "N2", "N3"))
})

test_that("load_cell_death_genes applies the same >= 3 genes filter to representative sets", {
  test_df <- data.frame(term = rep("BuiltInPath", 3),
                        gene = c("B1", "B2", "B3"),
                        stringsAsFactors = FALSE)
  rep_sets <- list(TinyRep = c("T1", "T2"),
                   OkRep   = c("O1", "O2", "O3"))

  res <- load_cell_death_genes(.test_df = test_df, .test_rep = rep_sets)
  expect_false("TinyRep" %in% names(res))
  expect_true("OkRep" %in% names(res))
})

test_that("load_cell_death_genes test mode (.test_df alone) skips the merge entirely", {
  test_df <- data.frame(term = rep("OnlyPath", 3),
                        gene = c("G1", "G2", "G3"),
                        stringsAsFactors = FALSE)
  res <- load_cell_death_genes(.test_df = test_df)
  expect_named(res, "OnlyPath")
})

test_that("load_cell_death_genes built-in result contains all recommended representative sets", {
  # Integration check against the real package data (no injection)
  result <- load_cell_death_genes()
  expected <- names(representative_genes)[sapply(representative_genes, length) >= 3]
  expect_true(all(expected %in% names(result)))
  expect_true(all(representative_pathways$Pathway %in% names(result)))
})
