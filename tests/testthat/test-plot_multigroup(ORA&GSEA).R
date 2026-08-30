library(testthat)
library(clusterProfiler)

# ========== Mock plotting functions ==========
local_mocked_bindings(
  gseaplot2 = function(x, geneSetID, title, color, base_size, ...) {
    ggplot2::ggplot()
  },
  gseaplot = function(x, geneSetID, by, color.line, base_size, title, ...) {
    ggplot2::ggplot()
  },
  .package = "enrichplot"
)

# ========== Test data: multi-group ORA comparison ==========
test_compare_df <- data.frame(
  Cluster = rep(c("Control", "Treatment"), each = 3),
  ID = rep(c("Ferroptosis", "Pyroptosis", "Apoptosis"), 2),
  Description = rep(c("Ferroptosis", "Pyroptosis", "Apoptosis"), 2),
  GeneRatio = rep(c("5/100", "4/100", "3/100"), 2),
  Count = rep(c(5, 4, 3), 2),
  p.adjust = c(0.001, 0.01, 0.05, 0.002, 0.03, 0.06),
  stringsAsFactors = FALSE
)

# Test data: targeted single-pathway comparison across groups
test_compare_single <- data.frame(
  Cluster = c("Control", "Treatment"),
  ID = c("Ferroptosis", "Ferroptosis"),
  Description = c("Ferroptosis", "Ferroptosis"),
  GeneRatio = c("5/100", "4/100"),
  Count = c(5, 4),
  p.adjust = c(0.001, 0.01),
  stringsAsFactors = FALSE
)

# ========== Test data: multi-group GSEA results ==========
test_gsea_g1 <- data.frame(
  ID = c("Ferroptosis_FerrDb V2(36305834)", "Pyroptosis_Pan Cancer_36703967"),
  Description = c("Ferroptosis_FerrDb V2(36305834)", "Pyroptosis_Pan Cancer_36703967"),
  NES = c(2.1, -1.8),
  p.adjust = c(0.001, 0.01),
  stringsAsFactors = FALSE
)
test_gsea_g2 <- data.frame(
  ID = c("Ferroptosis_FerrDb V2(36305834)", "Apoptosis_Pan Cancer_37968457"),
  Description = c("Ferroptosis_FerrDb V2(36305834)", "Apoptosis_Pan Cancer_37968457"),
  NES = c(1.5, 2.0),
  p.adjust = c(0.02, 0.005),
  stringsAsFactors = FALSE
)
test_gsea_list <- list(G1 = test_gsea_g1, G2 = test_gsea_g2)

# ---------- Test plot_death_compare_bubble ----------
test_that("plot_death_compare_bubble global mode returns faceted ggplot", {
  p <- plot_death_compare_bubble(test_compare_df, show_category = 2)
  expect_s3_class(p, "ggplot")
  # one row facet per group
  expect_s3_class(p$facet, "FacetGrid")

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))
  p2 <- plot_death_compare_bubble(test_compare_df, filename = tmp)
  expect_true(file.exists(tmp))
})

test_that("plot_death_compare_bubble return_data caps top pathways per group", {
  df_out <- plot_death_compare_bubble(test_compare_df, show_category = 2,
                                      return_data = TRUE)
  expect_true(is.data.frame(df_out))
  # at most 2 pathways per group
  expect_true(all(table(df_out$Cluster) <= 2))
  # sorted by p.adjust within each group
  ctrl <- df_out[df_out$Cluster == "Control", ]
  expect_true(all(diff(ctrl$p.adjust) >= 0))
  # original group order preserved as factor levels
  expect_equal(levels(df_out$Cluster), c("Control", "Treatment"))
})

test_that("plot_death_compare_bubble targeted single-pathway mode outputs INFO message", {
  expect_message(
    p <- plot_death_compare_bubble(test_compare_single),
    "INFO: Current multi-group comparison uses targeted single-pathway mode"
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_death_compare_bubble input error handling", {
  expect_error(plot_death_compare_bubble(NULL),
               "Please run celldeath_compare_enrich() first",
               fixed = TRUE)

  empty_df <- test_compare_df[0, ]
  expect_error(plot_death_compare_bubble(empty_df),
               "Multi-group enrichment result is empty",
               fixed = TRUE)
})

# ---------- Test plot_death_compare_heatmap ----------
test_that("plot_death_compare_heatmap computes -log10(p.adjust) and returns ggplot", {
  p <- plot_death_compare_heatmap(test_compare_df)
  expect_s3_class(p, "ggplot")

  df_out <- plot_death_compare_heatmap(test_compare_df, return_data = TRUE)
  expect_true(is.data.frame(df_out))
  expect_true("logp" %in% colnames(df_out))
  expect_equal(df_out$logp, -log10(df_out$p.adjust))
})

test_that("plot_death_compare_heatmap selects top pathways by lowest p.adjust", {
  df_out <- plot_death_compare_heatmap(test_compare_df, show_category = 2,
                                       return_data = TRUE)
  # only the 2 pathways with the best (cross-group minimum) p.adjust survive
  expect_equal(length(unique(as.character(df_out$Description))), 2)
  expect_true("Ferroptosis" %in% df_out$Description)   # best p.adjust overall
})

test_that("plot_death_compare_heatmap targeted single-pathway mode outputs INFO message", {
  expect_message(
    plot_death_compare_heatmap(test_compare_single),
    "INFO: Current multi-group comparison heatmap uses targeted single-pathway mode"
  )
})

test_that("plot_death_compare_heatmap input error handling and file saving", {
  expect_error(plot_death_compare_heatmap(NULL),
               "Please run celldeath_compare_enrich() first",
               fixed = TRUE)
  empty_df <- test_compare_df[0, ]
  expect_error(plot_death_compare_heatmap(empty_df),
               "Multi-group enrichment result is empty",
               fixed = TRUE)

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))
  plot_death_compare_heatmap(test_compare_df, filename = tmp)
  expect_true(file.exists(tmp))
})

# ---------- Test plot_death_gsea_curve_multiple ----------
test_that("plot_death_gsea_curve_multiple requires a named list", {
  expect_error(
    plot_death_gsea_curve_multiple(list(test_gsea_g1, test_gsea_g2), pathway = "Ferroptosis"),
    "must be a NAMED list"
  )
  expect_error(
    plot_death_gsea_curve_multiple(list(G1 = test_gsea_g1, test_gsea_g2), pathway = "Ferroptosis"),
    "must be a NAMED list"
  )
})

test_that("plot_death_gsea_curve_multiple exact pathway ID plots each group", {
  captured <- list()
  local_mocked_bindings(
    gseaplot2 = function(x, geneSetID, title, color, base_size, ...) {
      captured$geneSetID[[length(captured$geneSetID) + 1]] <<- geneSetID
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  res <- plot_death_gsea_curve_multiple(
    test_gsea_list, pathway = "Ferroptosis_FerrDb V2(36305834)"
  )
  expect_type(res, "list")
  expect_length(res, 2)
  expect_named(res, c("G1", "G2"))
  expect_true(all(sapply(res, function(p) inherits(p, "ggplot"))))
  # each group received the exact pathway ID
  expect_true(all(sapply(captured$geneSetID, function(id) "Ferroptosis_FerrDb V2(36305834)" %in% id)))
})

test_that("plot_death_gsea_curve_multiple death type resolves once to the recommended set", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) {
      list(Pathway = "Ferroptosis_FerrDb V2(36305834)", PMID = "36305834")
    },
    .package = "CDEnrich"
  )
  captured <- list()
  local_mocked_bindings(
    gseaplot2 = function(x, geneSetID, title, color, base_size, ...) {
      captured$geneSetID[[length(captured$geneSetID) + 1]] <<- geneSetID
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  res <- plot_death_gsea_curve_multiple(test_gsea_list, pathway = "Ferroptosis")
  expect_length(res, 2)
  # resolved recommended ID is what each group plots
  expect_true(all(sapply(captured$geneSetID, function(id) "Ferroptosis_FerrDb V2(36305834)" %in% id)))
})

test_that("plot_death_gsea_curve_multiple skips groups lacking the pathway with a warning", {
  expect_warning(
    res <- plot_death_gsea_curve_multiple(
      test_gsea_list, pathway = "Pyroptosis_Pan Cancer_36703967"
    ),
    "does not contain pathway"
  )
  # only G1 contains Pyroptosis
  expect_length(res, 1)
  expect_named(res, "G1")
})

test_that("plot_death_gsea_curve_multiple skips empty group results and errors when all empty", {
  empty_df <- test_gsea_g1[0, ]
  expect_warning(
    res <- plot_death_gsea_curve_multiple(
      list(G1 = empty_df, G2 = test_gsea_g2),
      pathway = "Ferroptosis_FerrDb V2(36305834)"
    ),
    "GSEA result is empty, skip plotting"
  )
  expect_length(res, 1)
  expect_named(res, "G2")

  expect_error(
    plot_death_gsea_curve_multiple(list(G1 = empty_df, G2 = empty_df),
                                   pathway = "Ferroptosis"),
    "All GSEA results are empty"
  )
})

test_that("plot_death_gsea_curve_multiple saves one file per group with filename_prefix", {
  prefix <- tempfile()
  files <- paste0(prefix, "_", c("G1", "G2"), ".png")
  on.exit(unlink(files), add = TRUE)

  res <- plot_death_gsea_curve_multiple(
    test_gsea_list, pathway = "Ferroptosis_FerrDb V2(36305834)",
    filename_prefix = prefix
  )
  expect_length(res, 2)
  expect_true(all(file.exists(files)))
})

# NEW: use_recommended = FALSE is now correctly forwarded to the single-group plot
test_that("plot_death_gsea_curve_multiple use_recommended=FALSE forwards all matched sets present per group", {
  # G1 gets a second Ferroptosis gene set; G2 keeps only the recommended one
  g1_extra <- rbind(
    test_gsea_g1,
    data.frame(
      ID = "Ferroptosis_ALI_37149072",
      Description = "Ferroptosis_ALI_37149072",
      NES = 1.7,
      p.adjust = 0.03,
      stringsAsFactors = FALSE
    )
  )
  local_mocked_bindings(
    get_representative_pathway = function(death_type) {
      list(Pathway = "Ferroptosis_FerrDb V2(36305834)", PMID = "36305834")
    },
    .package = "CDEnrich"
  )
  captured <- list()
  local_mocked_bindings(
    gseaplot2 = function(x, geneSetID, title, color, base_size, ...) {
      captured$geneSetID[[length(captured$geneSetID) + 1]] <<- geneSetID
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  expect_message(
    res <- plot_death_gsea_curve_multiple(
      list(G1 = g1_extra, G2 = test_gsea_g2),
      pathway = "Ferroptosis",
      use_recommended = FALSE
    ),
    "running against ALL"
  )
  expect_length(res, 2)
  # G1 contains both Ferroptosis gene sets -> both are plotted
  expect_equal(sort(captured$geneSetID[[1]]),
               c("Ferroptosis_ALI_37149072", "Ferroptosis_FerrDb V2(36305834)"))
  # G2 contains only the recommended one -> only it is plotted
  expect_equal(captured$geneSetID[[2]], "Ferroptosis_FerrDb V2(36305834)")
})

# ---------- Test plot_death_gsea_heatmap_multiple ----------
test_that("plot_death_gsea_heatmap_multiple NES mode returns ggplot and long data", {
  p <- plot_death_gsea_heatmap_multiple(test_gsea_list, statistic = "NES")
  expect_s3_class(p, "ggplot")

  df_out <- plot_death_gsea_heatmap_multiple(test_gsea_list, statistic = "NES",
                                             return_data = TRUE)
  expect_true(is.data.frame(df_out))
  expect_true(all(c("ID", "Group", "Value", "Description") %in% colnames(df_out)))
  # 3 unique pathways x 2 groups = 6 tiles (missing pathway kept as NA tile)
  expect_equal(nrow(df_out), 6)
  # Pyroptosis is absent from G2 -> NA value (grey tile)
  pyro_g2 <- df_out[df_out$ID == "Pyroptosis_Pan Cancer_36703967" & df_out$Group == "G2", ]
  expect_true(is.na(pyro_g2$Value))
})

test_that("plot_death_gsea_heatmap_multiple p.adjust mode applies -log10 transform", {
  df_out <- plot_death_gsea_heatmap_multiple(test_gsea_list, statistic = "p.adjust",
                                             return_data = TRUE)
  ferro_g1 <- df_out[df_out$ID == "Ferroptosis_FerrDb V2(36305834)" & df_out$Group == "G1", ]
  expect_equal(ferro_g1$Value, -log10(0.001))
})

test_that("plot_death_gsea_heatmap_multiple targeted single-pathway mode outputs INFO message", {
  single_list <- list(
    G1 = test_gsea_g1[1, ],
    G2 = test_gsea_g2[1, ]
  )
  expect_message(
    plot_death_gsea_heatmap_multiple(single_list),
    "contains only one pathway"
  )
})

test_that("plot_death_gsea_heatmap_multiple skips empty groups and validates input", {
  empty_df <- test_gsea_g1[0, ]
  expect_warning(
    plot_death_gsea_heatmap_multiple(list(G1 = empty_df, G2 = test_gsea_g2)),
    "has no result, skipped"
  )
  expect_warning(
    expect_error(
      plot_death_gsea_heatmap_multiple(list(G1 = empty_df)),
      "No valid GSEA results"
    ),
    "has no result, skipped"
  )
  expect_error(
    plot_death_gsea_heatmap_multiple(list(test_gsea_g1)),
    "must be a NAMED list"
  )
})

test_that("plot_death_gsea_heatmap_multiple saves plot to file", {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))
  p <- plot_death_gsea_heatmap_multiple(test_gsea_list, filename = tmp)
  expect_true(file.exists(tmp))
  expect_s3_class(p, "ggplot")
})
