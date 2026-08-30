library(testthat)
library(clusterProfiler)

# Mock enrichplot functions to avoid real plotting and computation
local_mocked_bindings(
  pairwise_termsim = function(x, ...) {
    structure(list(), class = "termsimResult")
  },
  emapplot = function(x, ...) {
    ggplot2::ggplot()
  },
  gseaplot2 = function(gsea_res, geneSetID, title, color, base_size, ...) {
    ggplot2::ggplot()
  },
  gseaplot = function(gsea_res, geneSetID, by, color.line, base_size, title, ...) {
    ggplot2::ggplot()
  },
  .package = "enrichplot"
)

# Test dataset: multiple pathways (global mode)
test_df_multi <- data.frame(
  ID = c("Ferroptosis", "Pyroptosis", "Apoptosis"),
  Description = c("Ferroptosis", "Pyroptosis", "Apoptosis"),
  GeneRatio = c("5/100", "4/100", "3/100"),
  Count = c(5, 4, 3),
  p.adjust = c(0.001, 0.01, 0.05),
  stringsAsFactors = FALSE
)

# Test dataset: targeted single-pathway mode
test_df_single <- data.frame(
  ID = "Ferroptosis",
  Description = "Ferroptosis",
  GeneRatio = "5/100",
  BgRatio = "486/14739",
  Count = 5,
  p.adjust = 0.001,
  stringsAsFactors = FALSE
)

# Construct valid fake enrichResult S4 object (simulate output from celldeath_enrich)
build_fake_enrichResult <- function(df){
  new("enrichResult",
      result = df,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      organism = "unknown",
      ontology = "custom",
      gene = c("GPX4","ACSL4","SLC7A11"),
      geneSets = list(Ferroptosis = c("GPX4","ACSL4","SLC7A11"))
  )
}
fake_enrich_multi  <- build_fake_enrichResult(test_df_multi)
fake_enrich_single <- build_fake_enrichResult(test_df_single)

# Test datasets for GSEA plot functions (result IDs use real gene set naming style)
test_gsea_multi <- data.frame(
  ID = c("Ferroptosis_FerrDb V2(36305834)", "Ferroptosis_ALI_37149072", "Apoptosis_Pan Cancer_37968457"),
  Description = c("Ferroptosis_FerrDb V2(36305834)", "Ferroptosis_ALI_37149072", "Apoptosis_Pan Cancer_37968457"),
  NES = c(2.1, 1.8, -1.5),
  p.adjust = c(0.001, 0.01, 0.05),
  stringsAsFactors = FALSE
)

test_gsea_single <- data.frame(
  ID = "Ferroptosis_FerrDb V2(36305834)",
  Description = "Ferroptosis_FerrDb V2(36305834)",
  NES = 2.1,
  p.adjust = 0.001,
  stringsAsFactors = FALSE
)

# ---------- Test plot_death_enrich_targeted ----------
test_that("plot_death_enrich_targeted returns compact plot and plotting data", {
  p <- plot_death_enrich_targeted(test_df_single)
  expect_s3_class(p, "ggplot")

  df_out <- plot_death_enrich_targeted(
    test_df_single,
    return_data = TRUE
  )
  expect_true(is.data.frame(df_out))
  expect_true("Label" %in% colnames(df_out))

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))

  p_save <- plot_death_enrich_targeted(
    test_df_single,
    filename = tmp
  )

  expect_true(file.exists(tmp))
  expect_s3_class(p_save, "ggplot")
})

test_that("plot_death_enrich_targeted validates single-row ORA input", {
  expect_error(
    plot_death_enrich_targeted(test_df_multi),
    "requires exactly one enrichment result"
  )

  incomplete_df <- test_df_single[
    ,
    setdiff(colnames(test_df_single), "BgRatio")
  ]

  expect_error(
    plot_death_enrich_targeted(incomplete_df),
    "Missing required columns: BgRatio"
  )
})

# ---------- Test plot_death_enrich_bubble ----------
test_that("plot_death_enrich_bubble global mode with multiple pathways", {
  p <- plot_death_enrich_bubble(test_df_multi, show_category = 2)
  expect_s3_class(p, "ggplot")

  df_out <- plot_death_enrich_bubble(test_df_multi, show_category = 2, return_data = TRUE)
  expect_true(is.data.frame(df_out))
  expect_equal(nrow(df_out), 2)

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))
  p2 <- plot_death_enrich_bubble(test_df_multi, filename = tmp)
  expect_true(file.exists(tmp))
  expect_s3_class(p2, "ggplot")
})

test_that("plot_death_enrich_bubble targeted branch delegates to the targeted summary plot", {
  # Input data frame
  expect_message(
    p <- plot_death_enrich_bubble(test_df_single),
    "Targeted single-pathway result detected"
  )
  expect_s3_class(p, "ggplot")

  # Input raw table extracted from an enrichResult object
  expect_message(
    p_obj <- plot_death_enrich_bubble(fake_enrich_single@result),
    "Targeted single-pathway result detected"
  )
  expect_s3_class(p_obj, "ggplot")
})

test_that("plot_death_enrich_bubble input error handling", {
  expect_error(plot_death_enrich_bubble(NULL),
               "Please run celldeath_enrich() first to obtain valid enrichment results!",
               fixed = TRUE
  )
  empty_df <- test_df_multi[0, ]
  expect_error(plot_death_enrich_bubble(empty_df), "Enrichment result is empty")
})

# ---------- Test plot_death_enrich_network ----------
test_that("plot_death_enrich_network global mode with multiple pathways", {
  p <- plot_death_enrich_network(test_df_multi, show_category = 2)
  expect_s3_class(p, "ggplot")

  expect_warning(
    plot_death_enrich_network(test_df_multi, color = "no_such_col"),
    "Parameter 'color' 'no_such_col' is not present in result columns"
  )

  # Additional test: file saving
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))
  p_save <- plot_death_enrich_network(test_df_multi, filename = tmp)
  expect_true(file.exists(tmp))
})

test_that("plot_death_enrich_network targeted single-pathway branch outputs message and switches plot type", {
  expect_message(
    p <- plot_death_enrich_network(test_df_single),
    "Targeted single-pathway result detected"
  )
  expect_s3_class(p, "ggplot")

  # Do not directly pass fake_enrich_single S4 object, extract internal raw table
  expect_message(
    p_obj <- plot_death_enrich_network(fake_enrich_single@result),
    "Targeted single-pathway result detected"
  )
  expect_s3_class(p_obj, "ggplot")
})

# ---------- Test plot_death_volcano ----------
test_that("plot_death_volcano multi-pathway mode supports data frame and enrichResult object", {
  p <- plot_death_volcano(test_df_multi, label_top = 2)
  expect_s3_class(p, "ggplot")

  p_obj <- plot_death_volcano(fake_enrich_multi@result, label_top = 2)
  expect_s3_class(p_obj, "ggplot")

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  p2 <- plot_death_volcano(test_df_multi, filename = tmp)
  expect_true(file.exists(tmp))
})

test_that("plot_death_volcano targeted single-pathway ORA mode delegates to the targeted summary plot", {
  expect_message(
    p <- plot_death_volcano(test_df_single),
    "Targeted single-pathway ORA result detected"
  )
  expect_s3_class(p, "ggplot")

  expect_message(
    p_obj <- plot_death_volcano(fake_enrich_single@result),
    "Targeted single-pathway ORA result detected"
  )
  expect_s3_class(p_obj, "ggplot")
})

test_that("plot_death_volcano targeted single-pathway GSEA keeps a one-point volcano with INFO message", {
  # Single-pathway GSEA input: no GeneRatio column, so the function must NOT
  # delegate to the ORA targeted plot; it draws a one-point volcano instead.
  expect_message(
    p <- plot_death_volcano(test_gsea_single),
    "Targeted single-pathway GSEA result detected"
  )
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 1)
  expect_true("NES" %in% colnames(p$data))
})

# ---------- Test plot_death_gsea_curve ----------
test_that("plot_death_gsea_curve exact pathway ID is used directly", {
  captured <- list()
  local_mocked_bindings(
    gseaplot2 = function(gsea_res, geneSetID, title, color, base_size, ...) {
      captured$geneSetID <<- geneSetID
      captured$title <<- title
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  p <- plot_death_gsea_curve(test_gsea_multi, pathway = "Apoptosis_Pan Cancer_37968457")
  expect_s3_class(p, "ggplot")
  expect_equal(captured$geneSetID, "Apoptosis_Pan Cancer_37968457")
  expect_equal(captured$title, "Apoptosis_Pan Cancer_37968457")
})

test_that("plot_death_gsea_curve death type resolves to the recommended gene set", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) {
      list(Pathway = "Ferroptosis_FerrDb V2(36305834)", PMID = "36305834")
    },
    .package = "CDEnrich"
  )
  captured <- list()
  local_mocked_bindings(
    gseaplot2 = function(gsea_res, geneSetID, title, color, base_size, ...) {
      captured$geneSetID <<- geneSetID
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  expect_message(
    p <- plot_death_gsea_curve(test_gsea_multi, pathway = "Ferroptosis"),
    "recommended representative gene set"
  )
  expect_s3_class(p, "ggplot")
  expect_equal(captured$geneSetID, "Ferroptosis_FerrDb V2(36305834)")
})

test_that("plot_death_gsea_curve use_recommended=FALSE plots all gene sets of the type", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) {
      list(Pathway = "Ferroptosis_FerrDb V2(36305834)", PMID = "36305834")
    },
    .package = "CDEnrich"
  )
  captured <- list()
  local_mocked_bindings(
    gseaplot2 = function(gsea_res, geneSetID, title, color, base_size, ...) {
      captured$geneSetID <<- geneSetID
      captured$title <<- title
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  expect_message(
    p <- plot_death_gsea_curve(test_gsea_multi, pathway = "Ferroptosis",
                               use_recommended = FALSE),
    "running against ALL"
  )
  expect_s3_class(p, "ggplot")
  # both Ferroptosis gene sets present in the result are plotted
  expect_equal(sort(captured$geneSetID),
               c("Ferroptosis_ALI_37149072", "Ferroptosis_FerrDb V2(36305834)"))
  # multi-pathway title joins IDs
  expect_true(grepl("vs", captured$title))
})

test_that("plot_death_gsea_curve handles invalid pathway and empty result", {
  expect_error(
    plot_death_gsea_curve(test_gsea_multi, pathway = "NotExistPathway"),
    "Pathway name does not exist"
  )
  empty_df <- test_gsea_multi[0, ]
  expect_error(
    plot_death_gsea_curve(empty_df, pathway = "Ferroptosis"),
    "GSEA result is empty"
  )
  expect_error(
    plot_death_gsea_curve(NULL, pathway = "Ferroptosis"),
    "GSEA result is empty"
  )
})

test_that("plot_death_gsea_curve delegates pathway resolution to .resolve_plot_pathway", {
  resolve_calls <- list()
  local_mocked_bindings(
    .resolve_plot_pathway = function(pathway, result_ids, use_recommended) {
      resolve_calls <<- append(resolve_calls,
                               list(list(pathway = pathway,
                                         result_ids = result_ids,
                                         use_recommended = use_recommended)))
      "Resolved_Pathway"
    },
    .package = "CDEnrich"
  )
  captured <- list()
  local_mocked_bindings(
    gseaplot2 = function(gsea_res, geneSetID, title, color, base_size, ...) {
      captured$geneSetID <<- geneSetID
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  plot_death_gsea_curve(test_gsea_multi, pathway = "Ferroptosis", use_recommended = FALSE)
  expect_length(resolve_calls, 1)
  expect_equal(resolve_calls[[1]]$pathway, "Ferroptosis")
  expect_false(resolve_calls[[1]]$use_recommended)
  expect_equal(sort(resolve_calls[[1]]$result_ids), sort(test_gsea_multi$ID))
  expect_equal(captured$geneSetID, "Resolved_Pathway")
})

# ---------- Test plot_death_gsea_ranked ----------
test_that("plot_death_gsea_ranked exact pathway ID with default and custom title", {
  captured <- list()
  local_mocked_bindings(
    gseaplot = function(gsea_res, geneSetID, by, color.line, base_size, title, ...) {
      captured$geneSetID <<- geneSetID
      captured$by <<- by
      captured$title <<- title
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  # default title falls back to the pathway name
  p <- plot_death_gsea_ranked(test_gsea_multi, pathway = "Apoptosis_Pan Cancer_37968457")
  expect_s3_class(p, "ggplot")
  expect_equal(captured$geneSetID, "Apoptosis_Pan Cancer_37968457")
  expect_equal(captured$by, "runningScore")
  expect_equal(captured$title, "Apoptosis_Pan Cancer_37968457")

  # custom title overrides the default
  p2 <- plot_death_gsea_ranked(test_gsea_multi, pathway = "Apoptosis_Pan Cancer_37968457",
                               title = "My Custom Title")
  expect_equal(captured$title, "My Custom Title")
})

test_that("plot_death_gsea_ranked death type resolves to the recommended gene set", {
  local_mocked_bindings(
    get_representative_pathway = function(death_type) {
      list(Pathway = "Ferroptosis_FerrDb V2(36305834)", PMID = "36305834")
    },
    .package = "CDEnrich"
  )
  captured <- list()
  local_mocked_bindings(
    gseaplot = function(gsea_res, geneSetID, by, color.line, base_size, title, ...) {
      captured$geneSetID <<- geneSetID
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  expect_message(
    p <- plot_death_gsea_ranked(test_gsea_single, pathway = "Ferroptosis"),
    "recommended representative gene set"
  )
  expect_s3_class(p, "ggplot")
  expect_equal(captured$geneSetID, "Ferroptosis_FerrDb V2(36305834)")
})

test_that("plot_death_gsea_ranked saves plot and forwards extra arguments", {
  captured <- list()
  local_mocked_bindings(
    gseaplot = function(gsea_res, geneSetID, by, color.line, base_size, title, ...) {
      captured$color.line <<- color.line
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))
  p <- plot_death_gsea_ranked(test_gsea_multi, pathway = "Apoptosis_Pan Cancer_37968457",
                              color = "red", filename = tmp)
  expect_true(file.exists(tmp))
  expect_equal(captured$color.line, "red")
})

test_that("plot_death_gsea_ranked handles invalid pathway and empty result", {
  expect_error(
    plot_death_gsea_ranked(test_gsea_multi, pathway = "NotExistPathway"),
    "Pathway name does not exist"
  )
  empty_df <- test_gsea_multi[0, ]
  expect_error(
    plot_death_gsea_ranked(empty_df, pathway = "Ferroptosis"),
    "GSEA result is empty"
  )
  expect_error(
    plot_death_gsea_ranked(NULL, pathway = "Ferroptosis"),
    "GSEA result is empty"
  )
})
