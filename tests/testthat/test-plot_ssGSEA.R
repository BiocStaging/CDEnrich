library(testthat)

# Mock pheatmap and grid.draw to avoid real plotting and graphics devices
local_mocked_bindings(
  grid.draw = function(x, ...) invisible(NULL),
  .package = "grid"
)

# Helper: construct a celldeath_ssgsea object directly (no real GSVA run needed)
make_ssgsea_obj <- function(score_mat, group) {
  structure(list(score = score_mat, group = group), class = "celldeath_ssgsea")
}

# Test data: multiple pathways x 4 samples (global mode)
score_mat_multi <- matrix(
  c(1, 2, 3, 4,
    5, 6, 7, 8,
    2, 3, 4, 5),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("Ferroptosis", "Pyroptosis", "Apoptosis"), paste0("S", 1:4))
)
group4 <- factor(c("Control", "Control", "Treatment", "Treatment"),
                 levels = c("Control", "Treatment"))
obj_multi <- make_ssgsea_obj(score_mat_multi, group4)
obj_multi_nogroup <- make_ssgsea_obj(score_mat_multi, NULL)

# Test data: single pathway (targeted mode)
score_mat_single <- matrix(c(1, 2, 3, 4), nrow = 1,
                           dimnames = list("Ferroptosis", paste0("S", 1:4)))
obj_single <- make_ssgsea_obj(score_mat_single, group4)

# ---------- Test plot_death_ssgsea_heatmap ----------
test_that("plot_death_ssgsea_heatmap input validation", {
  expect_error(plot_death_ssgsea_heatmap(data.frame()),
               "must be the output of celldeath_ssgsea",
               fixed = TRUE)

  # 修复：明确 nrow=0，ncol=4，0行4列空矩阵，dimnames才匹配
  empty_mat <- matrix(numeric(0), nrow = 0, ncol = 4,
                      dimnames = list(NULL, paste0("S", 1:4)))
  empty_obj <- make_ssgsea_obj(empty_mat, group4)

  expect_error(plot_death_ssgsea_heatmap(empty_obj),
               "score matrix is empty",
               fixed = TRUE)
})

test_that("plot_death_ssgsea_heatmap multi-pathway mode applies row z-score", {
  cap <- new.env()
  local_mocked_bindings(
    pheatmap = function(mat, ...) {
      cap$mat <<- mat
      cap$args <<- list(...)
      list(gtable = grid::nullGrob())
    },
    .package = "pheatmap"
  )

  p <- plot_death_ssgsea_heatmap(obj_multi, scale_rows = TRUE)
  # each row is z-scored across samples (mean ~ 0)
  expect_true(all(abs(rowMeans(cap$mat)) < 1e-10))
  expect_equal(nrow(cap$mat), 3)
  # multi-pathway defaults: row clustering on, row names shown
  expect_true(cap$args$cluster_rows)
  expect_true(cap$args$show_rownames)
  # group annotation attached
  expect_false(is.null(cap$args$annotation_col))
  expect_true("Group" %in% colnames(cap$args$annotation_col))
  # invisibly returns the pheatmap-like object
  expect_true("gtable" %in% names(p))
})

test_that("plot_death_ssgsea_heatmap scale_rows=FALSE keeps raw scores", {
  cap <- new.env()
  local_mocked_bindings(
    pheatmap = function(mat, ...) {
      cap$mat <<- mat
      list(gtable = grid::nullGrob())
    },
    .package = "pheatmap"
  )

  plot_death_ssgsea_heatmap(obj_multi, scale_rows = FALSE)
  expect_identical(unname(cap$mat), unname(score_mat_multi))
})

test_that("plot_death_ssgsea_heatmap targeted single-pathway adjustments", {
  cap <- new.env()
  local_mocked_bindings(
    pheatmap = function(mat, ...) {
      cap$mat <<- mat
      cap$args <<- list(...)
      list(gtable = grid::nullGrob())
    },
    .package = "pheatmap"
  )

  plot_death_ssgsea_heatmap(obj_single)
  # pathway name moved into the title, row label hidden
  expect_true(grepl("Ferroptosis", cap$args$main))
  expect_false(cap$args$show_rownames)
  # no clustering on a single row; samples ordered by group instead
  expect_false(cap$args$cluster_cols)
  expect_equal(colnames(cap$mat), c("S1", "S2", "S3", "S4"))
})

test_that("plot_death_ssgsea_heatmap orders samples by group in single-pathway mode", {
  shuffled_group <- factor(c("Treatment", "Control", "Treatment", "Control"),
                           levels = c("Control", "Treatment"))
  obj_shuffled <- make_ssgsea_obj(score_mat_single, shuffled_group)

  cap <- new.env()
  local_mocked_bindings(
    pheatmap = function(mat, ...) {
      cap$mat <<- mat
      list(gtable = grid::nullGrob())
    },
    .package = "pheatmap"
  )

  plot_death_ssgsea_heatmap(obj_shuffled)
  # Control samples (S2, S4) come before Treatment samples (S1, S3)
  expect_equal(colnames(cap$mat), c("S2", "S4", "S1", "S3"))
})

test_that("plot_death_ssgsea_heatmap works without group and saves to file", {
  cap <- new.env()
  local_mocked_bindings(
    pheatmap = function(mat, ...) {
      cap$args <<- list(...)
      list(gtable = grid::nullGrob())
    },
    .package = "pheatmap"
  )

  plot_death_ssgsea_heatmap(obj_multi_nogroup)
  expect_null(cap$args$annotation_col)
  expect_null(cap$args$annotation_colors)

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))
  plot_death_ssgsea_heatmap(obj_multi, filename = tmp)
  expect_true(file.exists(tmp))
})

# Show_category subsetting in heatmap
test_that("plot_death_ssgsea_heatmap show_category caps displayed pathways", {
  set.seed(1)
  mat_big <- matrix(rnorm(50 * 4), nrow = 50,
                    dimnames = list(paste0("Path", 1:50), paste0("S", 1:4)))
  obj_big <- make_ssgsea_obj(mat_big, group4)

  cap <- new.env()
  local_mocked_bindings(
    pheatmap = function(mat, ...) {
      cap$mat <<- mat
      list(gtable = grid::nullGrob())
    },
    .package = "pheatmap"
  )

  plot_death_ssgsea_heatmap(obj_big, show_category = 10)
  expect_equal(nrow(cap$mat), 10)
  expect_equal(ncol(cap$mat), 4)
})

# ---------- Test plot_death_ssgsea_boxplot ----------
test_that("plot_death_ssgsea_boxplot input validation", {
  expect_error(plot_death_ssgsea_boxplot(data.frame()),
               "must be the output of celldeath_ssgsea")
  expect_error(plot_death_ssgsea_boxplot(obj_multi_nogroup),
               "No group information found")
})

test_that("plot_death_ssgsea_boxplot multi-pathway mode returns faceted ggplot", {
  p <- plot_death_ssgsea_boxplot(obj_multi)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetWrap")

  df_out <- plot_death_ssgsea_boxplot(obj_multi, return_data = TRUE)
  expect_true(is.data.frame(df_out))
  expect_true(all(c("Sample", "Pathway", "Score", "Group") %in% colnames(df_out)))
  # 3 pathways x 4 samples in long format
  expect_equal(nrow(df_out), 12)
  # group labels carried into the long data frame
  expect_equal(sum(df_out$Group == "Control"), 6)

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))
  plot_death_ssgsea_boxplot(obj_multi, filename = tmp)
  expect_true(file.exists(tmp))
})

test_that("plot_death_ssgsea_boxplot show_category caps displayed pathways", {
  df_out <- plot_death_ssgsea_boxplot(obj_multi, show_category = 1,
                                      return_data = TRUE)
  expect_equal(length(unique(as.character(df_out$Pathway))), 1)
  expect_equal(nrow(df_out), 4)
})

test_that("plot_death_ssgsea_boxplot targeted single-pathway mode outputs INFO message", {
  expect_message(
    p <- plot_death_ssgsea_boxplot(obj_single),
    "INFO: Current result contains a single pathway"
  )
  expect_s3_class(p, "ggplot")
  # pathway name appended to the title
  expect_true(grepl("Ferroptosis", p$labels$title))
})

test_that("plot_death_ssgsea_boxplot add_points toggles the jitter layer", {
  p_with <- plot_death_ssgsea_boxplot(obj_multi, add_points = TRUE)
  p_without <- plot_death_ssgsea_boxplot(obj_multi, add_points = FALSE)
  expect_length(p_with$layers, 2)
  expect_length(p_without$layers, 1)
})

# Long pathway names are wrapped at underscores but underscores are kept
test_that("plot_death_ssgsea_boxplot wraps long pathway labels at underscores", {
  mat_long <- matrix(
    c(1, 2, 3, 4,
      5, 6, 7, 8),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("Necroptosis_Glioma_35719999", "Ferroptosis"),
                    paste0("S", 1:4))
  )
  obj_long <- make_ssgsea_obj(mat_long, group4)

  df_out <- plot_death_ssgsea_boxplot(obj_long, return_data = TRUE)
  lvls <- levels(df_out$Pathway)

  # Short name unchanged; long name wrapped after an underscore
  expect_true("Ferroptosis" %in% lvls)
  expect_true("Necroptosis_Glioma_\n35719999" %in% lvls)
  # Underscores are preserved (no replacement with spaces)
  expect_false(any(grepl("Necroptosis Glioma", lvls, fixed = TRUE)))
})

# ---------- Test plot_death_ssgsea_bar ----------
test_that("plot_death_ssgsea_bar input validation", {
  expect_error(plot_death_ssgsea_bar(data.frame()),
               "must be the output of celldeath_ssgsea")
  expect_error(plot_death_ssgsea_bar(obj_multi_nogroup),
               "requires group information")
})

test_that("plot_death_ssgsea_bar computes correct mean and SEM per group", {
  score_known <- matrix(c(1, 2, 3, 5, 6, 7), nrow = 1,
                        dimnames = list("Ferroptosis", paste0("S", 1:6)))
  group6 <- factor(rep(c("Control", "Treatment"), each = 3),
                   levels = c("Control", "Treatment"))
  obj_known <- make_ssgsea_obj(score_known, group6)

  suppressMessages(
    df_out <- plot_death_ssgsea_bar(obj_known, return_data = TRUE)
  )
  ctrl <- df_out[df_out$Group == "Control", ]
  treat <- df_out[df_out$Group == "Treatment", ]
  expect_equal(ctrl$Mean, 2)
  expect_equal(treat$Mean, 6)
  # SEM = sd / sqrt(n)
  expect_equal(ctrl$SEM, sd(c(1, 2, 3)) / sqrt(3))
  expect_equal(treat$SEM, sd(c(5, 6, 7)) / sqrt(3))
})

test_that("plot_death_ssgsea_bar orders pathways by overall mean score", {
  score_two <- matrix(
    c(1, 1, 1, 1,
      9, 9, 9, 9),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("LowPath", "HighPath"), paste0("S", 1:4))
  )
  obj_two <- make_ssgsea_obj(score_two, group4)
  df_out <- plot_death_ssgsea_bar(obj_two, return_data = TRUE)
  # ascending mean order: LowPath first, HighPath last
  expect_equal(levels(df_out$Pathway), c("LowPath", "HighPath"))
})

test_that("plot_death_ssgsea_bar returns flipped horizontal barplot and saves file", {
  p <- plot_death_ssgsea_bar(obj_multi)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$coordinates, "CoordFlip"))

  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))
  plot_death_ssgsea_bar(obj_multi, filename = tmp)
  expect_true(file.exists(tmp))
})

test_that("plot_death_ssgsea_bar targeted single-pathway mode outputs INFO message", {
  expect_message(
    plot_death_ssgsea_bar(obj_single),
    "INFO: Current result contains a single pathway"
  )
})

# Show_category subsetting in bar plot
test_that("plot_death_ssgsea_bar show_category caps displayed pathways", {
  set.seed(1)
  mat_big <- matrix(rnorm(30 * 4, mean = 1), nrow = 30,
                    dimnames = list(paste0("Path", 1:30), paste0("S", 1:4)))
  obj_big <- make_ssgsea_obj(mat_big, group4)

  df_out <- plot_death_ssgsea_bar(obj_big, show_category = 5,
                                  return_data = TRUE)
  expect_equal(length(unique(as.character(df_out$Pathway))), 5)
  # 5 pathways x 2 groups
  expect_equal(nrow(df_out), 10)
})
