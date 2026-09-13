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

# Helper: build a valid expression matrix (genes x samples) that overlaps
# with the curated cell death genes (so overlap checks pass).
build_expr <- function(n_sample = 4) {
  genes <- c("GPX4", "ACSL4", "SLC7A11", "FTH1", "TFRC",
             "GSDMD", "CASP1", "IL1B", "CASP3", "BAX", "BCL2",
             paste0("Filler", 1:8))
  mat <- matrix(rnorm(length(genes) * n_sample), nrow = length(genes),
                dimnames = list(genes, paste0("Sample", seq_len(n_sample))))
  mat
}

#    each test defines its own inline GSVA mock.

# ===================== Test celldeath_ssgsea =====================
test_that("celldeath_ssgsea passes correct parameters to ssgseaParam/gsva", {
  cap <- new.env()

  local_mocked_bindings(
    ssgseaParam = function(exprData, geneSets, minSize, normalize, ...) {
      cap$geneSets <<- geneSets
      cap$minSize <<- minSize
      cap$normalize <<- normalize
      cap$calls <- if(is.null(cap$calls)) 1 else cap$calls + 1
      list(exprData=exprData, geneSets=geneSets, minSize=minSize, normalize=normalize)
    },
    gsva = function(param, verbose, ...) {
      matrix(0.2, nrow=length(param$geneSets), ncol=ncol(param$exprData),
             dimnames = list(names(param$geneSets), colnames(param$exprData)))
    },
    .package = "GSVA"
  )

  expr <- build_expr()
  res <- celldeath_ssgsea(expr)

  expect_equal(sort(names(cap$geneSets)), c("Apoptosis", "Ferroptosis", "Pyroptosis"))
  expect_equal(cap$minSize, 3)
  expect_true(cap$normalize)
  expect_equal(cap$calls, 1)
})

test_that("celldeath_ssgsea term_select filters to a single pathway", {
  cap <- new.env()
  local_mocked_bindings(
    ssgseaParam = function(exprData, geneSets, minSize, normalize, ...) {
      cap$geneSets <<- geneSets
      cap$minSize <<- minSize
      cap$normalize <<- normalize
      cap$calls <- if(is.null(cap$calls)) 1 else cap$calls + 1
      list(exprData=exprData, geneSets=geneSets, minSize=minSize, normalize=normalize)
    },
    gsva = function(param, verbose, ...) {
      matrix(0.2, nrow=length(param$geneSets), ncol=ncol(param$exprData),
             dimnames = list(names(param$geneSets), colnames(param$exprData)))
    },
    .package = "GSVA"
  )

  expr <- build_expr()
  celldeath_ssgsea(expr, term_select = "Ferroptosis")
  expect_equal(names(cap$geneSets), "Ferroptosis")
})

# ---- term_select as an EXACT gene set name (used directly) ----
test_that("celldeath_ssgsea term_select as exact gene set name is used directly", {
  local_mocked_bindings(
    load_cell_death_genes = function() {
      list(
        Ferroptosis_Up_36747082   = c("GPX4", "ACSL4", "SLC7A11"),
        Ferroptosis_Down_36747082 = c("FTH1", "TFRC", "NCOA4"),
        Pyroptosis                = c("GSDMD", "CASP1", "IL1B")
      )
    },
    .package = "CDEnrich"
  )

  cap <- new.env()

  local_mocked_bindings(
    ssgseaParam = function(exprData, geneSets, minSize, normalize, ...) {
      cap$geneSets <<- geneSets
      cap$minSize <<- minSize
      cap$normalize <<- normalize
      cap$calls <- if(is.null(cap$calls)) 1 else cap$calls + 1
      list(exprData=exprData, geneSets=geneSets, minSize=minSize, normalize=normalize)
    },
    gsva = function(param, verbose, ...) {
      matrix(0.2,
             nrow = length(param$geneSets),
             ncol = ncol(param$exprData),
             dimnames = list(names(param$geneSets), colnames(param$exprData)))
    },
    .package = "GSVA"
  )

  expr <- build_expr()
  celldeath_ssgsea(expr, term_select = "Ferroptosis_Up_36747082")
  expect_equal(names(cap$geneSets), "Ferroptosis_Up_36747082")
})

# ---- use_recommended controls gene sets for a multi-set death type ----
test_that("celldeath_ssgsea use_recommended controls gene sets for a multi-set death type", {
  local_mocked_bindings(
    load_cell_death_genes = function() {
      list(
        Ferroptosis_Up_36747082   = c("GPX4", "ACSL4", "SLC7A11"),
        Ferroptosis_Down_36747082 = c("FTH1", "TFRC", "NCOA4"),
        Pyroptosis                = c("GSDMD", "CASP1", "IL1B")
      )
    },
    # added PMID to match the real function's return shape;
    # .resolve_term_select uses rec$PMID in its message.
    get_representative_pathway = function(death_type) {
      list(Pathway = "Ferroptosis_Up_36747082", PMID = "36747082")
    },
    .package = "CDEnrich"
  )

  cap <- new.env()
  local_mocked_bindings(
    ssgseaParam = function(exprData, geneSets, minSize, normalize, ...) {
      cap$geneSets <<- geneSets
      cap$minSize <<- minSize
      cap$normalize <<- normalize
      cap$calls <- if(is.null(cap$calls)) 1 else cap$calls + 1
      list(exprData=exprData, geneSets=geneSets, minSize=minSize, normalize=normalize)
    },
    gsva = function(param, verbose, ...) {
      matrix(0.2,
             nrow = length(param$geneSets),
             ncol = ncol(param$exprData),
             dimnames = list(names(param$geneSets), colnames(param$exprData)))
    },
    .package = "GSVA"
  )

  expr <- build_expr()

  # use_recommended = TRUE (default): only the recommended representative gene set
  celldeath_ssgsea(expr, term_select = "Ferroptosis", use_recommended = TRUE)
  expect_equal(sort(names(cap$geneSets)), "Ferroptosis_Up_36747082")

  # use_recommended = FALSE: ALL gene sets under Ferroptosis
  celldeath_ssgsea(expr, term_select = "Ferroptosis", use_recommended = FALSE)
  expect_equal(sort(names(cap$geneSets)),
               c("Ferroptosis_Down_36747082", "Ferroptosis_Up_36747082"))
})

# ---- .resolve_term_select is delegated the right arguments ----
test_that("celldeath_ssgsea delegates term_select/use_recommended to .resolve_term_select", {
  resolve_calls <- list()
  local_mocked_bindings(
    .resolve_term_select = function(term_select, death_gene_list_all, use_recommended) {
      resolve_calls <<- append(resolve_calls,
                               list(list(term_select = term_select,
                                         use_recommended = use_recommended,
                                         n_sets = length(death_gene_list_all))))
      # More genes than the default min_size = 3, so the set is not filtered out
      list(Resolved = c("GPX4", "ACSL4", "SLC7A11", "GSS"))
    },
    .package = "CDEnrich"
  )

  cap <- new.env()
  local_mocked_bindings(
    ssgseaParam = function(exprData, geneSets, minSize, normalize, ...) {
      cap$geneSets <<- geneSets
      cap$minSize <<- minSize
      cap$normalize <<- normalize
      cap$calls <- if(is.null(cap$calls)) 1 else cap$calls + 1
      list(exprData=exprData, geneSets=geneSets, minSize=minSize, normalize=normalize)
    },
    gsva = function(param, verbose, ...) {
      matrix(0.2,
             nrow = length(param$geneSets),
             ncol = ncol(param$exprData),
             dimnames = list(names(param$geneSets), colnames(param$exprData)))
    },
    .package = "GSVA"
  )

  expr <- build_expr()
  celldeath_ssgsea(expr, term_select = "Ferroptosis", use_recommended = FALSE)
  expect_length(resolve_calls, 1)
  expect_equal(resolve_calls[[1]]$term_select, "Ferroptosis")
  expect_false(resolve_calls[[1]]$use_recommended)
  expect_equal(resolve_calls[[1]]$n_sets, 3)               # full curated collection handed in
  expect_equal(names(cap$geneSets), "Resolved")            # resolved set flows into gsva
})

# ---- min_size filtering ----
test_that("celldeath_ssgsea drops gene sets smaller than min_size", {
  local_mocked_bindings(
    load_cell_death_genes = function() {
      list(
        BigSet   = c("GPX4", "ACSL4", "SLC7A11"),
        SmallSet = c("GSDMD", "CASP1")     # length 2
      )
    },
    .package = "CDEnrich"
  )

  cap <- new.env()
  local_mocked_bindings(
    ssgseaParam = function(exprData, geneSets, minSize, normalize, ...) {
      cap$geneSets <<- geneSets
      cap$minSize <<- minSize
      cap$normalize <<- normalize
      cap$calls <- if(is.null(cap$calls)) 1 else cap$calls + 1
      list(exprData=exprData, geneSets=geneSets, minSize=minSize, normalize=normalize)
    },
    gsva = function(param, verbose, ...) {
      matrix(0.2,
             nrow = length(param$geneSets),
             ncol = ncol(param$exprData),
             dimnames = list(names(param$geneSets), colnames(param$exprData)))
    },
    .package = "GSVA"
  )

  expr <- build_expr()
  celldeath_ssgsea(expr, min_size = 3)
  expect_equal(names(cap$geneSets), "BigSet")
})

test_that("celldeath_ssgsea errors when no pathway survives min_size filtering", {
  local_mocked_bindings(
    load_cell_death_genes = function() {
      list(SmallSet = c("GSDMD", "CASP1"))   # length 2 < min_size 3
    },
    .package = "CDEnrich"
  )
  expr <- build_expr()
  expect_error(
    celldeath_ssgsea(expr, min_size = 3),
    "No cell death pathway left after min_size filtering"
  )
})

# ---- overlap checks ----
test_that("celldeath_ssgsea errors when expression matrix has no gene overlap", {
  # matrix with only filler genes (no cell death genes)
  m <- matrix(rnorm(20), nrow = 5,
              dimnames = list(paste0("Filler", 1:5), paste0("Sample", 1:4)))
  expect_error(
    celldeath_ssgsea(m),
    "No overlap between expression matrix rownames and cell death genes"
  )
})

test_that("celldeath_ssgsea warns when overlap is smaller than min_size", {
  cap <- new.env()
  local_mocked_bindings(
    ssgseaParam = function(exprData, geneSets, minSize, normalize, ...) {
      cap$geneSets <<- geneSets
      list(exprData=exprData, geneSets=geneSets, minSize=minSize, normalize=normalize)
    },
    gsva = function(param, verbose, ...) {
      matrix(0.2,
             nrow = length(param$geneSets),
             ncol = ncol(param$exprData),
             dimnames = list(names(param$geneSets), colnames(param$exprData)))
    },
    .package = "GSVA"
  )

  genes <- c("GPX4", paste0("Filler", 1:8))
  m <- matrix(rnorm(length(genes) * 3), nrow = length(genes),
              dimnames = list(genes, paste0("Sample", 1:3)))
  expect_warning(
    celldeath_ssgsea(m),
    "cell death genes found"
  )
})

# ---- input validation on expr ----
test_that("celldeath_ssgsea input validation on expr", {
  # not a numeric matrix
  expect_error(celldeath_ssgsea(list()), "numeric matrix")
  expect_error(celldeath_ssgsea(matrix("a", 2, 2)), "numeric matrix")
  # numeric vector (not matrix)
  expect_error(celldeath_ssgsea(c(A = 1, B = 2, C = 3)), "numeric matrix")

  # missing rownames
  m1 <- matrix(rnorm(4), nrow = 2)
  colnames(m1) <- c("S1", "S2")
  expect_error(celldeath_ssgsea(m1), "gene symbols as rownames")

  # missing colnames
  m2 <- matrix(rnorm(4), nrow = 2)
  rownames(m2) <- c("G1", "G2")
  expect_error(celldeath_ssgsea(m2), "sample names as colnames")

  # duplicated rownames
  m3 <- matrix(rnorm(4), nrow = 2,
               dimnames = list(c("Gene1", "Gene1"), c("S1", "S2")))
  expect_error(celldeath_ssgsea(m3), "Duplicated gene symbols")

  # NA in matrix
  m4 <- matrix(rnorm(4), nrow = 2,
               dimnames = list(c("G1", "G2"), c("S1", "S2")))
  m4[1, 1] <- NA
  expect_error(celldeath_ssgsea(m4), "cannot contain missing values")
})

# ---- group validation (error cases only) ----
test_that("celldeath_ssgsea group label validation", {
  expr <- build_expr(n_sample = 3)

  # length mismatch
  expect_error(
    celldeath_ssgsea(expr, group = c("A", "B")),
    "length of group"
  )
  # NA / empty labels
  expect_error(
    celldeath_ssgsea(expr, group = c("A", NA, "B")),
    "Group labels cannot contain NA"
  )
  expect_error(
    celldeath_ssgsea(expr, group = c("A", "", "B")),
    "Group labels cannot contain NA"
  )
  # fewer than 2 groups
  expect_error(
    celldeath_ssgsea(expr, group = rep("A", 3)),
    "At least 2 distinct groups"
  )
})

test_that("celldeath_ssgsea warns when a group has fewer than 2 samples", {
  local_mocked_bindings(
    ssgseaParam = function(exprData, geneSets, minSize, normalize, ...) {
      list(exprData = exprData, geneSets = geneSets,
           minSize = minSize, normalize = normalize)
    },
    gsva = function(param, verbose, ...) {
      matrix(0.2,
             nrow = length(param$geneSets),
             ncol = ncol(param$exprData),
             dimnames = list(names(param$geneSets), colnames(param$exprData)))
    },
    .package = "GSVA"
  )

  expect_warning(
    celldeath_ssgsea(build_expr(n_sample = 3), group = c("A", "A", "B")),
    "fewer than 2 samples"
  )
})

# ---- result object structure ----
test_that("celldeath_ssgsea returns a celldeath_ssgsea object with score/group", {
  cap <- new.env()
  local_mocked_bindings(
    ssgseaParam = function(exprData, geneSets, minSize, normalize, ...) {
      cap$geneSets <<- geneSets
      list(exprData=exprData, geneSets=geneSets, minSize=minSize, normalize=normalize)
    },
    gsva = function(param, verbose, ...) {
      matrix(0.2,
             nrow = length(param$geneSets),
             ncol = ncol(param$exprData),
             dimnames = list(names(param$geneSets), colnames(param$exprData)))
    },
    .package = "GSVA"
  )

  expr <- build_expr(n_sample = 3)

  # without group
  res <- celldeath_ssgsea(expr)
  expect_s3_class(res, "celldeath_ssgsea")
  expect_true(is.matrix(res$score))
  expect_null(res$group)
  expect_equal(nrow(res$score), 3)   # 3 pathways
  expect_equal(ncol(res$score), 3)   # 3 samples

  # with group -> stored as factor
  grp <- factor(c("Control", "Treatment", "Treatment"))
  expect_warning(
    res_g <- celldeath_ssgsea(expr, group = grp),
    "fewer than 2 samples"
  )
  expect_s3_class(res_g, "celldeath_ssgsea")
  expect_true(is.factor(res_g$group))
  expect_equal(as.character(res_g$group), as.character(grp))
})

# ---- savefile ----
test_that("celldeath_ssgsea writes CSV when savefile=TRUE", {
  cap <- new.env()
  local_mocked_bindings(
    ssgseaParam = function(exprData, geneSets, minSize, normalize, ...) {
      cap$geneSets <<- geneSets
      list(exprData=exprData, geneSets=geneSets, minSize=minSize, normalize=normalize)
    },
    gsva = function(param, verbose, ...) {
      matrix(0.2,
             nrow = length(param$geneSets),
             ncol = ncol(param$exprData),
             dimnames = list(names(param$geneSets), colnames(param$exprData)))
    },
    .package = "GSVA"
  )

  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  expr <- build_expr()
  res <- celldeath_ssgsea(expr, savefile = TRUE, filename = tmp)
  expect_true(file.exists(tmp))
  csv <- read.csv(tmp)
  expect_true("Pathway" %in% names(csv))
})

# ===================== Test celldeath_ssgsea_diff =====================

# Build a celldeath_ssgsea object directly (no real GSVA run needed)
make_ssgsea_obj <- function(score_mat, group) {
  obj <- list(score = score_mat, group = group)
  class(obj) <- "celldeath_ssgsea"
  obj
}

test_that("celldeath_ssgsea_diff input validation", {
  # not a celldeath_ssgsea object
  expect_error(
    celldeath_ssgsea_diff(list(score = matrix(), group = factor("A"))),
    "must be the output of celldeath_ssgsea"
  )
  # missing group
  obj <- make_ssgsea_obj(matrix(1, nrow = 1, dimnames = list("P1", "S1")), NULL)
  expect_error(
    celldeath_ssgsea_diff(obj),
    "No group information found"
  )
})

# ---- test selection: 2 groups -> wilcox, >=3 groups -> kruskal ----
test_that("celldeath_ssgsea_diff uses wilcox.test for 2 groups", {
  captured <- list()
  local_mocked_bindings(
    wilcox.test  = function(x, y, ...) { captured$test <<- "wilcox";  list(statistic = 10, p.value = 0.05) },
    kruskal.test = function(x, g) { captured$test <<- "kruskal"; list(statistic = 8,  p.value = 0.03) },
    p.adjust     = function(p, method) { captured$method <<- method; p },
    .package = "stats"
  )
  score_mat <- matrix(c(1, 4, 2, 5, 3, 6), nrow = 2,
                      dimnames = list(c("PathA", "PathB"), c("S1", "S2", "S3")))
  group <- factor(c("Ctrl", "Ctrl", "Treat"))
  obj <- make_ssgsea_obj(score_mat, group)
  diff_df <- celldeath_ssgsea_diff(obj)
  expect_equal(captured$test, "wilcox")
  expect_true(all(c("Pathway", "Ctrl", "Treat", "Statistic", "P_value", "P_adjust") %in% names(diff_df)))
})

test_that("celldeath_ssgsea_diff uses kruskal.test for 3+ groups", {
  captured <- list()
  local_mocked_bindings(
    wilcox.test  = function(x, y, ...) { captured$test <<- "wilcox";  list(statistic = 10, p.value = 0.05) },
    kruskal.test = function(x, g) { captured$test <<- "kruskal"; list(statistic = 8,  p.value = 0.03) },
    p.adjust     = function(p, method) { captured$method <<- method; p },
    .package = "stats"
  )
  score_mat <- matrix(1:9, nrow = 3,
                      dimnames = list(c("PathA", "PathB", "PathC"),
                                      c("S1", "S2", "S3")))
  group <- factor(c("A", "B", "C"))
  obj <- make_ssgsea_obj(score_mat, group)
  celldeath_ssgsea_diff(obj)
  expect_equal(captured$test, "kruskal")
})

# ---- pAdjustMethod forwarded ----
test_that("celldeath_ssgsea_diff forwards pAdjustMethod to p.adjust", {
  captured <- list()
  local_mocked_bindings(
    wilcox.test  = function(x, y, ...) list(statistic = 1, p.value = 0.2),
    kruskal.test = function(x, g) list(statistic = 1, p.value = 0.2),
    p.adjust     = function(p, method) { captured$method <<- method; p },
    .package = "stats"
  )
  score_mat <- matrix(c(1,2,3,4,5,6,7,8), nrow = 2,
                      dimnames = list(c("P1", "P2"), c("S1", "S2", "S3", "S4")))
  group <- factor(c("G1", "G1", "G2", "G2"))
  obj <- make_ssgsea_obj(score_mat, group)
  # default BH
  celldeath_ssgsea_diff(obj)
  expect_equal(captured$method, "BH")
  # custom method
  celldeath_ssgsea_diff(obj, pAdjustMethod = "bonferroni")
  expect_equal(captured$method, "bonferroni")
})

# ---- output is sorted by p-value and has reset rownames ----
test_that("celldeath_ssgsea_diff output is sorted by p-value with reset rownames", {
  local_mocked_bindings(
    wilcox.test  = function(x, y, ...) list(statistic = 1, p.value = 0.2),
    kruskal.test = function(x, g) list(statistic = 1, p.value = 0.2),
    p.adjust     = function(p, method) p,
    .package = "stats"
  )
  score_mat <- matrix(rnorm(20), nrow = 5,
                      dimnames = list(paste0("Path", 1:5), paste0("S", 1:4)))
  group <- factor(c("G1", "G1", "G2", "G2"))

  obj <- list(
    score = score_mat,
    group = group
  )
  class(obj) <- "celldeath_ssgsea"

  diff_df <- celldeath_ssgsea_diff(obj)
  expect_true(is.data.frame(diff_df))
  expect_equal(rownames(diff_df), as.character(seq_len(nrow(diff_df))))
  expect_true(all(diff(diff_df$P_value) >= 0))
  expect_true("P_adjust" %in% names(diff_df))
})

# ---- savefile ----
test_that("celldeath_ssgsea_diff writes CSV when savefile=TRUE", {
  local_mocked_bindings(
    wilcox.test  = function(x, y, ...) list(statistic = 1, p.value = 0.2),
    kruskal.test = function(x, g) list(statistic = 1, p.value = 0.2),
    p.adjust     = function(p, method) p,
    .package = "stats"
  )
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  score_mat <- matrix(rnorm(8), nrow = 2,
                      dimnames = list(c("PathA", "PathB"), c("S1", "S2", "S3", "S4")))
  group <- factor(c("G1", "G1", "G2", "G2"))
  obj <- make_ssgsea_obj(score_mat, group)
  celldeath_ssgsea_diff(obj, savefile = TRUE, filename = tmp)
  expect_true(file.exists(tmp))
  csv <- read.csv(tmp)
  expect_true("Pathway" %in% names(csv))
})

# ===================== Test print.celldeath_ssgsea =====================

test_that("print.celldeath_ssgsea prints a summary and returns object invisibly", {
  score_mat <- matrix(1:6, nrow = 2, dimnames = list(c("P1", "P2"), c("S1", "S2", "S3")))
  group <- factor(c("Control", "Treatment", "Treatment"))
  obj <- make_ssgsea_obj(score_mat, group)
  expect_output(print(obj), "Cell death ssGSEA result")
  expect_output(print(obj), "Pathways scored:")
  expect_output(print(obj), "Samples scored")
  expect_output(print(obj), "Groups")
  expect_identical(print(obj), obj)
})
