# Tests for RMdifTree()
#
# Covers:
#   - Input validation (covariates, model, min_n_per_level)
#   - Auto model + effect-size selection (RM/MH vs PCM/pgamma)
#   - Output structures: kable (default), dataframe, tree, plot
#   - Flagged + Rescaled columns
#   - on_rescale modes (message / warning / stop)
#   - Stability assessment via stablelearner
#   - Covariate-name propagation (bare symbol, dif$age, dif[["age"]])
#
# All tests use small data so each fit completes in well under 1s.

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
need_tree_pkgs <- function() {
  skip_if_not_installed("psychotree")
  skip_if_not_installed("partykit")
  skip_if_not_installed("difR")
  skip_if_not_installed("iarm")
}

# Dichotomous data with strong gender DIF: items 1-2 favour F, items 3-4
# favour M. Guarantees the tree splits on gender at small minsize.
make_dif_dichotomous <- function(n = 400, seed = 1L) {
  set.seed(seed)
  gender <- factor(rep(c("F", "M"), length.out = n))
  age    <- runif(n, 18, 80)
  theta  <- stats::rnorm(n)
  dif_item <- function(easy_for) {
    b <- ifelse(gender == easy_for, -1, 1)
    stats::rbinom(n, 1, plogis(theta - b))
  }
  items <- data.frame(
    I1 = dif_item("F"),
    I2 = dif_item("F"),
    I3 = dif_item("M"),
    I4 = dif_item("M"),
    I5 = stats::rbinom(n, 1, plogis(theta)),
    I6 = stats::rbinom(n, 1, plogis(theta + 0.3))
  )
  list(items = items,
       covs  = data.frame(age = age, gender = gender))
}

# Polytomous data with a forced rescaling case: item I1 has no category-0
# responses when age > 50, so any age-based split rescales I1 in the old
# subgroup.
make_polytomous_with_rescale <- function(n = 300, seed = 7L) {
  set.seed(seed)
  age   <- runif(n, 18, 80)
  theta <- stats::rnorm(n)
  bake  <- function(theta, b) {
    pmax(0L, pmin(2L, as.integer(round(
      theta - b + stats::rnorm(length(theta), 0, 0.5) + 1
    ))))
  }
  items <- data.frame(
    I1 = bake(theta + (age - 50) * 0.05, 0),
    I2 = bake(theta, 0),
    I3 = bake(theta, 0.3),
    I4 = bake(theta, -0.3)
  )
  # Force no zeros for I1 in the older subgroup → guaranteed rescale.
  items$I1[age > 50 & items$I1 == 0L] <- 1L
  list(items = items,
       covs  = data.frame(age = age))
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMdifTree errors when covariates is missing", {
  d <- make_dif_dichotomous()
  expect_error(RMdifTree(d$items), regexp = "covariates")
})

test_that("RMdifTree errors when covariates row count mismatches data", {
  d <- make_dif_dichotomous()
  expect_error(
    RMdifTree(d$items, covariates = d$covs[seq_len(nrow(d$items) - 1L), ]),
    regexp = "same number of rows"
  )
})

test_that("RMdifTree errors when factor covariate has tiny levels", {
  d <- make_dif_dichotomous()
  d$covs$tiny <- factor(c(rep("a", nrow(d$items) - 2L),
                          rep("b", 2L)))
  expect_error(
    RMdifTree(d$items, covariates = d$covs, min_n_per_level = 20),
    regexp = "fewer than 20"
  )
})

test_that("RMdifTree errors when model = 'RM' on polytomous data", {
  d <- make_polytomous_with_rescale()
  expect_error(
    suppressMessages(
      RMdifTree(d$items, covariates = d$covs, model = "RM",
                minsize = 50)
    ),
    regexp = "dichotomous"
  )
})

# ---------------------------------------------------------------------
# Output structures
# ---------------------------------------------------------------------
test_that("RMdifTree default output is a knitr_kable", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  out <- RMdifTree(d$items, covariates = d$covs, minsize = 50)
  expect_s3_class(out, "knitr_kable")
  expect_equal(attr(out, "model"),       "RM")
  expect_equal(attr(out, "effect_size"), "MH")
})

test_that("RMdifTree output = 'dataframe' on dichotomous data picks RM/MH", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  df <- RMdifTree(d$items, covariates = d$covs,
                  minsize = 50, output = "dataframe")
  expect_s3_class(df, "data.frame")
  expect_true(all(c("NodeID", "Split", "Variable", "Direction",
                    "Item", "EffectSize", "SE", "Class",
                    "Flagged", "Rescaled",
                    "n_left", "n_right") %in% names(df)))
  expect_true(nrow(df) > 0L,
              info = "tree should have at least one split given strong DIF")
  expect_equal(attr(df, "model"),       "RM")
  expect_equal(attr(df, "effect_size"), "MH")
})

test_that("RMdifTree output = 'dataframe' on polytomous data picks PCM/pgamma", {
  need_tree_pkgs()
  d <- make_polytomous_with_rescale()
  set.seed(1)
  df <- suppressMessages(
    RMdifTree(d$items, covariates = d$covs,
              minsize = 50, output = "dataframe")
  )
  expect_s3_class(df, "data.frame")
  expect_equal(attr(df, "model"),       "PCM")
  expect_equal(attr(df, "effect_size"), "pgamma")
})

test_that("RMdifTree output = 'tree' returns an RMdifTree-classed partykit tree", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  tree <- RMdifTree(d$items, covariates = d$covs,
                    minsize = 50, output = "tree")
  expect_s3_class(tree, "RMdifTree")
  expect_s3_class(tree, "party")
  expect_true(!is.null(tree$info$effectsize))
  expect_true(is.list(tree$info$effectsize$nodes))
})

# ---------------------------------------------------------------------
# Flagged + Rescaled columns
# ---------------------------------------------------------------------
test_that("Flagged is logical and TRUE only for Class B/C", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  df <- suppressWarnings(
    RMdifTree(d$items, covariates = d$covs,
              minsize = 50, output = "dataframe")
  )
  expect_type(df$Flagged, "logical")
  # Class is one of A/B/C; rare NA permitted on sparse / unstable cells.
  expect_true(all(df$Class %in% c("A", "B", "C", NA)))
  # `Flagged` mirrors B/C (NA in Class → FALSE in Flagged).
  expected_flag <- !is.na(df$Class) & df$Class %in% c("B", "C")
  expect_equal(df$Flagged, expected_flag)
})

test_that("Rescaled column flags item × split cells affected by terminal-node rescaling", {
  need_tree_pkgs()
  d <- make_polytomous_with_rescale()
  set.seed(1)
  df <- suppressMessages(
    RMdifTree(d$items, covariates = d$covs,
              minsize = 50, output = "dataframe")
  )
  expect_type(df$Rescaled, "logical")
  rt <- attr(df, "rescaled_terminal_items")
  expect_s3_class(rt, "data.frame")
  expect_true("I1" %in% rt$Item)
  expect_true(any(df$Rescaled[df$Item == "I1"]))
})

# ---------------------------------------------------------------------
# on_rescale modes
# ---------------------------------------------------------------------
test_that("on_rescale = 'message' (default) emits message but does not error", {
  need_tree_pkgs()
  d <- make_polytomous_with_rescale()
  set.seed(1)
  expect_message(
    df <- RMdifTree(d$items, covariates = d$covs,
                    minsize = 50, output = "dataframe"),
    regexp = "no responses in their lowest category"
  )
  expect_s3_class(df, "data.frame")
})

test_that("on_rescale = 'warning' raises a warning instead", {
  need_tree_pkgs()
  d <- make_polytomous_with_rescale()
  set.seed(1)
  expect_warning(
    suppressMessages(
      RMdifTree(d$items, covariates = d$covs,
                minsize = 50, on_rescale = "warning",
                output = "dataframe")
    ),
    regexp = "no responses in their lowest category"
  )
})

test_that("on_rescale = 'stop' raises an error instead", {
  need_tree_pkgs()
  d <- make_polytomous_with_rescale()
  set.seed(1)
  expect_error(
    suppressMessages(
      RMdifTree(d$items, covariates = d$covs,
                minsize = 50, on_rescale = "stop",
                output = "dataframe")
    ),
    regexp = "no responses in their lowest category"
  )
})

# ---------------------------------------------------------------------
# Stability assessment
# ---------------------------------------------------------------------
test_that("stability = TRUE attaches a stability summary attribute", {
  need_tree_pkgs()
  skip_if_not_installed("stablelearner")
  d <- make_dif_dichotomous()
  set.seed(1)
  df <- RMdifTree(d$items, covariates = d$covs,
                  minsize = 50,
                  stability = TRUE, stability_B = 5L,
                  output = "dataframe")
  stab <- attr(df, "stability")
  expect_s3_class(stab, "data.frame")
  expect_true(all(c("Variable", "Type", "Selected_orig",
                    "Selection_freq_pct") %in% names(stab)))
  expect_setequal(stab$Variable, names(d$covs))
  # Pre-rendered kable also attached
  expect_s3_class(attr(df, "stability_kable"), "knitr_kable")
})

test_that("stability = TRUE does not clobber a pre-existing message sink", {
  need_tree_pkgs()
  skip_if_not_installed("stablelearner")
  d <- make_dif_dichotomous()
  # Mimic a front-end (e.g. RStudio) that has redirected the message
  # stream before the call. The function must leave that sink in place;
  # resetting it to the default (connection 2) silences console output.
  con <- file(tempfile(), open = "wt")
  withr::defer({
    if (sink.number(type = "message") != 2L) sink(NULL, type = "message")
    close(con)
  })
  sink(con, type = "message")
  before <- sink.number(type = "message")
  set.seed(1)
  invisible(RMdifTree(d$items, covariates = d$covs, minsize = 50,
                      stability = TRUE, stability_B = 5L,
                      output = "dataframe"))
  expect_equal(sink.number(type = "message"), before)
})

# ---------------------------------------------------------------------
# Covariate-name propagation
# ---------------------------------------------------------------------
test_that("Variable column reflects user's expression when covariates is passed as dif$col", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  df <- RMdifTree(d$items, d$covs$gender,
                  minsize = 50, output = "dataframe")
  expect_true("gender" %in% df$Variable)
})

test_that("Variable column reflects user's expression when covariates is a bare symbol", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  my_gender <- d$covs$gender
  set.seed(1)
  df <- RMdifTree(d$items, my_gender,
                  minsize = 50, output = "dataframe")
  expect_true("my_gender" %in% df$Variable)
})

test_that("Existing data.frame column names are preserved", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  df <- RMdifTree(d$items, covariates = d$covs,
                  minsize = 50, output = "dataframe")
  expect_true(all(df$Variable %in% names(d$covs)))
})

test_that("single covariate passed by numeric index (non-syntactic name) works", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  # d$covs[, 2] drops to a vector whose derived name is the non-syntactic
  # deparse "d$covs[, 2]"; the formula term must be matched as a literal
  # column name, not re-evaluated as an expression.
  df <- RMdifTree(d$items, covariates = d$covs[, 2],
                  minsize = 50, output = "dataframe")
  expect_s3_class(df, "data.frame")
  expect_true("d$covs[, 2]" %in% df$Variable)
})

# ---------------------------------------------------------------------
# Iterative purification, kable/plot output, pruning
# ---------------------------------------------------------------------

# Polytomous data with strong group DIF on I1/I2 so pctree() splits on the
# group and the partial-gamma effect sizes flag those items (exercising the
# purification loop).
make_dif_poly <- function(n = 400, seed = 11L) {
  set.seed(seed)
  grp   <- factor(rep(c("A", "B"), length.out = n))
  theta <- stats::rnorm(n)
  mk <- function(shift) {
    lin <- theta + ifelse(grp == "B", shift, 0)
    pmax(0L, pmin(2L, as.integer(round(lin + stats::rnorm(n, 0, 0.4) + 1))))
  }
  list(
    items = data.frame(I1 = mk(1.6), I2 = mk(1.6), I3 = mk(0),
                       I4 = mk(0), I5 = mk(0), I6 = mk(0)),
    covs  = data.frame(grp = grp)
  )
}

test_that("purification = 'iterative' (MH) runs on dichotomous DIF data", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  df <- RMdifTree(d$items, covariates = d$covs, minsize = 50,
                  purification = "iterative", output = "dataframe")
  expect_s3_class(df, "data.frame")
  expect_equal(attr(df, "effect_size"), "MH")
  expect_true(all(df$Class %in% c("A", "B", "C")))
  expect_true(any(df$Flagged))
})

test_that("purification = 'iterative' (partial gamma) with p_adj on polytomous data", {
  need_tree_pkgs()
  d <- make_dif_poly()
  set.seed(1)
  df <- RMdifTree(d$items, covariates = d$covs, minsize = 50,
                  purification = "iterative", p_adj = "fdr",
                  output = "dataframe")
  expect_s3_class(df, "data.frame")
  expect_equal(attr(df, "effect_size"), "pgamma")
  expect_true(all(df$Class %in% c("A", "B", "C")))
})

test_that("default kable output renders per-split effect sizes", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  kbl <- RMdifTree(d$items, covariates = d$covs, minsize = 50,
                   p_adj = "bonferroni")
  expect_s3_class(kbl, "knitr_kable")
})

test_that("output = 'plot' draws the tree without error", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(
    RMdifTree(d$items, covariates = d$covs, minsize = 50, output = "plot")
  )
})

test_that("prune_negligible = TRUE runs and returns a tree", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  set.seed(1)
  tr <- RMdifTree(d$items, covariates = d$covs, minsize = 50,
                  prune_negligible = TRUE, output = "tree")
  expect_s3_class(tr, "RMdifTree")
})

test_that("effect_size = 'MH' errors on polytomous data", {
  need_tree_pkgs()
  d <- make_dif_poly()
  expect_error(
    RMdifTree(d$items, covariates = d$covs, effect_size = "MH", minsize = 50),
    regexp = "dichotomous"
  )
})

test_that("invalid thresholds, alpha, and stability_B are rejected", {
  need_tree_pkgs()
  d <- make_dif_dichotomous()
  expect_error(
    RMdifTree(d$items, covariates = d$covs, thresholds = 0.5),
    regexp = "thresholds"
  )
  expect_error(
    RMdifTree(d$items, covariates = d$covs, alpha = 1.5),
    regexp = "alpha"
  )
  expect_error(
    RMdifTree(d$items, covariates = d$covs, stability_B = 0),
    regexp = "stability_B"
  )
})
