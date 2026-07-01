# Tests for RMdifLR()
#
# Covers:
#   - Input validation (data, dif_var, cutoff, conf, model)
#   - Auto-detect of model (PCM vs RM) by data shape
#   - Output structures: "dataframe", "kable", "ggplot"
#   - Flagged column / MaxDiff
#   - LR test attribute
#   - NA handling in dif_var
#   - sort = TRUE for kable

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
make_dichotomous <- function(n = 80, k = 5, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE),
                             nrow = n, ncol = k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_polytomous <- function(n = 80, k = 5, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:2, n * k, replace = TRUE),
                             nrow = n, ncol = k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMdifLR errors when data has non-zero minimum", {
  df <- make_dichotomous() + 1L  # min = 1
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  expect_error(RMdifLR(df, dif_var = grp), regexp = "scored starting at 0")
})

test_that("RMdifLR errors when dif_var length does not match nrow(data)", {
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df) - 1L))
  expect_error(RMdifLR(df, dif_var = grp),
               regexp = "same length as nrow")
})

test_that("RMdifLR errors when dif_var has fewer than 2 distinct levels", {
  df <- make_dichotomous()
  grp <- factor(rep("A", nrow(df)))
  expect_error(RMdifLR(df, dif_var = grp),
               regexp = "at least 2 distinct")
})

test_that("RMdifLR errors when model = 'RM' on polytomous data", {
  df <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  expect_error(
    RMdifLR(df, dif_var = grp, model = "RM"),
    regexp = "dichotomous"
  )
})

test_that("RMdifLR errors on invalid cutoff or conf", {
  df  <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  expect_error(RMdifLR(df, dif_var = grp, cutoff = -1),
               regexp = "non-negative")
  expect_error(RMdifLR(df, dif_var = grp, cutoff = c(0.1, 0.2)),
               regexp = "non-negative")
  expect_error(RMdifLR(df, dif_var = grp, conf = 0),
               regexp = "in \\(0, 1\\)")
  expect_error(RMdifLR(df, dif_var = grp, conf = 1),
               regexp = "in \\(0, 1\\)")
})

# ---------------------------------------------------------------------
# Output structure (dataframe)
# ---------------------------------------------------------------------
test_that("RMdifLR output = 'dataframe' on dichotomous data picks RM and returns one row per item", {
  skip_if_not_installed("eRm")
  df  <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))

  out <- RMdifLR(df, dif_var = grp, output = "dataframe")

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), ncol(df))                    # one row per item
  expect_true(all(c("Item", "A", "B", "All", "MaxDiff",
                    "SE_A", "SE_B", "SE_All") %in% names(out)))
  expect_equal(as.character(out$Item), names(df))

  lr <- attr(out, "lr_test")
  expect_type(lr, "list")
  expect_true(all(c("LR", "df", "p_value", "model") %in% names(lr)))
  expect_equal(lr$model, "RM")                          # auto-picked
})

test_that("RMdifLR level = 'threshold' on polytomous data returns one row per item × threshold", {
  skip_if_not_installed("eRm")
  df  <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))

  out <- RMdifLR(df, dif_var = grp, level = "threshold",
                 output = "dataframe")

  expect_s3_class(out, "data.frame")
  expect_true(all(c("Item", "Threshold", "A", "B", "All",
                    "MaxDiff") %in% names(out)))
  # 0/1/2 data → 2 thresholds per item
  expect_equal(nrow(out), ncol(df) * 2L)
  expect_equal(attr(out, "lr_test")$model, "PCM")       # auto-picked
})

# ---------------------------------------------------------------------
# Output structure (kable / ggplot)
# ---------------------------------------------------------------------
test_that("RMdifLR output = 'kable' returns a knitr_kable", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("knitr")
  df  <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  out <- RMdifLR(df, dif_var = grp, output = "kable")
  expect_s3_class(out, "knitr_kable")
})

test_that("RMdifLR defaults to kable (consistent with other kable/ggplot functions)", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("knitr")
  df  <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  expect_s3_class(RMdifLR(df, dif_var = grp), "knitr_kable")
})

test_that("RMdifLR output = 'ggplot' returns a ggplot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df  <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  out <- RMdifLR(df, dif_var = grp, output = "ggplot")
  expect_s3_class(out, "ggplot")
})

# ---------------------------------------------------------------------
# Flagged column / MaxDiff
# ---------------------------------------------------------------------
test_that("Flagged column is logical and absent when cutoff = NULL", {
  skip_if_not_installed("eRm")
  df  <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))

  with_cutoff <- RMdifLR(df, dif_var = grp, output = "dataframe",
                         cutoff = 0.5)
  expect_true("Flagged" %in% names(with_cutoff))
  expect_type(with_cutoff$Flagged, "logical")
  expect_true(all(with_cutoff$MaxDiff[with_cutoff$Flagged] > 0.5,
                  na.rm = TRUE))

  no_cutoff <- RMdifLR(df, dif_var = grp, output = "dataframe",
                       cutoff = NULL)
  expect_false("Flagged" %in% names(no_cutoff))
})

# ---------------------------------------------------------------------
# NA handling in dif_var
# ---------------------------------------------------------------------
test_that("RMdifLR drops NA rows in dif_var with a message", {
  skip_if_not_installed("eRm")
  df  <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  grp[c(1L, 5L, 9L)] <- NA

  expect_message(
    out <- RMdifLR(df, dif_var = grp, output = "dataframe"),
    regexp = "3 row\\(s\\) with NA"
  )
  expect_equal(attr(out, "lr_test")$n_persons, nrow(df) - 3L)
})

# ---------------------------------------------------------------------
# Sorting (kable)
# ---------------------------------------------------------------------
test_that("sort = TRUE in dataframe output sorts by MaxDiff descending", {
  skip_if_not_installed("eRm")
  df  <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  out <- RMdifLR(df, dif_var = grp, output = "dataframe", sort = TRUE)
  expect_true(all(diff(out$MaxDiff) <= 0))
})
