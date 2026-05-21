# Tests for RMciccPlot()

# Polytomous data is required to make class-interval and DIF modes work
# reliably at modest n. iarm::ICCplot errors when class_intervals is too
# large relative to the number of populated total-score levels, which
# happens easily with dichotomous data and small samples.
make_polytomous <- function(n = 200, k = 5, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMciccPlot errors when data has non-zero minimum", {
  df <- make_polytomous() + 1L
  expect_error(RMciccPlot(df), regexp = "scored starting at 0")
})

test_that("RMciccPlot errors when fewer than 2 items", {
  df <- make_polytomous()[, 1L, drop = FALSE]
  expect_error(RMciccPlot(df), regexp = "at least 2 items")
})

test_that("RMciccPlot errors when class_intervals < 2", {
  df <- make_polytomous()
  expect_error(RMciccPlot(df, class_intervals = 1),
               regexp = "class_intervals")
})

test_that("RMciccPlot errors when dif_var length mismatches nrow(data)", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df  <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df) - 1L))
  expect_error(RMciccPlot(df, dif_var = grp),
               regexp = "same length")
})

test_that("RMciccPlot errors when dif_var has fewer than 2 levels", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df  <- make_polytomous()
  grp <- factor(rep("A", nrow(df)))
  expect_error(RMciccPlot(df, dif_var = grp), regexp = "at least 2 distinct")
})

# ---------------------------------------------------------------------
# Output structures
# ---------------------------------------------------------------------
test_that("RMciccPlot default output is a patchwork/ggplot composite", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMciccPlot(df, class_intervals = 4L)
  expect_s3_class(p, "patchwork")
  expect_s3_class(p, "ggplot")
})

test_that("RMciccPlot output = 'list' returns one ggplot per item, named after the items", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  plist <- RMciccPlot(df, class_intervals = 4L, output = "list")
  expect_type(plist, "list")
  expect_equal(length(plist), ncol(df))
  expect_equal(names(plist), names(df))
  for (p in plist) expect_s3_class(p, "ggplot")
})

test_that("RMciccPlot DIF mode returns a composite", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df  <- make_polytomous(n = 300)
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  p <- RMciccPlot(df, class_intervals = 4L, dif_var = grp)
  expect_s3_class(p, "patchwork")
})

test_that("RMciccPlot method = 'score' returns a composite", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMciccPlot(df, method = "score")
  expect_s3_class(p, "patchwork")
})

# ---------------------------------------------------------------------
# NA handling
# ---------------------------------------------------------------------
test_that("RMciccPlot drops rows with NA in data or dif_var", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous(n = 250)
  df[1:5, 1] <- NA_integer_         # NA in items
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  grp[6:10] <- NA                   # NA in dif_var
  # Should not error: rows dropped jointly
  expect_s3_class(
    RMciccPlot(df, class_intervals = 4L, dif_var = grp),
    "patchwork"
  )
})
