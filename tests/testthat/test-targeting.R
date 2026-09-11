# Tests for RMtargeting()

make_polytomous <- function(n = 200, k = 6, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_dichotomous <- function(n = 200, k = 8, seed = 2L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMtargeting errors when data has non-zero minimum", {
  df <- make_dichotomous() + 1L
  expect_error(RMtargeting(df), regexp = "scored starting at 0")
})

# ---------------------------------------------------------------------
# Output structures
# ---------------------------------------------------------------------
test_that("RMtargeting default returns a patchwork/ggplot figure on polytomous data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  p  <- RMtargeting(df)
  # patchwork-composed plots also inherit ggplot
  expect_s3_class(p, "ggplot")
})

test_that("RMtargeting default returns a figure on dichotomous data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("eRm")
  df <- make_dichotomous()
  p  <- RMtargeting(df)
  expect_s3_class(p, "ggplot")
})

test_that("RMtargeting output = 'list' returns sub-plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  res <- RMtargeting(df, output = "list")
  expect_type(res, "list")
  expect_true(length(res) >= 1L)
})

test_that("RMtargeting accepts robust = TRUE and sort_items = 'location'", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  p <- RMtargeting(df, robust = TRUE, sort_items = "location",
                   ci_level = 0.84)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------
# Bottom panel selection
# ---------------------------------------------------------------------
test_that("RMtargeting panel = 'thresholds' returns a figure", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  df <- make_polytomous()
  p <- RMtargeting(df, panel = "thresholds")
  expect_s3_class(p, "ggplot")
})

test_that("RMtargeting category panel accepts labels and rejects wrong length", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  df <- make_polytomous()
  p <- RMtargeting(df, category_labels = c("None", "Some", "Most", "All"))
  expect_s3_class(p, "ggplot")
  expect_error(
    RMtargeting(df, category_labels = c("None", "Some")),
    regexp = "one label per response category"
  )
})

test_that("RMtargeting category panel works on dichotomous data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  df <- make_dichotomous()
  expect_s3_class(RMtargeting(df), "ggplot")
})

# ---------------------------------------------------------------------
# Modal-category boundaries
# ---------------------------------------------------------------------
test_that("modal boundaries reduce to the thresholds when they are ordered", {
  tau <- c(-2, -1, 0.6)
  mb <- .modal_boundaries(tau)
  expect_equal(mb$bounds, tau)
  expect_equal(mb$cats, 0:3)
  expect_length(mb$dropped, 0L)
})

test_that("a reversed threshold pair drops the category between them", {
  # tau2 < tau1, so category 1 is never the most likely response
  mb <- .modal_boundaries(c(-1.2, -1.6, 0.4))
  expect_equal(mb$dropped, 1L)
  expect_equal(mb$cats, c(0L, 2L, 3L))
  # the pooled boundary is the mean of the two reversed thresholds
  expect_equal(mb$bounds[1], mean(c(-1.2, -1.6)))
})

test_that("a run of three reversed thresholds drops two categories", {
  mb <- .modal_boundaries(c(-2.2, 0.6, 0.0, -0.6, 1.2, 2.2))
  expect_equal(mb$dropped, c(2L, 3L))
})

test_that("threshold reversals are found per item and NULL when ordered", {
  thr <- data.frame(
    Item = rep(c("a", "b"), each = 3L),
    k = rep(1:3, 2L),
    Location = c(-1, 0, 1, -1, 0.5, 0.2)
  )
  rev <- .threshold_reversals(thr, c("a", "b"))
  expect_equal(rev$Item, "b")
  expect_equal(rev$k, 2L)
  expect_equal(rev$gap, 0.3)
  expect_null(.threshold_reversals(thr[thr$Item == "a", ], "a"))
})

test_that("pale category fills are darkened for the error bars", {
  lum <- function(x) sum(grDevices::col2rgb(x)[, 1] * c(0.299, 0.587, 0.114))
  expect_lt(lum(.ci_contrast("#A8E1BC")), lum("#A8E1BC"))
  expect_equal(.ci_contrast("#382A54"), "#382A54")
})
