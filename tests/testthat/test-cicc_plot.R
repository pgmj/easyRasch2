# Tests for RMitemICCPlot()

# Polytomous data gives enough populated total-score levels for class
# intervals and DIF groups at modest n.
make_polytomous <- function(n = 200, k = 5, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation (reachable without the plotting packages)
# ---------------------------------------------------------------------
test_that("RMitemICCPlot errors when data has non-zero minimum", {
  df <- make_polytomous() + 1L
  expect_error(RMitemICCPlot(df), regexp = "scored starting at 0")
})

test_that("RMitemICCPlot errors when fewer than 2 items", {
  df <- make_polytomous()[, 1L, drop = FALSE]
  expect_error(RMitemICCPlot(df), regexp = "at least 2 items")
})

test_that("RMitemICCPlot errors when class_intervals < 2", {
  df <- make_polytomous()
  expect_error(RMitemICCPlot(df, class_intervals = 1),
               regexp = "class_intervals")
})

test_that("RMitemICCPlot errors when conf_level is out of range", {
  df <- make_polytomous()
  expect_error(RMitemICCPlot(df, conf_level = 1.2), regexp = "conf_level")
})

test_that("RMitemICCPlot errors when dif_var length mismatches nrow(data)", {
  df  <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df) - 1L))
  expect_error(RMitemICCPlot(df, dif_var = grp), regexp = "same length")
})

test_that("RMitemICCPlot errors when dif_var has fewer than 2 levels", {
  df  <- make_polytomous()
  grp <- factor(rep("A", nrow(df)))
  expect_error(RMitemICCPlot(df, dif_var = grp), regexp = "at least 2 distinct")
})

# ---------------------------------------------------------------------
# Output structures (no iarm needed for the non-DIF CML plot)
# ---------------------------------------------------------------------
test_that("RMitemICCPlot default output is a patchwork/ggplot composite", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMitemICCPlot(df, class_intervals = 4L)
  expect_s3_class(p, "patchwork")
  expect_s3_class(p, "ggplot")
})

test_that("RMitemICCPlot output = 'list' returns one ggplot per item, named", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  plist <- RMitemICCPlot(df, class_intervals = 4L, output = "list")
  expect_type(plist, "list")
  expect_equal(length(plist), ncol(df))
  expect_equal(names(plist), names(df))
  for (p in plist) expect_s3_class(p, "ggplot")
})

test_that("RMitemICCPlot method = 'score' and error_band return a composite", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  expect_s3_class(RMitemICCPlot(df, method = "score"), "patchwork")
  expect_s3_class(RMitemICCPlot(df, error_band = TRUE, ci = FALSE), "patchwork")
})

test_that("RMitemICCPlot `items` subsets the rendered panels", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  plist <- RMitemICCPlot(df, items = c("I2", "I4"), output = "list")
  expect_equal(names(plist), c("I2", "I4"))
})

# ---------------------------------------------------------------------
# DIF mode (needs iarm for the partial-gamma magnitude)
# ---------------------------------------------------------------------
test_that("RMitemICCPlot DIF mode returns a composite with a gamma table attr", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df  <- make_polytomous(n = 300)
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  p <- RMitemICCPlot(df, class_intervals = 4L, dif_var = grp)
  expect_s3_class(p, "patchwork")
  g <- attr(p, "dif_gamma")
  expect_s3_class(g, "data.frame")
  expect_true(all(c("Item", "gamma", "lower", "upper", "padj_bh") %in% names(g)))
})

# ---------------------------------------------------------------------
# NA handling
# ---------------------------------------------------------------------
test_that("RMitemICCPlot drops rows with NA in data or dif_var", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous(n = 300)
  df[1:5, 1] <- NA_integer_
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  grp[6:10] <- NA
  expect_s3_class(
    RMitemICCPlot(df, class_intervals = 4L, dif_var = grp), "patchwork"
  )
})

# ---------------------------------------------------------------------
# Grouping methods (quantile / width / score / manual, "cut" alias)
# ---------------------------------------------------------------------
test_that("all grouping methods run and 'cut' aliases 'quantile'", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  set.seed(42)
  df <- as.data.frame(matrix(sample(0:2, 150 * 4, replace = TRUE), nrow = 150))
  colnames(df) <- paste0("i", 1:4)

  for (m in c("quantile", "width", "score")) {
    p <- suppressWarnings(RMitemICCPlot(df, method = m))
    expect_s3_class(p, "patchwork")
  }
  # legacy alias gives the same grouping as "quantile"
  p_cut <- suppressWarnings(RMitemICCPlot(df, method = "cut", output = "list"))
  p_q   <- suppressWarnings(RMitemICCPlot(df, method = "quantile", output = "list"))
  expect_identical(p_cut[[1]]$data, p_q[[1]]$data)
})

test_that("manual grouping via score_breaks works and validates", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  set.seed(42)
  df <- as.data.frame(matrix(sample(0:2, 150 * 4, replace = TRUE), nrow = 150))
  colnames(df) <- paste0("i", 1:4)

  p <- suppressWarnings(
    RMitemICCPlot(df, method = "manual", score_breaks = c(3, 5, 7)))
  expect_s3_class(p, "patchwork")

  # validation
  expect_error(RMitemICCPlot(df, method = "manual"),
               "requires `score_breaks`")
  expect_error(RMitemICCPlot(df, method = "quantile", score_breaks = c(3)),
               "only used with method")
  expect_error(RMitemICCPlot(df, method = "manual", score_breaks = c(5, 3)),
               "strictly increasing")
  expect_error(RMitemICCPlot(df, method = "manual", score_breaks = c(3, 99)),
               "maximum possible total")
})

# ---------------------------------------------------------------------
# Class-interval caption (reports the grouping actually used)
# ---------------------------------------------------------------------
test_that(".cicc_interval_caption names each method and its group count", {
  rs <- c(rep(0:8, each = 10))
  q <- .cicc_class_bins(rs, "quantile", 4L)
  expect_match(.cicc_interval_caption("quantile", 4L, NULL, q),
               "^Class intervals: 4 groups of approximately equal size")

  w <- .cicc_class_bins(rs, "width", 4L)
  expect_match(.cicc_interval_caption("width", 4L, NULL, w),
               "^Class intervals: 4 equal-width groups")

  s <- .cicc_class_bins(rs, "score", 4L)
  expect_match(.cicc_interval_caption("score", 4L, NULL, s),
               "every observed total score is its own group \\(9 groups\\)")

  m <- .cicc_class_bins(rs, "manual", 4L, score_breaks = c(3L, 6L))
  expect_match(.cicc_interval_caption("manual", 4L, c(3L, 6L), m),
               "a new group starting at total scores 3, 6 \\(3 groups\\)")
  # singular when there is one break
  m1 <- .cicc_class_bins(rs, "manual", 4L, score_breaks = 4L)
  expect_match(.cicc_interval_caption("manual", 4L, 4L, m1),
               "at total score 4 ")
})

test_that(".cicc_interval_caption reports fewer quantile groups than requested", {
  # 9 distinct total scores cannot support 12 quantile bins
  rs <- rep(0:8, each = 10)
  b <- .cicc_class_bins(rs, "quantile", 12L)
  expect_lt(nlevels(b), 12L)
  expect_match(.cicc_interval_caption("quantile", 12L, NULL, b),
               "12 were requested, fewer were formed")
})

test_that(".cicc_interval_caption reports groups that no respondent falls into", {
  # equal-width bins over a range with only two occupied scores
  rs <- rep(c(4L, 5L), each = 50)
  b <- .cicc_class_bins(rs, "width", 6L)
  expect_equal(nlevels(b), 6L)
  cap <- .cicc_interval_caption("width", 6L, NULL, b)
  expect_match(cap, "6 equal-width groups")
  expect_match(cap, "2 of them contain respondents")
})

test_that(".cicc_class_bins flags the fall back to score level", {
  # Quantile grouping collapses when every respondent has the same total
  # score, since all the quantiles coincide.
  rs <- rep(4L, 50)
  b <- .cicc_class_bins(rs, "quantile", 4L)
  expect_true(isTRUE(attr(b, "fallback")))
  expect_match(.cicc_interval_caption("quantile", 4L, NULL, b),
               "too few distinct total scores")

  # Equal-width grouping pads the range by half a score either side, so it
  # still forms bins on a single distinct score and only collapses when one
  # interval is requested.
  expect_null(attr(.cicc_class_bins(rs, "width", 4L), "fallback"))
  expect_true(isTRUE(attr(.cicc_class_bins(rs, "width", 1L), "fallback")))

  # no flag when the grouping succeeded
  ok <- .cicc_class_bins(rep(0:8, each = 10), "width", 4L)
  expect_null(attr(ok, "fallback"))
})

test_that("the plot caption carries the class-interval clause", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  set.seed(11)
  n <- 250
  th <- stats::rnorm(n)
  d <- as.data.frame(sapply(1:8, function(j)
    stats::rbinom(n, 1, stats::plogis(th - (j - 4.5) / 3))))
  names(d) <- paste0("I", 1:8)
  p <- RMitemICCPlot(d, method = "width", class_intervals = 5, items = 1:2)
  # er2_caption() hard-wraps the note, so allow a line break inside it
  expect_match(p$patches$annotation$caption,
               "Class\\s+intervals:\\s+5\\s+equal-width")
})
