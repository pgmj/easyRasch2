# Tests for RMitemRestscoreCutoff(), the cutoff path of RMitemRestscore(),
# and RMitemRestscorePlot().

make_dich <- function(n = 200, k = 6, seed = 3L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_poly <- function(n = 200, k = 5, seed = 3L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:2, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# RMitemRestscoreCutoff()
# ---------------------------------------------------------------------
test_that("returns the documented structure", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  res <- RMitemRestscoreCutoff(make_dich(), iterations = 12, parallel = FALSE,
                               seed = 1, cutoff_method = "quantile")
  expect_named(res, c(
    "results", "item_cutoffs", "actual_iterations", "requested_iterations",
    "sample_n", "sample_n_total", "sample_has_na", "sample_summary",
    "item_names", "cutoff_method", "hdci_width", "dgp"
  ))
  expect_named(res$results,
               c("iteration", "Item", "Observed", "Expected", "Difference"))
  expect_named(res$item_cutoffs, c("Item", "diff_low", "diff_high"))
  expect_equal(res$actual_iterations, 12L)
  expect_equal(nrow(res$results), 12L * 6L)
  expect_equal(res$results$Difference,
               res$results$Observed - res$results$Expected)
  expect_true(all(res$item_cutoffs$diff_low < res$item_cutoffs$diff_high))
})

test_that("the expected gamma varies across iterations (thresholds refitted)", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  res <- RMitemRestscoreCutoff(make_dich(), iterations = 10, parallel = FALSE,
                               seed = 1, cutoff_method = "quantile")
  sds <- tapply(res$results$Expected, res$results$Item, stats::sd)
  expect_true(all(sds > 0))
})

test_that("the conditional DGP is the default and both DGPs run on polytomous data", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  res <- RMitemRestscoreCutoff(make_poly(), iterations = 8, parallel = FALSE,
                               seed = 1, cutoff_method = "quantile")
  expect_identical(res$dgp, "conditional")
  expect_equal(nrow(res$item_cutoffs), 5L)
  res_r <- RMitemRestscoreCutoff(make_poly(), iterations = 8, parallel = FALSE,
                                 seed = 1, dgp = "resample",
                                 cutoff_method = "quantile")
  expect_identical(res_r$dgp, "resample")
  expect_equal(nrow(res_r$item_cutoffs), 5L)
})

test_that("uses complete cases and records the raw sample size", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  d <- make_dich()
  d[1:5, 2] <- NA
  res <- RMitemRestscoreCutoff(d, iterations = 5, parallel = FALSE, seed = 1,
                               cutoff_method = "quantile")
  expect_equal(res$sample_n, 195L)
  expect_equal(res$sample_n_total, 200L)
  expect_true(res$sample_has_na)
})

test_that("the same seed reproduces the simulation", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  a <- RMitemRestscoreCutoff(make_dich(), iterations = 5, parallel = FALSE,
                             seed = 9, cutoff_method = "quantile")
  b <- RMitemRestscoreCutoff(make_dich(), iterations = 5, parallel = FALSE,
                             seed = 9, cutoff_method = "quantile")
  expect_identical(a$results, b$results)
})

test_that("parallel and sequential runs agree for the same seed", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("mirai")
  s <- RMitemRestscoreCutoff(make_dich(), iterations = 6, parallel = FALSE,
                             seed = 1, cutoff_method = "quantile")
  p <- suppressMessages(
    RMitemRestscoreCutoff(make_dich(), iterations = 6, parallel = TRUE,
                          n_cores = 2, seed = 1, cutoff_method = "quantile")
  )
  expect_identical(s$results, p$results)
})

test_that("errors when no complete cases remain", {
  skip_if_not_installed("iarm")
  d <- data.frame(a = c(0, NA, 1), b = c(NA, 1, NA))
  expect_error(
    RMitemRestscoreCutoff(d, iterations = 3, parallel = FALSE,
                          cutoff_method = "quantile"),
    "No complete cases"
  )
})

# ---------------------------------------------------------------------
# RMitemRestscore() with a cutoff
# ---------------------------------------------------------------------
test_that("the full cutoff object flags on bootstrap p-values by default", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  d <- make_dich()
  co <- RMitemRestscoreCutoff(d, iterations = 20, parallel = FALSE, seed = 1,
                              cutoff_method = "quantile")
  res <- suppressMessages(RMitemRestscore(d, cutoff = co,
                                          output = "dataframe"))
  expect_named(res, c(
    "Item", "Observed", "Expected", "Difference", "Diff_low", "Diff_high",
    "p_restscore", "padj_restscore", "Flagged", "Relative_location"
  ))
  expect_false("p_adjusted" %in% names(res))
  expect_true(all(res$p_restscore >= 1 / 21))
  expect_true(all(res$padj_restscore >= res$p_restscore))
  expect_true(all(res$Flagged %in% c("", "overfit", "underfit")))
})

test_that("p_value = FALSE and the bare item_cutoffs flag on the interval", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  d <- make_dich()
  co <- RMitemRestscoreCutoff(d, iterations = 10, parallel = FALSE, seed = 1,
                              cutoff_method = "quantile")
  a <- suppressMessages(RMitemRestscore(d, cutoff = co, p_value = FALSE,
                                        output = "dataframe"))
  b <- suppressMessages(RMitemRestscore(d, cutoff = co$item_cutoffs,
                                        output = "dataframe"))
  expect_identical(a, b)
  expect_false(any(c("p_restscore", "padj_restscore") %in% names(a)))
  expected_flag <- ifelse(a$Difference > a$Diff_high, "overfit",
                          ifelse(a$Difference < a$Diff_low, "underfit", ""))
  expect_identical(a$Flagged, expected_flag)
})

test_that("a planted underfitting item is flagged as underfit", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  set.seed(11)
  n <- 400
  th <- stats::rnorm(n)
  beta <- seq(-1, 1, length.out = 6)
  d <- as.data.frame(sapply(beta, function(b) {
    stats::rbinom(n, 1, stats::plogis(th - b))
  }))
  d[[6]] <- stats::rbinom(n, 1, 0.5) # unrelated to the trait
  colnames(d) <- paste0("I", 1:6)
  co <- RMitemRestscoreCutoff(d, iterations = 40, parallel = FALSE, seed = 1,
                              cutoff_method = "quantile")
  res <- suppressMessages(RMitemRestscore(d, cutoff = co,
                                          output = "dataframe"))
  expect_identical(res$Flagged[res$Item == "I6"], "underfit")
})

test_that("p_value = TRUE without the full object errors", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  d <- make_dich()
  co <- RMitemRestscoreCutoff(d, iterations = 5, parallel = FALSE, seed = 1,
                              cutoff_method = "quantile")
  expect_error(
    RMitemRestscore(d, cutoff = co$item_cutoffs, p_value = TRUE),
    "requires the full RMitemRestscoreCutoff"
  )
  expect_error(RMitemRestscore(d, p_value = TRUE),
               "requires the full RMitemRestscoreCutoff")
})

test_that("p_adj with a cutoff warns that it is ignored", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  d <- make_dich()
  co <- RMitemRestscoreCutoff(d, iterations = 5, parallel = FALSE, seed = 1,
                              cutoff_method = "quantile")
  expect_warning(
    suppressMessages(RMitemRestscore(d, cutoff = co, p_adj = "none",
                                     output = "dataframe")),
    "`p_adj` is ignored"
  )
})

test_that("mismatched item names error", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  d <- make_dich()
  co <- RMitemRestscoreCutoff(d, iterations = 5, parallel = FALSE, seed = 1,
                              cutoff_method = "quantile")
  d2 <- d
  colnames(d2)[1] <- "other"
  expect_error(suppressMessages(RMitemRestscore(d2, cutoff = co)),
               "do not match")
})

test_that("kable output has a caption naming the correction", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  d <- make_dich()
  co <- RMitemRestscoreCutoff(d, iterations = 10, parallel = FALSE, seed = 1,
                              cutoff_method = "quantile")
  kb <- suppressMessages(RMitemRestscore(d, cutoff = co))
  expect_s3_class(kb, "knitr_kable")
  expect_match(paste(kb, collapse = "\n"), "Westfall-Young", fixed = TRUE)
})

test_that("invalid p_value and alpha error", {
  skip_if_not_installed("iarm")
  d <- make_dich()
  expect_error(RMitemRestscore(d, p_value = "yes"), "must be NULL, TRUE")
  expect_error(RMitemRestscore(d, alpha = 2), "single number in")
})

# ---------------------------------------------------------------------
# RMitemRestscorePlot()
# ---------------------------------------------------------------------
test_that("RMitemRestscorePlot returns a ggplot with and without data", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  d <- make_dich()
  co <- RMitemRestscoreCutoff(d, iterations = 10, parallel = FALSE, seed = 1)
  expect_s3_class(RMitemRestscorePlot(co), "ggplot")
  expect_s3_class(RMitemRestscorePlot(co, d), "ggplot")
})

test_that("RMitemRestscorePlot rejects a non-restscore simfit", {
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  expect_error(RMitemRestscorePlot(list(results = data.frame(x = 1))),
               "not an RMitemRestscoreCutoff")
})
