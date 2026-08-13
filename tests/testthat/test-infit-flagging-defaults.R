# Tests for the 1.2.0 flagging defaults: `p_value = NULL` resolution, the
# descriptive interval width, and the two advisory notices.
#
# The notices use rlang's once-per-session frequency, so each test that
# expects one resets the counter first.

make_dich_fit <- function(n = 200, k = 6, seed = 5L) {
  set.seed(seed)
  th <- stats::rnorm(n)
  df <- as.data.frame(sapply(
    seq_len(k),
    function(j) stats::rbinom(n, 1, stats::plogis(th - (j - (k + 1) / 2) / 3))
  ))
  names(df) <- paste0("I", seq_len(k))
  df
}

quiet_cutoff <- function(df, ...) {
  suppressMessages(RMitemInfitCutoff(df, parallel = FALSE, seed = 1L, ...))
}

reset_notices <- function() {
  rlang::reset_message_verbosity("easyRasch2_band_flagging")
  rlang::reset_message_verbosity("easyRasch2_low_iterations")
}

# ---------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------
test_that("RMitemInfitCutoff defaults to a descriptive 95% interval", {
  expect_equal(formals(RMitemInfitCutoff)$hdci_width, 0.95)
})

test_that("RMitemInfit defaults p_value to NULL", {
  expect_null(formals(RMitemInfit)$p_value)
})

# ---------------------------------------------------------------------
# p_value = NULL resolution
# ---------------------------------------------------------------------
test_that("p_value = NULL uses p-values when the full cutoff object is given", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  co <- quiet_cutoff(df, iterations = 400L)
  reset_notices()
  res <- RMitemInfit(df, cutoff = co, output = "dataframe")
  expect_true(all(c("p_infit", "padj_infit") %in% names(res)))
  # and no band-flagging notice on this path
  reset_notices()
  expect_no_message(RMitemInfit(df, cutoff = co, output = "dataframe"))
})

test_that("p_value = NULL falls back to the band without simulations", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  co <- quiet_cutoff(df, iterations = 400L)
  reset_notices()
  res <- suppressMessages(
    RMitemInfit(df, cutoff = co$item_cutoffs, output = "dataframe")
  )
  expect_false(any(c("p_infit", "padj_infit") %in% names(res)))
  expect_true("Flagged" %in% names(res))
})

test_that("p_value = NULL with no cutoff computes nothing and says nothing", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  reset_notices()
  expect_no_message(res <- RMitemInfit(df, output = "dataframe"))
  expect_equal(names(res), c("Item", "Infit_MSQ", "Relative_location"))
})

test_that("explicit p_value = TRUE without simulations is still an error", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  co <- quiet_cutoff(df, iterations = 400L)
  expect_error(
    RMitemInfit(df, cutoff = co$item_cutoffs, p_value = TRUE),
    regexp = "requires the full"
  )
})

test_that("p_value rejects non-logical input", {
  df <- make_dich_fit()
  expect_error(RMitemInfit(df, p_value = "yes"), regexp = "must be NULL")
  expect_error(RMitemInfit(df, p_value = NA), regexp = "must be NULL")
  expect_error(RMitemInfit(df, p_value = c(TRUE, TRUE)), regexp = "must be NULL")
})

# ---------------------------------------------------------------------
# Advisory notices
# ---------------------------------------------------------------------
test_that("band flagging reports the family-wise error rate it implies", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit(k = 6L)
  co <- quiet_cutoff(df, iterations = 400L)
  reset_notices()
  # 1 - .95^6 = 26.5%, so the message should say 26%
  expect_message(
    RMitemInfit(df, cutoff = co, p_value = FALSE, output = "dataframe"),
    regexp = "family-wise error rate of about 26%"
  )
  # once per session
  expect_no_message(
    RMitemInfit(df, cutoff = co, p_value = FALSE, output = "dataframe")
  )
})

test_that("band flagging degrades gracefully when the width is unknown", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  co <- quiet_cutoff(df, iterations = 400L)
  reset_notices()
  expect_message(
    RMitemInfit(df, cutoff = co$item_cutoffs, output = "dataframe"),
    regexp = "1 - width\\^k"
  )
})

test_that("fewer than 400 iterations triggers the calibration notice", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  co <- quiet_cutoff(df, iterations = 200L)
  reset_notices()
  expect_message(
    RMitemInfit(df, cutoff = co, output = "dataframe"),
    regexp = "below the calibrated floor of 400"
  )
})

test_that("400 or more iterations is quiet on the console", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  co <- quiet_cutoff(df, iterations = 400L)
  reset_notices()
  expect_no_message(RMitemInfit(df, cutoff = co, output = "dataframe"))
})

# ---------------------------------------------------------------------
# Captions
# ---------------------------------------------------------------------
test_that("captions name the basis of the flag in both branches", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  co <- quiet_cutoff(df, iterations = 400L)

  reset_notices()
  pv <- paste(as.character(suppressMessages(
    RMitemInfit(df, cutoff = co)
  )), collapse = " ")
  expect_match(pv, "Flagged on the corrected p-value at alpha")
  # 400 is below 1000, so the reproducibility note appears
  expect_match(pv, "seed-dependent")

  reset_notices()
  band <- paste(as.character(suppressMessages(
    RMitemInfit(df, cutoff = co, p_value = FALSE)
  )), collapse = " ")
  expect_match(band, "tests all 6 items at once")
  expect_match(band, "Johansson, 2026")
})

test_that("the reproducibility note drops out at 1000 iterations", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  co <- quiet_cutoff(df, iterations = 1000L)
  reset_notices()
  cap <- paste(as.character(suppressMessages(
    RMitemInfit(df, cutoff = co)
  )), collapse = " ")
  expect_no_match(cap, "seed-dependent")
})

# ---------------------------------------------------------------------
# Sidak helpers
# ---------------------------------------------------------------------
test_that("the Sidak helpers invert each other", {
  for (k in c(5L, 9L, 20L)) {
    w <- .sidak_width(k, alpha = 0.05)
    expect_equal(.sidak_fwe(w, k), 0.05, tolerance = 1e-12)
  }
  expect_equal(.sidak_fwe(0.95, 9L), 1 - 0.95^9)
  expect_equal(.width_iterations(0.95), 400)
  expect_equal(.fmt_width(0.95), ".95")
  expect_equal(.fmt_width(0.99432), ".99432")
})

test_that(".attrition_clause explains lost iterations and stays silent otherwise", {
  msg <- .attrition_clause(actual = 391L, requested = 400L)
  expect_match(msg, "^ 9 of the 400 simulated datasets")
  expect_match(msg, "rest on 391")
  expect_match(msg, "Raise `iterations`")
  # nothing lost, or the requested count unknown (pre-1.2.0 cutoff objects)
  expect_null(.attrition_clause(400L, 400L))
  expect_null(.attrition_clause(400L, NULL))
  expect_null(.attrition_clause(NULL, 400L))
})

test_that("RMitemInfitCutoff records the requested iteration count", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  co <- quiet_cutoff(df, iterations = 400L)
  expect_equal(co$requested_iterations, 400L)
  expect_lte(co$actual_iterations, co$requested_iterations)
})

test_that(".band_error_clause returns NULL when the width is unusable", {
  expect_null(.band_error_clause(NULL, 9L))
  expect_null(.band_error_clause(1, 9L))
  expect_null(.band_error_clause(NA_real_, 9L))
  expect_match(.band_error_clause(0.95, 9L), "about 37%")
})

# ---------------------------------------------------------------------
# Plot follows the cutoff object's width
# ---------------------------------------------------------------------
test_that("RMitemInfitPlot interval follows simfit$hdci_width", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("iarm")
  df <- make_dich_fit()
  narrow <- quiet_cutoff(df, iterations = 400L, hdci_width = 0.95)
  wide <- quiet_cutoff(df, iterations = 400L, hdci_width = 0.999)

  seg <- function(p) {
    d <- p$layers[[2]]$data
    max(d$max_infit_msq) - min(d$min_infit_msq)
  }
  expect_lt(seg(RMitemInfitPlot(narrow, data = df)),
            seg(RMitemInfitPlot(wide, data = df)))
})
