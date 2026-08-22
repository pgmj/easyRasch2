# Tests for the 1.2.0 local dependence flagging defaults: `p_value = NULL`
# resolution in RMlocdepQ3() and RMlocdepGamma(), the descriptive interval
# width, the shared iteration default, the two advisory notices, and the
# false-discovery-rate reachability warning.
#
# The notices use rlang's once-per-session frequency, so each test that
# expects one resets the counter first.

make_poly <- function(n = 250, k = 7, seed = 5L) {
  set.seed(seed)
  th <- stats::rnorm(n)
  df <- as.data.frame(sapply(
    seq_len(k),
    function(j) stats::rbinom(n, 2, stats::plogis(th - (j - (k + 1) / 2) / 3))
  ))
  names(df) <- paste0("I", seq_len(k))
  df
}

quiet <- function(expr) suppressMessages(suppressWarnings(expr))

q3_cutoff <- function(df, ...) {
  quiet(RMlocdepQ3Cutoff(df, parallel = FALSE, seed = 1L, ...))
}
gamma_cutoff <- function(df, ...) {
  quiet(RMlocdepGammaCutoff(df, parallel = FALSE, seed = 1L, ...))
}

reset_notices <- function() {
  rlang::reset_message_verbosity("easyRasch2_band_flagging_locdep")
  rlang::reset_message_verbosity("easyRasch2_low_iterations_locdep")
}

# ---------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------
test_that("both cutoff functions default to a descriptive 95% interval", {
  expect_equal(formals(RMlocdepQ3Cutoff)$hdci_width, 0.95)
  expect_equal(formals(RMlocdepGammaCutoff)$hdci_width, 0.95)
})

test_that("all three cutoff functions share one iteration default", {
  expect_equal(formals(RMlocdepQ3Cutoff)$iterations, 400)
  expect_equal(formals(RMlocdepGammaCutoff)$iterations, 400)
  expect_equal(formals(RMitemInfitCutoff)$iterations, 400)
})

test_that("both local dependence functions default p_value to NULL", {
  expect_null(formals(RMlocdepQ3)$p_value)
  expect_null(formals(RMlocdepGamma)$p_value)
})

# ---------------------------------------------------------------------
# p_value = NULL resolution
# ---------------------------------------------------------------------
test_that("p_value = NULL uses p-values when the full cutoff object is given", {
  skip_on_cran()
  df <- make_poly()
  reset_notices()
  cu <- q3_cutoff(df, iterations = 40L)
  res <- quiet(RMlocdepQ3(df, cutoff = cu, output = "dataframe")$pairs)
  expect_true(all(c("p_q3", "padj_q3") %in% names(res)))

  reset_notices()
  cg <- gamma_cutoff(df, iterations = 40L)
  rg <- quiet(RMlocdepGamma(df, cutoff = cg, output = "dataframe")$direction1)
  expect_true(all(c("p_gamma", "padj_gamma") %in% names(rg)))
})

test_that("p_value = NULL falls back to the interval without simulations", {
  skip_on_cran()
  df <- make_poly()
  reset_notices()
  cg <- gamma_cutoff(df, iterations = 40L)
  # The bare $pair_cutoffs carries no simulated distributions.
  rg <- quiet(
    RMlocdepGamma(df, cutoff = cg$pair_cutoffs, output = "dataframe")$direction1
  )
  expect_false(any(c("p_gamma", "padj_gamma") %in% names(rg)))
  expect_true("flagged" %in% names(rg))
})

test_that("p_value = NULL with no cutoff leaves the asymptotic table alone", {
  skip_on_cran()
  df <- make_poly()
  reset_notices()
  rg <- quiet(RMlocdepGamma(df, output = "dataframe")$direction1)
  expect_false(any(c("p_gamma", "padj_gamma") %in% names(rg)))
  expect_true("padj_bh" %in% names(rg))
})

test_that("explicit p_value = TRUE without simulations is still an error", {
  skip_on_cran()
  df <- make_poly()
  expect_error(
    quiet(RMlocdepQ3(df, p_value = TRUE)),
    regexp = "requires the full"
  )
  expect_error(
    quiet(RMlocdepGamma(df, p_value = TRUE)),
    regexp = "requires the full"
  )
})

test_that("p_value rejects non-logical input", {
  skip_on_cran()
  df <- make_poly()
  expect_error(quiet(RMlocdepQ3(df, p_value = "yes")), regexp = "TRUE, FALSE")
  expect_error(
    quiet(RMlocdepGamma(df, p_value = c(TRUE, TRUE))),
    regexp = "TRUE, FALSE"
  )
})

# ---------------------------------------------------------------------
# Notices
# ---------------------------------------------------------------------
test_that("interval flagging reports the rate it implies over PAIRS", {
  skip_on_cran()
  df <- make_poly(k = 7L) # 21 pairs, not 7 items
  cu <- q3_cutoff(df, iterations = 40L)
  reset_notices()
  expect_message(
    suppressWarnings(
      RMlocdepQ3(df, cutoff = cu, p_value = FALSE, output = "dataframe")
    ),
    regexp = "21 item pairs"
  )
})

test_that("fewer than 400 iterations triggers the calibration notice", {
  skip_on_cran()
  df <- make_poly()
  cu <- q3_cutoff(df, iterations = 40L)
  reset_notices()
  expect_message(
    suppressWarnings(RMlocdepQ3(df, cutoff = cu, output = "dataframe")),
    regexp = "below the calibrated floor of 400"
  )
})

test_that("the notice names the cutoff function that produced the object", {
  skip_on_cran()
  df <- make_poly()
  cg <- gamma_cutoff(df, iterations = 40L)
  reset_notices()
  expect_message(
    suppressWarnings(RMlocdepGamma(df, cutoff = cg, output = "dataframe")),
    regexp = "RMlocdepGammaCutoff\\(\\)"
  )
})

# ---------------------------------------------------------------------
# Captions
# ---------------------------------------------------------------------
test_that("captions name the basis of the flag in both branches", {
  skip_on_cran()
  df <- make_poly()
  cu <- q3_cutoff(df, iterations = 40L)
  flat <- function(x) paste(as.character(x$pairs), collapse = " ")

  reset_notices()
  pv <- flat(quiet(RMlocdepQ3(df, cutoff = cu)))
  expect_match(pv, "bootstrap p-values", fixed = TRUE)
  expect_match(pv, "not the decision rule", fixed = TRUE)
  # 40 is below 400, so the calibration sentence appears
  expect_match(pv, "below the calibrated floor", fixed = TRUE)

  reset_notices()
  bd <- flat(quiet(RMlocdepQ3(df, cutoff = cu, p_value = FALSE)))
  expect_match(bd, "above the upper bound", fixed = TRUE)
  expect_match(bd, "family-wise error rate", fixed = TRUE)
  expect_match(bd, "21 item pairs", fixed = TRUE)
})

# ---------------------------------------------------------------------
# Cutoff objects
# ---------------------------------------------------------------------
test_that("both cutoff objects record the requested iteration count", {
  skip_on_cran()
  df <- make_poly()
  expect_equal(q3_cutoff(df, iterations = 40L)$requested_iterations, 40L)
  expect_equal(gamma_cutoff(df, iterations = 40L)$requested_iterations, 40L)
})

# ---------------------------------------------------------------------
# False discovery rate reachability
# ---------------------------------------------------------------------
test_that(".fdr_min_iterations matches p.adjust at the boundary", {
  # The closed form is m/(s(B+1)) for BH and c(m)*m/(s(B+1)) for BY, so with
  # one p-value at the Monte Carlo floor the smallest reachable adjusted value
  # crosses alpha exactly between B - 1 and B. Pinned here rather than left in
  # a comment, because the earlier `B >= m/alpha - 1` form was one short.
  smallest <- function(m, B, method) {
    min(stats::p.adjust(c(1 / (B + 1), rep(0.9, m - 1L)), method = method))
  }
  for (correction in c("fdr_bh", "fdr_by")) {
    method <- if (correction == "fdr_bh") "BH" else "BY"
    for (m in c(9L, 21L, 36L, 190L)) {
      for (alpha in c(0.05, 0.01)) {
        need <- .fdr_min_iterations(m, correction, alpha)
        expect_true(smallest(m, need, method) < alpha)
        expect_false(smallest(m, need - 1L, method) < alpha)
      }
    }
  }
})

test_that(".fdr_min_iterations is NULL for the family-wise and no correction", {
  expect_null(.fdr_min_iterations(36L, "fwer"))
  expect_null(.fdr_min_iterations(36L, "none"))
})

test_that("an FDR correction below its floor warns, and fwer does not", {
  skip_on_cran()
  df <- make_poly()
  cu <- q3_cutoff(df, iterations = 40L)
  reset_notices()
  expect_warning(
    suppressMessages(
      RMlocdepQ3(df, cutoff = cu, correction = "fdr_bh", output = "dataframe")
    ),
    regexp = "Benjamini-Hochberg threshold is out of reach"
  )
  reset_notices()
  expect_no_warning(
    suppressMessages(
      RMlocdepQ3(df, cutoff = cu, correction = "fwer", output = "dataframe")
    )
  )
})

test_that("the FDR warning is silent when there are enough iterations", {
  skip_on_cran()
  # 21 pairs at alpha = .05 needs 420, so a 500-iteration object clears it.
  df <- make_poly()
  cu <- q3_cutoff(df, iterations = 500L)
  reset_notices()
  expect_no_warning(
    suppressMessages(
      RMlocdepQ3(df, cutoff = cu, correction = "fdr_bh", output = "dataframe")
    )
  )
})

test_that(".capitalise leaves an empty string alone", {
  expect_equal(.capitalise(""), "")
  expect_equal(.capitalise("item pairs"), "Item pairs")
})
