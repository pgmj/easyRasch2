# Tests for the partial gamma DIF update in 1.3.1.9001: the vectorised
# coefficient, the `p_value = NULL` resolution in RMdifGamma(), the shared
# iteration and interval defaults, the two data-generating processes of
# RMdifGammaCutoff(), and the advisory notices.
#
# The notices use rlang's once-per-session frequency, so each test that
# expects one resets the counter first.

make_dif_data <- function(n1 = 200, n2 = 100, k = 8, mu2 = 1, seed = 11L) {
  set.seed(seed)
  th <- c(stats::rnorm(n1), stats::rnorm(n2, mu2))
  df <- as.data.frame(sapply(
    seq_len(k),
    function(j) stats::rbinom(n1 + n2, 2, stats::plogis(th - (j - (k + 1) / 2) / 3))
  ))
  names(df) <- paste0("I", seq_len(k))
  list(df = df, grp = rep(c("a", "b"), c(n1, n2)))
}

quiet <- function(expr) suppressMessages(suppressWarnings(expr))

dif_cutoff <- function(df, grp, ...) {
  quiet(RMdifGammaCutoff(df, grp, parallel = FALSE, seed = 1L, ...))
}

reset_notices <- function() {
  rlang::reset_message_verbosity("easyRasch2_band_flagging_dif")
  rlang::reset_message_verbosity("easyRasch2_low_iterations_dif")
}

iarm_gamma <- function(df, grp) {
  out <- NULL
  utils::capture.output(out <- iarm::partgam_DIF(as.data.frame(df), grp))
  as.numeric(out$gamma)
}

# ---------------------------------------------------------------------
# Vectorised partial gamma reproduces iarm::partgam_DIF()
# ---------------------------------------------------------------------
test_that(".partgam_dif_gamma equals iarm::partgam_DIF exactly", {
  skip_if_not_installed("iarm")
  d <- make_dif_data()
  grp <- d$grp
  cases <- list(
    character = grp,
    factor = factor(grp),
    reversed = factor(grp, levels = c("b", "a")),
    numeric = ifelse(grp == "a", 1L, 2L),
    unused_level = factor(grp, levels = c("a", "b", "c"))
  )
  for (nm in names(cases)) {
    fast <- .partgam_dif_gamma(d$df, cases[[nm]])$gamma
    expect_identical(fast, iarm_gamma(d$df, cases[[nm]]), info = nm)
  }
})

test_that(".partgam_dif_gamma handles NA, three groups and dichotomous items", {
  skip_if_not_installed("iarm")
  d <- make_dif_data()
  df_na <- d$df
  df_na[c(3, 40, 77), 2] <- NA
  grp_na <- d$grp
  grp_na[c(5, 90)] <- NA
  expect_identical(
    .partgam_dif_gamma(df_na, grp_na)$gamma,
    iarm_gamma(df_na, grp_na)
  )
  set.seed(2)
  g3 <- sample(c("x", "y", "z"), nrow(d$df), replace = TRUE)
  expect_identical(.partgam_dif_gamma(d$df, g3)$gamma, iarm_gamma(d$df, g3))
  dich <- as.data.frame(lapply(d$df, function(x) as.integer(x > 0)))
  expect_identical(
    .partgam_dif_gamma(dich, d$grp)$gamma,
    iarm_gamma(dich, d$grp)
  )
})

test_that("reversing the group levels negates every gamma", {
  d <- make_dif_data()
  ab <- .partgam_dif_gamma(d$df, factor(d$grp, levels = c("a", "b")))$gamma
  ba <- .partgam_dif_gamma(d$df, factor(d$grp, levels = c("b", "a")))$gamma
  expect_identical(ab, -ba)
})

test_that("a constant item gives NA for that item only", {
  d <- make_dif_data()
  df <- d$df
  df$I8 <- 0L
  g <- .partgam_dif_gamma(df, d$grp)$gamma
  expect_true(is.na(g[8]))
  expect_true(all(is.finite(g[-8])))
})

# ---------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------
test_that("RMdifGammaCutoff shares the item-fit defaults", {
  expect_equal(formals(RMdifGammaCutoff)$iterations, 400)
  expect_equal(formals(RMdifGammaCutoff)$hdci_width, 0.95)
  expect_equal(
    eval(formals(RMdifGammaCutoff)$dgp),
    c("conditional", "permutation")
  )
})

test_that("RMdifGamma defaults p_value to NULL", {
  expect_null(formals(RMdifGamma)$p_value)
})

# ---------------------------------------------------------------------
# p_value = NULL resolution
# ---------------------------------------------------------------------
test_that("p_value = NULL uses p-values when the full cutoff object is given", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  d <- make_dif_data()
  reset_notices()
  cu <- dif_cutoff(d$df, d$grp, iterations = 40L)
  res <- quiet(RMdifGamma(d$df, d$grp, cutoff = cu, output = "dataframe"))
  expect_true(all(c("p_gamma", "padj_gamma") %in% names(res)))
})

test_that("p_value = NULL falls back to the interval without simulations", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  d <- make_dif_data()
  cu <- dif_cutoff(d$df, d$grp, iterations = 40L)
  res <- quiet(
    RMdifGamma(d$df, d$grp, cutoff = cu$item_cutoffs, output = "dataframe")
  )
  expect_false(any(c("p_gamma", "padj_gamma") %in% names(res)))
  expect_true("flagged" %in% names(res))
})

test_that("p_value = NULL with no cutoff leaves the asymptotic table alone", {
  skip_if_not_installed("iarm")
  d <- make_dif_data()
  res <- RMdifGamma(d$df, d$grp, output = "dataframe")
  expect_false(any(c("p_gamma", "padj_gamma") %in% names(res)))
  expect_true("padj_bh" %in% names(res))
})

test_that("p_value rejects non-logical input", {
  skip_if_not_installed("iarm")
  d <- make_dif_data()
  expect_error(RMdifGamma(d$df, d$grp, p_value = "yes"), regexp = "TRUE, FALSE")
  expect_error(RMdifGamma(d$df, d$grp, p_value = NA), regexp = "TRUE, FALSE")
})

test_that("the observed gamma tested equals the gamma displayed", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  d <- make_dif_data()
  cu <- dif_cutoff(d$df, d$grp, iterations = 40L)
  res <- quiet(RMdifGamma(d$df, d$grp, cutoff = cu, output = "dataframe"))
  expect_equal(res$gamma, .partgam_dif_gamma(d$df, d$grp)$gamma)
})

# ---------------------------------------------------------------------
# Data-generating processes
# ---------------------------------------------------------------------
test_that("both dgp options run, are reproducible and record the dgp", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  d <- make_dif_data()
  for (dgp in c("conditional", "permutation")) {
    r1 <- dif_cutoff(d$df, d$grp, iterations = 30L, dgp = dgp)
    r2 <- dif_cutoff(d$df, d$grp, iterations = 30L, dgp = dgp)
    expect_equal(r1$results, r2$results, info = dgp)
    expect_identical(r1$dgp, dgp)
    expect_equal(r1$actual_iterations, 30L)
    expect_equal(r1$requested_iterations, 30L)
    expect_equal(r1$dif_group_sizes, c(200L, 100L))
  }
})

test_that("permutation keeps group sizes within every score stratum", {
  d <- make_dif_data()
  dm <- as.matrix(d$df)
  strata <- split(seq_len(nrow(dm)), rowSums(dm))
  dl <- list(
    dgp = "permutation",
    data_mat = dm,
    strata = strata[lengths(strata) > 1L],
    dif_factor = factor(d$grp),
    item_names = colnames(dm)
  )
  # Re-run the permutation step the iteration uses and check its invariant.
  set.seed(1)
  g <- dl$dif_factor
  for (s in dl$strata) g[s] <- g[s[sample.int(length(s))]]
  score <- rowSums(dm)
  expect_equal(unclass(table(score, g)), unclass(table(score, dl$dif_factor)),
    ignore_attr = TRUE)
  expect_false(identical(g, dl$dif_factor))
})

test_that("the conditional null keeps each respondent's total score", {
  skip_if_not_installed("psychotools")
  d <- make_dif_data()
  dm <- as.matrix(d$df)
  thr <- .fit_cml_thresholds(dm)
  dl <- list(
    thr_list = thr,
    cond_groups = .cond_groups(dm, thr),
    sample_n = nrow(dm),
    n_items = ncol(dm)
  )
  set.seed(1)
  sim <- as.matrix(.sim_cond_dataset(dl))
  expect_equal(rowSums(sim), rowSums(dm))
  expect_false(identical(unname(sim), unname(dm)))
})

test_that("the permutation null works when the model cannot be fitted", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  d <- make_dif_data()
  df <- d$df
  df$I8 <- 0L # a constant item breaks the CML fit
  res <- dif_cutoff(df, d$grp, iterations = 20L, dgp = "permutation")
  expect_equal(res$actual_iterations, 20L)
  expect_true(all(is.na(res$results$gamma[res$results$Item == "I8"])))
})

# ---------------------------------------------------------------------
# Notices
# ---------------------------------------------------------------------
test_that("interval flagging reports the rate it implies over items", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  d <- make_dif_data()
  cu <- dif_cutoff(d$df, d$grp, iterations = 40L)
  reset_notices()
  expect_message(
    RMdifGamma(d$df, d$grp, cutoff = cu, p_value = FALSE, output = "dataframe"),
    regexp = "over 8 items"
  )
})

test_that("p-values from fewer than 400 iterations trigger a notice", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  d <- make_dif_data()
  cu <- dif_cutoff(d$df, d$grp, iterations = 40L)
  reset_notices()
  expect_message(
    RMdifGamma(d$df, d$grp, cutoff = cu, output = "dataframe"),
    regexp = "below the calibrated floor"
  )
})

test_that("an unused factor level does not trigger a sample mismatch", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  d <- make_dif_data()
  grp <- factor(d$grp, levels = c("a", "b", "c"))
  cu <- dif_cutoff(d$df, grp, iterations = 20L)
  expect_equal(cu$dif_group_sizes, c(200L, 100L))
  expect_no_warning(
    suppressMessages(RMdifGamma(d$df, grp, cutoff = cu, output = "dataframe"))
  )
})

# ---------------------------------------------------------------------
# Asymptotic BH column
# ---------------------------------------------------------------------
test_that("RMdifGamma applies a real Benjamini-Hochberg adjustment", {
  skip_if_not_installed("iarm")
  d <- make_dif_data()
  raw <- NULL
  utils::capture.output(raw <- iarm::partgam_DIF(d$df, d$grp))
  res <- RMdifGamma(d$df, d$grp, output = "dataframe")
  expect_equal(
    res$padj_bh,
    stats::p.adjust(as.numeric(raw$pvalue), method = "BH")
  )
  expect_identical(res$Significance, .p_stars(res$padj_bh))
})

test_that("RMlocdepGamma adjusts over all tests in both directions", {
  skip_if_not_installed("iarm")
  d <- make_dif_data()
  raw <- NULL
  utils::capture.output(raw <- iarm::partgam_LD(d$df))
  res <- RMlocdepGamma(d$df, output = "dataframe")
  p_all <- c(as.numeric(raw[[1]][[5]]), as.numeric(raw[[2]][[5]]))
  bh <- stats::p.adjust(p_all, method = "BH")
  n1 <- nrow(res$direction1)
  expect_equal(res$direction1$padj_bh, bh[seq_len(n1)])
  expect_equal(res$direction2$padj_bh, bh[-seq_len(n1)])
})

test_that(".p_stars follows iarm's cutpoints", {
  expect_identical(
    .p_stars(c(0.0005, 0.005, 0.03, 0.07, 0.2, NA)),
    c("***", "**", "*", ".", "", "")
  )
})
