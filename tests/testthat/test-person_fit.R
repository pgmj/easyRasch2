# Tests for RMpersonFit()

# Genuine Rasch-null PCM data (so fit statistics behave under H0).
sim_pcm_null <- function(n = 300, J = 8, person_sd = 1.3, seed = 11L) {
  set.seed(seed)
  theta <- rnorm(n, 0, person_sd)
  thr   <- lapply(seq_len(J), function(j) sort(rnorm(2, 0, 0.9)))
  pcm_p <- function(th, t) {
    ce <- c(0, cumsum(th - t)); exp(ce - max(ce)) / sum(exp(ce - max(ce)))
  }
  m <- sapply(seq_len(J), function(j) {
    sapply(theta, function(th) sample(0:2, 1, prob = pcm_p(th, thr[[j]])))
  })
  df <- as.data.frame(m); colnames(df) <- paste0("I", seq_len(J))
  df
}

sim_rasch_null <- function(n = 300, J = 10, seed = 3L) {
  set.seed(seed)
  theta <- rnorm(n, 0, 1.3); beta <- seq(-2, 2, length.out = J)
  m <- sapply(seq_len(J), function(j) rbinom(n, 1, plogis(theta - beta[j])))
  df <- as.data.frame(m); colnames(df) <- paste0("I", seq_len(J))
  df
}

# ---------------------------------------------------------------------
# Input validation & structure
# ---------------------------------------------------------------------

test_that("RMpersonFit errors on non-zero minimum", {
  expect_error(RMpersonFit(sim_rasch_null() + 1L, iterations = 0),
               regexp = "scored starting at 0")
})

test_that("dataframe has expected columns with all statistics", {
  res <- RMpersonFit(sim_pcm_null(n = 120), iterations = 100, seed = 1,
                     output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("id", "n_answered", "sum_score", "infit_msq",
                    "outfit_msq", "lz", "p_infit", "p_outfit", "p_lz",
                    "flagged") %in% names(res)))
  expect_equal(nrow(res), 120L)
})

test_that("iterations = 0 returns statistics without p-values", {
  res <- RMpersonFit(sim_pcm_null(n = 100), iterations = 0, output = "dataframe")
  expect_false(any(grepl("^p_", names(res))))
  expect_false("flagged" %in% names(res))
  expect_true(all(c("infit_msq", "outfit_msq", "lz") %in% names(res)))
})

test_that("statistics argument subsets the output", {
  res <- RMpersonFit(sim_pcm_null(n = 100), statistics = "outfit",
                     iterations = 50, seed = 1, output = "dataframe")
  expect_true("outfit_msq" %in% names(res))
  expect_false(any(c("infit_msq", "lz") %in% names(res)))
})

test_that("p-values lie in [0, 1]", {
  res <- RMpersonFit(sim_pcm_null(n = 120), iterations = 100, seed = 2,
                     output = "dataframe")
  for (col in c("p_infit", "p_outfit", "p_lz")) {
    p <- res[[col]][!is.na(res[[col]])]
    expect_true(all(p >= 0 & p <= 1))
  }
})

# ---------------------------------------------------------------------
# Statistical behaviour
# ---------------------------------------------------------------------

test_that("conditional MSQ means are near 1 under the null", {
  res <- RMpersonFit(sim_pcm_null(n = 400), iterations = 0, output = "dataframe")
  expect_equal(mean(res$infit_msq,  na.rm = TRUE), 1, tolerance = 0.05)
  expect_equal(mean(res$outfit_msq, na.rm = TRUE), 1, tolerance = 0.05)
})

test_that("Type I error is near nominal under the null", {
  res <- RMpersonFit(sim_pcm_null(n = 400),
                     iterations = 400, seed = 4, output = "dataframe")
  expect_lt(mean(res$flagged, na.rm = TRUE), 0.12)
})

test_that("lz p-values are calibrated (not conservative) under the null", {
  # The lz resampling re-estimates the person location for each simulated
  # pattern, reproducing the ability-estimation effect (Snijders 2001;
  # Sinharay 2016) that otherwise makes the test conservative. With a fixed
  # location the rejection rate here would be roughly half nominal.
  res <- RMpersonFit(sim_rasch_null(n = 400, J = 10, seed = 3),
                     statistics = "lz", iterations = 300, seed = 7,
                     output = "dataframe")
  p <- res$p_lz[!is.na(res$p_lz)]
  expect_gt(mean(p), 0.45)
  expect_lt(mean(p), 0.55)
  expect_gt(mean(p < 0.10), 0.06)   # near nominal 0.10, not the conservative ~0.04
})

test_that("aberrant responders are detected (power)", {
  dat <- sim_rasch_null(n = 300, J = 10, seed = 3)
  ab  <- 1:15
  dat[ab, ] <- dat[ab, rev(seq_len(ncol(dat)))]   # reverse -> aberrant
  res <- RMpersonFit(dat, iterations = 300, seed = 2, output = "dataframe")
  expect_gt(mean(res$flagged[ab],  na.rm = TRUE), 0.6)
  expect_lt(mean(res$flagged[-ab], na.rm = TRUE), 0.12)
})

# ---------------------------------------------------------------------
# Missingness, extremes, schemes, outputs
# ---------------------------------------------------------------------

test_that("partial missingness is handled (rows kept, estimates finite)", {
  dat <- sim_pcm_null(n = 200)
  dat[cbind(sample(200, 40), sample(8, 40, replace = TRUE))] <- NA
  res <- RMpersonFit(dat, iterations = 0, output = "dataframe")
  expect_equal(nrow(res), 200L)
  expect_true(any(res$n_answered < 8L))
  finite_non_extreme <- res$outfit_msq[res$n_answered > 1L &
                                       res$sum_score > 0 &
                                       !is.na(res$outfit_msq)]
  expect_true(all(is.finite(finite_non_extreme)))
})

test_that("extreme scorers receive NA statistics", {
  dat <- sim_rasch_null(n = 200, J = 8, seed = 9)
  dat[1, ] <- 0L; dat[2, ] <- 1L
  res <- RMpersonFit(dat, iterations = 50, seed = 1, output = "dataframe")
  expect_true(is.na(res$outfit_msq[1]) && is.na(res$outfit_msq[2]))
})

test_that("all three statistics get valid p-values (per-statistic schemes)", {
  res <- RMpersonFit(sim_pcm_null(n = 150), iterations = 200, seed = 7,
                     output = "dataframe")
  for (col in c("p_infit", "p_outfit", "p_lz")) {
    p <- res[[col]][!is.na(res[[col]])]
    expect_true(length(p) > 0L && all(p >= 0 & p <= 1))
  }
})

test_that("MSQ-only run needs no person estimate but still yields p-values", {
  dat <- sim_pcm_null(n = 150)
  dat[cbind(sample(150, 25), sample(8, 25, replace = TRUE))] <- NA
  res <- RMpersonFit(dat, statistics = c("infit", "outfit"),
                     iterations = 200, seed = 1, output = "dataframe")
  ok <- res$n_answered > 1L & res$sum_score > 0 & !is.na(res$infit_msq)
  expect_true(all(is.finite(res$p_infit[ok])))   # conditional scheme, missingness handled
})

test_that("conditional sampler produces patterns with the fixed total", {
  thr <- list(0.2, -0.3, 0.1, 0.0, 0.4)
  sims <- easyRasch2:::.sim_conditional(thr, total = 2L, B = 300)
  expect_true(all(rowSums(sims) == 2L))
})

test_that("zstd adds non-inferential ZSTD columns", {
  res <- RMpersonFit(sim_pcm_null(n = 100), iterations = 50, zstd = TRUE,
                     seed = 1, output = "dataframe")
  expect_true(all(c("infit_zstd", "outfit_zstd") %in% names(res)))
})

test_that("kable output returns knitr_kable", {
  skip_if_not_installed("knitr")
  dat <- sim_pcm_null(n = 100)
  expect_s3_class(RMpersonFit(dat, iterations = 50), "knitr_kable")
})

test_that("ggplot output returns a named list of plots, one per statistic", {
  skip_if_not_installed("ggplot2")
  dat <- sim_pcm_null(n = 120)
  g <- RMpersonFit(dat, iterations = 100, seed = 1, output = "ggplot")
  expect_type(g, "list")
  expect_named(g, c("infit", "outfit", "lz"))
  expect_true(all(vapply(g, function(p) inherits(p, "ggplot"), logical(1))))
  # each plot is non-blank (has plotted points) -- the infit-blank regression
  expect_true(all(vapply(g, function(p) nrow(p$data) > 0L, logical(1))))
})

test_that("ggplot with a single statistic is a one-element list and not blank", {
  skip_if_not_installed("ggplot2")
  dat <- sim_pcm_null(n = 120)
  gi <- RMpersonFit(dat, statistics = "infit", iterations = 100, seed = 1,
                    output = "ggplot")
  expect_named(gi, "infit")
  expect_s3_class(gi$infit, "ggplot")
  expect_gt(nrow(gi$infit$data), 0L)            # not blank
  expect_match(gi$infit$labels$caption, "flagged")
})

test_that("MSQ plot captions split flagged into underfit/overfit; lz does not", {
  skip_if_not_installed("ggplot2")
  set.seed(2); n <- 400; J <- 8
  theta <- rnorm(n, 0, 1.3); beta <- seq(-1.5, 1.5, length.out = J)
  m <- sapply(seq_len(J), function(j) rbinom(n, 1, plogis(theta - beta[j])))
  m[1:8, ] <- m[1:8, J:1]                        # reversed -> underfit
  dat <- as.data.frame(m); colnames(dat) <- paste0("I", seq_len(J))
  g <- RMpersonFit(dat, iterations = 400, seed = 1, output = "ggplot")
  # (caption may wrap across lines, so check the words separately)
  expect_match(g$infit$labels$caption,  "underfit")
  expect_match(g$infit$labels$caption,  "overfit")
  expect_match(g$outfit$labels$caption, "underfit")
  expect_no_match(g$lz$labels$caption,  "underfit")   # one-sided: no split
})

test_that("flag = 'underfit' flags only the upper (MSQ > 1) direction", {
  set.seed(2); n <- 400; J <- 8
  theta <- rnorm(n, 0, 1.3); beta <- seq(-1.5, 1.5, length.out = J)
  m <- sapply(seq_len(J), function(j) rbinom(n, 1, plogis(theta - beta[j])))
  m[1:10, ] <- m[1:10, J:1]
  dat <- as.data.frame(m); colnames(dat) <- paste0("I", seq_len(J))
  res <- RMpersonFit(dat, statistics = "infit", iterations = 400, seed = 1,
                     flag = "underfit", output = "dataframe")
  flagged <- res[!is.na(res$flagged) & res$flagged, ]
  expect_true(all(flagged$infit_msq > 1))            # never flags overfit
})

test_that("flag = 'underfit' uses a one-sided MSQ p-value (differs from both)", {
  dat <- sim_pcm_null(n = 200)
  a <- RMpersonFit(dat, statistics = "infit", iterations = 400, seed = 1,
                   flag = "both",     output = "dataframe")
  b <- RMpersonFit(dat, statistics = "infit", iterations = 400, seed = 1,
                   flag = "underfit", output = "dataframe")
  expect_false(isTRUE(all.equal(a$p_infit, b$p_infit)))
})

test_that("plot status colours respondents by direction", {
  skip_if_not_installed("ggplot2")
  dat <- sim_pcm_null(n = 150)
  g <- RMpersonFit(dat, iterations = 300, seed = 1, output = "ggplot")
  expect_true("status" %in% names(g$infit$data))
  expect_true(all(levels(g$infit$data$status) %in%
                  c("Not flagged", "Underfit", "Overfit")))
  expect_true(all(levels(g$lz$data$status) %in% c("Not flagged", "Flagged")))
})

test_that("RMpersonFit Monte-Carlo p-values are never exactly 0", {
  dat <- sim_pcm_null()
  res <- RMpersonFit(dat, iterations = 50, seed = 1, output = "dataframe")
  p_cols <- grep("^p_", names(res), value = TRUE)
  expect_true(length(p_cols) > 0L)
  for (pc in p_cols) {
    p <- res[[pc]][!is.na(res[[pc]])]
    # (1 + count) / (B + 1) convention: strictly positive, never 0.
    expect_true(all(p > 0 & p <= 1))
  }
})
