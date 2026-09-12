# Tests for RMreliabilityCurve() and the shared information/latent-SD helpers

# Rasch-structured data, so information and reliability take realistic values
make_poly <- function(n = 300, k = 6, seed = 1L, theta_sd = 1.3) {
  set.seed(seed)
  deltas <- lapply(seq_len(k), function(i) {
    centre <- seq(-1.2, 1.2, length.out = k)[i]
    centre + c(-0.8, 0, 0.8)
  })
  df <- as.data.frame(
    sim_partial_score(deltas, stats::rnorm(n, 0, theta_sd))
  )
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# er2_caption() wraps with hard line breaks ("  \n"), so any assertion on
# caption text has to collapse whitespace first or it matches by luck.
flat_caption <- function(p) gsub("[[:space:]]+", " ", p$labels$caption)

# ---------------------------------------------------------------------
# Shared helpers extracted from .marginal_rxx()
# ---------------------------------------------------------------------
test_that(".test_information matches the inline sum it replaced", {
  thr_list <- list(c(-1, 0, 1), c(-0.5, 0.5), 0.2)
  theta <- c(-2, -0.5, 0, 1.5)

  inline <- vapply(theta, function(th) {
    sum(vapply(thr_list, function(thr) {
      cats <- 0:length(thr)
      P <- easyRasch2:::.pcm_cat_probs(th, thr)
      E <- sum(cats * P)
      sum((cats - E)^2 * P)
    }, numeric(1L)))
  }, numeric(1L))

  expect_identical(easyRasch2:::.test_information(thr_list, theta), inline)
})

test_that(".test_information is positive and peaks near the item locations", {
  thr_list <- list(0, 0, 0, 0)
  info <- easyRasch2:::.test_information(thr_list, c(-4, 0, 4))
  expect_true(all(info > 0))
  expect_gt(info[2L], info[1L])
  expect_gt(info[2L], info[3L])
})

test_that(".latent_sd reproduces the inline prior-SD estimate", {
  df <- make_poly()
  thr_list <- easyRasch2:::.fit_cml_thresholds(as.matrix(df))
  ge <- seq(-6, 6, length.out = 81L)
  inline <- easyRasch2:::.estimate_prior_sd(
    easyRasch2:::.grid_loglik(
      as.matrix(df), easyRasch2:::.logp_tables(thr_list, ge), ge
    ),
    ge, 0
  )
  expect_identical(easyRasch2:::.latent_sd(as.matrix(df), thr_list), inline)
})

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMreliabilityCurve errors when data has non-zero minimum", {
  df <- make_poly() + 1L
  expect_error(RMreliabilityCurve(df), regexp = "scored starting at 0")
})

test_that("RMreliabilityCurve rejects an out-of-range benchmark", {
  df <- make_poly()
  expect_error(RMreliabilityCurve(df, benchmark = 1.2), regexp = "between 0 and 1")
  expect_error(RMreliabilityCurve(df, benchmark = 0), regexp = "between 0 and 1")
})

test_that("RMreliabilityCurve rejects a reversed theta_range", {
  df <- make_poly()
  expect_error(RMreliabilityCurve(df, theta_range = c(2, 1)), regexp = "length 2")
})

test_that("RMreliabilityCurve rejects unknown items", {
  df <- make_poly()
  expect_error(RMreliabilityCurve(df, items = "nope"), regexp = "not found")
  expect_error(RMreliabilityCurve(df, items = c(1, 99)), regexp = "out-of-range")
})

test_that("RMreliabilityCurve validates item_params", {
  df <- make_poly()
  thr <- easyRasch2:::.fit_cml_thresholds(as.matrix(df))
  expect_error(RMreliabilityCurve(df, item_params = unname(thr)),
               regexp = "named by item")
  mismatch <- thr[1:2]
  names(mismatch) <- c("Z1", "Z2")
  expect_error(RMreliabilityCurve(df, item_params = mismatch),
               regexp = "does not match the number of columns")
  bad <- thr
  bad[[1L]] <- bad[[1L]][1L]           # too few thresholds for a 0-3 item
  expect_error(RMreliabilityCurve(df, item_params = bad),
               regexp = "responses up to category")
  expect_error(RMreliabilityCurve(df, item_params = "nope"),
               regexp = "named list of threshold vectors")
})

test_that("a named item_params subset restricts the scale, as in RMpersonParameters", {
  skip_on_cran()
  df <- make_poly()
  thr <- easyRasch2:::.fit_cml_thresholds(as.matrix(df))
  full <- RMreliabilityCurve(df, item_params = thr, output = "dataframe",
                             n_nodes = 21L, theta_range = c(-2, 2))
  part <- RMreliabilityCurve(df, item_params = thr[1:3], output = "dataframe",
                             n_nodes = 21L, theta_range = c(-2, 2))
  expect_true(all(part$information < full$information))
})

test_that("item_params accepts the RMitemParameters long data.frame", {
  skip_on_cran()
  df <- make_poly()
  long <- RMitemParameters(df, format = "long", se = FALSE,
                           output = "dataframe")
  a <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 21L)
  b <- RMreliabilityCurve(df, item_params = long, output = "dataframe",
                          n_nodes = 21L)
  expect_equal(a$information, b$information, tolerance = 1e-8)
})

# ---------------------------------------------------------------------
# Curve contents
# ---------------------------------------------------------------------
test_that("output = 'dataframe' returns the curve and its summary attributes", {
  skip_on_cran()
  df <- make_poly()
  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L)

  expect_s3_class(res, "data.frame")
  expect_identical(nrow(res), 41L)
  expect_named(res, c("theta", "information", "sem", "reliability"))
  expect_true(all(res$information > 0))
  expect_true(all(res$sem > 0))
  expect_true(all(res$reliability > 0 & res$reliability < 1))

  for (a in c("sigma", "marginal_ratio", "marginal_green", "sem_average")) {
    expect_true(is.finite(attr(res, a)))
  }
})

test_that("sem and reliability are the documented transforms of information", {
  skip_on_cran()
  df <- make_poly()
  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L)
  sigma <- attr(res, "sigma")

  expect_equal(res$sem, 1 / sqrt(res$information))
  # Ratio form, the same one .marginal_rxx() integrates
  expect_equal(res$reliability, sigma^2 / (sigma^2 + res$sem^2))
})

test_that("marginal_ratio matches RMreliability's marginal row", {
  skip_on_cran()
  df <- make_poly()
  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L)
  expect_identical(attr(res, "marginal_ratio"), easyRasch2:::.marginal_rxx(df))
})

test_that("both functions read the marginal off the same shared helper", {
  skip_on_cran()
  df <- make_poly()
  thr <- easyRasch2:::.fit_cml_thresholds(as.matrix(df))
  sigma <- easyRasch2:::.latent_sd(as.matrix(df), thr)
  marg <- easyRasch2:::.marginal_summaries(thr, sigma)
  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L)

  expect_identical(easyRasch2:::.marginal_rxx(df), marg$ratio)
  expect_identical(attr(res, "marginal_ratio"), marg$ratio)
  expect_identical(attr(res, "marginal_green"), marg$green)
  expect_identical(attr(res, "sem_average"), marg$sem_average)
})

test_that("the marginal does not move with the plotted grid", {
  skip_on_cran()
  df <- make_poly()
  vals <- vapply(list(
    list(n_nodes = 21L), list(n_nodes = 301L),
    list(theta_range = c(-2, 2)), list(theta_range = c(-8, 8))
  ), function(a) {
    attr(do.call(RMreliabilityCurve,
                 c(list(df, output = "dataframe"), a)), "marginal_ratio")
  }, numeric(1L))
  expect_equal(length(unique(vals)), 1L)
})

test_that("the RMreliability table row equals the curve attribute", {
  skip_on_cran()
  skip_if_not_installed("mirt")
  skip_if_not_installed("ggdist")
  df <- make_poly(n = 150L)
  tab <- RMreliability(df, draws = 30, rmu_iter = 3, seed = 1,
                       output = "dataframe")
  row <- tab[tab$metric == "Marginal (curve mean)", ]
  expect_equal(nrow(row), 1L)
  expect_equal(
    row$estimate,
    attr(RMreliabilityCurve(df, output = "dataframe"), "marginal_ratio")
  )
})

test_that("marginal_ratio is the weighted curve mean, not the ratio of averages", {
  skip_on_cran()
  df <- make_poly()
  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L)
  sigma <- attr(res, "sigma")

  # Rebuild both candidates on the integration grid the function uses
  thr <- easyRasch2:::.fit_cml_thresholds(as.matrix(df))
  wide <- seq(-6 * sigma, 6 * sigma, length.out = 161L)
  w <- stats::dnorm(wide, 0, sigma)
  w <- w / sum(w)
  sem2 <- 1 / easyRasch2:::.test_information(thr, wide)

  curve_mean <- sum(w * (sigma^2 / (sigma^2 + sem2)))
  ratio_of_averages <- sigma^2 / (sigma^2 + sum(w * sem2))

  expect_equal(attr(res, "marginal_ratio"), curve_mean)
  expect_false(isTRUE(all.equal(curve_mean, ratio_of_averages)))
  # sigma^2 / (sigma^2 + x) is convex in x, so averaging the curve gives the
  # larger value, and it is the one nearer observed-score reliability
  expect_gt(curve_mean, ratio_of_averages)
})

test_that("the reliability reference line matches the reported marginal_ratio", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  df <- make_poly()
  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L)
  p <- RMreliabilityCurve(df, statistic = "reliability", n_nodes = 41L)
  expect_match(
    flat_caption(p),
    sprintf("marginal reliability \\(%.3f\\)", attr(res, "marginal_ratio"))
  )
})

test_that("marginal_green keeps the superseded subtractive value", {
  skip_on_cran()
  df <- make_poly()
  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L)
  sigma <- attr(res, "sigma")
  expect_equal(
    attr(res, "marginal_green"),
    1 - attr(res, "sem_average")^2 / sigma^2
  )
  # and it is no longer what .marginal_rxx() returns
  expect_false(isTRUE(all.equal(
    attr(res, "marginal_green"),
    easyRasch2:::.marginal_rxx(df)
  )))
})

test_that(".marginal_rxx is bounded in (0, 1) where the old form went negative", {
  skip_on_cran()
  set.seed(3)
  n <- 500
  theta <- stats::rnorm(n, 0, 0.35)
  b <- seq(-1, 1, length.out = 6L)
  df <- as.data.frame(vapply(b, function(bi) {
    stats::rbinom(n, 1L, stats::plogis(theta - bi))
  }, numeric(n)))
  colnames(df) <- paste0("I", seq_along(b))

  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L)
  expect_lt(attr(res, "marginal_green"), 0)          # old form
  marg <- easyRasch2:::.marginal_rxx(df)             # new form
  expect_gt(marg, 0)
  expect_lt(marg, 1)
})

test_that("the ratio form is bounded where the Green form need not be", {
  # Few items and a narrow trait spread: the case Milanzi et al. (2015) flag,
  # where 1 - SEM^2 / sigma^2 goes below zero and the ratio form does not.
  set.seed(3)
  n <- 500
  theta <- stats::rnorm(n, 0, 0.35)
  b <- seq(-1, 1, length.out = 6L)
  df <- as.data.frame(vapply(b, function(bi) {
    stats::rbinom(n, 1L, stats::plogis(theta - bi))
  }, numeric(n)))
  colnames(df) <- paste0("I", seq_along(b))

  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L)
  expect_true(all(res$reliability > 0 & res$reliability < 1))
  expect_lt(attr(res, "marginal_green"), 0)
  expect_gt(attr(res, "marginal_ratio"), 0)
})

test_that("theta_range defaults to +/- 3 sigma and can be overridden", {
  skip_on_cran()
  df <- make_poly()
  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L)
  sigma <- attr(res, "sigma")
  expect_equal(range(res$theta), c(-3 * sigma, 3 * sigma))

  res2 <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 41L,
                             theta_range = c(-2, 2))
  expect_equal(range(res2$theta), c(-2, 2))
})

test_that("items subsets the scale and lowers information", {
  skip_on_cran()
  df <- make_poly()
  full <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 21L,
                             theta_range = c(-2, 2))
  short <- RMreliabilityCurve(df, items = c("I1", "I2", "I3"),
                              output = "dataframe", n_nodes = 21L,
                              theta_range = c(-2, 2))
  expect_true(all(short$information < full$information))
})

test_that("supplied item_params reproduce the estimated-threshold curve", {
  skip_on_cran()
  df <- make_poly()
  thr <- easyRasch2:::.fit_cml_thresholds(as.matrix(df))
  a <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 21L)
  b <- RMreliabilityCurve(df, item_params = thr, output = "dataframe",
                          n_nodes = 21L)
  expect_equal(a$information, b$information)
  expect_true(all(vapply(thr, is.numeric, logical(1L))))
})

# ---------------------------------------------------------------------
# Benchmark
# ---------------------------------------------------------------------
test_that("benchmark reports a contiguous range and a respondent percentage", {
  skip_on_cran()
  df <- make_poly()
  res <- RMreliabilityCurve(df, benchmark = 0.6, output = "dataframe",
                            n_nodes = 81L)
  rng <- attr(res, "benchmark_range")
  pct <- attr(res, "benchmark_percent")

  expect_s3_class(rng, "data.frame")
  expect_named(rng, c("xmin", "xmax"))
  expect_true(all(rng$xmin < rng$xmax))
  expect_true(pct >= 0 && pct <= 100)

  inside <- res$theta >= rng$xmin[1L] & res$theta <= rng$xmax[1L]
  expect_true(all(res$reliability[inside] >= 0.6))
})

test_that("an unreachable benchmark yields a NULL range and 0 percent", {
  skip_on_cran()
  df <- make_poly()
  res <- RMreliabilityCurve(df, benchmark = 0.999, output = "dataframe",
                            n_nodes = 41L)
  expect_null(attr(res, "benchmark_range"))
  expect_equal(attr(res, "benchmark_percent"), 0)
})

test_that(".curve_runs finds every contiguous run", {
  x <- 1:10
  expect_null(easyRasch2:::.curve_runs(x, rep(FALSE, 10)))
  one <- easyRasch2:::.curve_runs(x, c(F, F, T, T, T, F, F, F, F, F))
  expect_equal(one$xmin, 3)
  expect_equal(one$xmax, 5)
  two <- easyRasch2:::.curve_runs(x, c(T, T, F, F, T, T, T, F, F, F))
  expect_equal(nrow(two), 2L)
  expect_equal(two$xmax, c(2, 7))
})

# ---------------------------------------------------------------------
# Missing data and extreme scorers
# ---------------------------------------------------------------------
test_that("all-NA respondents are dropped and counted, not silently used", {
  skip_on_cran()
  df <- make_poly()
  df[1L, ] <- NA
  expect_message(
    res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 21L),
    regexp = "no responses dropped"
  )
  expect_true(is.finite(attr(res, "sigma")))
})

test_that("partial missingness is retained", {
  skip_on_cran()
  df <- make_poly()
  df[1:20, 1L] <- NA
  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 21L)
  expect_true(all(is.finite(res$sem)))
})

# ---------------------------------------------------------------------
# Output formats
# ---------------------------------------------------------------------
test_that("output = 'kable' returns a summary table, not the curve", {
  skip_on_cran()
  df <- make_poly()
  res <- RMreliabilityCurve(df, benchmark = 0.6, output = "kable")
  expect_s3_class(res, "knitr_kable")
  txt <- paste(as.character(res), collapse = "\n")
  expect_match(txt, "Marginal reliability \\(curve mean, as in RMreliability\\)")
  expect_match(txt, "Marginal reliability \\(Green/Lord, superseded\\)")
  expect_match(txt, "Theta range with reliability")
  expect_match(txt, "n = 300 respondents")
})

test_that("no benchmark means no benchmark rows in the summary", {
  skip_on_cran()
  df <- make_poly()
  txt <- gsub("[[:space:]]+", " ",
              paste(as.character(RMreliabilityCurve(df, output = "kable")),
                    collapse = " "))
  # attr() partial-matches by default, so an absent `benchmark` attribute used
  # to return `benchmark_percent` and print a row reading ">= NA"
  expect_false(grepl("Theta range with reliability", txt))
  expect_false(grepl("Respondents located", txt))
  expect_false(grepl("NA", txt))
})

test_that("output = 'ggplot' returns a plot for every statistic", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  df <- make_poly()
  for (st in c("sem", "information", "reliability")) {
    p <- RMreliabilityCurve(df, statistic = st, n_nodes = 21L)
    expect_s3_class(p, "ggplot")
  }
})

test_that("the plot caption reports n and the estimation source", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  df <- make_poly()
  cap <- flat_caption(RMreliabilityCurve(df, n_nodes = 21L))
  expect_match(cap, "n = 300 respondents")
  expect_match(cap, "CML thresholds")
  expect_match(cap, "property of the items")
})

test_that("x-axis breaks sit on whole logits with a half-logit minor", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  df <- make_poly()
  b <- ggplot2::ggplot_build(RMreliabilityCurve(df, n_nodes = 41L))
  major <- b$layout$panel_params[[1L]]$x$breaks
  minor <- b$layout$panel_params[[1L]]$x$minor_breaks
  major <- major[is.finite(major)]
  minor <- minor[is.finite(minor)]

  expect_true(all(major == round(major)))
  expect_true(all(diff(sort(major)) == 1))
  expect_true(all(diff(sort(minor)) == 0.5))
  # and a narrower requested range still lands on whole logits
  b2 <- ggplot2::ggplot_build(
    RMreliabilityCurve(df, theta_range = c(-2, 2), n_nodes = 41L)
  )
  m2 <- b2$layout$panel_params[[1L]]$x$breaks
  expect_equal(m2[is.finite(m2)], c(-2, -1, 0, 1, 2))
})

test_that("the caption names the statistic actually plotted", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  df <- make_poly()
  caps <- vapply(c("sem", "information", "reliability"), function(st) {
    flat_caption(RMreliabilityCurve(df, statistic = st, n_nodes = 21L))
  }, character(1L))

  expect_match(caps[["sem"]], "Conditional standard error of measurement")
  expect_match(caps[["information"]], "Test information I\\(theta\\), summed over")
  expect_match(caps[["reliability"]], "Conditional reliability, sigma\\^2")
  # an SEM plot must not open by announcing test information
  expect_false(grepl("Note\\.\\* Test information", caps[["sem"]]))
  # and the information branch must not stutter
  expect_false(grepl("Test information I\\(theta\\), test information",
                     caps[["information"]]))
})

test_that("the reference line is named in the units of its own axis", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  df <- make_poly()
  res <- RMreliabilityCurve(df, output = "dataframe", n_nodes = 21L)

  sem_cap <- flat_caption(RMreliabilityCurve(df, statistic = "sem",
                                             n_nodes = 21L))
  rel_cap <- flat_caption(RMreliabilityCurve(df, statistic = "reliability",
                                             n_nodes = 21L))
  inf_cap <- flat_caption(RMreliabilityCurve(df, statistic = "information",
                                             n_nodes = 21L))

  expect_match(sem_cap, sprintf("average SEM \\(%.2f logits\\)",
                                attr(res, "sem_average")))
  expect_match(rel_cap, sprintf("marginal reliability \\(%.3f\\)",
                                attr(res, "marginal_ratio")))
  expect_match(inf_cap, "average test information")
  # the old wording used one unitless phrase on every axis
  expect_false(grepl("flat marginal summary", sem_cap))
})

test_that("the estimator is named only when a respondent annotation uses it", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  df <- make_poly()
  expect_match(flat_caption(RMreliabilityCurve(df, n_nodes = 21L)),
               "Respondent locations by WLE")
  expect_match(
    flat_caption(RMreliabilityCurve(df, show_density = FALSE, benchmark = 0.6,
                                    n_nodes = 21L)),
    "Respondent locations by WLE"
  )
  expect_false(grepl(
    "Respondent locations",
    flat_caption(RMreliabilityCurve(df, show_density = FALSE, n_nodes = 21L))
  ))
})

test_that("mixing up statistic and method gives a pointed error", {
  df <- make_poly()
  expect_error(RMreliabilityCurve(df, method = "information"),
               regexp = "use statistic = \"information\"")
  expect_error(RMreliabilityCurve(df, method = "sem"),
               regexp = "person-location estimator")
  expect_error(RMreliabilityCurve(df, statistic = "WLE"),
               regexp = "use method = \"WLE\"")
  expect_error(RMreliabilityCurve(df, statistic = "EAP"),
               regexp = "quantity on the y-axis")
})

test_that("short scales get the asymptotic-SE caveat in the caption", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  short <- RMreliabilityCurve(make_poly(k = 5L), n_nodes = 21L)
  long <- RMreliabilityCurve(make_poly(k = 12L), n_nodes = 21L)
  expect_match(flat_caption(short), "asymptotic")
  expect_false(grepl("asymptotic", flat_caption(long)))
})

# ---------------------------------------------------------------------
# Bootstrap band
# ---------------------------------------------------------------------
test_that("boot = TRUE adds a band that brackets the point curve", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  df <- make_poly()
  res <- RMreliabilityCurve(df, boot = TRUE, boot_iter = 15, parallel = FALSE,
                            seed = 42, output = "dataframe", n_nodes = 21L)
  expect_true(all(c("sem_lower", "sem_upper", "reliability_lower",
                    "reliability_upper", "information_lower",
                    "information_upper") %in% names(res)))
  expect_true(all(res$sem_lower <= res$sem_upper))
  expect_identical(attr(res, "boot_iter"), 15L)
})

test_that("the bootstrap is reproducible under a fixed seed", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  df <- make_poly()
  args <- list(df, boot = TRUE, boot_iter = 10, parallel = FALSE,
               seed = 7, output = "dataframe", n_nodes = 21L)
  a <- do.call(RMreliabilityCurve, args)
  b <- do.call(RMreliabilityCurve, args)
  expect_identical(a$sem_lower, b$sem_lower)
  expect_identical(a$reliability_upper, b$reliability_upper)
})
