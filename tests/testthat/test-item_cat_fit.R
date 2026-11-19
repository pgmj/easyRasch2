# Tests for RMitemCatFit(), RMitemCatFitPlot() and their internal helpers.

# PCM data, items with four thresholds (five categories). `swap` recodes
# categories 1 and 2 of item 1 after simulation, a misfit with ordered
# threshold estimates.
make_pcm <- function(n = 300, k = 5, seed = 1L, swap = FALSE) {
  set.seed(seed)
  thr <- lapply(seq_len(k), function(i) c(-1.5, -0.5, 0.5, 1.5) + (i - (k + 1) / 2) * 0.3)
  X <- sim_partial_score(thr, stats::rnorm(n))
  if (swap) {
    x <- X[, 1L]
    X[x == 1L, 1L] <- 2L
    X[x == 2L, 1L] <- 1L
  }
  df <- as.data.frame(X)
  colnames(df) <- paste0("I", seq_len(k))
  df
}

reset_catfit_notices <- function() {
  for (cr in c("fwer", "fdr_bh", "fdr_by", "none")) {
    rlang::reset_message_verbosity(paste0("easyRasch2_catfit_low_iterations_", cr))
  }
}

# ---------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------
test_that(".cond_cat_probs matches brute-force enumeration", {
  thr <- list(c(-0.5, 0.7), c(0.2, -0.3), c(-1, 0, 1.2))
  P <- easyRasch2:::.cond_cat_probs(thr)
  m <- lengths(thr)
  pats <- as.matrix(expand.grid(lapply(m, function(mi) 0:mi)))
  w <- apply(pats, 1L, function(x) {
    prod(vapply(seq_along(thr), function(i) {
      exp(-sum(thr[[i]][seq_len(x[i])]))
    }, numeric(1L)))
  })
  r <- rowSums(pats)
  for (i in seq_along(thr)) {
    for (s in 0:sum(m)) {
      sel <- r == s
      brute <- vapply(0:m[i], function(j) {
        sum(w[sel & pats[, i] == j]) / sum(w[sel])
      }, numeric(1L))
      expect_equal(P[[i]][s + 1L, ], brute, tolerance = 1e-12)
    }
  }
})

test_that(".sim_given_totals keeps the totals and follows P(X = j | r)", {
  thr <- list(c(-0.5, 0.7), c(0.2, -0.3), c(-1, 0, 1.2))
  totals <- rep(3L, 20000L)
  set.seed(11)
  Xs <- easyRasch2:::.sim_given_totals(thr, totals)
  expect_identical(as.integer(rowSums(Xs)), totals)
  P <- easyRasch2:::.cond_cat_probs(thr)
  for (i in seq_along(thr)) {
    obs <- tabulate(Xs[, i] + 1L, nbins = length(thr[[i]]) + 1L) / length(totals)
    expect_equal(obs, P[[i]][4L, ], tolerance = 0.02)
  }
})

test_that(".pav agrees with stats::isoreg and pools ties", {
  set.seed(2)
  p <- sort(stats::runif(40))
  y <- stats::rbinom(40, 1, p)
  pv <- easyRasch2:::.pav(p, y)
  expect_equal(pv$cep, stats::isoreg(p, y)$yf, tolerance = 1e-12)
  expect_false(is.unsorted(pv$cep))
  # Tied forecasts are pooled before the regression.
  pv2 <- easyRasch2:::.pav(c(0.2, 0.2, 0.5, 0.5), c(1, 0, 0, 0))
  expect_equal(pv2$p, c(0.2, 0.5))
  expect_equal(pv2$n, c(2L, 2L))
  expect_equal(pv2$cep, c(0.25, 0.25))
})

test_that(".mcb_vector is non-negative, one value per event", {
  prep <- easyRasch2:::.catfit_prepare(make_pcm(n = 200, k = 4))
  mcb <- easyRasch2:::.mcb_vector(prep$X, prep$r, prep$P)
  expect_length(mcb, length(prep$event_item))
  expect_true(all(mcb >= -1e-12))
})

test_that(".catfit_pvalues: floors, ordering and the step-down", {
  set.seed(3)
  event_item <- rep(1:3, each = 2L)
  sim <- matrix(stats::rexp(6 * 99), 6)
  # Observed far above every simulated value for item 1, at the median for
  # the others.
  obs <- c(50, 50, apply(sim[3:6, ], 1L, stats::median))
  for (cr in c("fwer", "fdr_bh", "fdr_by", "none")) {
    pv <- easyRasch2:::.catfit_pvalues(obs, sim, event_item, cr)
    expect_equal(pv$p_item[1L], 1 / 100)
    expect_true(all(pv$p_item >= 1 / 100 & pv$p_item <= 1))
    expect_true(all(pv$p_event >= 1 / 100 & pv$p_event <= 1))
    expect_true(all(pv$padj >= pv$p_item - 1e-12))
  }
  pv <- easyRasch2:::.catfit_pvalues(obs, sim, event_item, "fwer")
  expect_equal(pv$padj[1L], 1 / 100)
  # An observed value equal to every simulated value gives p = 1.
  flat <- easyRasch2:::.catfit_pvalues(rep(1, 2), matrix(1, 2, 20), c(1L, 1L), "none")
  expect_equal(flat$p_item, 1)
  expect_equal(flat$p_event, c(1, 1))
})

test_that(".catfit_pvalues step-down equals an independent implementation", {
  set.seed(4)
  event_item <- rep(1:4, each = 3L)
  pool <- matrix(stats::rgamma(12 * 200, 2), 12)
  pv <- easyRasch2:::.catfit_pvalues(pool[, 1L], pool[, -1L], event_item, "fwer")
  Z <- t(apply(pool, 1L, function(t) (t - mean(t)) / stats::sd(t)))
  Tm <- vapply(1:4, function(i) apply(Z[event_item == i, ], 2L, max), numeric(200))
  p_item <- colMeans(sweep(Tm, 2L, Tm[1L, ], `>=`))
  ord <- order(Tm[1L, ], decreasing = TRUE)
  run <- rep(-Inf, 199)
  p <- numeric(4)
  for (j in 4:1) {
    run <- pmax(run, Tm[-1L, ord[j]])
    p[j] <- (1 + sum(run >= Tm[1L, ord[j]])) / 200
  }
  ref <- numeric(4)
  ref[ord] <- pmin(cummax(p), 1)
  expect_equal(pv$p_item, p_item)
  expect_equal(pv$padj, pmax(ref, p_item))
})

test_that(".score_group_factor forms at most the requested groups", {
  set.seed(5)
  r <- sample(1:15, 300, replace = TRUE)
  g <- easyRasch2:::.score_group_factor(r, 4L)
  expect_s3_class(g, "factor")
  expect_lte(nlevels(g), 4L)
  expect_length(g, 300L)
  expect_false(anyNA(g))
  # Groups do not overlap in score.
  # Unnamed, so the en-dash labels never pass through names() (C locale).
  lims <- t(vapply(levels(g), function(l) range(r[g == l]), numeric(2),
                   USE.NAMES = FALSE))
  expect_true(all(lims[-1L, 1L] > lims[-nrow(lims), 2L]))
})

# ---------------------------------------------------------------------
# RMitemCatFit()
# ---------------------------------------------------------------------
test_that("RMitemCatFit validates its arguments", {
  df <- make_pcm(n = 150, k = 4)
  expect_error(RMitemCatFit(df, iterations = 10), "integer >= 20")
  expect_error(RMitemCatFit(df, iterations = 20, item_labels = c("a", "b")),
               "same length")
  expect_error(RMitemCatFit(df, iterations = 20,
                            item_labels = c("a", "a", "b", "c")), "unique")
  expect_error(RMitemCatFit(df, iterations = 20, correction = "holm"))
  expect_error(RMitemCatFit(df, iterations = 20, output = "plot"))
})

test_that("RMitemCatFit stops on an empty category", {
  df <- make_pcm(n = 150, k = 4)
  df$I2[df$I2 == 2L] <- 1L
  expect_error(RMitemCatFit(df, iterations = 20),
               class = "easyRasch2_category_error")
})

test_that("RMitemCatFit and RMitemCatFitPlot stop on an extreme-only category", {
  df <- make_pcm(n = 150, k = 4)
  # Category 4 of I1 used only by respondents with a perfect total score.
  df$I1[df$I1 == 4L] <- 3L
  df[1:5, ] <- 4L
  err <- expect_error(RMitemCatFit(df, iterations = 20),
                      class = "easyRasch2_category_error")
  expect_identical(err$items, "I1")
  expect_match(conditionMessage(err), "I1 (category 4)", fixed = TRUE)
  expect_match(conditionMessage(err), "RMitemCatFit() stops here", fixed = TRUE)
  expect_error(RMitemCatFitPlot(df, type = "adjacent"),
               class = "easyRasch2_category_error")
})

test_that("RMitemCatFit dataframe output has the documented columns", {
  reset_catfit_notices()
  df <- make_pcm(n = 250, k = 4)
  res <- suppressMessages(RMitemCatFit(df, iterations = 40, seed = 1,
                                       output = "dataframe"))
  expect_s3_class(res, "data.frame")
  expect_identical(names(res), c("item", paste0("mcb_", 1:4),
                                 paste0("p_", 1:4), "p", "padj"))
  expect_identical(res$item, colnames(df))
  ps <- unlist(res[, c(paste0("p_", 1:4), "p", "padj")])
  expect_true(all(ps > 0 & ps <= 1))
  expect_true(all(res$padj >= res$p - 1e-12))
  expect_true(all(unlist(res[, paste0("mcb_", 1:4)]) >= -1e-12))
})

test_that("RMitemCatFit is reproducible with a seed", {
  df <- make_pcm(n = 200, k = 4)
  a <- suppressMessages(RMitemCatFit(df, iterations = 30, seed = 7, output = "dataframe"))
  b <- suppressMessages(RMitemCatFit(df, iterations = 30, seed = 7, output = "dataframe"))
  expect_identical(a, b)
})

test_that("RMitemCatFit flags an item with swapped categories", {
  df <- make_pcm(n = 500, k = 5, seed = 2L, swap = TRUE)
  res <- suppressMessages(RMitemCatFit(df, iterations = 60, seed = 1,
                                       output = "dataframe"))
  expect_equal(which.min(res$p), 1L)
  expect_lt(res$padj[1L], 0.05)
})

test_that("RMitemCatFit list output carries what RMitemCatFitPlot needs", {
  df <- make_pcm(n = 200, k = 4)
  # 30 iterations cannot reach the BH threshold, which the function says.
  expect_warning(
    fit <- suppressMessages(RMitemCatFit(df, iterations = 30, seed = 3,
                                         correction = "fdr_bh", output = "list")),
    "out of reach"
  )
  expect_s3_class(fit, "RMitemCatFit")
  expect_named(fit, c("table", "correction", "iterations", "iterations_used",
                      "n_used", "n_total", "qualifiers", "item_labels",
                      "prep", "sims"))
  expect_identical(fit$correction, "fdr_bh")
  expect_identical(fit$iterations, 30L)
  expect_lte(fit$iterations_used, 30L)
  expect_equal(ncol(fit$sims$mcb), fit$iterations_used)
  expect_equal(ncol(fit$sims$cep), 30L)
  expect_equal(ncol(fit$sims$rest), fit$iterations_used)
})

test_that("RMitemCatFit kable output has a caption naming the correction", {
  df <- make_pcm(n = 200, k = 4)
  tab <- suppressMessages(RMitemCatFit(df, iterations = 25, seed = 1))
  expect_s3_class(tab, "knitr_kable")
  txt <- paste(tab, collapse = "\n")
  expect_match(txt, "Item p")
  expect_match(txt, "Adjusted p")
  expect_match(txt, .correction_label("fwer"), fixed = TRUE)
})

test_that("RMitemCatFit uses complete cases and reports them", {
  df <- make_pcm(n = 250, k = 4)
  df$I3[c(1, 5, 9)] <- NA
  fit <- suppressMessages(RMitemCatFit(df, iterations = 25, seed = 1,
                                       output = "list"))
  expect_identical(fit$n_total, 250L)
  expect_true("complete cases" %in% fit$qualifiers)
  expect_lte(fit$n_used, 247L)
})

test_that("RMitemCatFit notes fewer than 400 iterations", {
  reset_catfit_notices()
  df <- make_pcm(n = 150, k = 4)
  expect_message(RMitemCatFit(df, iterations = 20, seed = 1, output = "dataframe"),
                 "below the recommended 400")
})

test_that("RMitemCatFit gives the same result in parallel", {
  skip_on_cran()
  skip_if_not_installed("mirai")
  df <- make_pcm(n = 200, k = 4)
  seq <- suppressMessages(RMitemCatFit(df, iterations = 24, seed = 5,
                                       output = "dataframe"))
  par <- suppressMessages(RMitemCatFit(df, iterations = 24, seed = 5,
                                       parallel = TRUE, n_cores = 2L,
                                       output = "dataframe"))
  expect_identical(par, seq)
})

# ---------------------------------------------------------------------
# RMitemCatFitPlot()
# ---------------------------------------------------------------------
test_that("RMitemCatFitPlot validates its arguments", {
  df <- make_pcm(n = 150, k = 4)
  expect_error(RMitemCatFitPlot(df, type = "bars"))
  expect_error(RMitemCatFitPlot(df, score_groups = 1), "integer >= 2")
  expect_error(RMitemCatFitPlot(df, score_groups = 2.5), "integer >= 2")
  expect_error(RMitemCatFitPlot(df, score_groups = Inf), "integer >= 2")
  expect_error(RMitemCatFitPlot(df, type = "restscore", iterations = 5),
               "integer >= 20")
  expect_error(RMitemCatFitPlot(df, conf_level = 1.2), "in \\(0, 1\\)")
  expect_error(RMitemCatFitPlot(df, fit = list()), "result of RMitemCatFit")
})

test_that("RMitemCatFitPlot returns a ggplot for every type", {
  skip_if_not_installed("ggplot2")
  df <- make_pcm(n = 200, k = 4)
  for (tp in c("residuals", "adjacent", "restscore", "calibration")) {
    p <- suppressMessages(RMitemCatFitPlot(df, type = tp, iterations = 20, seed = 1))
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("RMitemCatFitPlot dataframe output has the documented columns", {
  df <- make_pcm(n = 200, k = 4)
  cols <- list(
    residuals = c("item", "score_group", "category", "observed", "expected", "z"),
    adjacent = c("item", "pair", "score", "n_pair", "observed", "expected",
                 "lower", "upper"),
    restscore = c("item", "statistic", "category", "observed", "lower", "upper",
                "iterations_used"),
    calibration = c("item", "event", "predicted", "cep", "n", "lower", "upper")
  )
  for (tp in names(cols)) {
    d <- RMitemCatFitPlot(df, type = tp, iterations = 20, seed = 1,
                          output = "dataframe")
    expect_identical(names(d), cols[[tp]], info = tp)
    expect_identical(levels(d$item), colnames(df), info = tp)
  }
})

test_that("residuals: observed and expected counts agree per item and category", {
  # The category counts are sufficient for the CML estimates, so summed over
  # score groups the residual counts vanish.
  df <- make_pcm(n = 300, k = 4)
  d <- RMitemCatFitPlot(df, type = "residuals", output = "dataframe")
  tot <- stats::aggregate(cbind(observed, expected) ~ item + category, d, sum)
  expect_equal(tot$observed, tot$expected, tolerance = 1e-4)
})

test_that("adjacent: binomial bands are ordered and within [0, 1]", {
  # The band is a pair of binomial quantiles, so with few respondents and an
  # expected share near 0 or 1 it need not contain the expected share.
  df <- make_pcm(n = 200, k = 4)
  d <- RMitemCatFitPlot(df, type = "adjacent", output = "dataframe")
  expect_true(all(d$lower <= d$upper))
  expect_true(all(d$lower >= 0 & d$upper <= 1))
  expect_true(all(d$observed >= 0 & d$observed <= 1))
  expect_equal(d$lower, stats::qbinom(0.025, d$n_pair, d$expected) / d$n_pair)
})

test_that("calibration and rest-score plots reuse the simulations of a fit", {
  df <- make_pcm(n = 200, k = 4)
  fit <- suppressMessages(RMitemCatFit(df, iterations = 30, seed = 9,
                                       output = "list"))
  # Same seed and iterations: the per-iteration draws are the same, so the
  # calibration band from the fit equals the band simulated by the plot.
  a <- RMitemCatFitPlot(df, type = "calibration", fit = fit, output = "dataframe")
  b <- RMitemCatFitPlot(df, type = "calibration", iterations = 30, seed = 9,
                        output = "dataframe")
  expect_equal(a, b)
  f <- RMitemCatFitPlot(df, type = "restscore", fit = fit, output = "dataframe")
  expect_true(all(f$iterations_used == fit$iterations_used))
})

test_that("RMitemCatFitPlot shows each item's adjusted p with a fit", {
  skip_if_not_installed("ggplot2")
  df <- make_pcm(n = 200, k = 4)
  fit <- suppressMessages(RMitemCatFit(df, iterations = 30, seed = 9,
                                       output = "list"))
  p <- RMitemCatFitPlot(df, type = "residuals", fit = fit)
  strips <- levels(p$data$item)
  expect_length(strips, 4L)
  expect_true(all(grepl("\n", strips, fixed = TRUE)))
  expect_true(all(startsWith(strips, colnames(df))))
  expect_match(p$labels$caption, "RMitemCatFit")
})

test_that("RMitemCatFitPlot stops when the fit comes from other data", {
  df <- make_pcm(n = 200, k = 4)
  fit <- suppressMessages(RMitemCatFit(df, iterations = 25, seed = 1,
                                       output = "list"))
  expect_error(RMitemCatFitPlot(make_pcm(n = 200, k = 4, seed = 2L),
                                fit = fit),
               "different data")
})
