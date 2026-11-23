# Tests for RMdifImpact() and its engine (R/utils-dif-impact.R)
#
# Covers:
#   - Shift CML: equals pcmodel() without DIF items, analytic gradient,
#     recovery of a planted shift, dichotomous items
#   - Population checks: no DIF gives no projected impact, projection
#     scales with the target size
#   - Sign convention and `reference`, items by name or position
#   - Rates: definitions, agreement with the observed classification,
#     keep - split, the rule
#   - Per-item impacts, a category empty in one group, missing responses
#   - Second review (2026-11-19): a shift with no finite estimate stops by
#     name, default target moments carry their uncertainty into the rates,
#     very low or high thresholds stay on the grids, a category lost with
#     the incomplete rows is reported as such
#   - Output types and argument checks
#
# The numbers were also checked against the validated development pipeline
# (dev/dif_impact_ext_sim*.R) on one dataset: every estimate and SE agreed
# to 1e-6 or better (dev/dif-impact-design.md, "Package function").

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
sim_pcm_impact <- function(theta, tau) {
  sapply(tau, function(t) {
    m <- length(t)
    eta <- outer(theta, 0:m) -
      matrix(c(0, cumsum(t)), length(theta), m + 1L, byrow = TRUE)
    eta <- eta - apply(eta, 1L, max)
    p <- exp(eta)
    p <- p / rowSums(p)
    u <- stats::runif(length(theta))
    rowSums(u > t(apply(p, 1L, cumsum))[, -(m + 1L), drop = FALSE])
  })
}

# 8 items with 3 categories, n per group, focal mean 0.4, DIF in I3 (and I6).
make_impact_data <- function(n = 400, d3 = 0.8, d6 = 0.5, seed = 11L) {
  set.seed(seed)
  tau <- lapply(seq(-1.4, 1.4, length.out = 8), function(l) l + c(-0.7, 0.7))
  tau2 <- tau
  tau2[[3]] <- tau2[[3]] + d3
  tau2[[6]] <- tau2[[6]] + d6
  x <- rbind(sim_pcm_impact(stats::rnorm(n), tau),
             sim_pcm_impact(stats::rnorm(n, 0.4), tau2))
  x <- as.data.frame(x)
  names(x) <- paste0("I", seq_len(ncol(x)))
  list(data = x, group = rep(c("A", "B"), each = n), tau = tau)
}

dd <- make_impact_data()
res1 <- RMdifImpact(dd$data, dd$group, items = "I3", cut = 1,
                    target = list(n = c(2000, 2000)), output = "list")
res2 <- RMdifImpact(dd$data, dd$group, items = c("I3", "I6"), cut = 9,
                    cut_scale = "sum", output = "list")

# ---------------------------------------------------------------------
# Shift CML engine
# ---------------------------------------------------------------------
shift_fit <- function(X, focal, items) {
  m <- apply(X, 2L, max)
  prep <- .dimp_prep(m, items, .dimp_data_stats(X, focal, m))
  list(prep = prep, fit = .dimp_fit(prep, hessian = TRUE))
}

test_that("without DIF items the shift CML equals pcmodel()", {
  X <- as.matrix(dd$data)
  storage.mode(X) <- "integer"
  sf <- shift_fit(X, dd$group == "B", integer(0))
  thr <- .dimp_par_to_thr(sf$fit$par, sf$prep, colnames(X))$ref
  ref <- .fit_cml_thresholds(X)
  expect_equal(unname(unlist(thr)), unname(unlist(ref)), tolerance = 1e-4)
})

test_that("the analytic gradient matches numerical derivatives", {
  X <- as.matrix(dd$data)
  storage.mode(X) <- "integer"
  m <- apply(X, 2L, max)
  prep <- .dimp_prep(m, c(3L, 6L), .dimp_data_stats(X, dd$group == "B", m))
  set.seed(1)
  par <- stats::rnorm(prep$P - 1L + 2L, 0, 0.5)
  num <- vapply(seq_along(par), function(p) {
    e <- 1e-5 * (seq_along(par) == p)
    (.dimp_negll(par + e, prep) - .dimp_negll(par - e, prep)) / 2e-5
  }, numeric(1L))
  expect_equal(.dimp_grad(par, prep), num, tolerance = 1e-5)
})

test_that("a planted shift is recovered within its CI", {
  big <- make_impact_data(n = 1500, seed = 5L)
  r <- RMdifImpact(big$data, big$group, items = c("I3", "I6"),
                   table = "dif", output = "dataframe")
  expect_true(r$lower[1] < 0.8 && 0.8 < r$upper[1])
  expect_true(r$lower[2] < 0.5 && 0.5 < r$upper[2])
})

test_that("dichotomous items work", {
  set.seed(4)
  b <- seq(-1.5, 1.5, length.out = 10)
  b2 <- b
  b2[2] <- b2[2] + 0.8
  x <- rbind(sim_pcm_impact(stats::rnorm(400), as.list(b)),
             sim_pcm_impact(stats::rnorm(400), as.list(b2)))
  x <- as.data.frame(x)
  names(x) <- paste0("Q", 1:10)
  r <- RMdifImpact(x, rep(c("a", "b"), each = 400), items = "Q2",
                   output = "list")
  expect_equal(nrow(r$impact), 3L)
  expect_true(r$dif$DIF > 0)
  expect_true(all(r$impact$Impact < 0))
})

# ---------------------------------------------------------------------
# Population checks
# ---------------------------------------------------------------------
test_that("no DIF gives no projected impact, DIF gives a negative one", {
  thr <- list(ref = dd$tau, focal = dd$tau)
  thr$ref <- .center_thresholds(thr$ref)
  thr$focal <- thr$ref
  tg <- list(n = c(1000, 1000), mean = c(0, 0.4), sd = c(1, 1))
  expect_equal(unname(.dimp_project(thr, tg)$impact), c(0, 0, 0),
               tolerance = 1e-6)
  thr_d <- .dimp_with_shifts(thr, 3L, 0.8)
  expect_true(all(.dimp_project(thr_d, tg)$impact < 0))
})

test_that("the projected impact depends on the group ratio, its SE on n", {
  a <- RMdifImpact(dd$data, dd$group, items = "I3",
                   target = list(n = c(1000, 1000)), table = "projection",
                   output = "dataframe")
  b <- RMdifImpact(dd$data, dd$group, items = "I3",
                   target = list(n = c(4000, 4000)), table = "projection",
                   output = "dataframe")
  expect_equal(a$Impact, b$Impact, tolerance = 1e-6)
  expect_equal(a$SE_Delta / b$SE_Delta, rep(2, 3), tolerance = 1e-6)
  expect_equal(a$Impact_in_SE * 2, b$Impact_in_SE, tolerance = 1e-6)
})

test_that("target means and SDs default to the estimates", {
  s <- res1$settings
  expect_equal(c(s$Target_mean_ref, s$Target_mean_focal),
               c(s$Mean_ref, s$Mean_focal))
  expect_equal(c(s$Target_sd_ref, s$Target_sd_focal),
               c(s$SD_ref, s$SD_focal))
})

test_that("default target moments carry their uncertainty into the rates", {
  # With the target's means and SDs at the estimates, the target rates are
  # the sample rates, with the same SEs.
  r <- res1$rates
  s <- r[r$Population == "sample", ]
  t <- r[r$Population == "target", ]
  expect_equal(t$Estimate, s$Estimate)
  expect_equal(t$SE, s$SE)
})

# ---------------------------------------------------------------------
# Sign convention, reference, items
# ---------------------------------------------------------------------
test_that("DIF that makes an item harder for the focal group biases it down", {
  expect_true(res1$dif$DIF > 0)
  expect_true(all(res1$impact$Impact < 0))
  expect_true(all(res1$impact$upper < 0))
  expect_true(res1$impact$Delta_split[1] > res1$impact$Delta_keep[1])
  expect_equal(res1$settings$Comparison, "B minus A")
})

test_that("swapping the reference flips the signs", {
  r <- RMdifImpact(dd$data, dd$group, items = "I3", reference = "B",
                   output = "list")
  expect_equal(r$dif$DIF, -res1$dif$DIF, tolerance = 1e-4)
  # Latent and WLE means. The sum score under split equates the focal
  # group, which is now the other group, so it is not exactly symmetric.
  expect_equal(r$impact$Delta_keep[1:2], -res1$impact$Delta_keep[1:2],
               tolerance = 1e-4)
  expect_equal(r$impact$Delta_split[1:2], -res1$impact$Delta_split[1:2],
               tolerance = 1e-3)
})

test_that("items by position equal items by name", {
  a <- RMdifImpact(dd$data, dd$group, items = c(6, 3), output = "dataframe")
  b <- RMdifImpact(dd$data, dd$group, items = c("I3", "I6"),
                   output = "dataframe")
  expect_equal(a, b)
})

test_that("per-item impacts have the joint impact's sign and are smaller", {
  it <- res2$items
  expect_equal(nrow(it), 6L)
  expect_setequal(unique(it$Item), c("I3", "I6"))
  lat <- it[it$Estimator == "Latent mean", ]
  joint <- res2$impact$Impact[res2$impact$Estimator == "Latent mean"]
  # Not additive: making one item common moves Delta by more than the
  # linear gradient predicts, so the parts add up to somewhat more.
  expect_true(all(sign(lat$Impact) == sign(joint)))
  expect_true(all(abs(lat$Impact) < abs(joint)))
  expect_true(abs(sum(lat$Impact)) > abs(joint))
  expect_true(all(it$SE > 0))
  expect_null(res1$items)
})

# ---------------------------------------------------------------------
# Rates
# ---------------------------------------------------------------------
test_that("rates are shares and keep - split is their difference", {
  rt <- res1$rates
  expect_setequal(unique(rt$Population), c("sample", "target"))
  shares <- rt[rt$Strategy != "keep - split", ]
  expect_true(all(shares$Estimate >= 0 & shares$Estimate <= 1))
  for (pp in unique(rt$Population)) for (g in c("A", "B")) {
    for (rr in unique(rt$Rate)) {
      v <- function(st) {
        rt$Estimate[rt$Population == pp & rt$Group == g & rt$Rate == rr &
                      rt$Strategy == st]
      }
      expect_equal(v("keep - split"), v("keep") - v("split"))
    }
  }
})

test_that("model-based share at or above the cut matches the data", {
  x <- dd$data
  focal <- dd$group == "B"
  r <- rowSums(x)
  # Sum-score cut: keep classifies raw >= cut in both groups.
  rt <- res2$rates
  obs <- c(A = mean(r[!focal] >= 9), B = mean(r[focal] >= 9))
  est <- rt$Estimate[rt$Population == "sample" & rt$Strategy == "keep" &
                       rt$Rate == "at or above"]
  names(est) <- rt$Group[rt$Population == "sample" & rt$Strategy == "keep" &
                           rt$Rate == "at or above"]
  expect_equal(unname(est[c("A", "B")]), unname(obs), tolerance = 0.05)
})

test_that("the rule is reported and the reference sum-score rule is fixed", {
  ru <- res2$rule
  expect_equal(nrow(ru), 4L)
  expect_equal(ru$Lowest_raw_score[ru$Strategy == "keep"], c(9, 9))
  expect_equal(ru$Lowest_raw_score[ru$Strategy == "split" & ru$Group == "A"], 9)
  ref_diff <- res2$rates[res2$rates$Group == "A" &
                           res2$rates$Strategy == "keep - split", ]
  expect_equal(ref_diff$Estimate, rep(0, 3))
  expect_true("WLE" %in% names(res1$rule))
  expect_false("WLE" %in% names(res2$rule))
})

test_that("a tiny DIF in either direction leaves the sum-score rule alone", {
  thr0 <- .center_thresholds(dd$tau)
  for (d in c(-0.01, 0.01)) {
    thr <- .dimp_with_shifts(list(ref = thr0, focal = thr0), 3L, d)
    for (cut in 2:15) {
      pt <- .dimp_rate_parts(thr, 3L, cut)
      ru <- .dimp_rule(pt, NULL, "sum", cut)
      expect_identical(ru$split$focal, ru$keep$focal)
    }
  }
  # A shift of one whole raw point moves it.
  thr <- .dimp_with_shifts(list(ref = thr0, focal = thr0), 3L, 3)
  pt <- .dimp_rate_parts(thr, 3L, 8)
  ru <- .dimp_rule(pt, NULL, "sum", 8)
  expect_false(identical(ru$split$focal, ru$keep$focal))
})

test_that("keep and split classify alike when their parameters agree", {
  thr0 <- lapply(dd$tau, function(t) t - mean(unlist(dd$tau)))
  names(thr0) <- names(dd$data)
  thr <- list(ref = thr0, focal = thr0)
  pt <- .dimp_rate_parts(thr, 3L)
  wk <- .dimp_wle_by_score(thr0)
  # Cuts at the WLE of each raw score, where two solvers could disagree.
  for (cut in wk[2:16]) {
    ru <- .dimp_rule(pt, wk, "theta", cut)
    expect_identical(ru$split$ref, ru$keep$ref)
    expect_identical(ru$split$focal, ru$keep$focal)
  }
})

test_that("local bias at a cut is reported where the cut is", {
  lc <- res1$local
  expect_equal(nrow(lc), 4L)
  expect_equal(lc$Theta[3], 1)
  # The curve and the table share one grid: the tabled maximum is the
  # curve's, signed, with the same SE.
  cv <- res1$curve
  i <- which.max(abs(cv$b))
  expect_equal(lc$Estimate[1], cv$b[i])
  expect_equal(lc$Theta[1], cv$Theta[i])
  expect_equal(lc$SE[1], cv$SE[i], tolerance = 1e-6)
  expect_true(lc$Estimate[1] < 0)  # DIF makes I3 harder for the focal group
  r <- RMdifImpact(dd$data, dd$group, items = "I3", table = "local",
                   output = "dataframe")
  expect_equal(nrow(r), 2L)
})

# ---------------------------------------------------------------------
# Data handling
# ---------------------------------------------------------------------
test_that("a category empty in one group is allowed", {
  x <- dd$data
  focal <- dd$group == "B"
  x$I3[focal & x$I3 == 2] <- 1
  expect_error(suppressMessages(RMdifAnchor(x, dd$group)), "category")
  r <- RMdifImpact(x, dd$group, items = "I3", output = "dataframe")
  expect_equal(nrow(r), 3L)
  expect_true(all(is.finite(r$SE)))
})

test_that("a shift with no finite estimate stops by name", {
  x <- dd$data
  x$I3[dd$group == "B"] <- 0L
  expect_error(RMdifImpact(x, dd$group, items = "I3"),
               "I3 (lowest possible scores in group 'B')", fixed = TRUE,
               class = "easyRasch2_category_error")
  x$I3[dd$group == "B"] <- 2L
  expect_error(RMdifImpact(x, dd$group, items = "I3"),
               "I3 (highest possible scores in group 'B')", fixed = TRUE)
})

test_that("a category lost with the incomplete rows is reported as such", {
  x <- dd$data
  x$I2[x$I1 == 1] <- NA
  expect_error(RMdifImpact(x, dd$group, items = "I3"),
               "RMdifImpact() uses complete cases only", fixed = TRUE,
               class = "easyRasch2_category_error")
})

test_that("very low or high thresholds stay on the grids", {
  expect_equal(.dimp_span(list(-3, c(2, 4.9))), c(-10, 10))
  thr <- list(ref = list(-12, -11, -10, 0, 10, 11, 12, 1),
              focal = list(-12, -11.2, -10, 0, 10, 11, 12, 1))
  expect_equal(.dimp_span(thr$ref, thr$focal), c(-20, 20))
  pt <- .dimp_rate_parts(thr, 2L, cut_raw = 1)
  expect_true(is.finite(pt$theta_c_sum))
  expect_equal(.dimp_moments(thr$ref, pt$theta_c_sum)$E, 0.5,
               tolerance = 1e-3)
  # WLEs of the extreme scores are found, not clipped at the range ends
  w <- .dimp_wle_by_score(thr$ref)
  expect_true(all(abs(w) < 19))
})

test_that("respondents with missing responses are dropped with a message", {
  x <- dd$data
  x$I1[1:5] <- NA
  expect_message(
    r <- RMdifImpact(x, dd$group, items = "I3", output = "list"),
    "5 respondent"
  )
  expect_equal(r$settings$N, nrow(x) - 5L)
  expect_equal(r$settings$N_excluded, 5L)
  expect_equal(r$settings$N_incomplete, 5L)
})

test_that("the caption says complete cases only when incomplete rows went", {
  g <- dd$group
  g[1:3] <- NA
  k <- suppressMessages(RMdifImpact(dd$data, g, items = "I3"))
  expect_true(any(grepl("of 800 respondents", k)))
  expect_false(any(grepl("complete cases", k)))
  x <- dd$data
  x$I1[1:3] <- NA
  k <- suppressMessages(RMdifImpact(x, dd$group, items = "I3"))
  expect_true(any(grepl("complete cases", k)))
})

# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------
test_that("output types", {
  expect_named(res1, c("impact", "dif", "items", "local", "rates", "rule",
                       "projection", "curve", "settings", "kable"))
  expect_identical(res1$kable, RMdifImpact(dd$data, dd$group, items = "I3",
                   cut = 1, target = list(n = c(2000, 2000))))
  expect_s3_class(RMdifImpact(dd$data, dd$group, items = "I3"),
                  "knitr_kable")
  k <- RMdifImpact(dd$data, dd$group, items = "I3", cut = 1,
                   table = "rates")
  expect_s3_class(k, "knitr_kable")
  expect_true(any(grepl("theta = 1", k)))
  skip_if_not_installed("ggplot2")
  p <- RMdifImpact(dd$data, dd$group, items = "I3", cut = 1,
                   output = "ggplot")
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  # `table` does not apply to the plot.
  expect_s3_class(RMdifImpact(dd$data, dd$group, items = "I3",
                              table = "rates", output = "ggplot"), "ggplot")
})

test_that("nothing is attached as attributes", {
  df <- RMdifImpact(dd$data, dd$group, items = "I3", output = "dataframe")
  expect_setequal(names(attributes(df)), c("names", "class", "row.names"))
})

# ---------------------------------------------------------------------
# Argument checks
# ---------------------------------------------------------------------
test_that("RMdifImpact validates its arguments", {
  x <- dd$data
  g <- dd$group
  expect_error(RMdifImpact(x, g), "items")
  expect_error(RMdifImpact(x, rep("A", nrow(x)), items = "I3"), "exactly 2")
  expect_error(RMdifImpact(x, rep(c("A", "B", "C"), length.out = nrow(x)),
                           items = "I3"), "exactly 2")
  expect_error(RMdifImpact(x, g, items = "Z1"), "not found")
  expect_error(RMdifImpact(x, g, items = 99), "between 1 and")
  expect_error(RMdifImpact(x, g, items = c("I3", "I3")), "more than once")
  expect_error(RMdifImpact(x, g, items = names(x)), "common")
  expect_error(RMdifImpact(x, g, items = "I3", reference = "Z"), "reference")
  expect_error(RMdifImpact(x, g, items = "I3", cut = c(1, 2)), "cut")
  expect_error(RMdifImpact(x, g, items = "I3", cut = 2.5, cut_scale = "sum"),
               "whole number")
  expect_error(RMdifImpact(x, g, items = "I3", cut = 99, cut_scale = "sum"),
               "whole number")
  expect_error(RMdifImpact(x, g, items = "I3", conf = 1.2), "conf")
  expect_error(RMdifImpact(x, g, items = "I3", table = "rates"), "cut")
  expect_error(RMdifImpact(x, g, items = "I3", table = "projection"),
               "target")
  expect_error(RMdifImpact(x, g, items = "I3", table = "items"),
               "two or more")
  expect_error(RMdifImpact(x, g, items = "I3", target = list(N = 10)),
               "`n`")
  expect_error(RMdifImpact(x, g, items = "I3", target = list(n = 5)),
               "two finite")
  expect_error(RMdifImpact(x, g, items = "I3",
                           target = list(n = c(10, 10), sd = c(1, -1))),
               "positive")
  expect_error(RMdifImpact(x, g, items = "I3",
                           target = list(n = c(A = 10, Z = 10))),
               "names")
})

test_that("a target named by group is put in reference, focal order", {
  a <- RMdifImpact(dd$data, dd$group, items = "I3",
                   target = list(n = c(B = 1000, A = 3000)), output = "list")
  expect_equal(c(a$settings$Target_n_ref, a$settings$Target_n_focal),
               c(3000, 1000))
})
