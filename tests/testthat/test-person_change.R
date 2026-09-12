# Tests for RMpersonChange()

make_thr <- function(k = 6) {
  thr <- lapply(seq(-1.2, 1.2, length.out = k), function(b) b + c(-0.8, 0, 0.8))
  names(thr) <- paste0("I", seq_len(k))
  thr
}

# Two occasions with a known shift applied to every respondent
make_pair <- function(n = 120, k = 6, shift = 0.8, seed = 1L, theta_sd = 1.4) {
  set.seed(seed)
  thr <- make_thr(k)
  theta <- stats::rnorm(n, 0, theta_sd)
  t1 <- as.data.frame(sim_partial_score(thr, theta))
  t2 <- as.data.frame(sim_partial_score(thr, theta + shift))
  colnames(t1) <- colnames(t2) <- names(thr)
  list(t1 = t1, t2 = t2, thr = thr, theta = theta)
}

flat_caption <- function(p) gsub("[[:space:]]+", " ", p$labels$caption)

# ---------------------------------------------------------------------
# Argument handling
# ---------------------------------------------------------------------
test_that("values passed to the wrong argument name the right one", {
  d <- make_pair(n = 30)
  expect_error(RMpersonChange(d$t1, d$t2, method = "CML"),
               regexp = "belongs to `estimator`")
  expect_error(RMpersonChange(d$t1, d$t2, estimator = "WLE"),
               regexp = "belongs to `method`")
  expect_error(RMpersonChange(d$t1, d$t2, null = "two.sided"),
               regexp = "belongs to `direction`")
  expect_error(RMpersonChange(d$t1, d$t2, direction = "retest"),
               regexp = "belongs to `null`")
})

test_that("null = 'retest' requires retest_sd and explains why", {
  d <- make_pair(n = 30)
  err <- tryCatch(RMpersonChange(d$t1, d$t2, null = "retest"),
                  error = conditionMessage)
  expect_match(err, "per-occasion SD")
  expect_match(err, "RMretestSD")
  # it must say why these data cannot supply it
  expect_match(err, "that is the change being tested")
})

test_that("retest_sd is refused under the measurement null", {
  d <- make_pair(n = 30)
  expect_error(RMpersonChange(d$t1, d$t2, retest_sd = 0.3),
               regexp = "only applies when")
})

test_that("mismatched occasions are rejected", {
  d <- make_pair(n = 30)
  expect_error(RMpersonChange(d$t1, d$t2[1:10, ]),
               regexp = "same number of rows")
  expect_error(RMpersonChange(d$t1, d$t2[, 1:3]),
               regexp = "same items")
  swapped <- d$t2[, c(2, 1, 3, 4, 5, 6)]
  expect_error(RMpersonChange(d$t1, swapped),
               regexp = "same item names in the same order")
})

test_that("critical and alpha are validated", {
  d <- make_pair(n = 30)
  expect_error(RMpersonChange(d$t1, d$t2, critical = -1), regexp = "positive")
  expect_error(RMpersonChange(d$t1, d$t2, critical = "boot"),
               regexp = "\"simulate\" or a single positive number")
  expect_error(RMpersonChange(d$t1, d$t2, alpha = 0), regexp = "between 0 and 1")
})

test_that("id length is checked", {
  d <- make_pair(n = 30)
  expect_error(RMpersonChange(d$t1, d$t2, id = letters[1:3], critical = 1.96),
               regexp = "one entry per respondent")
})

# ---------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------
test_that("too few respondents to calibrate points at item_params", {
  d <- make_pair(n = 120)
  err <- tryCatch(RMpersonChange(d$t1[1:3, ], d$t2[1:3, ]),
                  error = conditionMessage)
  expect_match(err, "cannot be estimated from 3 respondent")
  expect_match(err, "single-subject")
})

test_that("a single respondent works with supplied item_params", {
  skip_on_cran()
  d <- make_pair(n = 120)
  res <- RMpersonChange(d$t1[1, ], d$t2[1, ], item_params = d$thr,
                        critical = "simulate", sim_iter = 50, parallel = FALSE, seed = 1)
  expect_equal(nrow(res), 1L)
  expect_true(is.finite(res$rci))
  expect_true(is.finite(res$p_value))
  expect_identical(attr(res, "anchor"), "supplied")
})

test_that("item_params overrides anchor, and says so only when anchor was set", {
  skip_on_cran()
  d <- make_pair(n = 60)
  expect_message(
    RMpersonChange(d$t1, d$t2, anchor = "t1", item_params = d$thr,
                   critical = 1.96),
    regexp = "`anchor` is ignored"
  )
  expect_silent(
    RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = 1.96)
  )
})

test_that("the three anchors give different calibrations but the same shape", {
  skip_on_cran()
  d <- make_pair(n = 150)
  out <- lapply(c("stack", "t1", "t2"), function(a) {
    RMpersonChange(d$t1, d$t2, anchor = a, critical = 1.96)
  })
  expect_true(all(vapply(out, nrow, integer(1L)) == 150L))
  expect_false(isTRUE(all.equal(out[[1L]]$theta_t1, out[[2L]]$theta_t1)))
  for (i in seq_along(out)) {
    expect_identical(attr(out[[i]], "anchor"), c("stack", "t1", "t2")[i])
  }
})

# ---------------------------------------------------------------------
# The statistic
# ---------------------------------------------------------------------
test_that("change, se_diff and rci are the documented transforms", {
  skip_on_cran()
  d <- make_pair(n = 120)
  res <- RMpersonChange(d$t1, d$t2, critical = 1.96)
  expect_equal(res$change, res$theta_t2 - res$theta_t1)
  expect_equal(res$se_diff, sqrt(res$se_t1^2 + res$se_t2^2))
  expect_equal(res$rci, res$change / res$se_diff)
  expect_identical(attr(res, "retest_sd"), 0)
})

test_that("retest_sd enters doubled, as a per-occasion SD", {
  skip_on_cran()
  d <- make_pair(n = 120)
  a <- RMpersonChange(d$t1, d$t2, critical = 1.96)
  b <- RMpersonChange(d$t1, d$t2, null = "retest", retest_sd = 0.4,
                      critical = 1.96)
  expect_equal(b$se_diff, sqrt(a$se_t1^2 + a$se_t2^2 + 2 * 0.4^2))
  expect_true(all(b$se_diff > a$se_diff))
  moved <- a$change != 0            # a zero change stays zero under any SE
  expect_true(all(abs(b$rci[moved]) < abs(a$rci[moved])))
})

test_that("a real shift is detected more often than no shift", {
  skip_on_cran()
  shifted <- make_pair(n = 200, shift = 1.2, seed = 3)
  null_pair <- make_pair(n = 200, shift = 0, seed = 3)
  a <- RMpersonChange(shifted$t1, shifted$t2, critical = 1.96)
  b <- RMpersonChange(null_pair$t1, null_pair$t2, critical = 1.96)
  n_flag <- function(x) sum(x$change_class != "none detected", na.rm = TRUE)
  expect_gt(n_flag(a), n_flag(b))
  # and the direction is the right one
  expect_gt(sum(a$change_class == "increase"),
            sum(a$change_class == "decrease"))
})

# ---------------------------------------------------------------------
# Critical values and classification
# ---------------------------------------------------------------------
test_that("numeric critical gives a symmetric cutoff and normal p-values", {
  skip_on_cran()
  d <- make_pair(n = 120)
  res <- RMpersonChange(d$t1, d$t2, critical = 1.96)
  expect_true(all(res$crit_lower == -1.96))
  expect_true(all(res$crit_upper == 1.96))
  expect_equal(res$p_value, 2 * stats::pnorm(-abs(res$rci)))
  expect_identical(attr(res, "critical"), 1.96)
})

test_that("the simulated null is symmetric when occasions are exchangeable", {
  skip_on_cran()
  d <- make_pair(n = 150)
  res <- RMpersonChange(d$t1, d$t2, critical = "simulate", sim_iter = 300, parallel = FALSE, seed = 2)
  # theta_1 and theta_2 are iid given theta under the null, so their
  # difference is symmetric about zero however skewed theta-hat is
  expect_equal(abs(res$crit_lower[1L]), res$crit_upper[1L], tolerance = 1e-8)
})

test_that("unequal missingness across occasions breaks that symmetry", {
  skip_on_cran()
  d <- make_pair(n = 200)
  d$t2[, 5:6] <- NA                  # occasion 2 carries less information
  res <- RMpersonChange(d$t1, d$t2, critical = "simulate", sim_iter = 400, parallel = FALSE, seed = 5)
  expect_gt(abs(abs(res$crit_lower[1L]) - res$crit_upper[1L]), 0.02)
})

test_that("the simulated critical value is tighter than the normal 1.96", {
  skip_on_cran()
  d <- make_pair(n = 200)
  res <- RMpersonChange(d$t1, d$t2, critical = "simulate", sim_iter = 400, parallel = FALSE, seed = 2)
  expect_lt(res$crit_upper[1L], 1.96)
  # and that difference is what makes simulating worth the cost
  fixed <- RMpersonChange(d$t1, d$t2, critical = 1.96)
  expect_gt(sum(res$change_class != "none detected"),
            sum(fixed$change_class != "none detected"))
})

test_that("conditional_crit gives per-respondent critical values", {
  skip_on_cran()
  d <- make_pair(n = 80)
  pooled <- RMpersonChange(d$t1, d$t2, critical = "simulate", sim_iter = 200, parallel = FALSE,
                           seed = 2)
  cond <- RMpersonChange(d$t1, d$t2, conditional_crit = TRUE, critical = "simulate", sim_iter = 200,
                         parallel = FALSE, seed = 2)
  expect_length(unique(pooled$crit_upper), 1L)
  expect_gt(length(unique(cond$crit_upper)), 1L)
})

test_that("one-sided tests put all alpha in one tail", {
  skip_on_cran()
  d <- make_pair(n = 120)
  up <- RMpersonChange(d$t1, d$t2, direction = "increase", critical = "simulate", sim_iter = 200,
                       parallel = FALSE, seed = 2)
  dn <- RMpersonChange(d$t1, d$t2, direction = "decrease", critical = "simulate", sim_iter = 200,
                       parallel = FALSE, seed = 2)
  expect_true(all(is.infinite(up$crit_lower) & up$crit_lower < 0))
  expect_true(all(is.finite(up$crit_upper)))
  expect_true(all(is.infinite(dn$crit_upper) & dn$crit_upper > 0))
  expect_true(all(dn$change_class != "increase"))
  # one-sided is more powerful in its own direction
  two <- RMpersonChange(d$t1, d$t2, critical = "simulate", sim_iter = 200, parallel = FALSE, seed = 2)
  expect_gte(sum(up$change_class == "increase"),
             sum(two$change_class == "increase"))
})

test_that("class levels name the direction of theta, not a clinical reading", {
  skip_on_cran()
  d <- make_pair(n = 60)
  res <- RMpersonChange(d$t1, d$t2, critical = 1.96)
  expect_identical(
    levels(res$change_class),
    c("decrease", "none detected", "increase")
  )
})

test_that("classification agrees with the critical values", {
  skip_on_cran()
  d <- make_pair(n = 150)
  res <- RMpersonChange(d$t1, d$t2, critical = "simulate", sim_iter = 200, parallel = FALSE, seed = 2)
  inc <- res$change_class == "increase"
  dec <- res$change_class == "decrease"
  expect_true(all(res$rci[inc] > res$crit_upper[inc]))
  expect_true(all(res$rci[dec] < res$crit_lower[dec]))
  none <- res$change_class == "none detected"
  expect_true(all(res$rci[none] <= res$crit_upper[none] &
                    res$rci[none] >= res$crit_lower[none]))
})

test_that("the simulated null has roughly the intended error rate", {
  skip_on_cran()
  # No true change, so flagging should sit near alpha
  d <- make_pair(n = 400, shift = 0, seed = 11)
  res <- RMpersonChange(d$t1, d$t2, alpha = 0.05, critical = "simulate", sim_iter = 400,
                        parallel = FALSE, seed = 11)
  rate <- mean(res$change_class != "none detected", na.rm = TRUE)
  expect_lt(rate, 0.12)
})

# ---------------------------------------------------------------------
# Exact null
# ---------------------------------------------------------------------
test_that("exact is the default and is recorded as such", {
  skip_on_cran()
  d <- make_pair(n = 80)
  res <- RMpersonChange(d$t1, d$t2, item_params = d$thr)
  expect_identical(attr(res, "critical"), "exact")
  expect_true(is.na(attr(res, "sim_iter")))
})

test_that("exact agrees with a long simulation on complete data", {
  skip_on_cran()
  d <- make_pair(n = 120)
  e <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "exact")
  s <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "simulate",
                      sim_iter = 3000, parallel = FALSE, seed = 2)
  expect_equal(e$crit_upper[1L], s$crit_upper[1L], tolerance = 0.05)
  expect_equal(e$crit_lower[1L], s$crit_lower[1L], tolerance = 0.05)
  expect_gt(stats::cor(e$p_value, s$p_value), 0.99)
})

test_that("the exact null is symmetric under exchangeable occasions", {
  skip_on_cran()
  d <- make_pair(n = 100)
  e <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "exact")
  # no Monte Carlo error, so this holds to machine precision
  expect_equal(abs(e$crit_lower[1L]), e$crit_upper[1L])
})

test_that("the exact critical value rises with test length toward 1.96", {
  skip_on_cran()
  crit <- vapply(c(5L, 10L, 20L), function(k) {
    d <- make_pair(n = 150, k = k, seed = 4L)
    RMpersonChange(d$t1, d$t2, item_params = d$thr,
                   critical = "exact")$crit_upper[1L]
  }, numeric(1L))
  expect_true(all(diff(crit) > 0))
  expect_true(all(crit < 1.96))
})

test_that("exact handles incomplete data by enumerating within patterns", {
  skip_on_cran()
  d <- make_pair(n = 150)
  d$t1[1:20, 1] <- NA
  d$t2[10:30, 3] <- NA
  e <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "exact")
  s <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "simulate",
                      sim_iter = 3000, parallel = FALSE, seed = 2)
  expect_true(all(is.finite(e$crit_upper)))
  expect_equal(e$crit_upper[1L], s$crit_upper[1L], tolerance = 0.06)
})

test_that("exact supports the retest null by integrating the deviations out", {
  skip_on_cran()
  d <- make_pair(n = 120)
  e <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "exact",
                      null = "retest", retest_sd = 0.3)
  s <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "simulate",
                      null = "retest", retest_sd = 0.3, sim_iter = 3000,
                      parallel = FALSE, seed = 2)
  expect_equal(e$crit_upper[1L], s$crit_upper[1L], tolerance = 0.05)
})

test_that("exact and simulate agree under EAP, with the prior held fixed", {
  skip_on_cran()
  d <- make_pair(n = 120)
  e <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "exact",
                      method = "EAP")
  s <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "simulate",
                      method = "EAP", sim_iter = 2000, parallel = FALSE,
                      seed = 2)
  # the simulated null must be scored under the same prior as the observation
  expect_equal(e$crit_upper[1L], s$crit_upper[1L], tolerance = 0.05)
})

test_that("exact conditional critical values vary and need no iterations", {
  skip_on_cran()
  d <- make_pair(n = 100)
  cond <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "exact",
                         conditional_crit = TRUE)
  pooled <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "exact")
  expect_gt(length(unique(round(cond$crit_upper, 6))), 1L)
  expect_length(unique(pooled$crit_upper), 1L)
})

test_that("exact one-sided tests keep the infinite bound", {
  skip_on_cran()
  d <- make_pair(n = 80)
  up <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "exact",
                       direction = "increase")
  expect_true(all(is.infinite(up$crit_lower) & up$crit_lower < 0))
  expect_true(all(is.finite(up$crit_upper)))
})

test_that("exact p-values are probabilities and track the classification", {
  skip_on_cran()
  d <- make_pair(n = 120)
  e <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "exact",
                      alpha = 0.05)
  expect_true(all(e$p_value >= 0 & e$p_value <= 1, na.rm = TRUE))
  flagged <- e$change_class != "none detected"
  expect_true(all(e$p_value[flagged] <= 0.05 + 1e-8))
})

test_that("critical is validated against the three accepted forms", {
  d <- make_pair(n = 30)
  expect_error(RMpersonChange(d$t1, d$t2, critical = "boot"),
               regexp = "\"exact\", \"simulate\" or a single positive")
})

test_that("the caption says the null was enumerated", {
  skip_on_cran()
  d <- make_pair(n = 60)
  k <- RMpersonChange(d$t1, d$t2, item_params = d$thr, critical = "exact",
                      output = "kable")
  expect_match(gsub("[[:space:]]+", " ", paste(as.character(k), collapse = " ")),
               "computed exactly by enumerating the null")
})

# ---------------------------------------------------------------------
# Tipping point
# ---------------------------------------------------------------------
test_that("retest_sd_tip is reported only for flagged respondents", {
  skip_on_cran()
  d <- make_pair(n = 150, shift = 1.2, seed = 4)
  res <- RMpersonChange(d$t1, d$t2, critical = 1.96)
  flagged <- res$change_class != "none detected"
  expect_true(all(is.finite(res$retest_sd_tip[flagged])))
  expect_true(all(is.na(res$retest_sd_tip[!flagged])))
  expect_true(all(res$retest_sd_tip[flagged] >= 0))
})

test_that("the tipping point is the SD that returns the RCI to the cutoff", {
  skip_on_cran()
  d <- make_pair(n = 150, shift = 1.2, seed = 4)
  res <- RMpersonChange(d$t1, d$t2, critical = 1.96)
  f <- which(res$change_class != "none detected" & res$retest_sd_tip > 0)
  skip_if(length(f) == 0L)
  i <- f[1L]
  se_at_tip <- sqrt(res$se_t1[i]^2 + res$se_t2[i]^2 + 2 * res$retest_sd_tip[i]^2)
  expect_equal(abs(res$change[i]) / se_at_tip, 1.96, tolerance = 1e-6)
})

test_that("a tipping SD of zero means no occasion noise is tolerated", {
  skip_on_cran()
  d <- make_pair(n = 200, shift = 1.2, seed = 4)
  res <- RMpersonChange(d$t1, d$t2, critical = 1.96)
  zero <- which(res$change_class != "none detected" & res$retest_sd_tip == 0)
  skip_if(length(zero) == 0L)
  # such a respondent sits essentially on the cutoff
  expect_true(all(abs(res$rci[zero]) - 1.96 < 0.05))
})

# ---------------------------------------------------------------------
# Missing data and extremes
# ---------------------------------------------------------------------
test_that("extreme scorers are flagged at each occasion", {
  skip_on_cran()
  d <- make_pair(n = 150, theta_sd = 2.5, seed = 5)
  res <- RMpersonChange(d$t1, d$t2, critical = 1.96)
  expect_type(res$extreme_t1, "logical")
  expect_true(any(res$extreme_t1 | res$extreme_t2))
  mx <- 3 * 6
  expect_true(all(res$sum_t1[res$extreme_t1] %in% c(0, mx)))
})

test_that("partial missingness is retained and simulated with", {
  skip_on_cran()
  d <- make_pair(n = 120)
  d$t1[1:15, 1] <- NA
  d$t2[5:20, 2] <- NA
  res <- RMpersonChange(d$t1, d$t2, critical = "simulate", sim_iter = 100, parallel = FALSE, seed = 6)
  expect_equal(nrow(res), 120L)
  expect_true(all(is.finite(res$se_diff)))
})

# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------
test_that("dataframe output has the documented columns and attributes", {
  skip_on_cran()
  d <- make_pair(n = 60)
  res <- RMpersonChange(d$t1, d$t2, critical = 1.96)
  expect_named(res, c(
    "id", "sum_t1", "sum_t2", "theta_t1", "se_t1", "theta_t2", "se_t2",
    "extreme_t1", "extreme_t2", "change", "se_diff", "rci", "p_value",
    "crit_lower", "crit_upper", "change_class", "retest_sd_tip"
  ))
  for (a in c("null", "retest_sd", "anchor", "alpha", "direction",
              "critical", "method", "conditional_crit")) {
    expect_false(is.null(attr(res, a)))
  }
})

test_that("kable and ggplot captions name the null in force", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  d <- make_pair(n = 60)

  k <- RMpersonChange(d$t1, d$t2, critical = 1.96, output = "kable")
  txt <- gsub("[[:space:]]+", " ", paste(as.character(k), collapse = " "))
  expect_match(txt, "no change beyond measurement error")
  expect_match(txt, "necessary condition for change")
  expect_false(grepl("uncorrected", txt))

  kr <- RMpersonChange(d$t1, d$t2, null = "retest", retest_sd = 0.3,
                       critical = 1.96, output = "kable")
  txtr <- gsub("[[:space:]]+", " ", paste(as.character(kr), collapse = " "))
  expect_match(txtr, "occasion-to-occasion fluctuation")

  p <- RMpersonChange(d$t1, d$t2, critical = 1.96, output = "ggplot")
  expect_s3_class(p, "ggplot")
  expect_match(flat_caption(p), "no change beyond measurement error")
})

test_that("the caption reports the critical values actually in force", {
  skip_on_cran()
  d <- make_pair(n = 60)
  flat <- function(k) gsub("[[:space:]]+", " ", paste(as.character(k),
                                                      collapse = " "))

  fixed <- RMpersonChange(d$t1, d$t2, critical = 1.96, output = "kable")
  expect_match(flat(fixed), "fixed at -1.96 and 1.96")

  res <- RMpersonChange(d$t1, d$t2, critical = "simulate", sim_iter = 100, parallel = FALSE, seed = 1)
  simmed <- RMpersonChange(d$t1, d$t2, critical = "simulate", sim_iter = 100, parallel = FALSE,
                           seed = 1, output = "kable")
  expect_match(flat(simmed), "simulated over 100 iterations")
  # the numbers themselves must be readable off the caption
  expect_match(flat(simmed), sprintf("%.2f and %.2f", res$crit_lower[1L],
                                     res$crit_upper[1L]))
})

test_that("a one-sided caption shows the infinite bound rather than dropping it", {
  skip_on_cran()
  d <- make_pair(n = 60)
  k <- RMpersonChange(d$t1, d$t2, direction = "increase", critical = 1.96,
                      output = "kable")
  expect_match(gsub("[[:space:]]+", " ", paste(as.character(k), collapse = " ")),
               "-Inf and 1.96")
})

test_that("per-respondent critical values reach the table, constant ones do not", {
  skip_on_cran()
  d <- make_pair(n = 60)
  pooled <- RMpersonChange(d$t1, d$t2, critical = "simulate", sim_iter = 100, parallel = FALSE,
                           seed = 1, output = "kable")
  cond <- RMpersonChange(d$t1, d$t2, conditional_crit = TRUE, critical = "simulate", sim_iter = 100,
                         parallel = FALSE, seed = 1, output = "kable")
  expect_false(grepl("Crit lower", paste(as.character(pooled), collapse = " ")))
  expect_match(paste(as.character(cond), collapse = " "), "Crit lower")
  expect_match(gsub("[[:space:]]+", " ",
                    paste(as.character(cond), collapse = " ")),
               "one pair per respondent")
})

# ---------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------
test_that("a fixed seed reproduces the simulated null", {
  skip_on_cran()
  d <- make_pair(n = 60)
  args <- list(d$t1, d$t2, sim_iter = 100, parallel = FALSE, seed = 99)
  a <- do.call(RMpersonChange, args)
  b <- do.call(RMpersonChange, args)
  expect_identical(a$crit_upper, b$crit_upper)
  expect_identical(a$p_value, b$p_value)
})
