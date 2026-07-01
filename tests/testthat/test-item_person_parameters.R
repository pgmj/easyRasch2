# Tests for RMitemParameters() and RMpersonParameters()

make_polytomous <- function(n = 250, k = 5, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:2, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_dichotomous <- function(n = 250, k = 6, seed = 2L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# Simulate data with genuine latent structure (for recovery-style checks).
make_pcm_signal <- function(n = 600, k = 6, person_sd = 1.4, seed = 7L) {
  set.seed(seed)
  theta <- rnorm(n, 0, person_sd)
  thr   <- lapply(seq_len(k), function(j) sort(rnorm(2, 0, 1)))
  pcm_p <- function(th, t) {
    ce <- c(0, cumsum(th - t)); exp(ce - max(ce)) / sum(exp(ce - max(ce)))
  }
  m <- sapply(seq_len(k), function(j) {
    sapply(theta, function(th) sample(0:2, 1, prob = pcm_p(th, thr[[j]])))
  })
  df <- as.data.frame(m); colnames(df) <- paste0("I", seq_len(k))
  list(data = df, theta = theta, person_sd = person_sd)
}

# =====================================================================
# RMitemParameters
# =====================================================================

test_that("RMitemParameters errors on non-zero minimum", {
  expect_error(RMitemParameters(make_dichotomous() + 1L),
               regexp = "scored starting at 0")
})

test_that("RMitemParameters errors on invalid ci_level", {
  skip_if_not_installed("eRm")
  expect_error(RMitemParameters(make_dichotomous(), ci_level = 1.2),
               regexp = "ci_level")
})

test_that("long format has expected columns and rows (polytomous)", {
  skip_if_not_installed("eRm")
  res <- RMitemParameters(make_polytomous(), output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_named(res, c("item", "threshold", "location",
                      "se", "ci_lower", "ci_upper"))
  # 5 items x 2 thresholds
  expect_equal(nrow(res), 10L)
  expect_true(all(res$ci_lower < res$location & res$location < res$ci_upper))
})

test_that("se = FALSE drops SE/CI columns", {
  skip_if_not_installed("eRm")
  res <- RMitemParameters(make_polytomous(), se = FALSE, output = "dataframe")
  expect_named(res, c("item", "threshold", "location"))
})

test_that("wide format returns one row per item with threshold columns", {
  skip_if_not_installed("eRm")
  res <- RMitemParameters(make_polytomous(), format = "wide",
                          se = FALSE, output = "dataframe")
  expect_equal(nrow(res), 5L)
  expect_true(all(c("item", "t1", "t2", "location") %in% names(res)))
})

test_that("wide format with SE appends se_t* columns", {
  skip_if_not_installed("eRm")
  res <- RMitemParameters(make_polytomous(), format = "wide",
                          se = TRUE, output = "dataframe")
  expect_true(all(c("se_t1", "se_t2") %in% names(res)))
})

test_that("dichotomous data yields one row per item", {
  skip_if_not_installed("eRm")
  res <- RMitemParameters(make_dichotomous(), output = "dataframe")
  expect_equal(nrow(res), 6L)
  expect_true(all(res$threshold == 1L))
})

test_that("kable output returns knitr_kable", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("knitr")
  expect_s3_class(RMitemParameters(make_polytomous()), "knitr_kable")
})

test_that("output = 'file' writes the dataframe to CSV and returns it invisibly", {
  skip_if_not_installed("eRm")
  path <- withr::local_tempfile(fileext = ".csv")
  df  <- RMitemParameters(make_polytomous(), output = "dataframe")
  ret <- suppressMessages(
    RMitemParameters(make_polytomous(), output = "file", filename = path)
  )
  expect_equal(ret, df)
  back <- read.csv(path)
  expect_equal(names(back), names(df))
  expect_equal(nrow(back), nrow(df))
  expect_equal(back$location, df$location, tolerance = 1e-6)
})

test_that("output = 'file' errors early when filename is missing", {
  skip_if_not_installed("eRm")
  expect_snapshot(
    RMitemParameters(make_polytomous(), output = "file"),
    error = TRUE
  )
})

test_that("wider ci_level gives wider intervals", {
  skip_if_not_installed("eRm")
  narrow <- RMitemParameters(make_polytomous(), ci_level = 0.90,
                             output = "dataframe")
  wide   <- RMitemParameters(make_polytomous(), ci_level = 0.99,
                             output = "dataframe")
  expect_true(all((wide$ci_upper - wide$ci_lower) >
                    (narrow$ci_upper - narrow$ci_lower)))
})

test_that("MML estimator agrees with CML on item locations", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("mirt")
  sim <- make_pcm_signal()
  cml <- RMitemParameters(sim$data, estimator = "CML", output = "dataframe")
  mml <- RMitemParameters(sim$data, estimator = "MML", output = "dataframe")
  ord <- match(paste(cml$item, cml$threshold),
               paste(mml$item, mml$threshold))
  expect_gt(cor(cml$location, mml$location[ord]), 0.98)
  expect_true(all(is.finite(mml$se)))
})

test_that("CML threshold SEs are recovered (delta method sane)", {
  skip_if_not_installed("eRm")
  res <- RMitemParameters(make_pcm_signal()$data, output = "dataframe")
  expect_true(all(res$se > 0 & res$se < 1))
})

# =====================================================================
# RMpersonParameters
# =====================================================================

test_that("RMpersonParameters errors on bad theta_range", {
  expect_error(RMpersonParameters(make_polytomous(), theta_range = c(5, -5)),
               regexp = "theta_range")
})

test_that("WLE dataframe has one row per person and expected columns", {
  skip_if_not_installed("eRm")
  df  <- make_polytomous(n = 200)
  res <- RMpersonParameters(df, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 200L)
  expect_named(res, c("theta", "sem", "sum_score", "n_answered", "extreme"))
})

test_that("WLE gives finite extrapolated estimates for extreme scores", {
  skip_if_not_installed("eRm")
  df <- make_dichotomous(n = 200, k = 6)
  df[1, ] <- 0L            # minimum score
  df[2, ] <- 1L            # maximum score
  res <- RMpersonParameters(df, output = "dataframe")
  expect_true(res$extreme[1] && res$extreme[2])
  # WLE is finite at the boundaries (Warm), well inside theta_range,
  # and a sensible step beyond the next-most-extreme estimate.
  expect_true(is.finite(res$theta[1]) && is.finite(res$theta[2]))
  expect_true(res$theta[1] > -10 && res$theta[2] < 10)
  second_low  <- min(res$theta[res$sum_score == 1])
  second_high <- max(res$theta[res$sum_score == 5])
  expect_lt(res$theta[1], second_low)    # zero score below the 1-score
  expect_gt(res$theta[2], second_high)   # perfect score above the 5-score
})

test_that("WLE satisfies Warm's weighted-likelihood score equation", {
  beta <- c(-1.5, -0.8, -0.2, 0.3, 0.9, 1.6)
  thr  <- as.list(beta)
  # Warm score: l'(theta) + J(theta) / (2 I(theta)); zero at the WLE.
  warm_score <- function(theta, x) {
    p   <- plogis(theta - beta)
    dll <- sum(x) - sum(p)
    I   <- sum(p * (1 - p))
    J   <- sum(p * (1 - p) * (1 - 2 * p))
    dll + J / (2 * I)
  }
  for (r in 0:length(beta)) {
    x  <- c(rep(1L, r), rep(0L, length(beta) - r))
    th <- easyRasch2:::.theta_wle(x, thr, c(-10, 10))[["theta"]]
    expect_lt(abs(warm_score(th, x)), 1e-4)   # solves the equation
  }
})

test_that("partial missingness is handled per-pattern", {
  skip_if_not_installed("eRm")
  df <- make_polytomous(n = 200)
  df[cbind(sample(200, 40), sample(5, 40, replace = TRUE))] <- NA
  res <- RMpersonParameters(df, output = "dataframe")
  expect_equal(nrow(res), 200L)
  expect_true(any(res$n_answered < 5L))
  # non-extreme persons get a finite estimate even with missingness
  expect_true(all(is.finite(res$theta[!res$extreme])))
})

test_that("EAP returns finite estimates at extreme scores", {
  skip_if_not_installed("eRm")
  df <- make_dichotomous(n = 200, k = 6)
  df[1, ] <- 0L; df[2, ] <- 1L
  res <- RMpersonParameters(df, method = "EAP", output = "dataframe")
  expect_true(all(is.finite(res$theta)))
  expect_true(all(is.finite(res$sem)))
})

test_that("EAP attaches the prior used; WLE does not", {
  skip_if_not_installed("eRm")
  df  <- make_pcm_signal()$data
  eap <- RMpersonParameters(df, method = "EAP", output = "dataframe")
  wle <- RMpersonParameters(df, method = "WLE", output = "dataframe")
  pr  <- attr(eap, "prior")
  expect_false(is.null(pr))
  expect_identical(unname(pr["estimated"]), 1)  # estimated
  expect_null(attr(wle, "prior"))
})

test_that("fixed prior_sd is respected and reported as not estimated", {
  skip_if_not_installed("eRm")
  df  <- make_pcm_signal()$data
  eap <- RMpersonParameters(df, method = "EAP", prior_sd = 1,
                            output = "dataframe")
  pr  <- attr(eap, "prior")
  expect_equal(unname(pr["sd"]), 1)
  expect_equal(unname(pr["estimated"]), 0)
})

test_that("marginal-ML prior SD recovers the simulated person SD", {
  skip_if_not_installed("eRm")
  sim <- make_pcm_signal(n = 800, k = 8, person_sd = 1.4)
  eap <- RMpersonParameters(sim$data, method = "EAP", output = "dataframe")
  est_sd <- unname(attr(eap, "prior")["sd"])
  expect_equal(est_sd, sim$person_sd, tolerance = 0.15)
})

test_that("invalid prior_sd is rejected", {
  skip_if_not_installed("eRm")
  expect_error(
    RMpersonParameters(make_polytomous(), method = "EAP", prior_sd = -1),
    regexp = "prior_sd"
  )
})

test_that("supplied item_params (list) bypasses model fitting", {
  skip_if_not_installed("eRm")
  df  <- make_dichotomous(n = 150, k = 4)
  ip  <- list(I1 = 0.2, I2 = -0.3, I3 = 0.1, I4 = 0.0)
  res <- RMpersonParameters(df, item_params = ip, output = "dataframe")
  expect_equal(nrow(res), 150L)
  expect_true(all(is.finite(res$theta[!res$extreme])))
})

test_that("supplied item_params (data.frame from RMitemParameters) works", {
  skip_if_not_installed("eRm")
  df <- make_polytomous(n = 150)
  ip <- RMitemParameters(df, output = "dataframe")[, c("item", "threshold",
                                                       "location")]
  res <- RMpersonParameters(df, item_params = ip, output = "dataframe")
  expect_equal(nrow(res), 150L)
})

test_that("kable and ggplot outputs return the right classes", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("knitr")
  df <- make_polytomous(n = 150)
  expect_s3_class(RMpersonParameters(df), "knitr_kable")
  skip_if_not_installed("ggplot2")
  expect_s3_class(RMpersonParameters(df, output = "ggplot"), "ggplot")
})

test_that("RMpersonParameters output = 'file' writes a CSV with one row per person", {
  skip_if_not_installed("eRm")
  df   <- make_polytomous(n = 150)
  path <- withr::local_tempfile(fileext = ".csv")
  ddf  <- RMpersonParameters(df, output = "dataframe")
  ret  <- suppressMessages(
    RMpersonParameters(df, output = "file", filename = path)
  )
  expect_equal(ret, ddf)
  back <- read.csv(path)
  expect_equal(nrow(back), nrow(df))
  expect_equal(back$theta, ddf$theta, tolerance = 1e-6)
})
