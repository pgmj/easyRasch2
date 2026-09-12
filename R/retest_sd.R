#' Occasion-to-occasion SD from a test-retest study
#'
#' Estimates the per-occasion SD, in logits, of a respondent's
#' occasion-specific deviation: the variation between two administrations that
#' is neither measurement error nor change in the construct. It is the
#' `retest_sd` argument of [RMpersonChange()] when `null = "retest"`.
#'
#' @param data_t1,data_t2 Data.frames or matrices of item responses from a
#'   test-retest study, same items in the same order, row `i` being the same
#'   respondent at both administrations.
#' @param anchor,item_params,method,estimator,theta_range As in
#'   [RMpersonChange()]. Use the same settings there as here, since the
#'   estimate is defined relative to the measurement-error variance the model
#'   implies.
#' @param sim_iter Integer. Iterations used to simulate the measurement-error
#'   component. Default `500`.
#' @param parallel Logical. Use `mirai` for the simulation if available.
#'   Default `TRUE`.
#' @param n_cores Integer or `NULL`. Parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first.
#' @param verbose Logical. Progress bar for the simulation. Default `FALSE`.
#' @param boot Logical. Bootstrap a confidence interval by resampling
#'   respondents. Default `TRUE`.
#' @param boot_iter Integer. Bootstrap iterations. Default `500`.
#' @param conf_int Numeric in (0, 1). HDCI width. Default `0.95`.
#' @param seed Integer or `NULL`. Random seed.
#' @param output Character. `"kable"` (default) or `"dataframe"`.
#'
#' @return A one-row data.frame (or its `knitr_kable`) with `variance`, `sd`,
#'   `lower`, `upper`, `var_change` (the observed variance of the change) and
#'   `var_error` (the simulated measurement-error variance of the change),
#'   plus `n`.
#'
#' @details
#' Writing \eqn{\hat\theta_{it} = \theta_i + u_{it} + e_{it}} with
#' \eqn{u_{it} \sim N(0, \sigma^2_{retest})} the occasion deviation and
#' \eqn{e_{it}} the measurement error, the observed change has variance
#' \eqn{2\sigma^2_{retest} + SE_{i1}^2 + SE_{i2}^2}. So
#'
#' \deqn{\hat\sigma^2_{retest} =
#'   \frac{\mathrm{Var}(d) - \mathrm{Var}_{err}(d)}{2}}
#'
#' **How the error term is obtained.** The obvious choice for
#' \eqn{\mathrm{Var}_{err}(d)} is \eqn{\overline{SE_1^2 + SE_2^2}} from the
#' information-based standard errors, and that is what Caronni et al. (2026)
#' and the classical formulations use. It is badly biased on short scales,
#' because the information standard error is asymptotic and overstates the
#' real sampling variance of \eqn{\hat\theta} when items are few, so the
#' subtraction removes too much. Simulation over a range of scale lengths, at a
#' true occasion variance of 0.09 and n = 2000:
#'
#' \tabular{rrr}{
#'   \strong{items} \tab \strong{asymptotic} \tab \strong{simulated} \cr
#'   6  \tab -0.023 \tab 0.072 \cr
#'   12 \tab  0.053 \tab 0.093 \cr
#'   20 \tab  0.086 \tab 0.099
#' }
#'
#' Across eight replications at n = 1000 the simulated version is close to
#' unbiased on a six-item scale (mean 0.002 against a truth of 0, and 0.081
#' against a truth of 0.09), with a replication SD near 0.015.
#'
#' This function therefore obtains \eqn{\mathrm{Var}_{err}(d)} by simulation:
#' respondents are held at their observed locations, two response vectors are
#' drawn with no occasion variance, and the variance of the resulting change is
#' averaged over `sim_iter` iterations. The observed missingness pattern is
#' carried into the simulated responses. `var_error` in the output is that
#' simulated quantity, not the mean of the squared standard errors.
#'
#' The retest interval has to be long enough for recall to fade and short
#' enough that no real change occurs. That is an assumption about the design
#' and the data cannot check it.
#'
#' **The estimate can be negative.** It is a variance obtained by subtraction,
#' so when the observed change is tighter than the measurement model says it
#' should be, the difference goes below zero. A negative value is reported as
#' such rather than floored, because flooring hides the signal. It points at
#' standard errors overstated by misfit or local dependence, at regression
#' toward a common value over the interval, or at sampling noise. `sd` is `NA`
#' whenever `variance` is negative.
#'
#' It is a variance estimated from `n` differences, so its own sampling error
#' is roughly \eqn{\sigma^2\sqrt{2/(n-1)}}. A few dozen respondents will not
#' pin it down, and the bootstrap interval is the honest way to see that. The
#' bootstrap resamples respondents and recomputes \eqn{\mathrm{Var}(d)} only,
#' holding the simulated error term at its full-sample value, since that term
#' is model-implied and its Monte Carlo error is small beside the sampling
#' error of the change variance.
#'
#' **Converting a published test-retest coefficient.** A published test-retest
#' ICC or Pearson correlation can be converted, but it is not itself the
#' quantity wanted here, and neither is a published test-retest SEM
#' (\eqn{SD\sqrt{1-ICC}}), which already contains measurement error. Since
#' \eqn{\hat\theta = \theta + u + e} gives \eqn{r = \mathrm{Var}(\theta) /
#' \mathrm{Var}(\hat\theta)},
#'
#' \deqn{\hat\sigma^2_{retest} =
#'   (1 - r)\,\mathrm{Var}(\hat\theta) - \mathrm{Var}_{err}(d)/2}
#'
#' with \eqn{\mathrm{Var}(\hat\theta)} the observed variance of the person
#' estimates and \eqn{\mathrm{Var}_{err}(d)} the simulated error term this
#' function reports as `var_error`. Checked against known truth at n = 3000,
#' true occasion variance 0.09:
#'
#' \tabular{rrrrr}{
#'   \strong{items} \tab \strong{r (logit)} \tab \strong{r (sum score)}
#'     \tab \strong{from logit r} \tab \strong{from sum-score r} \cr
#'   6  \tab .808 \tab .829 \tab 0.094 \tab 0.046 \cr
#'   12 \tab .873 \tab .890 \tab 0.084 \tab 0.048 \cr
#'   20 \tab .903 \tab .913 \tab 0.088 \tab 0.067
#' }
#'
#' Three conditions have to hold, and the first is usually violated:
#'
#' * **The coefficient must be on the logit metric.** The score-to-logit map is
#'   nonlinear, so a coefficient computed on sum scores runs consistently
#'   higher than the same coefficient on logits. After subtraction the residual
#'   is small, so that modest gap does real damage: on a six-item scale,
#'   substituting the sum-score value roughly halves the recovered occasion
#'   variance. The error has a known direction, understating occasion noise and
#'   leaving [RMpersonChange()] too permissive. Nearly all published
#'   coefficients are on sum scores.
#' * **The published sample must resemble yours.** A retest coefficient rises
#'   with the heterogeneity of the people measured, so it is a property of that
#'   sample rather than of the instrument. Combining someone else's `r` with
#'   your \eqn{\mathrm{Var}(\hat\theta)} assumes the structure transfers.
#' * **Pearson and an agreement ICC are not interchangeable.** Pearson ignores a
#'   systematic mean shift between occasions and so reports stability that
#'   practice or drift has removed. ICC(A,1) penalises it and is the closer
#'   choice here.
#'
#' The conversion also needs the simulated error term, which this function
#' produces only when given two occasions of raw data. Converting from a
#' published coefficient and a single administration is not implemented. When
#' that is the situation you are in, `retest_sd_tip` from [RMpersonChange()]
#' answers the practical question from the other end: it reports how much
#' occasion noise each flagged result would tolerate, which you can weigh
#' against a published coefficient without converting anything.
#'
#' @seealso [RMpersonChange()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' thr <- lapply(seq(-1.2, 1.2, length.out = 6), function(b) b + c(-0.8, 0, 0.8))
#' theta <- rnorm(200, 0, 1.4)
#' # a true occasion SD of 0.3 logits
#' r1 <- as.data.frame(easyRasch2:::sim_partial_score(thr, theta + rnorm(200, 0, 0.3)))
#' r2 <- as.data.frame(easyRasch2:::sim_partial_score(thr, theta + rnorm(200, 0, 0.3)))
#' colnames(r1) <- colnames(r2) <- paste0("I", 1:6)
#' RMretestSD(r1, r2, sim_iter = 100, parallel = FALSE, boot = FALSE)
#' }
RMretestSD <- function(
  data_t1,
  data_t2,
  anchor = "stack",
  item_params = NULL,
  method = "WLE",
  estimator = "CML",
  sim_iter = 500,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  boot = TRUE,
  boot_iter = 500,
  conf_int = 0.95,
  seed = NULL,
  theta_range = c(-10, 10),
  output = "kable"
) {
  anchor <- match.arg(anchor, c("stack", "t1", "t2"))
  method <- match.arg(method, c("WLE", "EAP"))
  estimator <- match.arg(estimator, c("CML", "MML"))
  output <- match.arg(output, c("kable", "dataframe"))

  if (
    !is.numeric(conf_int) || length(conf_int) != 1L ||
      conf_int <= 0 || conf_int >= 1
  ) {
    stop("`conf_int` must be a numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  validate_response_data(data_t1)
  validate_response_data(data_t2)
  data_t1 <- as.data.frame(data_t1)
  data_t2 <- as.data.frame(data_t2)
  .pc_check_paired(data_t1, data_t2)

  if (!is.null(item_params)) {
    thr_list <- .coerce_item_params(item_params)
    data_t1 <- .align_items(data_t1, thr_list)
    data_t2 <- .align_items(data_t2, thr_list)
  } else {
    cal_data <- switch(anchor,
      stack = rbind(data_t1, data_t2),
      t1 = data_t1,
      t2 = data_t2
    )
    thr_list <- if (estimator == "CML") {
      .center_thresholds(.rasch_fit_cml(cal_data, se = FALSE)$thr_list)
    } else {
      .center_thresholds(.rasch_fit_mml(cal_data, se = FALSE)$thr_list)
    }
  }

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max - 2L, 1L)
  }

  point <- .retest_sd_point(
    as.matrix(data_t1), as.matrix(data_t2), thr_list, method, theta_range,
    sim_iter = sim_iter, parallel = parallel, n_cores = n_cores,
    seed = seed, verbose = verbose
  )

  lower <- NA_real_
  upper <- NA_real_
  if (isTRUE(boot)) {
    if (!requireNamespace("ggdist", quietly = TRUE)) {
      stop(
        "Package 'ggdist' is required for boot = TRUE.\n",
        "Install it with: install.packages(\"ggdist\")",
        call. = FALSE
      )
    }
    set.seed(seed + 1L)
    # Only Var(d) is resampled; the simulated error term is held fixed (see
    # Details), which is what makes a bootstrap affordable here at all.
    d_obs <- point$change[point$keep]
    n <- length(d_obs)
    vals <- vapply(seq_len(boot_iter), function(i) {
      idx <- sample.int(n, n, replace = TRUE)
      (stats::var(d_obs[idx]) - point$var_error) / 2
    }, numeric(1L))
    vals <- vals[is.finite(vals)]
    if (length(vals) >= 2L) {
      int <- ggdist::hdci(vals, .width = conf_int)
      # The interval is on the variance, so it is carried to the SD only where
      # the bound is non-negative.
      lower <- if (int[1L, 1L] >= 0) sqrt(int[1L, 1L]) else NA_real_
      upper <- if (int[1L, 2L] >= 0) sqrt(int[1L, 2L]) else NA_real_
    }
  }

  res <- data.frame(
    n = nrow(data_t1),
    var_change = point$var_change,
    var_error = point$var_error,
    variance = point$variance,
    sd = point$sd,
    lower = lower,
    upper = upper,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (output == "dataframe") {
    return(res)
  }

  note <- if (is.na(point$sd)) {
    paste(
      "The estimated variance is negative, so no occasion SD is reported.",
      "The observed change is tighter than the measurement model implies,",
      "which points at standard errors overstated by misfit or local",
      "dependence, at regression toward a common value over the interval, or",
      "at sampling noise. It is reported unfloored."
    )
  } else {
    paste(
      "Per-occasion SD of occasion-specific deviation, in logits, for use as",
      "`retest_sd` in RMpersonChange(null = \"retest\")."
    )
  }

  knitr::kable(
    .round_display(
      res,
      c(var_change = 4, var_error = 4, variance = 4, sd = 3, lower = 3,
        upper = 3)
    ),
    format = "pipe",
    col.names = c(
      "n", "Var(change)", "Mean error var", "Occasion var", "Occasion SD",
      paste0("Lower (", round(conf_int * 100, 1), "% HDCI)"),
      paste0("Upper (", round(conf_int * 100, 1), "% HDCI)")
    ),
    caption = paste(
      note,
      "Estimated as (Var(change) - mean error variance) / 2.",
      .n_caption(nrow(data_t1), nrow(data_t1)),
      "measured twice."
    )
  )
}

#' Point estimate of the occasion variance
#'
#' @keywords internal
#' @noRd
.retest_sd_point <- function(
  mat_t1, mat_t2, thr_list, method, theta_range,
  sim_iter, parallel, n_cores, seed, verbose
) {
  e1 <- .estimate_thetas(mat_t1, thr_list, method = method,
                         theta_range = theta_range)
  e2 <- .estimate_thetas(mat_t2, thr_list, method = method,
                         theta_range = theta_range)

  d <- e2$theta - e1$theta
  keep <- is.finite(d)
  if (sum(keep) < 2L) {
    stop(
      "Fewer than two respondents have finite estimates at both occasions.",
      call. = FALSE
    )
  }

  var_change <- stats::var(d[keep])
  var_error <- .retest_error_var(
    theta_plug = ((e1$theta + e2$theta) / 2)[keep],
    thr_list = thr_list,
    method = method,
    theta_range = theta_range,
    na1 = is.na(mat_t1)[keep, , drop = FALSE],
    na2 = is.na(mat_t2)[keep, , drop = FALSE],
    sim_iter = sim_iter,
    parallel = parallel,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose
  )
  variance <- (var_change - var_error) / 2

  list(
    change = d,
    keep = keep,
    var_change = var_change,
    var_error = var_error,
    variance = variance,
    sd = if (variance >= 0) sqrt(variance) else NA_real_
  )
}

#' Simulated variance of the change under no occasion variance
#'
#' Reuses the `RMpersonChange()` null simulation, which returns RCIs; the
#' change itself is recovered by multiplying back through `se_diff`, so both
#' functions draw from one code path.
#'
#' @keywords internal
#' @noRd
.retest_error_var <- function(
  theta_plug, thr_list, method, theta_range, na1, na2,
  sim_iter, parallel, n_cores, seed, verbose, prior_sd = NULL
) {
  sim <- .pc_simulate(
    theta_null = theta_plug,
    thr_list = thr_list,
    sigma_retest = 0,
    method = method,
    theta_range = theta_range,
    na1 = na1,
    na2 = na2,
    sim_iter = sim_iter,
    parallel = parallel,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose,
    prior_sd = prior_sd,
    return_change = TRUE
  )
  mean(
    apply(sim, 2L, function(x) stats::var(x[is.finite(x)])),
    na.rm = TRUE
  )
}
