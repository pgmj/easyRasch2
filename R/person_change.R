#' Single-subject change between two occasions
#'
#' Tests, for each respondent, whether their person location moved between two
#' occasions by more than measurement error allows. The statistic is the Rasch
#' change index, \eqn{RCI = (\hat\theta_2 - \hat\theta_1)/SE_{diff}}, referred
#' by default to a simulated rather than a normal null.
#'
#' @param data_t1,data_t2 Data.frames or matrices of item responses at the two
#'   occasions. Same items, in the same order, one row per respondent, with row
#'   `i` of each being the same person. Items must be scored from 0.
#' @param id Optional vector of respondent identifiers, length `nrow(data_t1)`.
#'   Defaults to row numbers.
#' @param anchor Character. Where the item calibration comes from: `"stack"`
#'   (default, both occasions calibrated together), `"t1"`, or `"t2"`. Ignored
#'   when `item_params` is supplied.
#' @param item_params Optional pre-specified item parameters, either a named
#'   list of Andrich-threshold vectors or the long-format data.frame from
#'   [RMitemParameters()], as in [RMpersonParameters()]. **Required for a
#'   single respondent**, since item parameters cannot be estimated from one
#'   person. Supplied thresholds are used as given and are not re-centred.
#' @param method Character. Person-location estimator, `"WLE"` (default) or
#'   `"EAP"`.
#' @param estimator Character. How item parameters are estimated when
#'   `item_params` is `NULL`: `"CML"` (default) or `"MML"`.
#' @param null Character. Which null hypothesis is tested. `"measurement"`
#'   (default) treats the response process as the only source of variation.
#'   `"retest"` additionally treats occasion-to-occasion fluctuation as noise
#'   and requires `retest_sd`. See Details, and note that the two answer
#'   different questions.
#' @param retest_sd Numeric, or `NULL` (default). The **per-occasion** SD, in
#'   logits, of a respondent's occasion-specific deviation. Required when
#'   `null = "retest"` and not accepted otherwise. Estimate it with
#'   [RMretestSD()].
#' @param critical How the critical values are obtained. `"exact"` (default)
#'   enumerates the null distribution, `"simulate"` estimates it by Monte
#'   Carlo, and a positive number is used as a symmetric cutoff on the RCI
#'   referred to a normal null (for example `1.96`).
#' @param alpha Numeric in (0, 1). Error rate for a single respondent. Default
#'   `0.05`.
#' @param direction Character. `"two.sided"` (default), or `"increase"` /
#'   `"decrease"` for a one-sided test. These name the direction of
#'   \eqn{\theta}; whether an increase is an improvement depends on how the
#'   scale is oriented.
#' @param conditional_crit Logical. If `TRUE`, critical values are simulated
#'   separately for each respondent at their own location rather than pooled
#'   across the sample. Default `FALSE`. Ignored when `critical` is numeric.
#' @param sim_iter Integer. Simulation iterations when
#'   `critical = "simulate"`. Default `1000`.
#' @param parallel Logical. Use `mirai` for the simulation if available.
#'   Default `TRUE`.
#' @param n_cores Integer or `NULL`. Parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first; if neither is set, the
#'   simulation runs sequentially.
#' @param seed Integer or `NULL`. Random seed. See
#'   [easyRasch2-reproducibility].
#' @param theta_range Numeric length 2. Search range for the WLE root and
#'   bounds for the EAP grid. Default `c(-10, 10)`.
#' @param verbose Logical. Print a progress bar for the simulation. Default
#'   `FALSE`.
#' @param output Character. `"dataframe"` (default), `"kable"`, or `"ggplot"`.
#'
#' @return
#' * If `output = "dataframe"`: one row per respondent, with columns `id`,
#'   `sum_t1`, `sum_t2`, `theta_t1`, `se_t1`, `theta_t2`, `se_t2`,
#'   `extreme_t1`, `extreme_t2`, `change`, `se_diff`, `rci`, `p_value`,
#'   `crit_lower`, `crit_upper`, `change_class`, and `retest_sd_tip`.
#'   Attributes `null`, `retest_sd`, `anchor`, `alpha`, `direction`,
#'   `critical` and `sim_iter` record the analysis.
#' * If `output = "kable"`: the same content as a `knitr_kable`.
#' * If `output = "ggplot"`: occasion 1 against occasion 2, with the
#'   no-change band and points coloured by `change_class`.
#'
#' @details
#' **Which null is tested.** Following Zumbo (2026), the estimand is fixed by
#' what counts as systematic and what counts as residual, and that choice has
#' to be made before an estimator is picked. Two are available here:
#'
#' * `null = "measurement"`, the default. The respondent's \eqn{\theta} is
#'   fixed and identical at both occasions and only the response process is
#'   random, so \eqn{SE_{diff}^2 = SE_1^2 + SE_2^2}. This asks whether the two
#'   estimates differ by more than responding alone would produce.
#' * `null = "retest"`. Occasion-to-occasion fluctuation that is not change in
#'   the construct (state, mood, recall, practice) also counts as noise, so
#'   \eqn{SE_{diff}^2 = SE_1^2 + SE_2^2 + 2\sigma^2_{retest}}. The factor of 2
#'   is there because `retest_sd` is a per-occasion SD and both occasions
#'   carry one.
#'
#' The default is the narrower of the two, so it flags more change than a
#' retest-based null would. It is the only one computable from two response
#' vectors alone, which is why it is the default, but a significant result
#' under it is a necessary condition for change rather than evidence of it.
#' `retest_sd_tip` reports, for each flagged respondent, how large
#' \eqn{\sigma_{retest}} would have to be to overturn their result.
#'
#' Caronni et al. (2026) compute \eqn{1.96\sqrt{SE_A^2 + SE_B^2}} from the
#' Rasch score-to-measure map and call it a minimal detectable change. That is
#' the `null = "measurement"` quantity, whereas the classical MDC is built from
#' a test-retest ICC and is the `null = "retest"` quantity. The two are not
#' interchangeable, which is why the null in force is named in every caption
#' here.
#'
#' Regression to the mean is present under both nulls and is addressed by
#' neither. A regression-based reliable change index (Maassen, 2004) is out of
#' scope for a single-subject function, since it needs population parameters
#' the individual case does not supply.
#'
#' **Why 1.96 is the wrong cutoff.** The sum score is a sufficient statistic,
#' so \eqn{\hat\theta} takes only as many values as there are scores, and
#' \eqn{SE(\hat\theta)} is a deterministic function of the score rather than a
#' constant. The denominator of the RCI is therefore not independent of its
#' numerator: a respondent who moves to an extreme score produces a large
#' change and a large standard error together, and the ratio is damped. Both
#' tails of the null are pulled in, so the critical value sits below 1.96, and
#' the shortfall grows as the scale shortens. Enumerated values for one family
#' of items with four response categories: 1.62 at four items, 1.70 at six,
#' 1.79 at ten, 1.85 at twenty, 1.92 at forty. Response categories matter much
#' less than item count here, and almost all of their effect is the step from
#' dichotomous to polytomous: at six items the value moves from 1.55 with two
#' categories to 1.72 with three and only 1.78 with seven, while reliability
#' over the same range climbs from about .60 to about .91. Referring the RCI to
#' a normal null is conservative, and materially so on short scales.
#'
#' **How the null is obtained.** `critical = "exact"` enumerates it. Because
#' the score is sufficient, \eqn{\hat\theta(r)} and \eqn{SE(r)} are
#' deterministic lookups and \eqn{P(r \mid \theta)} follows from the
#' Lord-Wingersky recursion, so the null is a discrete distribution over score
#' pairs with no Monte Carlo error. Respondents are grouped by their pair of
#' answered-item sets, so incomplete data is enumerated within each pattern.
#' Under `null = "retest"` the occasion deviations are integrated out by
#' quadrature. `critical = "simulate"` estimates the same distribution by
#' drawing response vectors instead, which is slower and noisier but makes no
#' use of sufficiency.
#'
#' Quantiles are taken from the **signed** RCI, so `crit_lower` and
#' `crit_upper` are reported separately. Under the null the two occasions are
#' exchangeable, which makes the RCI symmetric about zero however skewed
#' \eqn{\hat\theta} itself is, and the two bounds then agree in magnitude.
#' They part company when exchangeability fails, most commonly when the
#' occasions differ in which items were answered.
#'
#' The pooled critical value is a property of the sample as well as the items,
#' because it averages over where the respondents sit. Targeting matters more
#' than the shape of the distribution: skewing the sample barely moves it,
#' while pushing the same sample off target until 40 percent of respondents
#' score at an extreme moved it from 1.70 to 1.52 in one check.
#'
#' **What the RCI measures, and what it does not.** The RCI answers whether a
#' change is distinguishable from the null. It is not a measure of how much
#' someone changed, and the two orderings genuinely differ. Its value mixes two
#' ingredients that it cannot separate: how far the respondent moved, and how
#' precisely each of their two positions was pinned down. On a scale of six
#' items with four response categories the standard error runs from 0.49 in the
#' middle to 1.45 at the boundary, a factor of three within one instrument. A respondent going from
#' the minimum to the maximum score moves 7.24 logits and scores RCI 3.53,
#' while one going from 3 to 15 moves 3.11 logits and scores 3.59. The first
#' traversed the whole scale, the second less than half of it, and the second
#' has the larger statistic. Both statements are correct about their own
#' question.
#'
#' In a group comparison every case shares one standard error, so ordering by
#' the test statistic matches ordering by effect. Here the scale factor changes
#' from respondent to respondent, and the output is read one person at a time,
#' which is where the misreading does damage. Read `change` for magnitude, in
#' logits, and `rci` only for the decision. Do not rank respondents by `rci`,
#' and do not compare RCIs across instruments or across regions of one scale as
#' though they were a common currency.
#'
#' Simulating fixes the reference distribution, not the estimand. A simulated
#' critical value under `null = "measurement"` is a well calibrated answer to
#' the measurement-error question, not to the retest question.
#'
#' **Calibration and its precondition.** `anchor = "stack"` calibrates on both
#' occasions row-bound together. It uses all the data and puts both occasions
#' on one metric by construction, at the cost of assuming the items behave the
#' same way on both occasions. When they do not, the compromise calibration
#' pulls the occasions toward each other and shrinks apparent change. Test that
#' precondition first with the package's DIF tools, using occasion as the
#' grouping variable on the stacked data. `anchor = "t1"` measures follow-up on
#' the baseline metric and keeps post-treatment data out of it, at the cost of
#' halving the calibration sample. `"t2"` is the mirror.
#'
#' Note that stacking puts each respondent in the data twice, so the item
#' parameters' own standard errors are optimistic. That does not propagate into
#' the person standard errors, which treat the thresholds as fixed.
#'
#' **`retest_sd_tip`.** For a respondent whose change is flagged, this is the
#' per-occasion retest SD that would bring their RCI back to the critical
#' value, \eqn{\sigma^2_{tip} = (\mathrm{change}^2/c^2 - SE_1^2 - SE_2^2)/2}.
#' It is `NA` for respondents not flagged, for whom the question does not
#' arise. The critical value is held fixed at the one in force, which is exact
#' when `critical` is numeric and an approximation when it is simulated, since
#' a simulated critical value drifts toward \eqn{\pm 1.96} as the added normal
#' component smooths the discreteness.
#'
#' @references
#' Caronni, A., et al. (2026). Improving single-subject change assessment:
#' deriving the minimal detectable change of questionnaires' ordinal scores
#' from the Rasch analysis measures.
#'
#' Jacobson, N. S., & Truax, P. (1991). Clinical significance: A statistical
#' approach to defining meaningful change in psychotherapy research.
#' *Journal of Consulting and Clinical Psychology, 59*(1), 12-19.
#' \doi{10.1037/0022-006X.59.1.12}
#'
#' Maassen, G. H. (2004). The standard error in the Jacobson and Truax
#' Reliable Change Index. *Journal of Clinical and Experimental
#' Neuropsychology, 26*(5), 643-657. \doi{10.1080/13803390409609791}
#'
#' Zumbo, B. D. (2026). Conditional standard error of measurement as an
#' estimand of individual score precision. *Psychometrika*.
#' \doi{10.1017/psy.2026.10141}
#'
#' @seealso [RMretestSD()], [RMpersonParameters()], [RMreliabilityCurve()],
#'   [RMscoreSE()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' thr <- lapply(seq(-1.2, 1.2, length.out = 6), function(b) b + c(-0.8, 0, 0.8))
#' theta <- rnorm(80, 0, 1.4)
#' t1 <- as.data.frame(easyRasch2:::sim_partial_score(thr, theta))
#' t2 <- as.data.frame(easyRasch2:::sim_partial_score(thr, theta + 0.8))
#' colnames(t1) <- colnames(t2) <- paste0("I", 1:6)
#'
#' RMpersonChange(t1, t2)
#'
#' # A single respondent needs an external calibration
#' names(thr) <- paste0("I", 1:6)
#' RMpersonChange(t1[1, ], t2[1, ], item_params = thr)
#'
#' # Monte Carlo alternative, for comparison
#' RMpersonChange(t1, t2, critical = "simulate", sim_iter = 200,
#'                parallel = FALSE, seed = 1)
#' }
RMpersonChange <- function(
  data_t1,
  data_t2,
  id = NULL,
  anchor = "stack",
  item_params = NULL,
  method = "WLE",
  estimator = "CML",
  null = "measurement",
  retest_sd = NULL,
  critical = "exact",
  alpha = 0.05,
  direction = "two.sided",
  conditional_crit = FALSE,
  sim_iter = 1000,
  parallel = TRUE,
  n_cores = NULL,
  seed = NULL,
  theta_range = c(-10, 10),
  verbose = FALSE,
  output = "dataframe"
) {
  anchor_set <- !missing(anchor)

  .pc_check_swapped_args(method, estimator, null, direction)

  anchor <- match.arg(anchor, c("stack", "t1", "t2"))
  method <- match.arg(method, c("WLE", "EAP"))
  estimator <- match.arg(estimator, c("CML", "MML"))
  null <- match.arg(null, c("measurement", "retest"))
  direction <- match.arg(direction, c("two.sided", "increase", "decrease"))
  output <- match.arg(output, c("dataframe", "kable", "ggplot"))

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  # --- The estimand, and the parameter it needs ------------------------------
  if (null == "retest") {
    if (is.null(retest_sd)) {
      stop(
        "`null = \"retest\"` needs `retest_sd`, the per-occasion SD in logits ",
        "of occasion-to-occasion fluctuation. It cannot be estimated from ",
        "these two occasions, since that is the change being tested. ",
        "Estimate it from a separate test-retest study with RMretestSD(), or ",
        "use null = \"measurement\" and read `retest_sd_tip` to see how much ",
        "occasion noise each result would tolerate.",
        call. = FALSE
      )
    }
    if (
      !is.numeric(retest_sd) || length(retest_sd) != 1L ||
        !is.finite(retest_sd) || retest_sd < 0
    ) {
      stop("`retest_sd` must be a single non-negative number.", call. = FALSE)
    }
  } else if (!is.null(retest_sd)) {
    stop(
      "`retest_sd` only applies when `null = \"retest\"`. Set null = ",
      "\"retest\" to test against occasion-to-occasion fluctuation as well ",
      "as measurement error.",
      call. = FALSE
    )
  }
  sigma_retest <- if (null == "retest") retest_sd else 0

  # --- Critical-value machinery ----------------------------------------------
  simulate_null <- identical(critical, "simulate")
  exact_null <- identical(critical, "exact")
  if (!simulate_null && !exact_null) {
    if (
      !is.numeric(critical) || length(critical) != 1L ||
        !is.finite(critical) || critical <= 0
    ) {
      stop(
        "`critical` must be \"exact\", \"simulate\" or a single positive ",
        "number.",
        call. = FALSE
      )
    }
  }
  if (simulate_null) {
    if (
      !is.numeric(sim_iter) || length(sim_iter) != 1L ||
        !is.finite(sim_iter) || sim_iter < 2
    ) {
      stop("`sim_iter` must be a single number of at least 2.", call. = FALSE)
    }
    sim_iter <- as.integer(sim_iter)
  }

  # --- Data ------------------------------------------------------------------
  validate_response_data(data_t1)
  validate_response_data(data_t2)
  data_t1 <- as.data.frame(data_t1)
  data_t2 <- as.data.frame(data_t2)
  .pc_check_paired(data_t1, data_t2)

  n_persons <- nrow(data_t1)
  if (is.null(id)) {
    id <- seq_len(n_persons)
  } else if (length(id) != n_persons) {
    stop(
      "`id` must have one entry per respondent (",
      n_persons,
      " expected, ",
      length(id),
      " supplied).",
      call. = FALSE
    )
  }

  # --- Item calibration ------------------------------------------------------
  if (!is.null(item_params)) {
    if (anchor_set) {
      message(
        "`item_params` supplied, so `anchor` is ignored and the supplied ",
        "item parameters are used."
      )
    }
    thr_list <- .coerce_item_params(item_params)
    data_t1 <- .align_items(data_t1, thr_list)
    data_t2 <- .align_items(data_t2, thr_list)
    anchor <- "supplied"
  } else {
    if (n_persons < 10L) {
      stop(
        "Item parameters cannot be estimated from ",
        n_persons,
        " respondent(s). Supply `item_params` from a published or previously ",
        "established calibration. This is the normal route for single-subject ",
        "use.",
        call. = FALSE
      )
    }
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

  mat_t1 <- as.matrix(data_t1)
  mat_t2 <- as.matrix(data_t2)
  .pc_check_categories(thr_list, mat_t1, mat_t2)

  # --- Person locations ------------------------------------------------------
  est1 <- .estimate_thetas(mat_t1, thr_list, method = method,
                           theta_range = theta_range)
  est2 <- .estimate_thetas(mat_t2, thr_list, method = method,
                           theta_range = theta_range)
  # The exact lookup must score its representative patterns under the same
  # prior the respondents were scored under, or the two would not line up.
  prior_sd <- if (method == "EAP") attr(est1, "prior_sd") else NULL

  steps <- vapply(thr_list, length, integer(1L))
  ext1 <- .pc_extreme(mat_t1, steps)
  ext2 <- .pc_extreme(mat_t2, steps)

  change <- est2$theta - est1$theta
  se_diff <- sqrt(est1$sem^2 + est2$sem^2 + 2 * sigma_retest^2)
  rci <- change / se_diff

  # --- Critical values and p-values ------------------------------------------
  probs <- .pc_tail_probs(direction, alpha)

  if (exact_null) {
    theta_null <- .pc_theta_null(est1, est2)
    ex <- .pc_exact_null(
      thr_list = thr_list,
      theta_null = theta_null,
      mat_t1 = mat_t1,
      mat_t2 = mat_t2,
      rci_obs = rci,
      sigma_retest = sigma_retest,
      method = method,
      theta_range = theta_range,
      prior_sd = prior_sd,
      probs = probs,
      direction = direction,
      conditional_crit = conditional_crit
    )
    crit <- ex$crit
    p_value <- ex$p_value
    actual_iter <- NA_integer_
  } else if (simulate_null) {
    theta_null <- .pc_theta_null(est1, est2)
    sim_mat <- .pc_simulate(
      theta_null = theta_null,
      thr_list = thr_list,
      sigma_retest = sigma_retest,
      method = method,
      theta_range = theta_range,
      na1 = is.na(mat_t1),
      na2 = is.na(mat_t2),
      prior_sd = prior_sd,
      sim_iter = sim_iter,
      parallel = parallel,
      n_cores = n_cores,
      seed = seed,
      verbose = verbose
    )
    crit <- .pc_critical(sim_mat, probs, conditional_crit, n_persons)
    p_value <- .pc_pvalues(rci, sim_mat, direction, conditional_crit)
    actual_iter <- ncol(sim_mat)
  } else {
    crit <- list(
      lower = rep(if (is.finite(probs$lo)) -critical else -Inf, n_persons),
      upper = rep(if (is.finite(probs$hi)) critical else Inf, n_persons)
    )
    p_value <- .pc_normal_p(rci, direction)
    actual_iter <- NA_integer_
  }

  change_class <- .pc_classify(rci, crit)
  tip <- .pc_tipping(change, est1$sem, est2$sem, crit, change_class)

  result <- data.frame(
    id = id,
    sum_t1 = rowSums(mat_t1, na.rm = TRUE),
    sum_t2 = rowSums(mat_t2, na.rm = TRUE),
    theta_t1 = est1$theta,
    se_t1 = est1$sem,
    theta_t2 = est2$theta,
    se_t2 = est2$sem,
    extreme_t1 = ext1,
    extreme_t2 = ext2,
    change = change,
    se_diff = se_diff,
    rci = rci,
    p_value = p_value,
    crit_lower = crit$lower,
    crit_upper = crit$upper,
    change_class = change_class,
    retest_sd_tip = tip,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  attr(result, "null") <- null
  attr(result, "retest_sd") <- sigma_retest
  attr(result, "anchor") <- anchor
  attr(result, "alpha") <- alpha
  attr(result, "direction") <- direction
  attr(result, "critical") <- if (simulate_null) {
    "simulate"
  } else if (exact_null) {
    "exact"
  } else {
    critical
  }
  attr(result, "sim_iter") <- actual_iter
  attr(result, "method") <- method
  attr(result, "conditional_crit") <- conditional_crit

  if (output == "dataframe") {
    return(result)
  }

  n_clause <- .n_caption(
    n_persons,
    n_persons,
    if (anyNA(mat_t1) || anyNA(mat_t2)) {
      "incomplete responses retained"
    } else {
      character()
    }
  )
  caption <- .pc_caption(result, length(thr_list), n_clause)

  if (output == "kable") {
    return(.pc_kable(result, caption))
  }

  .pc_plot(
    result, thr_list, sigma_retest, caption, method, theta_range, prior_sd,
    key1 = apply(!is.na(mat_t1), 1L, function(z) paste0(which(z), collapse = ",")),
    key2 = apply(!is.na(mat_t2), 1L, function(z) paste0(which(z), collapse = ","))
  )
}

#' No-change band from the score-to-theta lookup
#'
#' For each attainable baseline score, finds the follow-up locations that would
#' not be flagged. The result is one rectangle per baseline score, spanning to
#' the midpoints of its neighbours on the x axis.
#'
#' @keywords internal
#' @noRd
.pc_band <- function(fig, key1, key2, thr_list, sigma_retest, method,
                     theta_range, prior_sd, lims) {
  # The dominant pair of answered-item sets; with complete data there is one.
  key <- paste(key1, key2, sep = "|")
  dom <- names(sort(table(key), decreasing = TRUE))[1L]
  parts <- strsplit(dom, "|", fixed = TRUE)[[1L]]
  a1 <- as.integer(strsplit(parts[1L], ",", fixed = TRUE)[[1L]])
  a2 <- as.integer(strsplit(parts[2L], ",", fixed = TRUE)[[1L]])

  l1 <- .pc_score_lookup(thr_list, a1, method, theta_range, prior_sd)
  l2 <- .pc_score_lookup(thr_list, a2, method, theta_range, prior_sd)

  cl <- stats::median(fig$crit_lower, na.rm = TRUE)
  cu <- stats::median(fig$crit_upper, na.rm = TRUE)

  x <- l1$theta
  o <- order(x)
  x <- x[o]
  edges <- c(lims[1L], (x[-1L] + x[-length(x)]) / 2, lims[2L])

  rows <- lapply(seq_along(x), function(i) {
    idx <- o[i]
    rci <- (l2$theta - l1$theta[idx]) /
      sqrt(l1$sem[idx]^2 + l2$sem^2 + 2 * sigma_retest^2)
    keep <- (!is.finite(cu) | rci <= cu) & (!is.finite(cl) | rci >= cl)
    if (!any(keep)) {
      return(NULL)
    }
    data.frame(
      xmin = edges[i], xmax = edges[i + 1L],
      lower = min(l2$theta[keep]), upper = max(l2$theta[keep]),
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
  out$lower <- pmax(out$lower, lims[1L])
  out$upper <- pmin(out$upper, lims[2L])
  out
}

# ---------------------------------------------------------------------------
# Internal: validation
# ---------------------------------------------------------------------------

#' Catch values passed to the wrong argument
#'
#' `method`, `estimator`, `null` and `direction` all take short character
#' values and are easy to mix up. Name the right argument rather than letting
#' `match.arg()` print a bare list of choices.
#'
#' @keywords internal
#' @noRd
.pc_check_swapped_args <- function(method, estimator, null, direction) {
  owners <- list(
    method = c("WLE", "EAP"),
    estimator = c("CML", "MML"),
    null = c("measurement", "retest"),
    direction = c("two.sided", "increase", "decrease")
  )
  supplied <- list(
    method = method,
    estimator = estimator,
    null = null,
    direction = direction
  )
  for (arg in names(supplied)) {
    val <- supplied[[arg]]
    if (!is.character(val) || length(val) != 1L) next
    if (val %in% owners[[arg]]) next
    hit <- names(owners)[vapply(owners, function(v) val %in% v, logical(1L))]
    if (length(hit) == 1L) {
      stop(
        "`", arg, "` was given \"", val, "\", which belongs to `", hit,
        "`. Use ", hit, " = \"", val, "\".",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' Check that the two occasions describe the same items and the same people
#'
#' @keywords internal
#' @noRd
.pc_check_paired <- function(data_t1, data_t2) {
  if (nrow(data_t1) != nrow(data_t2)) {
    stop(
      "`data_t1` and `data_t2` must have the same number of rows, with row i ",
      "being the same respondent at both occasions (",
      nrow(data_t1), " and ", nrow(data_t2), " supplied).",
      call. = FALSE
    )
  }
  if (nrow(data_t1) == 0L) {
    stop("`data_t1` has no rows.", call. = FALSE)
  }
  if (ncol(data_t1) != ncol(data_t2)) {
    stop(
      "`data_t1` and `data_t2` must contain the same items (",
      ncol(data_t1), " and ", ncol(data_t2), " columns supplied).",
      call. = FALSE
    )
  }
  if (!identical(colnames(data_t1), colnames(data_t2))) {
    stop(
      "`data_t1` and `data_t2` must have the same item names in the same ",
      "order.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Check that the thresholds cover the categories observed at either occasion
#'
#' @keywords internal
#' @noRd
.pc_check_categories <- function(thr_list, mat_t1, mat_t2) {
  for (i in seq_along(thr_list)) {
    thr <- thr_list[[i]]
    if (!is.numeric(thr) || length(thr) < 1L || !all(is.finite(thr))) {
      stop(
        "`item_params` entry ", i, " must be a finite numeric vector of ",
        "Andrich thresholds.",
        call. = FALSE
      )
    }
    max_cat <- suppressWarnings(max(
      c(mat_t1[, i], mat_t2[, i]),
      na.rm = TRUE
    ))
    if (is.finite(max_cat) && length(thr) < max_cat) {
      stop(
        "Item ", i, " has responses up to category ", max_cat,
        " but only ", length(thr), " threshold(s).",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' Flag respondents at the minimum or maximum score they could have reached
#'
#' @keywords internal
#' @noRd
.pc_extreme <- function(mat, steps) {
  answered <- !is.na(mat)
  n_answered <- rowSums(answered)
  rs <- rowSums(mat, na.rm = TRUE)
  max_score <- as.numeric(answered %*% steps)
  n_answered > 0L & (rs == 0 | rs == max_score)
}

# ---------------------------------------------------------------------------
# Internal: the test
# ---------------------------------------------------------------------------

#' Tail probabilities for the requested direction
#'
#' @keywords internal
#' @noRd
.pc_tail_probs <- function(direction, alpha) {
  switch(direction,
    two.sided = list(lo = alpha / 2, hi = 1 - alpha / 2),
    increase = list(lo = -Inf, hi = 1 - alpha),
    decrease = list(lo = alpha, hi = Inf)
  )
}

#' Plug-in location for the no-change null
#'
#' Under the null the respondent occupies one location at both occasions, so
#' the two estimates are combined by inverse-variance weighting. Any
#' reasonable plug-in serves, since the null distribution of the RCI varies
#' slowly with theta.
#'
#' @keywords internal
#' @noRd
.pc_theta_null <- function(est1, est2) {
  w1 <- 1 / est1$sem^2
  w2 <- 1 / est2$sem^2
  out <- (est1$theta * w1 + est2$theta * w2) / (w1 + w2)
  # Fall back to whichever estimate is usable when a weight is not finite
  bad <- !is.finite(out)
  out[bad] <- rowMeans(
    cbind(est1$theta, est2$theta)[bad, , drop = FALSE],
    na.rm = TRUE
  )
  out[!is.finite(out)] <- 0
  out
}

#' Critical values from the simulated null
#'
#' Quantiles of the *signed* RCI. With exchangeable occasions the null is
#' symmetric about zero and the two bounds agree, but they do not when the
#' occasions differ in missingness, so both are carried rather than assuming
#' symmetry and reading a quantile of the absolute RCI.
#'
#' @keywords internal
#' @noRd
.pc_critical <- function(sim_mat, probs, conditional_crit, n_persons) {
  q <- function(x, p) {
    if (!is.finite(p)) {
      return(if (p < 0) -Inf else Inf)
    }
    unname(stats::quantile(x, probs = p, na.rm = TRUE, type = 7))
  }
  if (isTRUE(conditional_crit)) {
    lower <- vapply(seq_len(nrow(sim_mat)), function(i) q(sim_mat[i, ], probs$lo), numeric(1L))
    upper <- vapply(seq_len(nrow(sim_mat)), function(i) q(sim_mat[i, ], probs$hi), numeric(1L))
  } else {
    pooled <- as.numeric(sim_mat)
    lower <- rep(q(pooled, probs$lo), n_persons)
    upper <- rep(q(pooled, probs$hi), n_persons)
  }
  list(lower = lower, upper = upper)
}

#' Monte-Carlo p-values against the simulated null
#'
#' Uses the house `(1 + count) / (B + 1)` convention, so a p-value is never
#' exactly zero. Ties matter here, because a discrete sum score gives the RCI
#' finitely many values, so counting is done with `left.open = TRUE`.
#'
#' @keywords internal
#' @noRd
.pc_pvalues <- function(rci, sim_mat, direction, conditional_crit) {
  count_ge <- function(sorted_null, x) {
    length(sorted_null) - findInterval(x, sorted_null, left.open = TRUE)
  }
  count_le <- function(sorted_null, x) findInterval(x, sorted_null)

  one <- function(null_vec, x) {
    null_vec <- sort(null_vec[is.finite(null_vec)])
    B <- length(null_vec)
    if (B == 0L || !is.finite(x)) {
      return(NA_real_)
    }
    cnt <- switch(direction,
      two.sided = count_ge(sort(abs(null_vec)), abs(x)),
      increase = count_ge(null_vec, x),
      decrease = count_le(null_vec, x)
    )
    (1 + cnt) / (B + 1)
  }

  if (isTRUE(conditional_crit)) {
    return(vapply(
      seq_along(rci),
      function(i) one(sim_mat[i, ], rci[i]),
      numeric(1L)
    ))
  }

  pooled <- sort(as.numeric(sim_mat)[is.finite(as.numeric(sim_mat))])
  B <- length(pooled)
  if (B == 0L) {
    return(rep(NA_real_, length(rci)))
  }
  pooled_abs <- sort(abs(pooled))
  cnt <- switch(direction,
    two.sided = count_ge(pooled_abs, abs(rci)),
    increase = count_ge(pooled, rci),
    decrease = count_le(pooled, rci)
  )
  p <- (1 + cnt) / (B + 1)
  p[!is.finite(rci)] <- NA_real_
  p
}

#' Normal-reference p-values, used when `critical` is numeric
#'
#' @keywords internal
#' @noRd
.pc_normal_p <- function(rci, direction) {
  switch(direction,
    two.sided = 2 * stats::pnorm(-abs(rci)),
    increase = stats::pnorm(rci, lower.tail = FALSE),
    decrease = stats::pnorm(rci)
  )
}

#' Classify each respondent against the critical values
#'
#' Levels name the direction of theta, not a clinical reading, because whether
#' an increase is an improvement depends on how the scale is oriented.
#'
#' @keywords internal
#' @noRd
.pc_classify <- function(rci, crit) {
  out <- rep("none detected", length(rci))
  out[is.finite(rci) & rci > crit$upper] <- "increase"
  out[is.finite(rci) & rci < crit$lower] <- "decrease"
  out[!is.finite(rci)] <- NA_character_
  factor(out, levels = c("decrease", "none detected", "increase"))
}

#' Retest SD that would overturn a flagged result
#'
#' @keywords internal
#' @noRd
.pc_tipping <- function(change, se1, se2, crit, change_class) {
  cval <- rep(NA_real_, length(change))
  cval[!is.na(change_class) & change_class == "increase"] <-
    crit$upper[!is.na(change_class) & change_class == "increase"]
  cval[!is.na(change_class) & change_class == "decrease"] <-
    abs(crit$lower[!is.na(change_class) & change_class == "decrease"])

  v <- (change^2 / cval^2 - se1^2 - se2^2) / 2
  out <- rep(NA_real_, length(change))
  ok <- is.finite(v) & v > 0
  out[ok] <- sqrt(v[ok])
  # Flagged but the algebra gives a non-positive variance: no occasion noise
  # at all is tolerated.
  out[is.finite(v) & v <= 0] <- 0
  out
}

# ---------------------------------------------------------------------------
# Internal: null simulation
# ---------------------------------------------------------------------------

#' Simulate the no-change null
#'
#' Item parameters are held fixed at the observed calibration, so no iteration
#' can fail to converge. Each respondent is placed at `theta_null`, two
#' response vectors are drawn, and the RCI is recomputed. Returns a
#' respondents-by-iterations matrix, which serves both the pooled and the
#' per-respondent critical values from one simulation.
#'
#' `return_change = TRUE` returns the raw change instead of the RCI, which is
#' what [RMretestSD()] needs for its simulated error-variance term.
#'
#' @keywords internal
#' @noRd
.pc_simulate <- function(
  theta_null,
  thr_list,
  sigma_retest,
  method,
  theta_range,
  na1,
  na2,
  sim_iter,
  parallel,
  n_cores,
  seed,
  verbose,
  prior_sd = NULL,
  return_change = FALSE
) {
  use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)
  if (parallel && !use_parallel) {
    message(
      "Install 'mirai' for parallel processing: install.packages(\"mirai\")"
    )
    message("Running sequentially...")
  }
  if (use_parallel) {
    if (is.null(n_cores)) {
      n_cores <- getOption("mc.cores")
    }
    if (is.null(n_cores)) {
      warning(
        "For parallel processing, specify n_cores or set options(mc.cores = N).\n",
        "Falling back to sequential.",
        call. = FALSE
      )
      use_parallel <- FALSE
    } else {
      n_cores <- min(n_cores, sim_iter)
    }
  }

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max - 2L, 1L)
  }
  set.seed(seed)
  sim_seeds <- sample.int(.Machine$integer.max, sim_iter)

  sim_args <- list(
    theta_null = theta_null,
    thr_list = thr_list,
    sigma_retest = sigma_retest,
    method = method,
    theta_range = theta_range,
    na1 = if (any(na1)) na1 else NULL,
    na2 = if (any(na2)) na2 else NULL,
    prior_sd = prior_sd,
    return_change = return_change
  )

  results <- if (use_parallel) {
    .pc_sim_parallel(sim_iter, sim_seeds, sim_args, n_cores, verbose)
  } else {
    .pc_sim_sequential(sim_iter, sim_seeds, sim_args, verbose)
  }

  ok <- vapply(results, is.numeric, logical(1L))
  if (sum(ok) < 2L) {
    stop(
      "Fewer than 2 simulation iterations succeeded, so the null could not ",
      "be calibrated. First message: ",
      if (any(!ok)) as.character(results[!ok][[1L]]) else "none",
      call. = FALSE
    )
  }
  if (any(!ok)) {
    warning(
      sum(!ok), " of ", sim_iter, " simulation iterations failed and were ",
      "dropped.",
      call. = FALSE
    )
  }

  matrix(
    unlist(results[ok], use.names = FALSE),
    nrow = length(theta_null),
    ncol = sum(ok)
  )
}

#' One iteration of the no-change null
#'
#' @keywords internal
#' @noRd
run_single_change_sim <- function(seed, sim_args) {
  # The RNG kind is pinned, not just the seed: mirai daemons start under
  # L'Ecuyer-CMRG while the calling session uses the Mersenne-Twister
  # default, so seeding alone would make the parallel and sequential paths
  # draw different streams from the same `seed`.
  set.seed(
    seed,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )

  tryCatch(
    {
      th <- sim_args$theta_null
      n <- length(th)
      s <- sim_args$sigma_retest
      # Under null = "retest" each occasion carries its own deviation, which
      # is why retest_sd enters the analytic SE_diff doubled.
      th1 <- if (s > 0) th + stats::rnorm(n, 0, s) else th
      th2 <- if (s > 0) th + stats::rnorm(n, 0, s) else th

      r1 <- .pc_sim_responses(sim_args$thr_list, th1, sim_args$na1)
      r2 <- .pc_sim_responses(sim_args$thr_list, th2, sim_args$na2)

      # The prior is held at the one the observed respondents were scored
      # under. Re-estimating it from each null dataset would score the null on
      # a different scale from the observation it is being compared with.
      e1 <- .estimate_thetas(r1, sim_args$thr_list, method = sim_args$method,
                             theta_range = sim_args$theta_range,
                             prior_sd = sim_args$prior_sd)
      e2 <- .estimate_thetas(r2, sim_args$thr_list, method = sim_args$method,
                             theta_range = sim_args$theta_range,
                             prior_sd = sim_args$prior_sd)

      change <- e2$theta - e1$theta
      if (isTRUE(sim_args$return_change)) {
        change
      } else {
        change / sqrt(e1$sem^2 + e2$sem^2 + 2 * s^2)
      }
    },
    error = function(e) as.character(conditionMessage(e))
  )
}

#' Draw one occasion's responses, carrying the observed missingness pattern
#'
#' `sim_partial_score()` drops to a vector for a single respondent, so the
#' result is reshaped explicitly. That matters because a single respondent is
#' the headline use case.
#'
#' @keywords internal
#' @noRd
.pc_sim_responses <- function(thr_list, theta, na_mask) {
  out <- matrix(
    unlist(sim_partial_score(thr_list, theta), use.names = FALSE),
    nrow = length(theta),
    ncol = length(thr_list)
  )
  if (!is.null(na_mask)) {
    out[na_mask] <- NA_integer_
  }
  out
}

#' @keywords internal
#' @noRd
.pc_sim_parallel <- function(
  sim_iter,
  sim_seeds,
  sim_args,
  n_cores,
  verbose = FALSE
) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = sim_iter, style = 3)
  }

  tasks <- lapply(seq_len(sim_iter), function(i) {
    mirai::mirai(
      {
        run_single_change_sim(seed, sim_args)
      },
      seed = sim_seeds[i],
      sim_args = sim_args,
      run_single_change_sim = run_single_change_sim,
      .pc_sim_responses = .pc_sim_responses,
      sim_partial_score = sim_partial_score,
      sim_poly_item = sim_poly_item,
      .estimate_thetas = .estimate_thetas,
      .theta_wle = .theta_wle,
      .theta_eap = .theta_eap,
      .pcm_cat_probs = .pcm_cat_probs,
      .logp_tables = .logp_tables,
      .grid_loglik = .grid_loglik,
      .estimate_prior_sd = .estimate_prior_sd
    )
  })

  results <- vector("list", sim_iter)
  for (i in seq_len(sim_iter)) {
    res <- mirai::call_mirai(tasks[[i]])$data
    results[[i]] <- if (inherits(res, "errorValue")) "mirai_error" else res
    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}

#' @keywords internal
#' @noRd
.pc_sim_sequential <- function(sim_iter, sim_seeds, sim_args, verbose = FALSE) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = sim_iter, style = 3)
  }

  results <- vector("list", sim_iter)
  for (i in seq_len(sim_iter)) {
    results[[i]] <- run_single_change_sim(sim_seeds[i], sim_args)
    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}

# ---------------------------------------------------------------------------
# Internal: output
# ---------------------------------------------------------------------------

#' Caption naming the null in force
#'
#' The null is named on every reporting surface, because the measurement and
#' retest nulls answer different questions and the same RCI machinery serves
#' both.
#'
#' @keywords internal
#' @noRd
.pc_caption <- function(result, n_items, n_clause) {
  null <- attr(result, "null")
  sigma <- attr(result, "retest_sd")
  crit <- attr(result, "critical")
  alpha <- attr(result, "alpha")
  direction <- attr(result, "direction")
  anchor <- attr(result, "anchor")

  null_txt <- if (null == "measurement") {
    paste0(
      "Null: no change beyond measurement error, ",
      "SE_diff = sqrt(SE_1^2 + SE_2^2). This is the narrower of the two ",
      "nulls available and does not allow for occasion-to-occasion ",
      "fluctuation, so a flagged change is a necessary condition for change ",
      "rather than evidence of it."
    )
  } else {
    sprintf(
      paste0(
        "Null: no change beyond measurement error and occasion-to-occasion ",
        "fluctuation, SE_diff = sqrt(SE_1^2 + SE_2^2 + 2 x %.3f^2)."
      ),
      sigma
    )
  }

  crit_txt <- if (identical(crit, "exact")) {
    if (isTRUE(attr(result, "conditional_crit"))) {
      sprintf(
        paste0(
          "Critical values computed exactly by enumerating the null, one pair ",
          "per respondent (see crit_lower and crit_upper); medians %s and %s."
        ),
        .pc_fmt_crit(stats::median(result$crit_lower, na.rm = TRUE)),
        .pc_fmt_crit(stats::median(result$crit_upper, na.rm = TRUE))
      )
    } else {
      sprintf(
        paste0(
          "Critical values computed exactly by enumerating the null, pooled ",
          "across respondents: %s and %s."
        ),
        .pc_fmt_crit(result$crit_lower[1L]),
        .pc_fmt_crit(result$crit_upper[1L])
      )
    }
  } else if (identical(crit, "simulate")) {
    if (isTRUE(attr(result, "conditional_crit"))) {
      sprintf(
        paste0(
          "Critical values simulated over %d iterations, one pair per ",
          "respondent (see crit_lower and crit_upper); medians %s and %s."
        ),
        attr(result, "sim_iter"),
        .pc_fmt_crit(stats::median(result$crit_lower, na.rm = TRUE)),
        .pc_fmt_crit(stats::median(result$crit_upper, na.rm = TRUE))
      )
    } else {
      sprintf(
        paste0(
          "Critical values simulated over %d iterations, pooled across ",
          "respondents: %s and %s."
        ),
        attr(result, "sim_iter"),
        .pc_fmt_crit(result$crit_lower[1L]),
        .pc_fmt_crit(result$crit_upper[1L])
      )
    }
  } else {
    sprintf(
      "Critical values fixed at %s and %s, referred to a normal null.",
      .pc_fmt_crit(result$crit_lower[1L]),
      .pc_fmt_crit(result$crit_upper[1L])
    )
  }

  cal_txt <- if (anchor == "supplied") {
    "Item parameters supplied, not estimated from these data."
  } else {
    sprintf("Item parameters from %d items, calibrated on %s.", n_items,
      switch(anchor,
        stack = "both occasions stacked",
        t1 = "occasion 1 only",
        t2 = "occasion 2 only"
      )
    )
  }

  flagged <- sum(!is.na(result$change_class) &
                   result$change_class != "none detected")

  paste(
    null_txt,
    crit_txt,
    cal_txt,
    sprintf(
      "%s test at alpha = %.3g per respondent.",
      switch(direction,
        two.sided = "Two-sided",
        increase = "One-sided (increase)",
        decrease = "One-sided (decrease)"
      ),
      alpha
    ),
    sprintf("%d of %d respondents flagged.", flagged, nrow(result)),
    paste0(n_clause, "."),
    sep = " "
  )
}

#' Format a critical value, which is infinite on the unused side of a
#' one-sided test
#'
#' @keywords internal
#' @noRd
.pc_fmt_crit <- function(x) {
  if (!is.finite(x)) {
    return(if (x < 0) "-Inf" else "Inf")
  }
  sprintf("%.2f", x)
}

#' @keywords internal
#' @noRd
.pc_kable <- function(result, caption) {
  display <- .round_display(
    result,
    c(theta_t1 = 2, se_t1 = 3, theta_t2 = 2, se_t2 = 3, change = 2,
      se_diff = 3, rci = 2, p_value = 4, crit_lower = 2, crit_upper = 2,
      retest_sd_tip = 3)
  )
  cols <- c(
    "id", "sum_t1", "sum_t2", "theta_t1", "se_t1", "theta_t2", "se_t2",
    "change", "se_diff", "rci", "p_value"
  )
  names_out <- c(
    "ID", "Sum T1", "Sum T2", "Theta T1", "SE T1", "Theta T2", "SE T2",
    "Change", "SE diff", "RCI", "p"
  )
  # The critical values are constant unless they were simulated per
  # respondent, so they are reported once in the caption and shown per row
  # only when they actually vary.
  if (isTRUE(attr(result, "conditional_crit"))) {
    cols <- c(cols, "crit_lower", "crit_upper")
    names_out <- c(names_out, "Crit lower", "Crit upper")
  }
  cols <- c(cols, "change_class", "retest_sd_tip")
  names_out <- c(names_out, "Class", "Tipping SD")

  knitr::kable(
    display[, cols, drop = FALSE],
    format = "pipe",
    col.names = names_out,
    caption = caption
  )
}

#' Occasion 1 against occasion 2, with the no-change band
#'
#' The band is built from the same score-to-theta lookup the test uses, so a
#' point lies outside it if and only if that respondent was flagged. An earlier
#' version drew a smooth band of half-width crit * sqrt(2/I(theta)), placing
#' both occasions at the baseline location; that disagreed with the
#' classification for about a tenth of respondents, worst where a large change
#' lands near a boundary and the follow-up standard error is much the larger of
#' the two.
#'
#' Because attainable scores are discrete the band is a step function. It is
#' widest away from the centre, where the standard error is large but the scale
#' has not run out of room, and narrows again at the boundary itself because
#' there are no further attainable scores to move to. With mixed missingness
#' the band is drawn for the most common pair of answered-item sets, and under
#' `conditional_crit = TRUE` it uses the median critical values, so in those
#' two cases it is indicative rather than exact.
#'
#' @keywords internal
#' @noRd
.pc_plot <- function(result, thr_list, sigma_retest, caption,
                     method = "WLE", theta_range = c(-10, 10),
                     prior_sd = NULL, key1 = NULL, key2 = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for output = \"ggplot\". ",
      "Install it with: install.packages(\"ggplot2\")",
      call. = FALSE
    )
  }

  finite <- is.finite(result$theta_t1) & is.finite(result$theta_t2)
  fig <- result[finite, , drop = FALSE]
  key1 <- key1[finite]
  key2 <- key2[finite]
  if (nrow(fig) == 0L) {
    stop("No respondent has finite locations at both occasions.", call. = FALSE)
  }

  rng <- range(c(fig$theta_t1, fig$theta_t2), na.rm = TRUE)
  pad <- 0.05 * diff(rng)
  lims <- c(rng[1L] - pad, rng[2L] + pad)

  band <- .pc_band(
    fig, key1, key2, thr_list, sigma_retest, method, theta_range, prior_sd,
    lims
  )

  ggplot2::ggplot(fig, ggplot2::aes(x = .data$theta_t1, y = .data$theta_t2)) +
    ggplot2::geom_rect(
      data = band,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$lower, ymax = .data$upper),
      inherit.aes = FALSE,
      fill = "grey70",
      alpha = 0.35
    ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", colour = "grey30", linewidth = 0.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$change_class),
      size = 2, alpha = 0.85
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "decrease" = "#D95F02",
        "none detected" = "grey45",
        "increase" = "#1B9E77"
      ),
      drop = FALSE,
      name = "Change in theta"
    ) +
    # The band flares where information collapses, so the panel is held to the
    # range of the data rather than to the range of the band.
    ggplot2::coord_equal(xlim = lims, ylim = lims) +
    ggplot2::labs(
      x = "Person location, occasion 1 (logits)",
      y = "Person location, occasion 2 (logits)",
      caption = er2_caption(caption)
    ) +
    ggplot2::theme_bw() +
    er2_axis_margins() +
    er2_plot_caption()
}
