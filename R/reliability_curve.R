#' Conditional measurement precision across the latent scale
#'
#' Plots (or tabulates) how precisely a scale measures at each point of the
#' latent continuum, rather than collapsing precision to the single number
#' reported by [RMreliability()]. The curve is the conditional standard error
#' of measurement \eqn{SEM(\theta) = 1/\sqrt{I(\theta)}}, the test information
#' \eqn{I(\theta)} it comes from, or the conditional reliability derived from
#' either.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed.
#' @param statistic Character. The quantity on the y-axis: `"sem"` (default,
#'   the conditional standard error in logits), `"information"` (test
#'   information), or `"reliability"` (see Details for the formula used, and
#'   why). This is the argument that chooses what the figure shows, not
#'   `method`.
#' @param method Character. Person-location estimator, `"WLE"` (default) or
#'   `"EAP"`. **This does not change the curve**, which is always built from
#'   test information. It governs only the respondent locations behind the
#'   density overlay and the `benchmark` percentage, and is ignored when
#'   `show_density = FALSE` and `benchmark = NULL`.
#' @param benchmark Numeric in (0, 1) or `NULL` (default). When supplied, the
#'   region of the scale whose conditional reliability reaches `benchmark` is
#'   shaded, and the percentage of respondents located inside it is reported.
#'   There is deliberately no default value; see Details.
#' @param reference Character. Horizontal line showing the flat summary that a
#'   single coefficient implies: `"marginal"` (default) or `"none"`.
#' @param items Optional character vector of column names, or numeric column
#'   indices, selecting a subset of items. Useful for comparing a short form
#'   against the full scale.
#' @param item_params Optional pre-specified item parameters, used in place of
#'   the thresholds estimated from `data`. Either a named list of Andrich
#'   threshold vectors or a long-format data.frame from
#'   [RMitemParameters()], the same two forms [RMpersonParameters()] accepts.
#'   Supply this to draw the curve from anchored or previously published
#'   parameters. As in [RMpersonParameters()], supplied thresholds are used as
#'   given and are **not** re-centred, so their origin defines the origin of
#'   the theta axis.
#' @param boot Logical. If `TRUE`, add a bootstrap confidence band to the curve
#'   by resampling respondents. Default `FALSE`.
#' @param boot_iter Integer. Bootstrap iterations when `boot = TRUE`. Default
#'   `200`.
#' @param conf_int Numeric in (0, 1). HDCI width for the bootstrap band.
#'   Default `0.95`.
#' @param parallel Logical. Use `mirai` for the bootstrap if available.
#'   Default `TRUE`.
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first; if neither is set, the bootstrap
#'   falls back to sequential.
#' @param seed Integer or `NULL`. Random seed for the bootstrap. See
#'   [easyRasch2-reproducibility].
#' @param show_density Logical. Draw the respondent location distribution as a
#'   background band. Default `TRUE`.
#' @param theta_range Numeric length-2 vector, or `NULL` (default) for
#'   \eqn{\pm 3\sigma}.
#' @param n_nodes Integer. Number of points on the theta grid. Default `161`.
#' @param verbose Logical. Print a progress bar for the bootstrap. Default
#'   `FALSE`.
#' @param output Character. `"ggplot"` (default), `"dataframe"` for the curve,
#'   or `"kable"` for a summary table.
#'
#' @return
#' * If `output = "ggplot"`: a `ggplot` object.
#' * If `output = "dataframe"`: a data.frame with one row per grid point and
#'   columns `theta`, `information`, `sem`, `reliability`, plus
#'   `<statistic>_lower` / `<statistic>_upper` band columns when
#'   `boot = TRUE`. Attributes carry the summary quantities: `sigma`,
#'   `marginal_ratio` (the latent-density-weighted mean of the reliability
#'   curve, matching [RMreliability()]), `marginal_green` (the superseded
#'   subtractive coefficient), `sem_average` (root mean error variance),
#'   `benchmark`,
#'   `benchmark_range`, `benchmark_percent` and `n_not_estimable`.
#' * If `output = "kable"`: a `knitr_kable` summary table of those quantities.
#'
#' @details
#' **What the curve is.** Test information \eqn{I(\theta) = \sum_i
#' \mathrm{Var}_i(\mathrm{score} \mid \theta)} is summed from the CML item
#' thresholds (`psychotools`), and the conditional standard error follows as
#' \eqn{1/\sqrt{I(\theta)}}. Both are properties of the item set alone. The
#' respondent distribution is drawn behind the curve because precision only
#' matters where people actually are, but it does not enter the curve.
#'
#' **Conditional reliability, and which formula.** Two definitions circulate:
#'
#' \deqn{\rho(\theta) = \frac{\sigma^2}{\sigma^2 + SEM(\theta)^2}
#'   \qquad \mathrm{and} \qquad
#'   \rho(\theta) = 1 - \frac{SEM(\theta)^2}{\sigma^2}}
#'
#' with \eqn{\sigma} the latent SD estimated by marginal maximum likelihood.
#' This function uses the first, the ratio form, for three reasons. It is
#' bounded in (0, 1), whereas the subtractive form returns negative values
#' whenever \eqn{SEM(\theta) > \sigma}, which is common at the floor of a
#' skewed scale and at both tails of a short one. It treats \eqn{\sigma^2} as
#' true-score variance, which is what the MML estimate of the latent SD is.
#' And it tracks the reliability of the observed scores far more closely.
#' Milanzi et al. (2015) report the subtractive coefficient falling to
#' \eqn{-0.120} against a true value of \eqn{0.480} for a 1PL with
#' \eqn{\sigma^2 = 0.25}. Replicating their design (`dev/milanzi_check.R`,
#' 1000 replications per cell) gives a mean absolute error against the exact
#' expected-sum-score reliability of \eqn{0.22} for the subtractive form and
#' \eqn{0.018} for the ratio form, with the subtractive form negative in 99.7%
#' of replications of their lowest-variance cell. For polytomous Rasch data at
#' 5 to 20 items the corresponding errors are \eqn{0.18} and \eqn{0.015}.
#'
#' The advantage is concentrated where it matters rather than uniform. Where
#' information is high relative to the trait variance the two forms are close
#' and either can be nearer the truth. The gap comes from the low-information
#' cases, where the subtractive form does not merely lose accuracy but leaves
#' the (0, 1) interval altogether.
#'
#' `marginal_ratio` is the latent-density-weighted mean of the reliability
#' curve rather than the ratio formed from the averaged error variance. The two
#' differ by Jensen's inequality, and the curve mean is the more accurate of
#' the pair (mean absolute error 0.015 against 0.018 for binary data, 0.010
#' against 0.015 for polytomous).
#'
#' [RMreliability()] reports the same ratio-form coefficient in its
#' "Marginal (curve mean)" row, so the scalar and the curve agree. The
#' superseded subtractive value is still returned as `marginal_green` for
#' comparison with easyRasch2 1.2.0 and earlier, and with `mirt::marginal_rxx()`
#' and similar software.
#'
#' **Limits.** The information-based standard error is asymptotic, and is only
#' an approximation for short scales; Milanzi et al. (2015) find the worst
#' behaviour with six items. The curve describes the precision of the person
#' location estimate on the logit scale, which is not the same quantity as the
#' reliability of a raw sum score, and reliability computed on a latent scale
#' is consistently higher than its manifest counterpart. Use [RMscoreSE()] for
#' the raw-score view.
#'
#' **No default benchmark.** `benchmark` is `NULL` unless asked for, in keeping
#' with the package's treatment of fixed rules of thumb. McNeish and Dumas
#' (2025), whose respondent-weighted summary this borrows, are explicit that
#' their own interpretive percentages are heuristic and should not be used as
#' cutoffs.
#'
#' @references
#' Green, B. F., Bock, R. D., Humphreys, L. G., Linn, R. L., & Reckase, M. D.
#' (1984). Technical Guidelines for Assessing Computerized Adaptive Tests.
#' *Journal of Educational Measurement, 21*(4), 347-360.
#' \doi{10.1111/j.1745-3984.1984.tb01039.x}
#'
#' McNeish, D., & Dumas, D. (2025). Reliability representativeness: How well
#' does coefficient alpha summarize reliability across the score distribution?
#' *Behavior Research Methods, 57*(3), 93.
#' \doi{10.3758/s13428-025-02611-8}
#'
#' Milanzi, E., Molenberghs, G., Alonso, A., Verbeke, G., & De Boeck, P.
#' (2015). Reliability measures in item response theory: Manifest versus latent
#' correlation functions. *British Journal of Mathematical and Statistical
#' Psychology, 68*(1), 43-64. \doi{10.1111/bmsp.12033}
#'
#' @seealso [RMreliability()], [RMscoreSE()], [RMtargeting()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # Conditional SEM across the scale
#'   RMreliabilityCurve(phq9[, 1:9])
#'
#'   # Conditional reliability, with the region reaching .8 shaded
#'   RMreliabilityCurve(phq9[, 1:9], statistic = "reliability",
#'                      benchmark = 0.8)
#' }
#'
#' # Summary quantities
#' RMreliabilityCurve(phq9[, 1:9], benchmark = 0.8, output = "kable")
#' }
RMreliabilityCurve <- function(
  data,
  statistic = "sem",
  method = "WLE",
  benchmark = NULL,
  reference = "marginal",
  items = NULL,
  item_params = NULL,
  boot = FALSE,
  boot_iter = 200,
  conf_int = 0.95,
  parallel = TRUE,
  n_cores = NULL,
  seed = NULL,
  show_density = TRUE,
  theta_range = NULL,
  n_nodes = 161L,
  verbose = FALSE,
  output = "ggplot"
) {
  # `statistic` picks the quantity on the y-axis and `method` the person-location
  # estimator. Reaching for the wrong one is easy, so say which is which rather
  # than letting match.arg() report a bare list of choices.
  stat_choices <- c("sem", "information", "reliability")
  method_choices <- c("WLE", "EAP")
  if (is.character(method) && length(method) == 1L && method %in% stat_choices) {
    stop(
      "`method` selects the person-location estimator (\"WLE\" or \"EAP\"). ",
      "To put ", method, " on the y-axis, use statistic = \"", method, "\".",
      call. = FALSE
    )
  }
  if (
    is.character(statistic) && length(statistic) == 1L &&
      statistic %in% method_choices
  ) {
    stop(
      "`statistic` selects the quantity on the y-axis (\"sem\", ",
      "\"information\" or \"reliability\"). To choose the person-location ",
      "estimator, use method = \"", statistic, "\".",
      call. = FALSE
    )
  }
  statistic <- match.arg(statistic, stat_choices)
  method <- match.arg(method, method_choices)
  reference <- match.arg(reference, c("marginal", "none"))
  output <- match.arg(output, c("ggplot", "dataframe", "kable"))

  validate_response_data(data)

  if (!is.null(benchmark)) {
    if (
      !is.numeric(benchmark) ||
        length(benchmark) != 1L ||
        !is.finite(benchmark) ||
        benchmark <= 0 ||
        benchmark >= 1
    ) {
      stop(
        "`benchmark` must be a single number strictly between 0 and 1, ",
        "or NULL.",
        call. = FALSE
      )
    }
  }

  if (
    !is.numeric(conf_int) || length(conf_int) != 1L ||
      conf_int <= 0 || conf_int >= 1
  ) {
    stop("`conf_int` must be a numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(n_nodes) || length(n_nodes) != 1L ||
      !is.finite(n_nodes) || n_nodes < 11
  ) {
    stop("`n_nodes` must be a single number of at least 11.", call. = FALSE)
  }
  n_nodes <- as.integer(n_nodes)

  if (!is.null(theta_range)) {
    if (
      !is.numeric(theta_range) || length(theta_range) != 2L ||
        !all(is.finite(theta_range)) || theta_range[1L] >= theta_range[2L]
    ) {
      stop(
        "`theta_range` must be a numeric vector of length 2 with ",
        "theta_range[1] < theta_range[2].",
        call. = FALSE
      )
    }
  }

  data <- as.data.frame(data)
  if (!is.null(items)) {
    data <- .curve_select_items(data, items)
  }

  # Supplied parameters decide which items are on the scale, so they are
  # coerced and aligned before the sample bookkeeping below.
  fixed_params <- !is.null(item_params)
  thr_list <- NULL
  if (fixed_params) {
    thr_list <- .coerce_item_params(item_params)
    data <- .align_items(data, thr_list)
  }

  if (ncol(data) < 2L) {
    stop("`data` must contain at least two items.", call. = FALSE)
  }

  n_total <- nrow(data)
  has_na <- anyNA(data)
  data <- .drop_empty_respondents(data)
  n_used <- nrow(data)
  if (n_used < 2L) {
    stop("Fewer than two respondents with any responses.", call. = FALSE)
  }
  data_mat <- as.matrix(data)
  n_items <- ncol(data_mat)

  # --- Item thresholds -------------------------------------------------------
  if (fixed_params) {
    .curve_check_categories(thr_list, data_mat)
  } else {
    thr_list <- .fit_cml_thresholds(data_mat)
  }

  # --- Latent SD -------------------------------------------------------------
  sigma <- .latent_sd(data_mat, thr_list)
  if (!is.finite(sigma) || sigma <= 0) {
    stop(
      "The latent SD could not be estimated, so the reliability scaling is ",
      "undefined. Check for items with no variation.",
      call. = FALSE
    )
  }

  # --- The curve -------------------------------------------------------------
  if (is.null(theta_range)) {
    theta_range <- c(-3 * sigma, 3 * sigma)
  }
  grid <- seq(theta_range[1L], theta_range[2L], length.out = n_nodes)
  curve <- .curve_stats(thr_list, grid, sigma)

  # --- Marginal summaries ----------------------------------------------------
  # Same helper `.marginal_rxx()` uses, so `marginal_ratio` is the identical
  # quantity to the "Marginal (curve mean)" row of RMreliability() by
  # construction rather than by coincidence. Independent of `n_nodes` and
  # `theta_range`, which govern the plotted grid only.
  marg <- .marginal_summaries(thr_list, sigma)
  sem_bar <- marg$sem_average
  marginal_ratio <- marg$ratio
  marginal_green <- marg$green

  # --- Respondent locations --------------------------------------------------
  est <- .estimate_thetas(
    data_mat,
    thr_list,
    method = method,
    theta_range = c(-10, 10)
  )
  theta_hat <- est$theta[is.finite(est$theta)]
  n_not_estimable <- n_used - length(theta_hat)

  # --- Benchmark -------------------------------------------------------------
  benchmark_range <- NULL
  benchmark_percent <- NA_real_
  if (!is.null(benchmark)) {
    if (length(theta_hat) > 0L) {
      rxx_person <- .curve_stats(thr_list, theta_hat, sigma)$reliability
      benchmark_percent <- 100 * mean(rxx_person >= benchmark)
    }
    benchmark_range <- .curve_runs(grid, curve$reliability >= benchmark)
  }

  # --- Bootstrap band --------------------------------------------------------
  band <- NULL
  actual_boot <- NA_integer_
  if (isTRUE(boot)) {
    boot_out <- .curve_bootstrap(
      data = data,
      thr_list = thr_list,
      fixed_params = fixed_params,
      grid = grid,
      boot_iter = boot_iter,
      conf_int = conf_int,
      parallel = parallel,
      n_cores = n_cores,
      seed = seed,
      verbose = verbose
    )
    band <- boot_out$band
    actual_boot <- boot_out$n_ok
  }

  # --- Assemble --------------------------------------------------------------
  curve_df <- data.frame(
    theta = grid,
    information = curve$information,
    sem = curve$sem,
    reliability = curve$reliability,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  if (!is.null(band)) {
    curve_df$information_lower <- band$information[, 1L]
    curve_df$information_upper <- band$information[, 2L]
    curve_df$sem_lower <- band$sem[, 1L]
    curve_df$sem_upper <- band$sem[, 2L]
    curve_df$reliability_lower <- band$reliability[, 1L]
    curve_df$reliability_upper <- band$reliability[, 2L]
  }

  attr(curve_df, "sigma") <- sigma
  attr(curve_df, "marginal_ratio") <- marginal_ratio
  attr(curve_df, "marginal_green") <- marginal_green
  attr(curve_df, "sem_average") <- sem_bar
  attr(curve_df, "benchmark") <- benchmark
  attr(curve_df, "benchmark_range") <- benchmark_range
  attr(curve_df, "benchmark_percent") <- benchmark_percent
  attr(curve_df, "n_not_estimable") <- n_not_estimable
  attr(curve_df, "boot_iter") <- actual_boot

  if (output == "dataframe") {
    return(curve_df)
  }

  # --- Caption pieces --------------------------------------------------------
  qualifiers <- if (has_na) "incomplete responses retained" else character()
  n_clause <- .n_caption(n_used, n_total, qualifiers)

  if (output == "kable") {
    return(.curve_kable(
      curve_df,
      n_items = n_items,
      n_clause = n_clause,
      method = method,
      fixed_params = fixed_params
    ))
  }

  # --- ggplot ----------------------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for output = \"ggplot\". ",
      "Install it with: install.packages(\"ggplot2\")",
      call. = FALSE
    )
  }

  .curve_plot(
    curve_df,
    statistic = statistic,
    reference = reference,
    sem_bar = sem_bar,
    marginal_ratio = marginal_ratio,
    sigma = sigma,
    benchmark = benchmark,
    benchmark_range = benchmark_range,
    benchmark_percent = benchmark_percent,
    theta_hat = if (isTRUE(show_density)) theta_hat else numeric(0),
    n_items = n_items,
    n_clause = n_clause,
    method = method,
    fixed_params = fixed_params,
    conf_int = conf_int,
    actual_boot = actual_boot
  )
}

# ---------------------------------------------------------------------------
# Internal: curve construction
# ---------------------------------------------------------------------------

#' Information, conditional SEM and conditional reliability at given thetas
#'
#' The reliability is the bounded ratio form
#' \eqn{\sigma^2/(\sigma^2 + SEM(\theta)^2)}, the same form `.marginal_rxx()`
#' integrates. See the Details section of [RMreliabilityCurve()].
#'
#' @param thr_list List of Andrich threshold vectors.
#' @param theta Numeric vector of latent locations.
#' @param sigma Latent SD.
#' @return A list with `information`, `sem` and `reliability`.
#' @keywords internal
#' @noRd
.curve_stats <- function(thr_list, theta, sigma) {
  info <- .test_information(thr_list, theta)
  sem <- 1 / sqrt(info)
  list(
    information = info,
    sem = sem,
    reliability = sigma^2 / (sigma^2 + sem^2)
  )
}

#' Contiguous runs of a logical vector, as x-axis intervals
#'
#' @param x Numeric grid.
#' @param ok Logical vector, same length as `x`.
#' @return data.frame with `xmin` and `xmax`, or `NULL` when nothing qualifies.
#' @keywords internal
#' @noRd
.curve_runs <- function(x, ok) {
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) {
    return(NULL)
  }
  r <- rle(ok)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  keep <- r$values
  data.frame(
    xmin = x[starts[keep]],
    xmax = x[ends[keep]],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Resolve the `items` argument to a column subset
#'
#' @param data A data.frame of item responses.
#' @param items Character column names or numeric indices.
#' @return `data` restricted to the selected columns.
#' @keywords internal
#' @noRd
.curve_select_items <- function(data, items) {
  if (is.character(items)) {
    missing_items <- setdiff(items, names(data))
    if (length(missing_items) > 0L) {
      stop(
        "`items` not found in `data`: ",
        paste(missing_items, collapse = ", "),
        call. = FALSE
      )
    }
    return(data[, items, drop = FALSE])
  }
  if (is.numeric(items)) {
    if (any(items < 1) || any(items > ncol(data))) {
      stop("`items` contains out-of-range column indices.", call. = FALSE)
    }
    return(data[, items, drop = FALSE])
  }
  stop(
    "`items` must be a character vector of column names or a numeric ",
    "vector of column indices.",
    call. = FALSE
  )
}

#' Check that supplied thresholds cover the observed response categories
#'
#' `.coerce_item_params()` and `.align_items()` (person_parameters.R) handle
#' the shape and the item matching. This adds the one check the curve needs:
#' an item whose thresholds stop short of its observed top category would make
#' the information sum undefined there.
#'
#' @param thr_list Coerced list of Andrich threshold vectors.
#' @param data_mat Response matrix, aligned with `thr_list`.
#' @return Invisibly `TRUE`; otherwise stops.
#' @keywords internal
#' @noRd
.curve_check_categories <- function(thr_list, data_mat) {
  for (i in seq_along(thr_list)) {
    thr <- thr_list[[i]]
    if (!is.numeric(thr) || length(thr) < 1L || !all(is.finite(thr))) {
      stop(
        "`item_params` entry ", i, " must be a finite numeric vector of ",
        "Andrich thresholds.",
        call. = FALSE
      )
    }
    max_cat <- suppressWarnings(max(data_mat[, i], na.rm = TRUE))
    if (is.finite(max_cat) && length(thr) < max_cat) {
      stop(
        "`item_params` entry ", i, " has ", length(thr), " threshold(s) but ",
        "item ", i, " has responses up to category ", max_cat, ".",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Internal: bootstrap band
# ---------------------------------------------------------------------------

#' Bootstrap confidence band for the conditional precision curve
#'
#' Resamples respondents, recomputes the curve on the same theta grid, and
#' takes a pointwise HDCI. When `item_params` were supplied by the user the
#' thresholds are held fixed and only the latent SD is resampled.
#'
#' @return A list with `band` (a list of two-column matrices, one per
#'   statistic) and `n_ok` (successful iterations).
#' @keywords internal
#' @noRd
.curve_bootstrap <- function(
  data,
  thr_list,
  fixed_params,
  grid,
  boot_iter,
  conf_int,
  parallel,
  n_cores,
  seed,
  verbose
) {
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required for boot = TRUE.\n",
      "Install it with: install.packages(\"ggdist\")",
      call. = FALSE
    )
  }

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
      n_cores <- min(n_cores, boot_iter)
    }
  }

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max - 2L, 1L)
  }
  set.seed(seed)
  boot_seeds <- sample.int(.Machine$integer.max, boot_iter)

  boot_args <- list(
    data = data,
    thr_list = thr_list,
    fixed_params = fixed_params,
    grid = grid
  )

  results <- if (use_parallel) {
    .curve_boot_parallel(boot_iter, boot_seeds, boot_args, n_cores, verbose)
  } else {
    .curve_boot_sequential(boot_iter, boot_seeds, boot_args, verbose)
  }

  ok <- vapply(results, is.list, logical(1L))
  successful <- results[ok]
  n_ok <- length(successful)

  if (n_ok < 2L) {
    warning(
      "Fewer than 2 bootstrap iterations succeeded; band not drawn.",
      call. = FALSE
    )
    return(list(band = NULL, n_ok = n_ok))
  }

  band <- lapply(c("information", "sem", "reliability"), function(nm) {
    mat <- vapply(successful, function(x) x[[nm]], numeric(length(grid)))
    t(vapply(
      seq_len(nrow(mat)),
      function(i) {
        int <- ggdist::hdci(mat[i, ], .width = conf_int)
        c(int[1L, 1L], int[1L, 2L])
      },
      numeric(2L)
    ))
  })
  names(band) <- c("information", "sem", "reliability")

  list(band = band, n_ok = n_ok)
}

#' Single bootstrap iteration for the precision curve
#'
#' @keywords internal
#' @noRd
run_single_curve_boot <- function(seed, data_list) {
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
  idx <- sample.int(nrow(data_list$data), nrow(data_list$data), replace = TRUE)
  dat_b <- data_list$data[idx, , drop = FALSE]

  tryCatch(
    {
      thr_b <- if (isTRUE(data_list$fixed_params)) {
        data_list$thr_list
      } else {
        .fit_cml_thresholds(dat_b)
      }
      sigma_b <- .latent_sd(dat_b, thr_b)
      if (!is.finite(sigma_b) || sigma_b <= 0) {
        return("Latent SD not estimable for this resample")
      }
      .curve_stats(thr_b, data_list$grid, sigma_b)
    },
    error = function(e) as.character(conditionMessage(e))
  )
}

#' @keywords internal
#' @noRd
.curve_boot_parallel <- function(
  boot_iter,
  boot_seeds,
  boot_args,
  n_cores,
  verbose = FALSE
) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = boot_iter, style = 3)
  }

  tasks <- lapply(seq_len(boot_iter), function(i) {
    mirai::mirai(
      {
        run_single_curve_boot(seed, data_list)
      },
      seed = boot_seeds[i],
      data_list = boot_args,
      run_single_curve_boot = run_single_curve_boot,
      .curve_stats = .curve_stats,
      .test_information = .test_information,
      .latent_sd = .latent_sd,
      .fit_cml_thresholds = .fit_cml_thresholds,
      .pcm_cat_probs = .pcm_cat_probs,
      .center_thresholds = .center_thresholds,
      .estimate_prior_sd = .estimate_prior_sd,
      .logp_tables = .logp_tables,
      .grid_loglik = .grid_loglik
    )
  })

  results <- vector("list", boot_iter)
  for (i in seq_len(boot_iter)) {
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
.curve_boot_sequential <- function(
  boot_iter,
  boot_seeds,
  boot_args,
  verbose = FALSE
) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = boot_iter, style = 3)
  }

  results <- vector("list", boot_iter)
  for (i in seq_len(boot_iter)) {
    results[[i]] <- run_single_curve_boot(boot_seeds[i], boot_args)
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

#' Shared caption/footnote clause describing how the curve was built
#'
#' @keywords internal
#' @noRd
.curve_source_clause <- function(n_items, method, fixed_params) {
  paste0(
    "Test information summed over ",
    n_items,
    " items from ",
    if (fixed_params) "supplied" else "CML",
    " thresholds, respondent locations by ",
    method,
    "."
  )
}

#' Summary table of the conditional precision curve
#'
#' @keywords internal
#' @noRd
.curve_kable <- function(curve_df, n_items, n_clause, method, fixed_params) {
  # `attr()` partial-matches by default. When `benchmark` is NULL its attribute
  # is absent, and a non-exact lookup would silently return `benchmark_percent`
  # instead, printing a benchmark row that was never asked for. Every read here
  # is therefore exact.
  sigma <- attr(curve_df, "sigma", exact = TRUE)
  benchmark <- attr(curve_df, "benchmark", exact = TRUE)
  bench_range <- attr(curve_df, "benchmark_range", exact = TRUE)
  bench_pct <- attr(curve_df, "benchmark_percent", exact = TRUE)
  n_ne <- attr(curve_df, "n_not_estimable", exact = TRUE)

  best <- which.min(curve_df$sem)

  quantity <- c(
    "Latent SD (sigma)",
    "Marginal reliability (curve mean, as in RMreliability)",
    "Marginal reliability (Green/Lord, superseded)",
    "Average SEM (logits)",
    "Minimum SEM (logits)",
    "Theta at minimum SEM"
  )
  value <- c(
    sprintf("%.3f", sigma),
    sprintf("%.3f", attr(curve_df, "marginal_ratio", exact = TRUE)),
    sprintf("%.3f", attr(curve_df, "marginal_green", exact = TRUE)),
    sprintf("%.3f", attr(curve_df, "sem_average", exact = TRUE)),
    sprintf("%.3f", curve_df$sem[best]),
    sprintf("%.2f", curve_df$theta[best])
  )

  if (!is.null(benchmark)) {
    range_txt <- if (is.null(bench_range)) {
      "none"
    } else {
      paste(
        sprintf("%.2f to %.2f", bench_range$xmin, bench_range$xmax),
        collapse = "; "
      )
    }
    quantity <- c(
      quantity,
      sprintf("Theta range with reliability >= %.2f", benchmark),
      "Respondents located in that range"
    )
    value <- c(
      value,
      range_txt,
      if (is.na(bench_pct)) "not estimable" else sprintf("%.1f%%", bench_pct)
    )
  }

  if (n_ne > 0L) {
    quantity <- c(quantity, "Respondents with no estimable location")
    value <- c(value, as.character(n_ne))
  }

  knitr::kable(
    data.frame(
      quantity = quantity,
      value = value,
      stringsAsFactors = FALSE,
      row.names = NULL
    ),
    format = "pipe",
    col.names = c("Quantity", "Value"),
    caption = paste0(
      "Conditional measurement precision. ",
      .curve_source_clause(n_items, method, fixed_params),
      " Conditional reliability is the ratio form, ",
      "sigma^2 / (sigma^2 + SEM^2). ",
      n_clause,
      "."
    )
  )
}

#' Build the conditional precision figure
#'
#' @keywords internal
#' @noRd
.curve_plot <- function(
  curve_df,
  statistic,
  reference,
  sem_bar,
  marginal_ratio,
  sigma,
  benchmark,
  benchmark_range,
  benchmark_percent,
  theta_hat,
  n_items,
  n_clause,
  method,
  fixed_params,
  conf_int,
  actual_boot
) {
  y_lab <- switch(statistic,
    sem = "Conditional SEM (logits)",
    information = "Test information",
    reliability = "Conditional reliability"
  )

  fig <- data.frame(
    theta = curve_df$theta,
    value = curve_df[[statistic]],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  has_band <- paste0(statistic, "_lower") %in% names(curve_df)
  if (has_band) {
    fig$lower <- curve_df[[paste0(statistic, "_lower")]]
    fig$upper <- curve_df[[paste0(statistic, "_upper")]]
  }

  # The flat summary a single coefficient implies, on this axis. Each axis
  # gets its own natural summary: the root mean error variance for SEM and
  # information, and the mean of the reliability curve for reliability, which
  # is the `marginal_ratio` the table and attributes report.
  ref_value <- switch(statistic,
    sem = sem_bar,
    information = 1 / sem_bar^2,
    reliability = marginal_ratio
  )

  # The benchmark is always stated on the reliability scale, so convert it.
  bench_value <- if (is.null(benchmark)) {
    NULL
  } else {
    sem_crit <- sigma * sqrt((1 - benchmark) / benchmark)
    switch(statistic,
      sem = sem_crit,
      information = 1 / sem_crit^2,
      reliability = benchmark
    )
  }

  y_vals <- c(fig$value, if (has_band) c(fig$lower, fig$upper), ref_value)
  y_vals <- y_vals[is.finite(y_vals)]
  y_min <- if (statistic == "reliability") 0 else 0
  y_max <- if (statistic == "reliability") 1 else max(y_vals, na.rm = TRUE)

  p <- ggplot2::ggplot(fig, ggplot2::aes(x = .data$theta, y = .data$value))

  # Respondent distribution, drawn first so the curve sits on top of it.
  if (length(theta_hat) > 1L) {
    dens <- stats::density(
      theta_hat,
      from = min(fig$theta),
      to = max(fig$theta),
      n = nrow(fig)
    )
    dens_df <- data.frame(
      theta = dens$x,
      value = y_min + (dens$y / max(dens$y)) * 0.22 * (y_max - y_min),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    p <- p +
      ggplot2::geom_area(
        data = dens_df,
        fill = "grey70",
        colour = NA,
        alpha = 0.45
      )
  }

  if (!is.null(benchmark) && !is.null(benchmark_range)) {
    p <- p +
      ggplot2::geom_rect(
        data = benchmark_range,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax),
        ymin = -Inf,
        ymax = Inf,
        inherit.aes = FALSE,
        fill = "#1B9E77",
        alpha = 0.12
      )
  }

  if (has_band) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        fill = "grey40",
        alpha = 0.25
      )
  }

  if (reference == "marginal") {
    p <- p +
      ggplot2::geom_hline(
        yintercept = ref_value,
        linetype = "dashed",
        colour = "grey30",
        linewidth = 0.4
      )
  }

  if (!is.null(bench_value)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = bench_value,
        linetype = "dotted",
        colour = "#1B9E77",
        linewidth = 0.5
      )
  }

  p +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
    ggplot2::labs(
      x = "Person location (logits)",
      y = y_lab,
      caption = er2_caption(.curve_caption(
        statistic = statistic,
        reference = reference,
        ref_value = ref_value,
        sigma = sigma,
        benchmark = benchmark,
        benchmark_percent = benchmark_percent,
        has_density = length(theta_hat) > 1L,
        has_band = has_band,
        conf_int = conf_int,
        actual_boot = actual_boot,
        n_items = n_items,
        n_clause = n_clause,
        method = method,
        fixed_params = fixed_params
      ))
    ) +
    ggplot2::theme_bw() +
    er2_axis_margins() +
    er2_plot_caption()
}

#' Caption text for the conditional precision figure
#'
#' Leads with the quantity actually on the y-axis, and names the reference
#' line in that quantity's own units. Both were previously generic, which made
#' an SEM plot read as though information were plotted and made the average
#' SEM in logits look like a reliability coefficient.
#'
#' @keywords internal
#' @noRd
.curve_caption <- function(
  statistic,
  reference,
  ref_value,
  sigma,
  benchmark,
  benchmark_percent,
  has_density,
  has_band,
  conf_int,
  actual_boot,
  n_items,
  n_clause,
  method,
  fixed_params
) {
  source_txt <- paste0(
    "summed over ", n_items, " items using ",
    if (fixed_params) "supplied" else "CML", " thresholds"
  )

  parts <- switch(statistic,
    sem = paste0(
      "Conditional standard error of measurement, 1 / sqrt(I(theta)), from ",
      "test information ", source_txt, "."
    ),
    information = paste0(
      "Test information I(theta), ", source_txt, "."
    ),
    reliability = paste0(
      "Conditional reliability, sigma^2 / (sigma^2 + SEM(theta)^2) with ",
      sprintf("sigma = %.2f", sigma), ", from test information ", source_txt, "."
    )
  )

  parts <- c(
    parts,
    "The curve is a property of the items, not of the sample."
  )

  # The estimator only touches the respondent-derived annotations, so it is
  # named only when one of them is present.
  if (has_density || !is.null(benchmark)) {
    parts <- c(parts, paste0("Respondent locations by ", method, "."))
  }

  if (has_density) {
    parts <- c(parts, "Shaded distribution shows where respondents fall.")
  }

  if (reference == "marginal") {
    ref_txt <- switch(statistic,
      sem = sprintf("the average SEM (%.2f logits)", ref_value),
      information = sprintf("the average test information (%.2f)", ref_value),
      reliability = sprintf("the marginal reliability (%.3f)", ref_value)
    )
    parts <- c(
      parts,
      paste0(
        "Dashed line is ", ref_txt,
        ", the flat summary a single coefficient reports."
      )
    )
  }

  if (!is.null(benchmark)) {
    pct <- if (is.na(benchmark_percent)) {
      "an unknown share of"
    } else {
      sprintf("%.0f%% of", benchmark_percent)
    }
    sem_crit <- sigma * sqrt((1 - benchmark) / benchmark)
    crit_txt <- switch(statistic,
      reliability = sprintf("reliability >= %.2f", benchmark),
      sem = sprintf("SEM <= %.2f logits (reliability %.2f)", sem_crit, benchmark),
      information = sprintf(
        "information >= %.2f (reliability %.2f)",
        1 / sem_crit^2,
        benchmark
      )
    )
    parts <- c(
      parts,
      sprintf(
        "Dotted line and shaded band mark %s, reached by %s respondents.",
        crit_txt,
        pct
      )
    )
  }

  if (has_band && !is.na(actual_boot)) {
    parts <- c(
      parts,
      sprintf(
        "Band is a pointwise %.0f%% HDCI over %d bootstrap resamples.",
        100 * conf_int,
        actual_boot
      )
    )
  }

  if (n_items < 10L) {
    parts <- c(
      parts,
      paste(
        "Information-based standard errors are asymptotic and only",
        "approximate for a scale this short."
      )
    )
  }

  paste(c(parts, paste0(n_clause, ".")), collapse = " ")
}
