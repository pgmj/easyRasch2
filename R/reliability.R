#' Cronbach's alpha
#'
#' Closed-form alpha from item-score and total-score variances. Complete
#' cases only (rows with any `NA` are dropped).
#'
#' @param data A data.frame or matrix of item responses.
#'
#' @return Numeric scalar; `NA_real_` if fewer than two items, fewer than
#'   two complete cases, or zero total-score variance.
#'
#' @keywords internal
#' @noRd
cronbach_alpha <- function(data) {
  d <- stats::na.omit(as.data.frame(data))
  k <- ncol(d)
  if (k < 2L || nrow(d) < 2L) {
    return(NA_real_)
  }
  total_var <- stats::var(rowSums(d))
  if (!is.finite(total_var) || total_var == 0) {
    return(NA_real_)
  }
  item_vars <- vapply(d, stats::var, numeric(1L))
  (k / (k - 1L)) * (1 - sum(item_vars) / total_var)
}

# ---------------------------------------------------------------------------

#' Relative measurement uncertainty (RMU)
#'
#' Reliability estimated from a matrix of posterior or plausible-value
#' draws (Bignardi, Kievit & Bürkner, 2026). The columns of `input_draws`
#' are split at random into two halves, and the reliability is the Pearson
#' correlation across persons between paired columns of the two halves,
#' summarised over the pairs as a mean with an HDCI. [RMreliability()] uses
#' it with plausible values from \pkg{mirt}.
#'
#' Adapted (with permission, GPL-2/3) from `rmu()`, formerly `reliability()`,
#' in the 'gbtoolbox' package, \url{https://github.com/giac01/gbtoolbox}.
#'
#' @param input_draws Numeric matrix or data.frame of draws, one row per
#'   person and one column per draw. At least two columns, ideally many.
#' @param level Numeric in (0, 1). Width of the HDCI returned. Default `0.95`.
#' @param verbose Logical. Print summary information about the input. Default
#'   `FALSE`.
#'
#' @return A 1-row data.frame with columns `rmu_estimate`, `hdci_lowerbound`,
#'   `hdci_upperbound`, plus the `.width`/`.point`/`.interval` metadata
#'   columns added by `ggdist::mean_hdci()`.
#'
#' @details
#' A column pair in which either column has zero variance gets a
#' correlation of 0, without a message, since the correlation is undefined
#' there.
#'
#' The function takes no `seed`. The split of the columns is random and
#' follows the session's random number stream (see
#' [easyRasch2-reproducibility]). Needs the \pkg{ggdist} package (in
#' Suggests).
#'
#' @references
#' Bignardi, G., Kievit, R., & Bürkner, P. C. (2026). A general method for
#' estimating reliability using Bayesian Measurement Uncertainty. *PsyArXiv*.
#' \url{https://osf.io/preprints/psyarxiv/h54k8}
#'
#' @seealso [RMreliability()]
#'
#' @export
RMUreliability <- function(input_draws, level = 0.95, verbose = FALSE) {
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required for RMUreliability().\n",
      "Install it with: install.packages(\"ggdist\")",
      call. = FALSE
    )
  }
  if (!is.numeric(level) || level <= 0 || level >= 1) {
    stop(
      "`level` must be a numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  input_draws <- as.matrix(input_draws)
  if (ncol(input_draws) < 2L) {
    stop("`input_draws` must have at least 2 columns.", call. = FALSE)
  }

  sds <- apply(input_draws, 2L, stats::sd, na.rm = TRUE)
  zero_sd_cols <- which(sds == 0)
  if (length(zero_sd_cols) > 0L) {
    warning(
      sprintf(
        "Found %d column(s) with zero standard deviation (columns: %s).",
        length(zero_sd_cols),
        paste(zero_sd_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  na_count <- sum(is.na(input_draws))
  if (na_count > 0L) {
    warning(
      sprintf("Found %d NA value(s) in `input_draws`.", na_count),
      call. = FALSE
    )
  }

  if (verbose) {
    message("Subjects: ", nrow(input_draws), "; draws: ", ncol(input_draws))
  }

  col_select <- sample.int(ncol(input_draws), replace = FALSE)
  half <- floor(length(col_select) / 2)
  cols_a <- col_select[seq_len(half)]
  cols_b <- col_select[(half + 1L):(2L * half)]

  draws_a <- input_draws[, cols_a, drop = FALSE]
  draws_b <- input_draws[, cols_b, drop = FALSE]

  rel_post <- vapply(
    seq_len(ncol(draws_a)),
    function(i) {
      x <- draws_a[, i]
      y <- draws_b[, i]
      if (
        stats::var(x, na.rm = TRUE) == 0 ||
          stats::var(y, na.rm = TRUE) == 0
      ) {
        return(0)
      }
      stats::cor(x, y, method = "pearson", use = "complete.obs")
    },
    numeric(1L)
  )

  hdci <- ggdist::mean_hdci(rel_post, .width = level)
  colnames(hdci)[1:3] <- c("rmu_estimate", "hdci_lowerbound", "hdci_upperbound")
  hdci
}

# ---------------------------------------------------------------------------

#' Reliability coefficients for a Rasch model
#'
#' Computes four reliability coefficients for a Rasch or partial credit
#' model: Cronbach's alpha, the person separation index (PSI, the
#' WLE-based separation reliability), the marginal reliability (the mean of
#' the conditional reliability curve over the estimated latent
#' distribution), and the relative measurement uncertainty (RMU,
#' [RMUreliability()]) from plausible values drawn with `mirt::fscores()`.
#'
#' With `boot = TRUE`, alpha, PSI and the marginal reliability get
#' confidence intervals from a nonparametric bootstrap that resamples
#' respondents and recomputes the coefficients without refitting a `mirt`
#' model. The RMU interval is the HDCI of the correlations across
#' plausible-value draws, averaged over `rmu_iter` random splits of the
#' draws.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. See Details for how missing values are handled.
#' @param conf_int Numeric in (0, 1). HDCI width of the bootstrap intervals
#'   and of RMU. Default `0.95`.
#' @param draws Integer. Number of plausible values drawn from the `mirt`
#'   model for RMU. Default `1000`. More draws give a more stable RMU, at a
#'   cost that grows about linearly.
#' @param rmu_iter Integer. How many times [RMUreliability()] is repeated on
#'   the same draws, each with a new random split of the columns. The
#'   estimates are averaged over the repetitions. Default `50`.
#' @param estim Character. Estimator used by `mirt::fscores()` as the
#'   starting point for the plausible values: `"WLE"` (default), `"EAP"`,
#'   `"MAP"` or `"ML"`. The plausible values themselves are drawn by
#'   Metropolis-Hastings. PSI and the marginal reliability do not use it.
#' @param boot Logical. Bootstrap confidence intervals for alpha, PSI and
#'   the marginal reliability. Default `FALSE`.
#' @param boot_iter Integer. Number of bootstrap samples with
#'   `boot = TRUE`. Default `200`.
#' @param parallel Logical. Run the bootstrap in parallel via \pkg{mirai}
#'   when it is installed. Default `TRUE`.
#' @param n_cores Integer or `NULL`. Number of parallel workers. When
#'   `NULL`, `getOption("mc.cores")` is used. If that is unset too, the
#'   bootstrap runs sequentially.
#' @param seed Integer or `NULL`. Random seed for reproducibility. See
#'   [easyRasch2-reproducibility] for what this guarantees and how it
#'   interacts with `parallel`.
#' @param verbose Logical. Print progress messages and a progress bar for
#'   the bootstrap. Default `FALSE`.
#' @param theta_range Numeric vector of length 2. Theta limits passed to
#'   `mirt::fscores()`. Default `c(-10, 10)`.
#' @param output Character. `"kable"` (default) for a `knitr::kable()`
#'   table, or `"dataframe"`.
#'
#' @return
#' * `output = "kable"`: a `knitr_kable` table with one row per coefficient.
#' * `output = "dataframe"`: a data.frame with the columns `metric`
#'   (`"Cronbach's alpha"`, `"PSI"`, `"Marginal (curve mean)"` and
#'   `"RMU (<estim>)"`), `estimate`, `lower`, `upper` and `notes`.
#'
#' @details
#' \strong{Data.} Respondents with no responses are dropped. Cronbach's
#' alpha uses complete cases. PSI and the marginal reliability use every
#' respondent's answered items, and PSI leaves out respondents with the
#' lowest or highest possible score (see below). At least one complete case
#' is needed.
#'
#' \strong{Marginal reliability} is the mean of the conditional reliability
#' curve over the latent density, \eqn{\int \sigma^2/(\sigma^2 +
#' 1/I(\theta))\, g(\theta)\,d\theta}. The test information
#' \eqn{I(\theta)} is summed from the CML item parameters, and \eqn{g} is the
#' normal latent density \eqn{N(\mu, \sigma^2)} with both parameters
#' estimated by marginal maximum likelihood with the items fixed.
#' Integrating over the estimated distribution, rather than the
#' \eqn{N(0,1)} that `mirt::marginal_rxx()` assumes, keeps it correct on the
#' Rasch logit scale, where \eqn{\sigma} is often well above 1 and an
#' \eqn{N(0,1)} assumption underestimates reliability.
#'
#' Up to version 1.3.0 \eqn{\mu} was held at 0. Centring the item thresholds
#' fixes the item mean, not where the respondents are, so a normal density
#' pinned to 0 could reach an off-target sample only by widening.
#' \eqn{\sigma} then absorbed the mistargeting, and marginal reliability rose
#' as targeting got worse. See the note under Value in
#' [RMreliabilityCurve()] for the size of the effect.
#'
#' The marginal reliability is the model-based counterpart of the
#' sample-based PSI. PSI divides by the observed spread of the WLE
#' estimates, the marginal reliability uses the fitted latent density. A
#' large gap between them flags an off-target or non-normal sample, and
#' since version 1.3.1 both fall when targeting gets worse.
#'
#' \strong{How it relates to published coefficients.} The reported value is
#' none of the coefficients named in the literature.
#'
#' *It is the mean of a curve, not a ratio of averages.* The classical
#' latent-scale reliability is one variance ratio,
#' \eqn{\sigma^2/(\sigma^2 + \overline{SEM^2})} (Milanzi et al., 2015,
#' section 3.2, with a constant error variance). Here the ratio is averaged
#' over the latent density instead. Since \eqn{\sigma^2/(\sigma^2 + x)} is
#' convex in \eqn{x}, the mean of the curve is the larger of the two
#' (Jensen's inequality). The choice was made on accuracy: averaging the
#' curve tracked the exact expected-sum-score reliability more closely than
#' the ratio of averages in an internal replication (mean absolute error
#' 0.015 against 0.018 for binary data, 0.010 against 0.015 for
#' polytomous).
#'
#' *It is a latent-scale quantity checked against a manifest criterion.*
#' Milanzi et al. (2015) argue that interest usually lies in the
#' reliability of observed scores, and that latent coefficients run
#' consistently higher. This one was therefore validated against the exact
#' expected-sum-score reliability, \eqn{Var(\mu)/(Var(\mu) +
#' Var(\epsilon))}, and tracked it to within about 0.02 in the conditions
#' tested. It remains a latent-scale coefficient that approximates a
#' manifest one.
#'
#' *Milanzi et al. recommend something this package does not implement.*
#' They favour Taylor-series manifest reliability measures, which
#' approximate \eqn{Var(\mu)} and \eqn{Var(\epsilon)} on the observed-score
#' scale directly. In the same internal replication those were the most
#' accurate (0.012). They are cited here for their criticism of the
#' subtractive coefficient, which this package acted on, not for their
#' remedy.
#'
#' Up to version 1.2.0 this row used the subtractive coefficient
#' \eqn{1 - \overline{1/I(\theta)}/\sigma^2} (Green et al., 1984, Milanzi et
#' al.'s equation 12, attributed there to Lord, 1980). It can fall below
#' zero when the average error variance exceeds the trait variance, and
#' much of the gap between PSI and that coefficient came from the different
#' formulas rather than the sample. It is still returned as the
#' `marginal_green` attribute of [RMreliabilityCurve()] for comparison.
#'
#' \strong{PSI} is the WLE-based separation reliability,
#' \eqn{1 - \overline{SEM^2} / \mathrm{Var}(\hat\theta)}, from the CML item
#' thresholds (`psychotools`) and the WLE person locations and their SEMs.
#' Respondents with the lowest or highest possible raw score are left out,
#' since their boundary estimates would inflate the person variance and
#' overstate reliability. Earlier versions used `eRm::SepRel()` with MLE,
#' and the values can differ, most for scales with many extreme scores,
#' such as dichotomous items.
#'
#' \strong{RMU} follows Bignardi, Kievit and Bürkner (2026), here with
#' plausible values from `mirt` rather than fully Bayesian posterior draws
#' (see Mislevy, 1991, for plausible values).
#'
#' Bootstrap samples that cannot be fitted are dropped, with a message
#' giving how many and why. The usual reason is a sparse lowest or middle
#' response category that is missing from the sample.
#'
#' @references
#' Bignardi, G., Kievit, R., & Bürkner, P. C. (2026). A general method for
#' estimating reliability using Bayesian Measurement Uncertainty. *PsyArXiv*.
#' \url{https://osf.io/preprints/psyarxiv/h54k8}
#'
#' Green, B. F., Bock, R. D., Humphreys, L. G., Linn, R. L., & Reckase, M. D.
#' (1984). Technical Guidelines for Assessing Computerized Adaptive Tests.
#' *Journal of Educational Measurement, 21*(4), 347–360.
#' \doi{10.1111/j.1745-3984.1984.tb01039.x}
#'
#' Milanzi, E., Molenberghs, G., Alonso, A., Verbeke, G., & De Boeck, P.
#' (2015). Reliability measures in item response theory: Manifest versus latent
#' correlation functions. *British Journal of Mathematical and Statistical
#' Psychology, 68*(1), 43-64. \doi{10.1111/bmsp.12033}
#'
#' Mislevy, R. J. (1991). Randomization-Based Inference about Latent Variables
#' from Complex Samples. *Psychometrika, 56*(2), 177-196.
#' \doi{10.1007/BF02294457}
#'
#' @seealso [RMreliabilityCurve()] for precision along the scale,
#'   [RMUreliability()], [RMpersonParameters()].
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggdist", quietly = TRUE) &&
#'     requireNamespace("eRm", quietly = TRUE)) {
#'   set.seed(1)
#'   RMreliability(eRm::raschdat1[, 1:20], draws = 1000)
#'
#'   # Bootstrap intervals for alpha, PSI and the marginal reliability
#'   # (use 200 or more bootstrap samples in real analyses)
#'   RMreliability(eRm::raschdat1[, 1:20], draws = 1000,
#'                 boot = TRUE, boot_iter = 25, parallel = FALSE, seed = 42)
#' }
#' }
RMreliability <- function(
  data,
  conf_int = 0.95,
  draws = 1000,
  rmu_iter = 50,
  estim = "WLE",
  boot = FALSE,
  boot_iter = 200,
  parallel = TRUE,
  n_cores = NULL,
  seed = NULL,
  verbose = FALSE,
  theta_range = c(-10, 10),
  output = "kable"
) {
  output <- match.arg(output, c("kable", "dataframe"))
  estim <- match.arg(estim, c("WLE", "EAP", "MAP", "ML"))

  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required for RMreliability().\n",
      "Install it with: install.packages(\"ggdist\")",
      call. = FALSE
    )
  }

  validate_response_data(data)

  if (nrow(stats::na.omit(data)) == 0L) {
    stop("No complete cases in `data`.", call. = FALSE)
  }

  # Respondents with no responses at all contribute nothing and break the CML
  # fit behind the PSI (psychotools errors on all-NA rows); drop them, keeping
  # the raw total for the caption.
  n_total <- nrow(as.data.frame(data))
  data <- .drop_empty_respondents(data)

  data_mat <- as.matrix(data)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L
  n_persons <- nrow(data)
  n_items <- ncol(data)

  # rgl workaround for any iarm-related fallout
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # A NULL `seed` is resolved here rather than left to inherit the session
  # stream. The plausible-value sampler below advances that stream
  # nondeterministically (see the comment at the RMU block), so without a
  # seed of its own everything after it would differ between two calls made
  # under the same set.seed().
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max - 2L, 1L)
  }

  # --- Full-sample fits ------------------------------------------------------
  # .mml_fit() runs the shared category check before mirt, which renumbers
  # items with an empty lowest or interior category (the CML paths check in
  # .fit_cml_thresholds()).
  mirt_fit <- .mml_fit(data, accelerate = "squarem")
  # --- Marginal reliability point estimate (native CML test info) ------------
  marg_rel <- .marginal_rxx(data)

  # --- PSI point estimate (WLE-based separation reliability) -----------------
  psi <- .wle_psi(data)

  # --- Cronbach's alpha point estimate ---------------------------------------
  alpha <- cronbach_alpha(data)

  # --- Plausible values + RMU ------------------------------------------------
  set.seed(seed)
  pvs <- mirt::fscores(
    mirt_fit,
    method = estim,
    theta_lim = theta_range,
    plausible.draws = draws,
    plausible.type = "MH",
    verbose = FALSE
  )
  rmu_input <- do.call(cbind, lapply(pvs, as.numeric))

  # mirt's MH plausible-value sampler leaves the R RNG in a nondeterministic
  # state even under set.seed() (the draws themselves are reproducible, the
  # stream advancement is not), so the RMU column splits below would differ
  # between identical calls. Re-seed to make the whole result reproducible;
  # + 2L keeps the stream distinct from the bootstrap's seed + 1L.
  set.seed(seed + 2L)

  rmu_iter_results <- do.call(
    rbind,
    lapply(seq_len(rmu_iter), function(i) {
      RMUreliability(rmu_input, level = conf_int)[, 1:3, drop = FALSE]
    })
  )
  rmu_summary <- list(
    estimate = mean(rmu_iter_results$rmu_estimate),
    lower = mean(rmu_iter_results$hdci_lowerbound),
    upper = mean(rmu_iter_results$hdci_upperbound)
  )

  # --- Bootstrap (optional) --------------------------------------------------
  alpha_lower <- NA_real_
  alpha_upper <- NA_real_
  psi_lower <- NA_real_
  psi_upper <- NA_real_
  marg_lower <- NA_real_
  marg_upper <- NA_real_
  actual_boot <- NA_integer_

  if (isTRUE(boot)) {
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

    set.seed(seed + 1L)
    boot_seeds <- sample.int(.Machine$integer.max, boot_iter)

    boot_args <- list(
      data = data,
      is_polytomous = is_polytomous,
      estim = estim,
      theta_range = theta_range
    )

    if (use_parallel) {
      boot_results <- run_reliability_boot_parallel(
        boot_iter,
        boot_seeds,
        boot_args,
        n_cores,
        verbose
      )
    } else {
      boot_results <- run_reliability_boot_sequential(
        boot_iter,
        boot_seeds,
        boot_args,
        verbose
      )
    }

    ok <- vapply(boot_results, is.list, logical(1L))
    successful <- boot_results[ok]
    actual_boot <- length(successful)
    .notify_dropped_resamples(boot_results, what = "interval")

    if (actual_boot < 2L) {
      warning(
        "Fewer than 2 bootstrap iterations succeeded; CI not reported.",
        call. = FALSE
      )
    } else {
      alpha_vec <- vapply(successful, function(x) x$alpha, numeric(1L))
      psi_vec <- vapply(successful, function(x) x$psi, numeric(1L))
      marg_vec <- vapply(successful, function(x) x$marginal, numeric(1L))

      alpha_int <- ggdist::hdci(alpha_vec, .width = conf_int)
      psi_int <- ggdist::hdci(psi_vec, .width = conf_int)
      marg_int <- ggdist::hdci(marg_vec, .width = conf_int)

      alpha_lower <- alpha_int[1L, 1L]
      alpha_upper <- alpha_int[1L, 2L]
      psi_lower <- psi_int[1L, 1L]
      psi_upper <- psi_int[1L, 2L]
      marg_lower <- marg_int[1L, 1L]
      marg_upper <- marg_int[1L, 2L]
    }
  }

  # --- Assemble output table -------------------------------------------------
  conf_pct <- round(conf_int * 100, 1)
  boot_note <- if (isTRUE(boot)) {
    if (!is.na(actual_boot) && actual_boot >= 2L) {
      paste0(actual_boot, " bootstrap resamples")
    } else {
      "bootstrap failed"
    }
  } else {
    "no bootstrap"
  }
  rmu_note <- paste0(draws, " PVs, ", rmu_iter, " RMU iterations")

  result_df <- data.frame(
    metric = c(
      "Cronbach's alpha",
      "PSI",
      "Marginal (curve mean)",
      paste0("RMU (", estim, ")")
    ),
    estimate = c(alpha, psi, marg_rel, rmu_summary$estimate),
    lower = c(alpha_lower, psi_lower, marg_lower, rmu_summary$lower),
    upper = c(alpha_upper, psi_upper, marg_upper, rmu_summary$upper),
    notes = c(boot_note, boot_note, boot_note, rmu_note),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (output == "dataframe") {
    return(result_df)
  }

  # Kable display rounding (the dataframe output above stays unrounded)
  result_df <- .round_display(result_df, c(estimate = 3, lower = 3, upper = 3))

  knitr::kable(
    result_df,
    format = "pipe",
    col.names = c(
      "Metric",
      "Estimate",
      paste0("Lower (", conf_pct, "% HDCI)"),
      paste0("Upper (", conf_pct, "% HDCI)"),
      "Notes"
    ),
    caption = paste0(
      "Reliability for ",
      n_items,
      " items, ",
      .n_caption(n_persons, n_total),
      ". PSI is the WLE-based separation reliability and excludes min/max ",
      "scoring respondents. Marginal is the latent-density-weighted mean of ",
      "the conditional reliability curve, sigma^2 / (sigma^2 + SEM(theta)^2)."
    )
  )
}

# ---------------------------------------------------------------------------
# Internal: single bootstrap iteration
# ---------------------------------------------------------------------------

#' WLE-based Person Separation Index (separation reliability)
#'
#' PSI = `1 - mean(SEM^2) / Var(theta)` on Warm's WLE person locations and their
#' analytic SEMs (CML item thresholds via `psychotools`). Extreme (minimum and
#' maximum) scorers are excluded: their boundary WLE estimates would inflate the
#' person variance and overstate reliability, so the standard
#' separation-reliability convention drops them.
#'
#' @param data A data.frame or matrix of item responses (items from 0).
#' @param thr_list Optional pre-fitted CML thresholds (else fitted here).
#' @return The PSI (numeric), or `NA` if fewer than two non-extreme persons.
#' @keywords internal
#' @noRd
.wle_psi <- function(data, thr_list = NULL) {
  data_mat <- as.matrix(data)
  if (is.null(thr_list)) {
    thr_list <- .fit_cml_thresholds(data_mat)
  }
  est <- .estimate_thetas(data_mat, thr_list, method = "WLE")
  rs <- rowSums(data_mat, na.rm = TRUE)
  max_p <- vapply(
    seq_len(nrow(data_mat)),
    function(i) {
      ans <- !is.na(data_mat[i, ])
      sum(vapply(thr_list[ans], length, integer(1L)))
    },
    integer(1L)
  )
  keep <- is.finite(est$theta) & is.finite(est$sem) & rs > 0L & rs < max_p
  if (sum(keep) < 2L) {
    return(NA_real_)
  }
  1 - mean(est$sem[keep]^2) / stats::var(est$theta[keep])
}

#' Native marginal reliability (CML test information over an assumed normal)
#'
#' The latent-density-weighted mean of the conditional reliability curve,
#' \eqn{\int \sigma^2 / (\sigma^2 + 1/I(\theta))\, g(\theta)\, d\theta}, with the
#' test information \eqn{I(\theta) = \sum_i \mathrm{Var}_i(\mathrm{score}\mid\theta)}
#' summed from the CML item parameters and \eqn{g(\theta) = N(0, \sigma^2)} whose
#' SD is estimated by marginal maximum likelihood. Integrating over the
#' *estimated* latent variance (rather than the \eqn{N(0,1)} assumed by
#' `mirt::marginal_rxx()`) keeps it correct on the Rasch logit scale where
#' \eqn{\sigma \neq 1}.
#'
#' This is the bounded ratio form. It replaced Green's subtractive
#' \eqn{1 - \overline{1/I(\theta)}/\sigma^2} in the development version after the
#' replication in `dev/milanzi_check.R`; see `dev/TODO-reliability-form.md` for
#' the evidence and the decision. Being bounded in (0, 1) it needs no floor,
#' where the subtractive form did.
#'
#' @param data Response matrix/data.frame (items from 0).
#' @param thr_list Optional pre-fitted CML thresholds.
#' @param n_nodes Number of quadrature nodes.
#' @return Marginal reliability (numeric), or `NA` if the latent distribution
#'   is not estimable.
#' @keywords internal
#' @noRd
.marginal_rxx <- function(data, thr_list = NULL, n_nodes = 161L) {
  data_mat <- as.matrix(data)
  if (is.null(thr_list)) {
    thr_list <- .fit_cml_thresholds(data_mat)
  }
  lat <- .latent_moments(data_mat, thr_list)
  if (!is.finite(lat$sd) || lat$sd <= 0) {
    return(NA_real_)
  }
  .marginal_summaries(thr_list, lat$sd, mu = lat$mean, n_nodes = n_nodes)$ratio
}

#' Run a single reliability bootstrap iteration
#'
#' Resamples respondents with replacement and reads Cronbach's alpha
#' (closed-form), the WLE-based PSI (`.wle_psi()`), and the native marginal
#' reliability (`.marginal_rxx()`) off the resample -- no model fit by `mirt`.
#'
#' @keywords internal
#' @noRd
run_single_reliability_boot <- function(seed, data_list) {
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
      # All three indices are computed natively (CML/WLE) -- no mirt fit.
      thr_b <- .fit_cml_thresholds(dat_b)
      alpha_b <- cronbach_alpha(dat_b) # closed-form
      psi_b <- .wle_psi(dat_b, thr_list = thr_b) # WLE separation
      if (!is.finite(psi_b)) {
        return("PSI not estimable for this resample")
      }
      marg_b <- .marginal_rxx(dat_b, thr_list = thr_b) # CML marginal

      list(alpha = alpha_b, psi = psi_b, marginal = marg_b)
    },
    # Item names travel as an attribute (it survives mirai serialization)
    # so .notify_dropped_resamples() need not parse the message.
    error = function(e) structure(conditionMessage(e), items = e$items)
  )
}

# ---------------------------------------------------------------------------
# Internal: parallel runner (mirai)
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
run_reliability_boot_parallel <- function(
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
        run_single_reliability_boot(seed, data_list)
      },
      seed = boot_seeds[i],
      data_list = boot_args,
      run_single_reliability_boot = run_single_reliability_boot,
      cronbach_alpha = cronbach_alpha,
      # Native engine helpers used by .wle_psi() / .marginal_rxx() in the daemon.
      .wle_psi = .wle_psi,
      .marginal_rxx = .marginal_rxx,
      .test_information = .test_information,
      .latent_moments = .latent_moments,
      .marginal_summaries = .marginal_summaries,
      .fit_cml_thresholds = .fit_cml_thresholds,
      .cml_threshold_summary = .cml_threshold_summary,
      .cml_fit = .cml_fit,
      .check_response_categories = .check_response_categories,
      .estimate_thetas = .estimate_thetas,
      .theta_wle = .theta_wle,
      .pcm_cat_probs = .pcm_cat_probs,
      .center_thresholds = .center_thresholds,
      .estimate_prior_sd = .estimate_prior_sd,
      .estimate_prior_moments = .estimate_prior_moments,
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

# ---------------------------------------------------------------------------
# Internal: sequential runner
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
run_reliability_boot_sequential <- function(
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
    results[[i]] <- run_single_reliability_boot(boot_seeds[i], boot_args)
    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}
