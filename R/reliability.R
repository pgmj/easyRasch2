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

#' Relative Measurement Uncertainty (RMU)
#'
#' Bayesian-style reliability estimate (Bignardi, Kievit & Bürkner, 2025)
#' computed from a matrix of posterior or plausible-value draws. The columns of
#' `input_draws` are split at random into two halves; reliability is the
#' Pearson correlation across persons of paired columns from the two halves,
#' summarised across pairs as a posterior mean with HDCI.
#'
#' Adapted (with permission, GPL-2/3) from
#' \url{https://github.com/giac01/gbtoolbox/blob/main/R/reliability.R}.
#'
#' @param input_draws Numeric matrix or data.frame of draws. Rows are
#'   subjects; columns are draws. Must have at least two columns; ideally many.
#' @param level Numeric in (0, 1). Width of the HDCI returned. Default `0.95`.
#' @param verbose Logical. Print summary information about the input. Default
#'   `FALSE`.
#'
#' @return A 1-row data.frame with columns `rmu_estimate`, `hdci_lowerbound`,
#'   `hdci_upperbound`, plus the `.width`/`.point`/`.interval` metadata
#'   columns added by `ggdist::mean_hdci()`.
#'
#' @details
#' The function silently returns 0 for any column pair where either side has
#' zero variance (the correlation is undefined there).
#'
#' Requires the `ggdist` package (Suggests).
#'
#' @references
#' Bignardi, G., Kievit, R., & Bürkner, P. C. (2025). A general method for
#' estimating reliability using Bayesian Measurement Uncertainty. *PsyArXiv*.
#' \doi{10.31234/osf.io/h54k8_v1}
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

#' Reliability metrics for a Rasch model
#'
#' Computes three reliability indices for a Rasch / partial credit model:
#' the Person Separation Index (PSI) -- WLE-based separation reliability -- the
#' marginal reliability (native, CML test information integrated over the
#' estimated normal latent density), and Relative Measurement Uncertainty (RMU)
#' via [RMUreliability()] applied to plausible values from `mirt::fscores()`.
#'
#' Confidence intervals for **PSI** and **Marginal** reliability are obtained
#' by non-parametric bootstrap (resampling respondents; all three indices are
#' recomputed natively per resample, no model is refitted by `mirt`). The RMU
#' interval is the HDCI of correlations across plausible-value draws, averaged
#' over `rmu_iter` random splits of the draws.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers).
#' @param conf_int Numeric in (0, 1). HDCI width for both bootstrap CIs and
#'   RMU. Default `0.95`.
#' @param draws Integer. Number of plausible-value draws drawn from the mirt
#'   model for the RMU calculation. Default `1000`. More gives a more stable
#'   RMU; computational cost is mostly linear.
#' @param rmu_iter Integer. Number of times [RMUreliability()] is repeated on
#'   the same set of plausible-value draws (each repetition uses a fresh
#'   random column split). Estimates are averaged across repetitions to
#'   stabilise against split-induced variability. Default `50`.
#' @param estim Character. Theta estimator used by `mirt::fscores()` for the
#'   RMU plausible-value seed. One of `"WLE"` (default), `"EAP"`, `"MAP"`,
#'   `"ML"`. Plausible draws themselves are produced by Metropolis-Hastings.
#'   (PSI and marginal reliability are computed natively and do not use this.)
#' @param boot Logical. If `TRUE`, run a non-parametric bootstrap to obtain
#'   CIs for PSI and Marginal reliability. Default `FALSE`.
#' @param boot_iter Integer. Number of bootstrap iterations when
#'   `boot = TRUE`. Default `200`.
#' @param parallel Logical. Use parallel processing via `mirai` for the
#'   bootstrap if available. Default `TRUE`.
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first; if neither is set,
#'   bootstrapping falls back to sequential.
#' @param seed Integer or `NULL`. Master random seed for reproducibility.
#' @param verbose Logical. Print progress messages and a progress bar for
#'   the bootstrap. Default `FALSE`.
#' @param theta_range Numeric length-2 vector. Theta limits passed to
#'   `mirt::fscores()`. Default `c(-10, 10)`.
#' @param output Character. `"kable"` (default) for a formatted
#'   `knitr::kable()` table, or `"dataframe"` for the underlying data.frame.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object with one row per metric.
#' * If `output = "dataframe"`: a data.frame with columns `metric`,
#'   `estimate`, `lower`, `upper`, `notes`.
#'
#' @details
#' Marginal reliability is the native Green (1984) coefficient,
#' \eqn{1 - \overline{1/I(\theta)}/\sigma^2}, where the test information
#' \eqn{I(\theta)} is summed from the CML item parameters and the average error
#' variance is taken over the estimated normal latent density \eqn{N(0,
#' \sigma^2)} (\eqn{\sigma} from marginal ML). Unlike `mirt::marginal_rxx()`,
#' which assumes \eqn{N(0,1)}, this integrates over the estimated latent
#' variance, so it is correct on the Rasch logit scale (where \eqn{\sigma}
#' is typically well above 1, and the \eqn{N(0,1)} assumption
#' underestimates reliability). It is the model-based complement to the
#' sample-based PSI; a large gap between the two flags an off-target or
#' non-normal sample.
#'
#' PSI is the WLE-based separation reliability,
#' \eqn{1 - \overline{SEM^2} / \mathrm{Var}(\hat\theta)}, computed from CML item
#' thresholds (`psychotools`) and Warm's WLE person locations / analytic SEMs.
#' Respondents with extreme (min/max) raw scores are excluded -- their boundary
#' estimates would inflate the person variance and overstate reliability.
#' (Earlier versions used `eRm::SepRel()` with MLE; the values can differ, most
#' noticeably for scales with many extreme scorers, e.g. dichotomous items.)
#'
#' RMU is from Bignardi, Kievit, & Bürkner (2025), modified here to use mirt
#' plausible values rather than fully Bayesian posterior draws (see Mislevy,
#' 1991, for the plausible-values framework).
#'
#' Bootstrap iterations that fail to converge are silently dropped.
#'
#' @references
#' Bignardi, G., Kievit, R., & Bürkner, P. C. (2025). A general method for
#' estimating reliability using Bayesian Measurement Uncertainty. *PsyArXiv*.
#' \doi{10.31234/osf.io/h54k8_v1}
#'
#' Green, B. F., Bock, R. D., Humphreys, L. G., Linn, R. L., & Reckase, M. D.
#' (1984). Technical Guidelines for Assessing Computerized Adaptive Tests.
#' *Journal of Educational Measurement, 21*(4), 347–360.
#' \doi{10.1111/j.1745-3984.1984.tb01039.x}
#'
#' Mislevy, R. J. (1991). Randomization-Based Inference about Latent Variables
#' from Complex Samples. *Psychometrika, 56*(2), 177-196.
#' \doi{10.1007/BF02294457}
#'
#' Adams, R. J. (2005). Reliability as a measurement design effect.
#' *Studies in Educational Evaluation, 31*(2), 162-172.
#' \doi{10.1016/j.stueduc.2005.05.008}
#'
#' @seealso [RMUreliability()]
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
#'   # Bootstrap CI for PSI and Marginal
#'   # (use more bootstrap iterations, e.g. 200+, in real analyses)
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

  # --- Full-sample fits ------------------------------------------------------
  mirt_fit <- mirt::mirt(
    data = data,
    model = 1,
    itemtype = "Rasch",
    verbose = FALSE,
    accelerate = "squarem"
  )
  # --- Marginal reliability point estimate (native CML test info) ------------
  marg_rel <- .marginal_rxx(data)

  # --- PSI point estimate (WLE-based separation reliability) -----------------
  psi <- .wle_psi(data)

  # --- Cronbach's alpha point estimate ---------------------------------------
  alpha <- cronbach_alpha(data)

  # --- Plausible values + RMU ------------------------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }
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
  if (!is.null(seed)) {
    set.seed(seed + 2L)
  }

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

    if (!is.null(seed)) {
      set.seed(seed + 1L)
    }
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
      "Marginal",
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
      "scoring respondents."
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
#' Green's (1984) marginal reliability,
#' \eqn{1 - \int [1/I(\theta)]\, g(\theta)\, d\theta / \sigma^2}, with the test
#' information \eqn{I(\theta) = \sum_i \mathrm{Var}_i(\mathrm{score}\mid\theta)}
#' summed from the CML item parameters and the expectation taken over a normal
#' latent density \eqn{g(\theta) = N(0, \sigma^2)} whose SD is estimated by
#' marginal maximum likelihood. Unlike `mirt::marginal_rxx()` (which assumes
#' \eqn{N(0,1)}), this integrates over the *estimated* latent variance, so it is
#' correct on the Rasch logit scale where \eqn{\sigma \neq 1}. Floored at 0.
#'
#' @param data Response matrix/data.frame (items from 0).
#' @param thr_list Optional pre-fitted CML thresholds.
#' @param n_nodes Number of quadrature nodes.
#' @return Marginal reliability (numeric), or `NA` if the latent SD is not
#'   estimable.
#' @keywords internal
#' @noRd
.marginal_rxx <- function(data, thr_list = NULL, n_nodes = 161L) {
  data_mat <- as.matrix(data)
  if (is.null(thr_list)) {
    thr_list <- .fit_cml_thresholds(data_mat)
  }
  ge <- seq(-6, 6, length.out = 81L)
  sigma <- .estimate_prior_sd(
    .grid_loglik(data_mat, .logp_tables(thr_list, ge), ge),
    ge,
    0
  )
  if (!is.finite(sigma) || sigma <= 0) {
    return(NA_real_)
  }
  g <- seq(-6 * sigma, 6 * sigma, length.out = n_nodes)
  I <- vapply(
    g,
    function(th) {
      sum(vapply(
        thr_list,
        function(thr) {
          cats <- 0:length(thr)
          P <- .pcm_cat_probs(th, thr)
          E <- sum(cats * P)
          sum((cats - E)^2 * P)
        },
        numeric(1L)
      ))
    },
    numeric(1L)
  )
  w <- stats::dnorm(g, 0, sigma)
  w <- w / sum(w)
  max(1 - sum(w * (1 / I)) / sigma^2, 0)
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
  set.seed(seed)
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
    error = function(e) as.character(conditionMessage(e))
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
      .fit_cml_thresholds = .fit_cml_thresholds,
      .estimate_thetas = .estimate_thetas,
      .theta_wle = .theta_wle,
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
