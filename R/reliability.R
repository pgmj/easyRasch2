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
#' \doi{10.31234/osf.io/h54k8}
#'
#' @seealso [RMreliability()]
#'
#' @export
RMUreliability <- function(input_draws, level = 0.95, verbose = FALSE) {

  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required for RMUreliability().\n",
         "Install it with: install.packages(\"ggdist\")",
         call. = FALSE)
  }
  if (!is.numeric(level) || level <= 0 || level >= 1) {
    stop("`level` must be a numeric value strictly between 0 and 1.",
         call. = FALSE)
  }

  input_draws <- as.matrix(input_draws)
  if (ncol(input_draws) < 2L) {
    stop("`input_draws` must have at least 2 columns.", call. = FALSE)
  }

  sds <- apply(input_draws, 2L, stats::sd, na.rm = TRUE)
  zero_sd_cols <- which(sds == 0)
  if (length(zero_sd_cols) > 0L) {
    warning(sprintf(
      "Found %d column(s) with zero standard deviation (columns: %s).",
      length(zero_sd_cols), paste(zero_sd_cols, collapse = ", ")
    ), call. = FALSE)
  }

  na_count <- sum(is.na(input_draws))
  if (na_count > 0L) {
    warning(sprintf("Found %d NA value(s) in `input_draws`.", na_count),
            call. = FALSE)
  }

  if (verbose) {
    message("Subjects: ", nrow(input_draws),
            "; draws: ", ncol(input_draws))
  }

  col_select <- sample.int(ncol(input_draws), replace = FALSE)
  half       <- floor(length(col_select) / 2)
  cols_a     <- col_select[seq_len(half)]
  cols_b     <- col_select[(half + 1L):(2L * half)]

  draws_a <- input_draws[, cols_a, drop = FALSE]
  draws_b <- input_draws[, cols_b, drop = FALSE]

  rel_post <- vapply(seq_len(ncol(draws_a)), function(i) {
    x <- draws_a[, i]
    y <- draws_b[, i]
    if (stats::var(x, na.rm = TRUE) == 0 ||
        stats::var(y, na.rm = TRUE) == 0) {
      return(0)
    }
    stats::cor(x, y, method = "pearson", use = "complete.obs")
  }, numeric(1L))

  hdci <- ggdist::mean_hdci(rel_post, .width = level)
  colnames(hdci)[1:3] <- c("rmu_estimate", "hdci_lowerbound", "hdci_upperbound")
  hdci
}

# ---------------------------------------------------------------------------

#' Reliability metrics for a Rasch model
#'
#' Computes three reliability indices for a Rasch / partial credit model:
#' the Person Separation Index (PSI) via `eRm::SepRel()`, empirical
#' reliability via `mirt::empirical_rxx()`, and Relative Measurement
#' Uncertainty (RMU) via [RMUreliability()] applied to plausible values
#' from `mirt::fscores()`.
#'
#' Confidence intervals for **PSI** and **Empirical** reliability are obtained
#' by non-parametric bootstrap (refitting both the eRm and mirt models on
#' each resample). The RMU interval is the HDCI of correlations across
#' plausible-value draws, averaged over `rmu_iter` random splits of the draws.
#'
#' This is a CRAN-friendly replacement for `easyRasch::RIreliability()`. The
#' TAM-based plausible-values path has been removed; mirt is the only MML
#' backend.
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
#'   empirical reliability and the plausible-value seed. One of `"WLE"`
#'   (default), `"EAP"`, `"MAP"`, `"ML"`. Plausible draws themselves are
#'   produced by Metropolis-Hastings.
#' @param boot Logical. If `TRUE`, run a non-parametric bootstrap to obtain
#'   CIs for PSI and Empirical reliability. Default `FALSE`.
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
#' For **dichotomous** data (max = 1), `mirt::mirt(..., itemtype = "Rasch")`
#' and `eRm::RM()` are fitted. For **polytomous** data, the same `itemtype`
#' is used in mirt and `eRm::PCM()` in eRm.
#'
#' PSI is computed by `eRm::SepRel()` on the person-parameter object from
#' `eRm::person.parameter()`. Note that PSI excludes respondents with extreme
#' (min/max) raw scores by construction.
#'
#' RMU is from Bignardi, Kievit, & Bürkner (2025), modified here to use mirt
#' plausible values rather than fully Bayesian posterior draws (see Mislevy,
#' 1991, for the plausible-values framework).
#'
#' Bootstrap iterations that fail to converge in either mirt or eRm are
#' silently dropped.
#'
#' @references
#' Bignardi, G., Kievit, R., & Bürkner, P. C. (2025). A general method for
#' estimating reliability using Bayesian Measurement Uncertainty. *PsyArXiv*.
#' \doi{10.31234/osf.io/h54k8}
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
#' \dontrun{
#' set.seed(1)
#' RMreliability(eRm::raschdat1[, 1:20], draws = 1000)
#'
#' # With bootstrap CI for PSI and Empirical (slower)
#' RMreliability(eRm::raschdat1[, 1:20], draws = 1000,
#'               boot = TRUE, boot_iter = 200, parallel = FALSE, seed = 42)
#' }
RMreliability <- function(data,
                          conf_int    = 0.95,
                          draws       = 1000,
                          rmu_iter    = 50,
                          estim       = "WLE",
                          boot        = FALSE,
                          boot_iter   = 200,
                          parallel    = TRUE,
                          n_cores     = NULL,
                          seed        = NULL,
                          verbose     = FALSE,
                          theta_range = c(-10, 10),
                          output      = "kable") {

  output <- match.arg(output, c("kable", "dataframe"))
  estim  <- match.arg(estim,  c("WLE", "EAP", "MAP", "ML"))

  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required for RMreliability().\n",
         "Install it with: install.packages(\"ggdist\")",
         call. = FALSE)
  }

  validate_response_data(data)

  if (nrow(stats::na.omit(data)) == 0L) {
    stop("No complete cases in `data`.", call. = FALSE)
  }

  data_mat       <- as.matrix(data)
  is_polytomous  <- max(data_mat, na.rm = TRUE) > 1L
  n_persons      <- nrow(data)
  n_items        <- ncol(data)

  # rgl workaround for any iarm-related fallout
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # --- Full-sample fits ------------------------------------------------------
  mirt_fit <- mirt::mirt(
    data       = data,
    model      = 1,
    itemtype   = "Rasch",
    verbose    = FALSE,
    accelerate = "squarem"
  )
  erm_fit <- if (is_polytomous) eRm::PCM(data) else eRm::RM(data)

  # --- Empirical reliability point estimate ----------------------------------
  emp_rel <- as.numeric(
    mirt::empirical_rxx(
      mirt::fscores(mirt_fit,
                    method         = estim,
                    theta_lim      = theta_range,
                    full.scores.SE = TRUE,
                    verbose        = FALSE)
    )
  )

  # --- PSI point estimate ----------------------------------------------------
  psi <- as.numeric(eRm::SepRel(eRm::person.parameter(erm_fit))$sep.rel)

  # --- Plausible values + RMU ------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  pvs <- mirt::fscores(
    mirt_fit,
    method          = estim,
    theta_lim       = theta_range,
    plausible.draws = draws,
    plausible.type  = "MH",
    verbose         = FALSE
  )
  rmu_input <- do.call(cbind, lapply(pvs, as.numeric))

  rmu_iter_results <- do.call(
    rbind,
    lapply(seq_len(rmu_iter), function(i) {
      RMUreliability(rmu_input, level = conf_int)[, 1:3, drop = FALSE]
    })
  )
  rmu_summary <- list(
    estimate = mean(rmu_iter_results$rmu_estimate),
    lower    = mean(rmu_iter_results$hdci_lowerbound),
    upper    = mean(rmu_iter_results$hdci_upperbound)
  )

  # --- Bootstrap (optional) --------------------------------------------------
  psi_lower <- NA_real_
  psi_upper <- NA_real_
  emp_lower <- NA_real_
  emp_upper <- NA_real_
  actual_boot <- NA_integer_

  if (isTRUE(boot)) {
    use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)
    if (parallel && !use_parallel) {
      message("Install 'mirai' for parallel processing: install.packages(\"mirai\")")
      message("Running sequentially...")
    }
    if (use_parallel) {
      if (is.null(n_cores)) n_cores <- getOption("mc.cores")
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

    if (!is.null(seed)) set.seed(seed + 1L)
    boot_seeds <- sample.int(.Machine$integer.max, boot_iter)

    boot_args <- list(
      data          = data,
      is_polytomous = is_polytomous,
      estim         = estim,
      theta_range   = theta_range
    )

    if (use_parallel) {
      boot_results <- run_reliability_boot_parallel(
        boot_iter, boot_seeds, boot_args, n_cores, verbose
      )
    } else {
      boot_results <- run_reliability_boot_sequential(
        boot_iter, boot_seeds, boot_args, verbose
      )
    }

    ok          <- vapply(boot_results, is.list, logical(1L))
    successful  <- boot_results[ok]
    actual_boot <- length(successful)

    if (actual_boot < 2L) {
      warning("Fewer than 2 bootstrap iterations succeeded; CI not reported.",
              call. = FALSE)
    } else {
      psi_vec <- vapply(successful, function(x) x$psi,       numeric(1L))
      emp_vec <- vapply(successful, function(x) x$empirical, numeric(1L))

      psi_int <- ggdist::hdci(psi_vec, .width = conf_int)
      emp_int <- ggdist::hdci(emp_vec, .width = conf_int)

      psi_lower <- psi_int[1L, 1L]
      psi_upper <- psi_int[1L, 2L]
      emp_lower <- emp_int[1L, 1L]
      emp_upper <- emp_int[1L, 2L]
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
    metric   = c("PSI",
                 paste0("Empirical (", estim, ")"),
                 paste0("RMU (", estim, ")")),
    estimate = round(c(psi, emp_rel, rmu_summary$estimate), 3),
    lower    = round(c(psi_lower, emp_lower, rmu_summary$lower), 3),
    upper    = round(c(psi_upper, emp_upper, rmu_summary$upper), 3),
    notes    = c(boot_note, boot_note, rmu_note),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  if (output == "dataframe") {
    return(result_df)
  }

  knitr::kable(
    result_df,
    format    = "pipe",
    col.names = c("Metric", "Estimate",
                  paste0("Lower (", conf_pct, "% HDCI)"),
                  paste0("Upper (", conf_pct, "% HDCI)"),
                  "Notes"),
    caption   = paste0(
      "Reliability for ", n_items, " items, n = ", n_persons,
      ". PSI excludes min/max scoring respondents."
    )
  )
}

# ---------------------------------------------------------------------------
# Internal: single bootstrap iteration
# ---------------------------------------------------------------------------

#' Run a single reliability bootstrap iteration
#'
#' @keywords internal
#' @noRd
run_single_reliability_boot <- function(seed, data_list) {
  set.seed(seed)
  idx     <- sample.int(nrow(data_list$data), nrow(data_list$data), replace = TRUE)
  dat_b   <- data_list$data[idx, , drop = FALSE]

  tryCatch({
    mirt_fit <- mirt::mirt(
      data       = dat_b,
      model      = 1,
      itemtype   = "Rasch",
      verbose    = FALSE,
      accelerate = "squarem"
    )
    emp <- as.numeric(
      mirt::empirical_rxx(
        mirt::fscores(mirt_fit,
                      method         = data_list$estim,
                      theta_lim      = data_list$theta_range,
                      full.scores.SE = TRUE,
                      verbose        = FALSE)
      )
    )

    erm_fit <- if (data_list$is_polytomous) eRm::PCM(dat_b) else eRm::RM(dat_b)
    psi <- as.numeric(eRm::SepRel(eRm::person.parameter(erm_fit))$sep.rel)

    list(empirical = emp, psi = psi)
  }, error = function(e) as.character(conditionMessage(e)))
}

# ---------------------------------------------------------------------------
# Internal: parallel runner (mirai)
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
run_reliability_boot_parallel <- function(boot_iter, boot_seeds, boot_args,
                                          n_cores, verbose = FALSE) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = boot_iter, style = 3)
  }

  tasks <- lapply(seq_len(boot_iter), function(i) {
    mirai::mirai(
      { run_single_reliability_boot(seed, data_list) },
      seed                          = boot_seeds[i],
      data_list                     = boot_args,
      run_single_reliability_boot   = run_single_reliability_boot
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
run_reliability_boot_sequential <- function(boot_iter, boot_seeds, boot_args,
                                            verbose = FALSE) {
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = boot_iter, style = 3)

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
