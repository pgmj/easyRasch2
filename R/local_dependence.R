#' Q3 Residual Correlations for Local Dependence Assessment
#'
#' Computes Yen's Q3 residual correlations between item pairs using a Rasch
#' model fitted via Marginal Maximum Likelihood (MML) in `mirt`. When `cutoff`
#' is provided, high correlations (above the dynamic cut-off) are noted in the
#' table caption.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed.
#' @param cutoff Optional. A single numeric value added to the mean
#'   off-diagonal Q3 correlation to produce the dynamic cut-off threshold.
#'   When `NULL` (default), the raw Q3 matrix is returned without any cut-off.
#'   Use [RMlocdepQ3cutoff()] to derive a simulation-based cutoff.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying numeric data.frame.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object showing the lower triangle
#'   of the Q3 correlation matrix. When `cutoff` is provided, the caption
#'   includes the dynamic cut-off value.
#' * If `output = "dataframe"`: a data.frame (rounded to 2 decimal places) with
#'   the lower triangle of the Q3 correlation matrix; the upper triangle and
#'   diagonal are set to `NA`.
#'
#' @details
#' The Q3 statistic (Yen, 1984) is the correlation between residuals of pairs
#' of items after accounting for the latent trait. Under local independence,
#' Q3 values are expected to be around \eqn{-1/(k-1)} where \eqn{k} is the
#' number of items. When `cutoff` is supplied, the dynamic cut-off is the mean
#' of all off-diagonal Q3 values plus `cutoff`, following the approach of
#' Christensen et al. (2017).
#'
#' `mirt` is used for model fitting here because Q3 requires model-based
#' expected responses, which are most readily available from MML estimation.
#'
#' @seealso [RMlocdepQ3cutoff()] for simulation-based cutoff determination.
#'
#' @references
#' Yen, W. M. (1984). Effects of local item dependence on the fit and
#' equating performance of the three-parameter logistic model.
#' *Applied Psychological Measurement*, 8(2), 125–145.
#'
#' Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
#' for Yen's Q3: Identification of local dependence in the Rasch model.
#' *Applied Psychological Measurement*, 41(3), 178–194.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate binary item response data (10 items, 200 persons)
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#' )
#' colnames(sim_data) <- paste0("Item", 1:10)
#'
#' # Raw Q3 matrix (no cutoff)
#' RMlocdepQ3(sim_data)
#'
#' # Display Q3 table with a simulation-based cutoff
#' cutoff_res <- RMlocdepQ3cutoff(sim_data)
#' RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
#'
#' # Get the underlying data.frame
#' q3_df <- RMlocdepQ3(sim_data, output = "dataframe")
#' }
RMlocdepQ3 <- function(data, cutoff = NULL, output = "kable") {
  # --- Input validation -------------------------------------------------------
  if (!is.null(cutoff)) {
    if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
      stop("`cutoff` must be a single numeric value or NULL.", call. = FALSE)
    }
  }

  output <- match.arg(output, c("kable", "dataframe"))

  validate_response_data(data)

  # --- Model fitting ----------------------------------------------------------
  mirt_model <- mirt::mirt(
    data,
    model = 1,
    itemtype = "Rasch",
    verbose = FALSE,
    accelerate = "squarem"
  )

  # --- Compute Q3 residual correlations ---------------------------------------
  resid_mat <- mirt::residuals(mirt_model, type = "Q3", digits = 2, verbose = FALSE)

  diag(resid_mat) <- NA
  resid_df <- as.data.frame(resid_mat)

  # Round values
  resid_df[] <- lapply(resid_df, function(x) round(x, 2))

  # Blank out upper triangle and diagonal (set to NA)
  resid_df[upper.tri(resid_df)] <- NA
  diag(resid_df) <- NA

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(resid_df)
  }

  # "kable" output
  if (!is.null(cutoff)) {
    mean_resid <- mean(as.matrix(resid_df), na.rm = TRUE)
    dyn_cutoff <- mean_resid + cutoff
    caption_text <- paste0(
      "Dynamic cut-off: ", round(dyn_cutoff, 3),
      " (mean Q3 = ", round(mean_resid, 3),
      " + ", round(cutoff, 3), ").",
      " Correlations exceeding the cut-off may indicate local dependence."
    )
  } else {
    caption_text <- "Q3 residual correlations (lower triangle). No cut-off applied."
  }

  # Replace NA with "" for display (knitr::kable doesn't render NA nicely)
  resid_display <- as.data.frame(
    lapply(resid_df, function(x) {
      ifelse(is.na(x), "", as.character(x))
    }),
    stringsAsFactors = FALSE
  )
  rownames(resid_display) <- rownames(resid_df)

  knitr::kable(
    resid_display,
    caption = caption_text,
    format = "simple"
  )
}

#' Simulation-Based Q3 Cutoff for Local Dependence Assessment
#'
#' Uses parametric simulation to determine an appropriate cutoff for the
#' difference between the maximum and mean Q3 residual correlations under
#' the Rasch model. The resulting percentile values can be used as the
#' `cutoff` argument to [RMlocdepQ3()].
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed.
#' @param iterations Integer. Number of simulation iterations. Default `500`.
#' @param n_cores Integer or `NULL`. Number of parallel workers. If `NULL`
#'   (default), uses `getOption("mc.cores")` when that option is set, otherwise
#'   runs sequentially. Requires the **mirai** package for parallel execution.
#' @param seed Integer. Random seed for reproducibility. Default `123`.
#' @param progress Logical. Show a text progress bar? Default `TRUE`.
#'
#' @return A named list with:
#' \describe{
#'   \item{`results`}{data.frame with columns `mean`, `max`, and `diff`
#'     (max minus mean Q3) for each successful iteration.}
#'   \item{`actual_iterations`}{Number of successful iterations.}
#'   \item{`sample_n`}{Number of rows (persons) in `data`.}
#'   \item{`sample_summary`}{Summary of the estimated theta distribution.}
#'   \item{`max_diff`}{Maximum of the `diff` column.}
#'   \item{`sd_diff`}{Standard deviation of the `diff` column.}
#'   \item{`p95`}{95th percentile of `diff`.}
#'   \item{`p99`}{99th percentile of `diff`.}
#'   \item{`p995`}{99.5th percentile of `diff`.}
#'   \item{`p999`}{99.9th percentile of `diff`.}
#'   \item{`suggested_cutoff`}{The 99th percentile (recommended default).}
#' }
#'
#' @details
#' Item parameters are estimated via Conditional Maximum Likelihood (CML)
#' using `eRm::RM()` (dichotomous data) or `eRm::PCM()` (polytomous data).
#' Person parameters (thetas) are estimated with `eRm::person.parameter()`.
#' New datasets are simulated by re-sampling thetas with replacement and
#' generating responses from the estimated item parameters:
#' \itemize{
#'   \item Dichotomous data: `psychotools::rrm()`.
#'   \item Polytomous data: internal `sim_partial_score()`.
#' }
#' For each simulated dataset a Rasch model is fitted with `mirt::mirt()` and
#' the Q3 residual correlation matrix is computed. The `diff` statistic
#' (max Q3 − mean Q3) is collected across iterations; its percentile
#' distribution provides the cutoff values.
#'
#' Parallel execution via **mirai** is optional. When `n_cores > 1` and
#' **mirai** is installed, tasks are dispatched across `n_cores` daemons.
#' Set `options(mc.cores = 4)` to default to 4 cores.
#'
#' @seealso [RMlocdepQ3()] for using the resulting cutoff.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 300 * 8, replace = TRUE), nrow = 300, ncol = 8)
#' )
#' colnames(sim_data) <- paste0("Item", 1:8)
#'
#' cutoff_res <- RMlocdepQ3cutoff(sim_data, iterations = 200, seed = 1)
#' cutoff_res$suggested_cutoff
#'
#' # Use the result with RMlocdepQ3()
#' RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
#' }
RMlocdepQ3cutoff <- function(data, iterations = 500, n_cores = NULL,
                              seed = 123, progress = TRUE) {
  # --- Input validation -------------------------------------------------------
  validate_response_data(data)

  if (!is.numeric(iterations) || length(iterations) != 1L ||
      iterations < 1L || iterations != floor(iterations)) {
    stop("`iterations` must be a positive integer.", call. = FALSE)
  }
  iterations <- as.integer(iterations)

  if (!is.null(n_cores)) {
    if (!is.numeric(n_cores) || length(n_cores) != 1L ||
        n_cores < 1L || n_cores != floor(n_cores)) {
      stop("`n_cores` must be a positive integer or NULL.", call. = FALSE)
    }
    n_cores <- as.integer(n_cores)
  } else {
    mc_opt <- getOption("mc.cores")
    if (!is.null(mc_opt) && is.numeric(mc_opt) && mc_opt >= 1L) {
      n_cores <- as.integer(mc_opt)
    }
  }

  use_parallel <- !is.null(n_cores) && n_cores > 1L &&
    requireNamespace("mirai", quietly = TRUE)

  if (!is.null(n_cores) && n_cores > 1L && !use_parallel) {
    message("Install the 'mirai' package for parallel processing. Running sequentially.")
  }

  # --- Data characteristics ---------------------------------------------------
  data_mat <- as.matrix(data)
  sample_n  <- nrow(data_mat)
  data_max  <- max(data_mat, na.rm = TRUE)

  # --- Estimate item parameters and thetas (CML) ------------------------------
  if (data_max == 1L) {
    # Dichotomous
    erm_out        <- eRm::RM(data)
    item_locations <- erm_out$betapar * -1
    names(item_locations) <- colnames(data)
    thetas <- eRm::person.parameter(erm_out)[["theta.table"]][["Person Parameter"]]
    mode_flag <- "dichotomous"
  } else {
    # Polytomous
    erm_out   <- eRm::PCM(data)
    thresholds_mat <- extract_item_thresholds(data, pcm_fit = erm_out)
    thetas    <- eRm::person.parameter(erm_out)[["theta.table"]][["Person Parameter"]]
    n_items   <- ncol(data)
    itemlist  <- vector("list", n_items)
    for (i in seq_len(n_items)) {
      itemlist[[i]] <- list(stats::na.omit(thresholds_mat[i, ]))
    }
    mode_flag <- "polytomous"
  }

  item_names <- colnames(data)
  message(
    "RMlocdepQ3cutoff: starting ", iterations, " iterations",
    if (use_parallel) paste0(" on ", n_cores, " cores") else " sequentially",
    "."
  )

  # Pre-compute data_max_cols once (used in polytomous validation)
  data_max_cols <- if (mode_flag == "polytomous") {
    apply(data_mat, 2L, max, na.rm = TRUE)
  } else {
    NULL
  }

  # --- Per-iteration worker function ------------------------------------------
  # NOTE: This closure captures only plain R objects plus package::fun() calls
  # so that it serialises safely into mirai daemon processes.
  run_one <- function(seed_i) {
    set.seed(seed_i)
    input_thetas <- sample(thetas, size = sample_n, replace = TRUE)

    if (mode_flag == "dichotomous") {
      test_data <- as.data.frame(
        psychotools::rrm(input_thetas, item_locations, return_setting = FALSE)
      )
      # Discard if any item has fewer than 8 responses in either category
      n_resp <- colSums(as.matrix(test_data), na.rm = TRUE)
      if (min(n_resp, na.rm = TRUE) < 8L ||
          min(sample_n - n_resp, na.rm = TRUE) < 8L) {
        return("Insufficient responses in generated data.")
      }
    } else {
      # Inline polytomous simulation. Mirrors sim_poly_item() in
      # R/utils-simulation.R but is inlined here so that run_one() is a
      # self-contained closure that serialises safely for mirai daemons.
      deltas_unlisted <- lapply(itemlist, function(x) unlist(x))
      sim_item_inline <- function(deltas, thetas_v) {
        k <- length(deltas) + 1L
        n <- length(thetas_v)
        Y <- outer(thetas_v, deltas, "-")
        cs <- t(rbind(rep(0, n), apply(Y, 1L, cumsum)))
        ecs <- exp(cs)
        nrm <- rowSums(ecs)
        z <- ecs / nrm
        vapply(seq_len(n), function(i) {
          sample(0L:(k - 1L), 1L, prob = z[i, ])
        }, integer(1L))
      }
      test_data <- as.data.frame(
        sapply(deltas_unlisted, sim_item_inline, thetas_v = input_thetas)
      )
      colnames(test_data) <- item_names
      # Check all expected categories are present for each item
      for (j in seq_len(n_items)) {
        obs_cats <- sort(unique(test_data[[j]]))
        exp_cats <- seq(0L, data_max_cols[j], by = 1L)
        if (!all(exp_cats %in% obs_cats)) {
          return("Missing response categories in generated data.")
        }
      }
    }

    colnames(test_data) <- item_names
    mirt_fit <- tryCatch(
      mirt::mirt(test_data, model = 1L, itemtype = "Rasch",
                 verbose = FALSE, accelerate = "squarem"),
      error = function(e) NULL
    )
    if (is.null(mirt_fit)) return("mirt fitting failed.")

    resid <- mirt::residuals(mirt_fit, type = "Q3", digits = 2L, verbose = FALSE)
    diag(resid) <- NA
    data.frame(
      mean = mean(resid, na.rm = TRUE),
      max  = max(resid, na.rm = TRUE)
    )
  }

  # --- Generate per-iteration seeds -------------------------------------------
  set.seed(seed)
  iter_seeds <- sample.int(.Machine$integer.max, iterations)

  # --- Run iterations ---------------------------------------------------------
  if (use_parallel) {
    residcor <- .run_cutoff_parallel(run_one, iter_seeds, n_cores, progress)
  } else {
    residcor <- .run_cutoff_sequential(run_one, iter_seeds, progress)
  }

  # --- Post-processing --------------------------------------------------------
  is_bad <- vapply(residcor, is.character, logical(1L))
  actual_iterations <- iterations - sum(is_bad)

  good_results <- residcor[!is_bad]
  if (length(good_results) == 0L) {
    stop("All iterations failed. Check your data.", call. = FALSE)
  }

  results <- do.call(rbind, good_results)
  results$diff <- results$max - results$mean

  out <- list()
  out$results           <- results
  out$actual_iterations <- actual_iterations
  out$sample_n          <- sample_n
  out$sample_summary    <- summary(thetas)
  out$max_diff          <- max(results$diff)
  out$sd_diff           <- stats::sd(results$diff)
  out$p95               <- stats::quantile(results$diff, 0.95)
  out$p99               <- stats::quantile(results$diff, 0.99)
  out$p995              <- stats::quantile(results$diff, 0.995)
  out$p999              <- stats::quantile(results$diff, 0.999)
  out$suggested_cutoff  <- out$p99

  message(
    "RMlocdepQ3cutoff: finished. ",
    actual_iterations, " of ", iterations, " iterations succeeded.",
    " Suggested cutoff (p99): ", round(out$suggested_cutoff, 4), "."
  )

  out
}

# ---------------------------------------------------------------------------
# Internal parallel / sequential runners
# ---------------------------------------------------------------------------

#' @noRd
.run_cutoff_parallel <- function(run_one, iter_seeds, n_cores, progress) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  n <- length(iter_seeds)
  tasks <- vector("list", n)
  for (i in seq_len(n)) {
    tasks[[i]] <- mirai::mirai(
      { run_one(seed_i) },
      seed_i  = iter_seeds[[i]],
      run_one = run_one
    )
  }

  results <- vector("list", n)
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    on.exit({ close(pb) }, add = TRUE)
  }
  for (i in seq_len(n)) {
    res <- mirai::call_mirai(tasks[[i]])$data
    results[[i]] <- if (inherits(res, "errorValue")) "mirai error" else res
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  if (progress) cat("\n")
  results
}

#' @noRd
.run_cutoff_sequential <- function(run_one, iter_seeds, progress) {
  n <- length(iter_seeds)
  results <- vector("list", n)
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  for (i in seq_len(n)) {
    results[[i]] <- run_one(iter_seeds[[i]])
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  if (progress) cat("\n")
  results
}
