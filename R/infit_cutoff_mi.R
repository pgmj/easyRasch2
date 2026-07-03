#' Simulation-Based Infit MSQ Cutoff Determination for Multiply Imputed Data
#'
#' Extends \code{\link{RMitemInfitCutoff}} to work with multiply imputed datasets
#' produced by the `mice` package. Runs the parametric bootstrap simulation on
#' each imputed dataset and stacks the resulting distributions, so that the
#' final cutoff intervals reflect both sampling variability and imputation
#' uncertainty.
#'
#' @param mids_object A `mids` object (multiply imputed dataset) as returned by
#'   `mice::mice()`. Each completed dataset must contain only the item response
#'   columns to be analysed (i.e., no ID or grouping variables). Items must be
#'   scored starting at 0 (non-negative integers).
#' @param iterations Integer. Total number of simulation iterations to run
#'   across all imputations. These are distributed approximately evenly across
#'   the `m` imputed datasets (default 500).
#' @param parallel Logical. Use parallel processing via `mirai` within each
#'   imputed dataset (default `TRUE`). Passed to \code{\link{RMitemInfitCutoff}}.
#' @param n_cores Integer or `NULL`. Number of parallel workers. Passed to
#'   \code{\link{RMitemInfitCutoff}}.
#' @param verbose Logical. Show progress messages (default `FALSE`).
#' @param seed Integer or `NULL`. Master random seed for reproducibility. A
#'   unique per-imputation seed is derived from this value.
#' @param cutoff_method Character string specifying how cutoff intervals are
#'   computed from the stacked distribution. Either `"hdci"` (default) for the
#'   Highest Density Interval via `ggdist::hdci()`, or `"quantile"` for the
#'   2.5th/97.5th percentiles via `stats::quantile()`.
#' @param hdci_width Numeric. Width of the HDCI when `cutoff_method = "hdci"`.
#'   Default is `0.999` (99.9% HDCI). Ignored when
#'   `cutoff_method = "quantile"`.
#'
#' @return A list with the same structure as \code{\link{RMitemInfitCutoff}}, so
#'   that the result can be passed directly to \code{\link{RMitemInfit}},
#'   \code{\link{RMitemInfitMI}}, and \code{\link{RMitemInfitPlot}}:
#' \describe{
#'   \item{`results`}{data.frame with columns `iteration`, `imputation`,
#'     `Item`, `InfitMSQ`, `OutfitMSQ` — the stacked simulation results from
#'     all imputed datasets.}
#'   \item{`item_cutoffs`}{data.frame with per-item cutoff summaries: `Item`,
#'     `infit_low`, `infit_high`, `outfit_low`, `outfit_high`. Computed from
#'     the stacked distribution.}
#'   \item{`actual_iterations`}{Total number of successful iterations across
#'     all imputations.}
#'   \item{`sample_n`}{Number of rows (respondents) per imputed dataset.}
#'   \item{`sample_summary`}{Summary statistics of estimated person parameters
#'     from the first imputed dataset.}
#'   \item{`item_names`}{Character vector of item names.}
#'   \item{`cutoff_method`}{The method used to compute cutoffs.}
#'   \item{`hdci_width`}{The HDCI width used.}
#'   \item{`n_imputations`}{Number of imputed datasets used.}
#'   \item{`iterations_per_imputation`}{Integer vector of requested iterations
#'     per imputed dataset.}
#'   \item{`actual_iterations_per_imputation`}{Integer vector of successful
#'     iterations per imputed dataset.}
#' }
#'
#' @details
#' The function completes each of the `m` imputed datasets via
#' `mice::complete()`, then calls \code{\link{RMitemInfitCutoff}} on each one. The
#' total number of iterations is split approximately evenly across imputations
#' (i.e., each imputed dataset receives `ceiling(iterations / m)` or
#' `floor(iterations / m)` iterations). The per-imputation simulation results
#' are stacked into a single distribution from which cutoff intervals are
#' computed, naturally incorporating imputation uncertainty.
#'
#' Imputed datasets that cause model convergence failures are dropped with a
#' warning. If all imputations fail, the function stops with an error.
#'
#' The `mice` package must be installed (it is in Suggests, not Imports).
#'
#' @seealso \code{\link{RMitemInfitCutoff}}, \code{\link{RMitemInfitMI}},
#'   \code{\link{RMitemInfitPlot}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("mice", quietly = TRUE) &&
#'     requireNamespace("iarm", quietly = TRUE) &&
#'     requireNamespace("ggdist", quietly = TRUE)) {
#'   # Create example data with ~10% MCAR missingness
#'   set.seed(42)
#'   mat <- matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
#'   mat[sample(length(mat), round(0.10 * length(mat)))] <- NA
#'   sim_data <- as.data.frame(mat)
#'   colnames(sim_data) <- paste0("Item", 1:8)
#'
#'   # mice's ordinal method (`polr`) requires the items to be ordered
#'   # factors, so code them as such before imputing. RMitemInfitCutoffMI()
#'   # converts the completed factors back to numeric internally.
#'   sim_data[] <- lapply(sim_data, function(x) factor(x, ordered = TRUE))
#'
#'   # Impute (use more imputations, e.g. m = 5+, in real analyses)
#'   imp <- mice::mice(sim_data, m = 2, method = "polr", seed = 123,
#'                     printFlag = FALSE)
#'
#'   # Compute simulation-based cutoffs across imputations
#'   # (use more iterations, e.g. 250+, in real analyses)
#'   cutoff_mi <- RMitemInfitCutoffMI(imp, iterations = 50, parallel = FALSE,
#'                                 seed = 42)
#'   cutoff_mi$item_cutoffs
#'
#'   # Use with RMitemInfitMI()
#'   RMitemInfitMI(imp, cutoff = cutoff_mi)
#' }
#' }
RMitemInfitCutoffMI <- function(
  mids_object,
  iterations = 500,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.999
) {
  # --- Check mice -----------------------------------------------------------

  if (!requireNamespace("mice", quietly = TRUE)) {
    stop(
      "Package 'mice' is required for RMitemInfitCutoffMI() but is not installed.\n",
      "Install it with: install.packages(\"mice\")",
      call. = FALSE
    )
  }

  if (!inherits(mids_object, "mids")) {
    stop(
      "`mids_object` must be a 'mids' object as returned by mice::mice().",
      call. = FALSE
    )
  }

  cutoff_method <- match.arg(cutoff_method, c("hdci", "quantile"))

  if (cutoff_method == "hdci" && !requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required when cutoff_method = \"hdci\" but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")\n",
      "Alternatively, use cutoff_method = \"quantile\" to avoid this dependency.",
      call. = FALSE
    )
  }

  m <- mids_object$m

  # --- Distribute iterations across imputations -------------------------------
  base_iter <- iterations %/% m
  remainder <- iterations %% m
  iters_per_imp <- rep(base_iter, m)
  if (remainder > 0L) {
    iters_per_imp[seq_len(remainder)] <- iters_per_imp[seq_len(remainder)] + 1L
  }

  # --- Derive per-imputation seeds -------------------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }
  imp_seeds <- sample.int(.Machine$integer.max, m)

  # --- Run RMitemInfitCutoff on each imputed dataset ------------------------------
  all_results <- vector("list", m)
  actual_per_imp <- integer(m)
  sample_summary_first <- NULL
  sample_n <- NULL
  item_names <- NULL
  n_failed <- 0L

  for (i in seq_len(m)) {
    completed_data <- .mi_completed_to_numeric(
      mice::complete(mids_object, action = i)
    )

    if (verbose) {
      message(sprintf(
        "Processing imputation %d/%d (%d iterations)...",
        i,
        m,
        iters_per_imp[i]
      ))
    }

    result_i <- tryCatch(
      {
        RMitemInfitCutoff(
          data = completed_data,
          iterations = iters_per_imp[i],
          parallel = parallel,
          n_cores = n_cores,
          verbose = verbose,
          seed = imp_seeds[i],
          cutoff_method = "quantile",
          hdci_width = hdci_width
        )
      },
      error = function(e) {
        warning(
          sprintf(
            "RMitemInfitCutoff() failed for imputation %d: %s",
            i,
            conditionMessage(e)
          ),
          call. = FALSE
        )
        NULL
      }
    )

    if (is.null(result_i)) {
      n_failed <- n_failed + 1L
      next
    }

    # Tag rows with imputation number
    res_df <- result_i$results
    res_df$imputation <- i
    all_results[[i]] <- res_df
    actual_per_imp[i] <- result_i$actual_iterations

    # Capture metadata from first successful imputation
    if (is.null(item_names)) {
      item_names <- result_i$item_names
      sample_n <- result_i$sample_n
      sample_summary_first <- result_i$sample_summary
    }
  }

  if (n_failed == m) {
    stop(
      "RMitemInfitCutoff() failed for all ",
      m,
      " imputed datasets. ",
      "Check your data.",
      call. = FALSE
    )
  }

  if (n_failed > 0L) {
    warning(
      sprintf(
        "%d of %d imputed datasets failed and were excluded.",
        n_failed,
        m
      ),
      call. = FALSE
    )
  }

  # --- Stack results ----------------------------------------------------------
  stacked_df <- do.call(
    rbind,
    all_results[!vapply(all_results, is.null, logical(1L))]
  )
  rownames(stacked_df) <- NULL

  # Re-number iterations sequentially across imputations
  stacked_df$iteration <- seq_len(nrow(stacked_df) %/% length(item_names))
  # Actually, iteration numbering per row isn't meaningful after stacking;
  # keep original per-imputation iteration numbers and add imputation column.
  # The $iteration column from each imputation is already present; we just

  # leave it as-is.

  total_actual <- sum(actual_per_imp)

  # --- Compute cutoffs from stacked distribution ------------------------------
  unique_items <- unique(stacked_df$Item)
  item_cutoffs <- do.call(
    rbind,
    lapply(unique_items, function(item) {
      sub <- stacked_df[stacked_df$Item == item, ]
      if (cutoff_method == "hdci") {
        infit_interval <- ggdist::hdci(sub$InfitMSQ, .width = hdci_width)
        outfit_interval <- ggdist::hdci(sub$OutfitMSQ, .width = hdci_width)
        data.frame(
          Item = item,
          infit_low = infit_interval[1L, 1L],
          infit_high = infit_interval[1L, 2L],
          outfit_low = outfit_interval[1L, 1L],
          outfit_high = outfit_interval[1L, 2L],
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      } else {
        data.frame(
          Item = item,
          infit_low = stats::quantile(sub$InfitMSQ, 0.025, na.rm = TRUE),
          infit_high = stats::quantile(sub$InfitMSQ, 0.975, na.rm = TRUE),
          outfit_low = stats::quantile(sub$OutfitMSQ, 0.025, na.rm = TRUE),
          outfit_high = stats::quantile(sub$OutfitMSQ, 0.975, na.rm = TRUE),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
    })
  )
  rownames(item_cutoffs) <- NULL

  # --- Return -----------------------------------------------------------------
  list(
    results = stacked_df,
    item_cutoffs = item_cutoffs,
    actual_iterations = total_actual,
    sample_n = sample_n,
    sample_summary = sample_summary_first,
    item_names = item_names,
    cutoff_method = cutoff_method,
    hdci_width = hdci_width,
    n_imputations = m - n_failed,
    iterations_per_imputation = iters_per_imp,
    actual_iterations_per_imputation = actual_per_imp
  )
}
