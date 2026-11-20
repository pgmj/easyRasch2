#' Simulation-based infit MSQ cutoffs for multiply imputed data
#'
#' The counterpart of [RMitemInfitCutoff()] for data imputed with the
#' \pkg{mice} package. Runs the simulations in each imputed dataset and
#' stacks them, so that the interval per item reflects both sampling
#' variability and the uncertainty from imputation. The result is used by
#' [RMitemInfitMI()].
#'
#' @param mids_object A `mids` object from `mice::mice()`. Each completed
#'   dataset must contain only the item response columns, with no ID or
#'   grouping variables. Items are scored from 0. Items imputed as
#'   factors, as the `polr` method needs, are converted back to numbers.
#'   Their levels must be the response codes (`"0"`, `"1"`, ...).
#' @param iterations Integer. Total number of simulated datasets across all
#'   imputations, split about evenly between the `m` imputed datasets.
#'   Default `500`. See the section on interval width.
#' @param parallel Logical. Run the simulations within each imputed
#'   dataset in parallel via \pkg{mirai}. Default `TRUE`. Passed to
#'   [RMitemInfitCutoff()].
#' @param n_cores Integer or `NULL`. Number of parallel workers. Passed to
#'   \code{\link{RMitemInfitCutoff}}.
#' @param verbose Logical. Show progress messages (default `FALSE`).
#' @param seed Integer or `NULL`. Random seed for reproducibility, from
#'   which a seed for each imputation is derived. See
#'   [easyRasch2-reproducibility] for what this guarantees and how it
#'   interacts with `parallel`.
#' @param cutoff_method Character. How the interval per item is computed
#'   from the stacked simulations: `"hdci"` (default), the highest-density
#'   continuous interval from `ggdist::hdci()` (needs \pkg{ggdist}), or
#'   `"quantile"`, the 2.5th to 97.5th percentiles.
#' @param hdci_width Numeric. Width of the interval with
#'   `cutoff_method = "hdci"`. Default `0.999`, see the section on interval
#'   width. Ignored with `"quantile"`.
#'
#' @return A list like that of [RMitemInfitCutoff()], which
#'   [RMitemInfitMI()] and [RMitemInfitPlot()] accept, with:
#' \describe{
#'   \item{`results`}{data.frame with the columns `iteration`,
#'     `imputation`, `Item`, `InfitMSQ` and `OutfitMSQ`, the stacked
#'     simulations from all imputed datasets.}
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
#' Each of the `m` imputed datasets is completed with `mice::complete()` and
#' passed to [RMitemInfitCutoff()], with `ceiling(iterations / m)` or
#' `floor(iterations / m)` iterations each. The simulations from all
#' imputations are stacked, and the interval per item is computed from the
#' stack.
#'
#' Imputed datasets in which the model cannot be fitted are dropped with a
#' warning. If all imputations fail, the function stops with an error.
#'
#' The \pkg{mice} package must be installed (it is in Suggests).
#'
#' @section Interval width under multiple imputation:
#' `hdci_width` defaults to `0.999` here, not to the `0.95` of
#' [RMitemInfitCutoff()]. [RMitemInfitMI()] has no corrected p-value yet, so
#' under imputation the interval is still the decision rule, and it keeps a
#' width chosen for that. Its family-wise error rate is `1 - 0.999^k` over
#' `k` items. An interval this wide needs roughly 5000 iterations in all
#' before it reaches its stated width (Johansson, 2026). Combining bootstrap
#' p-values across imputations is planned but needs its own calibration,
#' since Johansson (2026) covers complete data only.
#'
#' @references
#' Johansson, M. (2026). Simulation-based cutoffs for conditional item fit in
#' Rasch models: Iterations, multiplicity correction, and decision stability.
#' *PsyArXiv*. \doi{10.31234/osf.io/7pqz4_v2}
#'
#' @seealso [RMitemInfitMI()] for the pooled test, [RMitemInfitPlot()] to
#'   plot the simulated values, [RMitemInfitCutoff()] for complete data.
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
#'   # Simulation-based cutoffs across imputations (use far more
#'   # iterations in real analyses, see the section on interval width)
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

  .mi_validate_structure(mids_object) # once, not per imputation

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
