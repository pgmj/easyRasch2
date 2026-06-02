#' Bootstrap Item-Restscore Misfit Detection
#'
#' Non-parametric bootstrap of item-restscore fit using
#' `iarm::item_restscore()`. For each iteration, a sample of size `samplesize`
#' is drawn from `data` with replacement, the appropriate Rasch model is
#' refitted, and item-restscore results are classified as `"overfit"`,
#' `"underfit"`, or `"no misfit"` based on the BH-adjusted p-value (< .05) and
#' the sign of `expected - observed`. The function returns the percentage of
#' iterations in which each item is flagged.
#'
#' Useful with large samples, where the asymptotic test underlying
#' \code{\link{RMitemRestscore}} can flag items that are not practically
#' misfitting; bootstrapping gives a more nuanced view of the probability of
#' an item actually being misfit.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers).
#' @param iterations Integer. Number of bootstrap samples (default 200).
#' @param samplesize Integer. Size of each bootstrap sample (default 600). Must
#'   not exceed `nrow(data)`.
#' @param parallel Logical. Use parallel processing via `mirai` if available
#'   (default `TRUE`).
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first. If neither is set and
#'   `parallel = TRUE`, a warning is issued and execution falls back to
#'   sequential processing.
#' @param cutoff Numeric. Items flagged in fewer than this percentage of
#'   iterations are excluded from the kable output (default 5). Has no effect
#'   when `output = "dataframe"` or `output = "raw"`.
#' @param verbose Logical. Show a progress bar (default `FALSE`).
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param output Either `"kable"` (default) for a formatted `knitr::kable()`
#'   table, `"dataframe"` for the per-item summary data.frame, or `"raw"` for
#'   the per-iteration long data.frame (useful for custom plotting).
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object listing items flagged in
#'   more than `cutoff`% of iterations, with columns Item, Item-restscore
#'   result, % of iterations, Conditional MSQ infit, and Relative item
#'   location, and a caption noting iteration count, bootstrap size, and
#'   number of complete cases.
#' * If `output = "dataframe"`: a data.frame with one row per item ×
#'   classification combination (`Item`, `item_restscore`, `n`, `percent`,
#'   `Infit_MSQ`, `Relative_location`), including `"no misfit"` rows.
#' * If `output = "raw"`: a long data.frame with one row per item ×
#'   successful iteration (`iteration`, `Item`, `item_restscore`, `diff`,
#'   `diff_abs`), where `diff = expected - observed`.
#'
#' @details
#' For **dichotomous** data (maximum score = 1), the full-sample model is
#' fitted via `eRm::RM()` and bootstrap iterations refit via `eRm::RM()` with
#' `se = FALSE` for speed.
#'
#' For **polytomous** data (maximum score > 1), the full-sample model is
#' fitted via `eRm::PCM()` (so item locations can be obtained from
#' `eRm::thresholds()`) and bootstrap iterations refit via
#' `psychotools::pcmodel(..., hessian = FALSE)` for speed.
#' `iarm::item_restscore()` accepts both model classes.
#'
#' Conditional infit MSQ (computed once on the full sample via
#' `iarm::out_infit()`) and relative item locations are reported alongside the
#' bootstrap percentages for context.
#'
#' Iterations that fail (e.g., due to convergence issues on a degenerate
#' bootstrap sample) are silently discarded; the caption / `actual_iterations`
#' reflects only successful runs.
#'
#' Parallel processing is provided by the `mirai` package (optional). Install
#' it with `install.packages("mirai")` to enable parallelisation.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' @references
#' Kreiner, S. (2011). A Note on Item-Restscore Association in Rasch Models.
#' *Applied Psychological Measurement, 35*(7), 557-561.
#' \doi{10.1177/0146621611410227}
#'
#' @seealso \code{\link{RMitemRestscore}}, \code{\link{RMitemInfit}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 400 * 8, replace = TRUE), nrow = 400, ncol = 8)
#' )
#' colnames(sim_data) <- paste0("Item", 1:8)
#'
#' # Few iterations for a fast example; use 100+ in real analyses
#' # Default kable output (only items flagged > cutoff%)
#' RMitemRestscoreBoot(sim_data, iterations = 50, samplesize = 300,
#'                 parallel = FALSE, seed = 1)
#'
#' # Per-item summary data.frame (all classifications, including "no misfit")
#' summary_df <- RMitemRestscoreBoot(sim_data, iterations = 50, samplesize = 300,
#'                               parallel = FALSE, seed = 1,
#'                               output = "dataframe")
#'
#' # Per-iteration long data for custom plotting
#' raw_df <- RMitemRestscoreBoot(sim_data, iterations = 50, samplesize = 300,
#'                           parallel = FALSE, seed = 1, output = "raw")
#'
#' # Distribution of (expected - observed) across iterations, per item
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(raw_df, aes(x = Item, y = diff)) +
#'     geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
#'     geom_violin(fill = "grey90", colour = NA) +
#'     geom_jitter(aes(colour = item_restscore),
#'                 width = 0.15, alpha = 0.5, size = 0.8) +
#'     scale_colour_manual(values = c("overfit"   = "#377eb8",
#'                                    "underfit"  = "#e41a1c",
#'                                    "no misfit" = "grey60")) +
#'     labs(y = "Expected - observed restscore correlation",
#'          colour = NULL) +
#'     theme_minimal() +
#'     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#' }
#' }
RMitemRestscoreBoot <- function(data,
                            iterations = 200,
                            samplesize = 600,
                            parallel   = TRUE,
                            n_cores    = NULL,
                            cutoff     = 5,
                            verbose    = FALSE,
                            seed       = NULL,
                            output     = "kable") {

  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMitemRestscoreBoot() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe", "raw"))

  validate_response_data(data)

  if (samplesize > nrow(data)) {
    stop(
      paste0("`samplesize` (", samplesize,
             ") cannot be larger than the number of rows in `data` (",
             nrow(data), ")."),
      call. = FALSE
    )
  }

  # rgl workaround for iarm dependency chain (vcdExtra -> rgl)
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  data_mat      <- as.matrix(data)
  item_names    <- colnames(data_mat)
  n_items       <- ncol(data_mat)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  # --- Full-sample model: locations + conditional infit ----------------------
  if (is_polytomous) {
    erm_full <- eRm::PCM(data)
    thresh_table <- eRm::thresholds(erm_full)$threshtable[[1L]]
    if ("Location" %in% colnames(thresh_table)) {
      item_avg_locations <- thresh_table[, "Location"]
    } else {
      item_avg_locations <- rowMeans(thresh_table, na.rm = TRUE)
    }
  } else {
    erm_full <- eRm::RM(data)
    item_avg_locations <- stats::coef(erm_full, "beta") * -1
  }
  pp <- eRm::person.parameter(erm_full)
  person_avg_location <- mean(pp$theta.table[["Person Parameter"]], na.rm = TRUE)
  relative_item_avg_locations <- item_avg_locations - person_avg_location

  cfit       <- iarm::out_infit(erm_full)
  n_complete <- nrow(stats::na.omit(data))

  # --- Parallel setup --------------------------------------------------------
  use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)
  if (parallel && !use_parallel) {
    message("Install 'mirai' package for parallel processing: install.packages(\"mirai\")")
    message("Running sequentially...")
  }
  if (use_parallel) {
    if (is.null(n_cores)) {
      n_cores <- getOption("mc.cores")
    }
    if (is.null(n_cores)) {
      warning(
        paste0(
          "For parallel processing, specify n_cores or set options(mc.cores = N).\n",
          "(Use `parallel::detectCores()` to see how many cores are available.)\n",
          "Falling back to sequential (single core) processing."
        ),
        call. = FALSE
      )
      use_parallel <- FALSE
    } else {
      n_cores <- min(n_cores, iterations)
    }
  }

  # --- Per-iteration seeds for reproducibility -------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }
  boot_seeds <- sample.int(.Machine$integer.max, iterations)

  boot_data_list <- list(
    data          = data,
    samplesize    = samplesize,
    is_polytomous = is_polytomous,
    item_names    = item_names
  )

  # --- Run bootstrap ---------------------------------------------------------
  if (use_parallel) {
    results_raw <- run_boot_restscore_parallel(
      iterations, boot_seeds, boot_data_list, n_cores, verbose
    )
  } else {
    results_raw <- run_boot_restscore_sequential(
      iterations, boot_seeds, boot_data_list, verbose
    )
  }

  ok <- vapply(results_raw, is.data.frame, logical(1L))
  successful <- results_raw[ok]

  if (length(successful) == 0L) {
    stop("All bootstrap iterations failed. Check your data.", call. = FALSE)
  }
  actual_iterations <- length(successful)

  # Tag with iteration index and stack
  iter_dfs <- lapply(seq_along(successful), function(i) {
    df <- successful[[i]]
    df$iteration <- i
    df[, c("iteration", "Item", "item_restscore", "diff", "diff_abs")]
  })
  fit_all <- do.call(rbind, iter_dfs)
  rownames(fit_all) <- NULL

  # --- Return raw per-iteration long data ------------------------------------
  if (output == "raw") {
    return(fit_all)
  }

  # --- Per-item classification counts ----------------------------------------
  classes <- c("overfit", "underfit", "no misfit")
  combos <- expand.grid(
    Item           = item_names,
    item_restscore = classes,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  counts <- as.data.frame(
    table(Item = factor(fit_all$Item, levels = item_names),
          item_restscore = factor(fit_all$item_restscore, levels = classes)),
    responseName = "n",
    stringsAsFactors = FALSE
  )
  counts$Item           <- as.character(counts$Item)
  counts$item_restscore <- as.character(counts$item_restscore)

  fit_tbl <- merge(combos, counts,
                   by = c("Item", "item_restscore"),
                   all.x = TRUE, sort = FALSE)
  fit_tbl$n[is.na(fit_tbl$n)] <- 0L
  per_item_total <- tapply(fit_tbl$n, fit_tbl$Item, sum)
  fit_tbl$percent <- round(fit_tbl$n * 100 / per_item_total[fit_tbl$Item], 1)

  cfit_df <- data.frame(
    Item              = item_names,
    Infit_MSQ         = round(cfit$Infit, 2),
    Relative_location = round(relative_item_avg_locations, 2),
    stringsAsFactors  = FALSE,
    row.names         = NULL
  )

  result_df <- merge(fit_tbl, cfit_df, by = "Item", sort = FALSE)
  result_df <- result_df[order(match(result_df$Item, item_names),
                               result_df$item_restscore), ]
  rownames(result_df) <- NULL

  # --- Return per-item summary dataframe -------------------------------------
  if (output == "dataframe") {
    return(result_df)
  }

  # --- Kable: only flagged items above cutoff --------------------------------
  show <- result_df[result_df$item_restscore != "no misfit" &
                      result_df$percent > cutoff, , drop = FALSE]

  if (nrow(show) == 0L) {
    message(paste0("No item indicates misfit in more than ", cutoff,
                   "% of iterations."))
    return(invisible(NULL))
  }

  show <- show[order(show$item_restscore, -show$percent), ]
  rownames(show) <- NULL
  show <- show[, c("Item", "item_restscore", "percent",
                   "Infit_MSQ", "Relative_location")]

  knitr::kable(
    show,
    format    = "pipe",
    col.names = c("Item", "Item-restscore result", "% of iterations",
                  "Conditional MSQ infit", "Relative item location"),
    caption   = paste0(
      "Results based on ", actual_iterations,
      " successful bootstrap iterations with n = ", samplesize,
      " and ", n_items, " items. ",
      "Conditional mean-square infit based on complete responders only (n = ",
      n_complete, ")."
    )
  )
}

# ---------------------------------------------------------------------------
# Internal: single bootstrap iteration
# ---------------------------------------------------------------------------

#' Run a single item-restscore bootstrap iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMitemRestscoreBoot()] containing `data`,
#'   `samplesize`, `is_polytomous`, `item_names`.
#' @return A data.frame with columns `Item`, `item_restscore`, `diff`,
#'   `diff_abs`, or a character string on failure.
#' @keywords internal
run_single_boot_restscore <- function(seed, data_list) {
  set.seed(seed)
  idx <- sample.int(nrow(data_list$data), data_list$samplesize, replace = TRUE)
  d   <- data_list$data[idx, , drop = FALSE]

  tryCatch({
    if (data_list$is_polytomous) {
      model_fit <- psychotools::pcmodel(d, hessian = FALSE)
    } else {
      model_fit <- eRm::RM(d, se = FALSE)
    }

    i1 <- as.data.frame(iarm::item_restscore(model_fit))
    res_mat <- i1[[1L]]
    n_items <- length(data_list$item_names)

    observed <- as.numeric(res_mat[seq_len(n_items), 1L])
    expected <- as.numeric(res_mat[seq_len(n_items), 2L])
    p_adj    <- as.numeric(res_mat[seq_len(n_items), 5L])
    diff_val <- expected - observed

    cls <- ifelse(p_adj < 0.05 & diff_val < 0, "overfit",
                  ifelse(p_adj < 0.05 & diff_val > 0, "underfit", "no misfit"))

    data.frame(
      Item           = data_list$item_names,
      item_restscore = cls,
      diff           = diff_val,
      diff_abs       = abs(diff_val),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }, error = function(e) as.character(conditionMessage(e)))
}

# ---------------------------------------------------------------------------
# Internal: parallel runner (mirai)
# ---------------------------------------------------------------------------

#' Run item-restscore bootstrap iterations in parallel using mirai
#'
#' @param iterations Number of iterations.
#' @param boot_seeds Integer vector of per-iteration seeds.
#' @param boot_data_list List of data passed to each worker.
#' @param n_cores Number of mirai daemons.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_boot_restscore_parallel <- function(iterations, boot_seeds, boot_data_list,
                                        n_cores, verbose = FALSE) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  tasks <- lapply(seq_len(iterations), function(i) {
    mirai::mirai(
      { run_single_boot_restscore(seed, data_list) },
      seed                      = boot_seeds[i],
      data_list                 = boot_data_list,
      run_single_boot_restscore = run_single_boot_restscore
    )
  })

  results <- vector("list", iterations)
  for (i in seq_len(iterations)) {
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

#' Run item-restscore bootstrap iterations sequentially
#'
#' @param iterations Number of iterations.
#' @param boot_seeds Integer vector of per-iteration seeds.
#' @param boot_data_list List of data passed to each worker.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_boot_restscore_sequential <- function(iterations, boot_seeds, boot_data_list,
                                          verbose = FALSE) {
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)

  results <- vector("list", iterations)
  for (i in seq_len(iterations)) {
    results[[i]] <- run_single_boot_restscore(boot_seeds[i], boot_data_list)
    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}
