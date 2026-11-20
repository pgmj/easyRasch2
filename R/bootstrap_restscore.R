#' Bootstrap flag rates of the asymptotic item-restscore test
#'
#' Draws bootstrap samples of size `samplesize` from `data`, runs the
#' asymptotic item-restscore test of `iarm::item_restscore()` in each, and
#' reports how often each item is flagged as overfit or underfit. This
#' describes how stable the asymptotic flags are across samples of a given
#' size. It is not a calibrated test of item fit. For that, use
#' [RMitemRestscoreCutoff()] with [RMitemRestscore()].
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0.
#' @param iterations Integer. Number of bootstrap samples. Default `200`.
#' @param samplesize Integer. Size of each bootstrap sample. Default `600`.
#'   Must not exceed `nrow(data)`.
#' @param parallel Logical. Run the bootstrap in parallel via \pkg{mirai}
#'   when it is installed. Default `TRUE`.
#' @param n_cores Integer or `NULL`. Number of parallel workers. When
#'   `NULL`, `getOption("mc.cores")` is used. If that is unset too, a
#'   warning is given and the bootstrap runs sequentially.
#' @param cutoff Numeric. The kable shows only items flagged in more than
#'   `cutoff` percent of the samples. Default `5`. Not used by the other
#'   outputs.
#' @param verbose Logical. Show a progress bar. Default `FALSE`.
#' @param seed Integer or `NULL`. Random seed for reproducibility. See
#'   [easyRasch2-reproducibility] for what this guarantees and how it
#'   interacts with `parallel`.
#' @param output Character. `"kable"` (default) for a `knitr::kable()`
#'   table, `"dataframe"` for the summary per item, or `"raw"` for one row
#'   per item and bootstrap sample (for custom plots).
#'
#' @return
#' * `output = "kable"`: a `knitr_kable` table of the items flagged in more
#'   than `cutoff` percent of the samples, with the columns Item,
#'   Item-restscore result, % of iterations, Conditional MSQ infit and
#'   Relative item location. The caption gives the number of samples, the
#'   sample size and the number of complete cases.
#' * `output = "dataframe"`: a data.frame with one row per item and result
#'   (`"overfit"`, `"underfit"`, `"no misfit"`), with the columns `Item`,
#'   `item_restscore`, `n`, `n_tested`, `percent`, `Infit_MSQ` and
#'   `Relative_location`. `percent` is `n` out of `n_tested`, the samples
#'   in which the item could be tested, and `NA` when it never could.
#' * `output = "raw"`: a data.frame with one row per item and bootstrap
#'   sample, with the columns `iteration`, `Item`, `item_restscore`, `diff`
#'   and `diff_abs`. Note that `diff` is expected minus observed gamma, the
#'   opposite sign of `Difference` in [RMitemRestscore()]. Rows for an item
#'   that could not be tested in that sample have `NA` in
#'   `item_restscore`, `diff` and `diff_abs`.
#'
#' @details
#' **What the percentage means.** Each sample is tested with the asymptotic
#' test of `iarm::item_restscore()` with Benjamini-Hochberg adjustment, and
#' an item is flagged when its adjusted p is below .05, as overfit or
#' underfit by the sign of the difference. The result inherits that test's
#' calibration at `samplesize`. With the Rasch model true, it flags too many
#' items as overfit in small or mistargeted samples of dichotomous items,
#' too few items of any kind in large polytomous samples, and too few
#' underfitting items at every sample size (see [RMitemRestscore()] and
#' \url{https://pgmj.github.io/rasch_book/chapters/restscore-asymptotic-null.html}).
#' The percentage is therefore how often that test flags an item in
#' samples of size `samplesize`, not the probability that the item misfits,
#' and it changes with `samplesize`. A low underfit percentage is weak
#' evidence of fit.
#'
#' **Model.** The full sample is fitted by CML with
#' `psychotools::pcmodel()` (a dichotomous item is a 2-category partial
#' credit model). The conditional infit MSQ (`iarm::out_infit()`, on
#' complete cases) and the relative item location (the mean of the item's
#' CML thresholds minus the mean WLE person location) are computed once on
#' the full sample and shown for context. Each bootstrap sample is refitted
#' by CML.
#'
#' **Samples that fail.** A sample in which the model cannot be fitted is
#' left out, and the caption and the counts include only the samples that
#' could be used. A polytomous item that is constant in a sample makes that
#' sample fail. A dichotomous item that is constant, or a polytomous item
#' with an empty middle category, cannot be tested in that sample. It is
#' then left out of that item's percentages, the other items keep the
#' sample, and the caption lists the items with fewer samples. When a
#' constant dichotomous item is left out, the Benjamini-Hochberg adjustment
#' in that sample covers only the items tested there, since iarm drops the
#' item before testing. An item with an empty lowest category is shifted to
#' start at 0 by psychotools and is still tested.
#'
#' Parallel processing needs the optional \pkg{mirai} package. The
#' \pkg{iarm} package must be installed (it is in Suggests).
#'
#' @references
#' Kreiner, S. (2011). A Note on Item-Restscore Association in Rasch Models.
#' *Applied Psychological Measurement, 35*(7), 557-561.
#' \doi{10.1177/0146621611410227}
#'
#' @seealso [RMitemRestscore()] and [RMitemRestscoreCutoff()] for a
#'   calibrated test, [RMitemInfit()].
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("iarm", quietly = TRUE)) {
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 400 * 8, replace = TRUE), nrow = 400, ncol = 8)
#' )
#' colnames(sim_data) <- paste0("Item", 1:8)
#'
#' # Few iterations for a fast example. Use 100 or more in real analyses.
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
#'     geom_jitter(aes(colour = item_restscore, shape = item_restscore),
#'                 width = 0.15, alpha = 0.5, size = 0.8) +
#'     scale_colour_manual(values = c("overfit"   = "#0072B2",
#'                                    "underfit"  = "#D55E00",
#'                                    "no misfit" = "grey60")) +
#'     labs(y = "Expected minus observed item-restscore gamma",
#'          colour = NULL, shape = NULL) +
#'     theme_bw() +
#'     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#' }
#' }
#' }
RMitemRestscoreBoot <- function(
  data,
  iterations = 200,
  samplesize = 600,
  parallel = TRUE,
  n_cores = NULL,
  cutoff = 5,
  verbose = FALSE,
  seed = NULL,
  output = "kable"
) {
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMitemRestscoreBoot() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe", "raw"))

  validate_response_data(data)

  # Respondents with no responses at all contribute nothing, break the CML
  # fit (psychotools errors on all-NA rows), and would poison the bootstrap
  # resampling pool; drop them, keeping the raw totals for the caption.
  n_total <- nrow(as.data.frame(data))
  has_na <- anyNA(data)
  data <- .drop_empty_respondents(data)

  if (samplesize > nrow(data)) {
    stop(
      paste0(
        "`samplesize` (",
        samplesize,
        ") cannot be larger than the number of rows in `data` (",
        nrow(data),
        ")."
      ),
      call. = FALSE
    )
  }

  # rgl workaround for iarm dependency chain (vcdExtra -> rgl)
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  data_mat <- as.matrix(data)
  item_names <- colnames(data_mat)
  n_items <- ncol(data_mat)

  # --- Full-sample model: locations + conditional infit ----------------------
  # CML item parameters (psychotools; a dichotomous item is a 2-category PCM)
  # and WLE person locations, consistent with RMitemRestscore() and the rest of
  # the package. The item-restscore and infit statistics from iarm are
  # conditional and engine-invariant; only the relative-location reference
  # shifts slightly (WLE vs eRm MLE person mean).
  # iarm's statistics refit on complete cases, so check them too.
  .check_complete_case_categories(data)
  fit_full <- .cml_fit(data) # shared category check, then pcmodel()

  thr_list <- .center_thresholds(lapply(
    psychotools::threshpar(fit_full),
    as.numeric
  ))
  item_avg_locations <- vapply(thr_list, mean, numeric(1L))
  person_avg_location <- mean(
    .estimate_thetas(data_mat, thr_list, method = "WLE")$theta,
    na.rm = TRUE
  )
  relative_item_avg_locations <- item_avg_locations - person_avg_location

  cfit <- iarm::out_infit(fit_full)
  n_complete <- nrow(stats::na.omit(data))

  # --- Parallel setup --------------------------------------------------------
  use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)
  if (parallel && !use_parallel) {
    message(
      "Install 'mirai' package for parallel processing: install.packages(\"mirai\")"
    )
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
    data = data,
    samplesize = samplesize,
    item_names = item_names
  )

  # --- Run bootstrap ---------------------------------------------------------
  if (use_parallel) {
    results_raw <- run_boot_restscore_parallel(
      iterations,
      boot_seeds,
      boot_data_list,
      n_cores,
      verbose
    )
  } else {
    results_raw <- run_boot_restscore_sequential(
      iterations,
      boot_seeds,
      boot_data_list,
      verbose
    )
  }

  ok <- vapply(results_raw, is.data.frame, logical(1L))
  successful <- results_raw[ok]

  if (length(successful) == 0L) {
    stop("All bootstrap iterations failed. Check your data.", call. = FALSE)
  }
  actual_iterations <- length(successful)

  # Shown once at least one bootstrap iteration succeeded, so a call that
  # stops earlier (data checks, all iterations failed) does not spend the
  # once-per-session warning, and before any output, so it also reaches the
  # dataframe/raw outputs and the no-flags path, where a low underfit rate
  # reads as evidence of fit.
  rlang::warn(
    c(
      paste0(
        "RMitemRestscoreBoot() flags items with the asymptotic ",
        "item-restscore test, which is miscalibrated under the Rasch model."
      ),
      i = paste0(
        "It flags too many items as overfit in small or mistargeted ",
        "dichotomous samples and too few as underfit at any sample size."
      ),
      i = "The percentages are flag rates at `samplesize`, not probabilities of misfit.",
      i = "For a calibrated test, use RMitemRestscoreCutoff() with RMitemRestscore()."
    ),
    class = "easyRasch2_asymptotic_restscore_warning",
    .frequency = "once",
    .frequency_id = "easyRasch2_restscoreboot_asymptotic"
  )

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
    Item = item_names,
    item_restscore = classes,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  counts <- as.data.frame(
    table(
      Item = factor(fit_all$Item, levels = item_names),
      item_restscore = factor(fit_all$item_restscore, levels = classes)
    ),
    responseName = "n",
    stringsAsFactors = FALSE
  )
  counts$Item <- as.character(counts$Item)
  counts$item_restscore <- as.character(counts$item_restscore)

  fit_tbl <- merge(
    combos,
    counts,
    by = c("Item", "item_restscore"),
    all.x = TRUE,
    sort = FALSE
  )
  fit_tbl$n[is.na(fit_tbl$n)] <- 0L
  per_item_total <- tapply(fit_tbl$n, fit_tbl$Item, sum)
  # Each item's base is the resamples in which it could be tested: a
  # dichotomous item that is constant, or a polytomous item with an empty
  # middle category, gets NA there and is left out. (A constant polytomous
  # item makes iarm fail, so the whole resample is dropped; an empty lowest
  # category is shifted by psychotools and still tested.) An item never
  # testable gets NA, not NaN.
  fit_tbl$n_tested <- as.integer(per_item_total[fit_tbl$Item])
  # Unrounded values; the kable path rounds a display copy before rendering.
  fit_tbl$percent <- ifelse(
    fit_tbl$n_tested > 0L,
    fit_tbl$n * 100 / fit_tbl$n_tested,
    NA_real_
  )
  short <- per_item_total[item_names] < actual_iterations
  short_note <- if (any(short)) {
    paste0(
      " Percentages are of the resamples in which an item could be tested,",
      " which was fewer for: ",
      paste0(item_names[short], " (", per_item_total[item_names][short], ")",
             collapse = ", "),
      "."
    )
  } else {
    ""
  }

  cfit_df <- data.frame(
    Item = item_names,
    Infit_MSQ = as.numeric(cfit$Infit),
    Relative_location = as.numeric(relative_item_avg_locations),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  result_df <- merge(fit_tbl, cfit_df, by = "Item", sort = FALSE)
  result_df <- result_df[
    order(match(result_df$Item, item_names), result_df$item_restscore),
  ]
  rownames(result_df) <- NULL

  # --- Return per-item summary dataframe -------------------------------------
  if (output == "dataframe") {
    return(result_df)
  }

  # Kable display rounding (the dataframe output above stays unrounded); done
  # before the percent filter below so the kable keeps its previous behaviour
  # (filtering on the displayed, rounded percentage).
  result_df <- .round_display(result_df, c(
    percent = 1, Infit_MSQ = 2, Relative_location = 2
  ))

  # --- Kable: only flagged items above cutoff --------------------------------
  show <- result_df[
    result_df$item_restscore != "no misfit" &
      !is.na(result_df$percent) &
      result_df$percent > cutoff,
    ,
    drop = FALSE
  ]

  if (nrow(show) == 0L) {
    message(paste0(
      "No item indicates misfit in more than ",
      cutoff,
      "% of iterations."
    ))
    return(invisible(NULL))
  }

  show <- show[order(show$item_restscore, -show$percent), ]
  rownames(show) <- NULL
  show <- show[, c(
    "Item",
    "item_restscore",
    "percent",
    "Infit_MSQ",
    "Relative_location"
  )]

  knitr::kable(
    show,
    format = "pipe",
    col.names = c(
      "Item",
      "Item-restscore result",
      "% of iterations",
      "Conditional MSQ infit",
      "Relative item location"
    ),
    caption = paste0(
      "Item-restscore flag rates from ",
      actual_iterations,
      " successful bootstrap iterations (resample size ",
      samplesize,
      ", ",
      n_items,
      " items). ",
      .n_caption(
        n_complete,
        n_total,
        if (has_na) "complete cases" else character()
      ),
      ". Conditional MSQ infit and relative item location are from the",
      " full sample. Items are flagged by the asymptotic test",
      " (BH-adjusted p < .05), which is miscalibrated under the Rasch model,",
      " so the percentages are how often that test flags an item at this",
      " resample size, not the probability that it misfits. Use",
      " RMitemRestscoreCutoff() with RMitemRestscore() for simulation-based",
      " p-values.",
      short_note
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
#'   `samplesize`, `item_names`.
#' @return A data.frame with columns `Item`, `item_restscore`, `diff`,
#'   `diff_abs`, or a character string on failure.
#' @noRd
run_single_boot_restscore <- function(seed, data_list) {
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
  idx <- sample.int(nrow(data_list$data), data_list$samplesize, replace = TRUE)
  d <- data_list$data[idx, , drop = FALSE]

  tryCatch(
    {
      # CML refit via psychotools::pcmodel() for both dichotomous (2-category PCM)
      # and polytomous data; accepted by iarm::item_restscore() and faster than
      # eRm. hessian = FALSE skips the (unused) SE computation. check = FALSE:
      # the full data were checked, and in a resample an item that loses a
      # category is reshaped by psychotools and handled below (NA or shifted).
      model_fit <- .cml_fit(d, hessian = FALSE, check = FALSE)

      i1 <- as.data.frame(iarm::item_restscore(model_fit))
      res_mat <- i1[[1L]]

      # Match rows by item name, not position: when a dichotomous item is
      # constant in this resample, psychotools drops it and iarm returns one
      # row fewer (a constant polytomous item makes iarm error instead).
      # That item gets NA here (left out of its own percentages), and the
      # other items keep this resample's results.
      idx <- if (is.null(rownames(res_mat))) {
        if (nrow(res_mat) != length(data_list$item_names)) {
          stop("iarm returned unnamed rows for a reduced item set.")
        }
        seq_len(nrow(res_mat))
      } else {
        match(data_list$item_names, rownames(res_mat))
      }
      observed <- as.numeric(res_mat[idx, 1L])
      expected <- as.numeric(res_mat[idx, 2L])
      p_adj <- as.numeric(res_mat[idx, 5L])
      diff_val <- expected - observed

      cls <- ifelse(
        p_adj < 0.05 & diff_val < 0,
        "overfit",
        ifelse(p_adj < 0.05 & diff_val > 0, "underfit", "no misfit")
      )

      data.frame(
        Item = data_list$item_names,
        item_restscore = cls,
        diff = diff_val,
        diff_abs = abs(diff_val),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    },
    error = function(e) as.character(conditionMessage(e))
  )
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
#' @noRd
run_boot_restscore_parallel <- function(
  iterations,
  boot_seeds,
  boot_data_list,
  n_cores,
  verbose = FALSE
) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  tasks <- lapply(seq_len(iterations), function(i) {
    mirai::mirai(
      {
        run_single_boot_restscore(seed, data_list)
      },
      seed = boot_seeds[i],
      data_list = boot_data_list,
      run_single_boot_restscore = run_single_boot_restscore,
      .cml_fit = .cml_fit
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
#' @noRd
run_boot_restscore_sequential <- function(
  iterations,
  boot_seeds,
  boot_data_list,
  verbose = FALSE
) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

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
