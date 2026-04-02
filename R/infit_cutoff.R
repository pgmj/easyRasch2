#' Simulation-Based Infit MSQ Cutoff Determination
#'
#' Uses parametric bootstrap simulation to determine appropriate cutoff values
#' for [RMiteminfit()]. Under a correctly fitting Rasch model, infit MSQ
#' statistics have a known distribution; this function simulates that
#' distribution and returns per-item empirical cutoffs.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Only complete cases (rows without
#'   any `NA`) are used.
#' @param iterations Integer. Number of simulation iterations (default 250).
#' @param parallel Logical. Use parallel processing via `mirai` if available
#'   (default `TRUE`).
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first. If neither is set and
#'   `parallel = TRUE`, a warning is issued and execution falls back to
#'   sequential (single core) processing.
#' @param verbose Logical. Show a progress bar (default `FALSE`).
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param cutoff_method Character string specifying how cutoff intervals are
#'   computed. Either `"hdi"` (default) for the Highest Density Interval via
#'   `ggdist::hdi()`, or `"quantile"` for the 2.5th/97.5th percentiles via
#'   `stats::quantile()`.
#' @param hdi_width Numeric. Width of the HDI when `cutoff_method = "hdi"`.
#'   Default is `0.999` (99.9% HDI). Ignored when `cutoff_method = "quantile"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`results`}{data.frame with columns `iteration`, `Item`,
#'     `InfitMSQ`, `OutfitMSQ` (one row per item per successful iteration).}
#'   \item{`item_cutoffs`}{data.frame with per-item cutoff summaries: `Item`,
#'     `infit_low`, `infit_high`, `outfit_low`, `outfit_high`. Bounds are
#'     computed using the method specified by `cutoff_method`.}
#'   \item{`actual_iterations`}{Number of successful iterations.}
#'   \item{`sample_n`}{Number of complete cases used.}
#'   \item{`sample_summary`}{Summary statistics of estimated person parameters.}
#'   \item{`item_names`}{Character vector of item names from data.}
#'   \item{`cutoff_method`}{The method used to compute cutoffs (`"hdi"` or
#'     `"quantile"`).}
#'   \item{`hdi_width`}{The HDI width used (only meaningful when
#'     `cutoff_method = "hdi"`).}
#' }
#'
#' @details
#' For each simulation iteration, person parameters (thetas) are resampled
#' with replacement from CML estimates, response data are simulated under the
#' Rasch model, the model is refitted, and conditional infit and outfit MSQ
#' statistics are computed via `iarm::out_infit()`. The distribution of these
#' statistics across iterations provides empirical critical values per item.
#' Failed iterations (e.g., due to convergence issues or degenerate data) are
#' silently discarded.
#'
#' Supports both **dichotomous** data (via `eRm::RM()` and
#' `psychotools::rrm()`) and **polytomous** data (via `eRm::PCM()` and an
#' internal partial credit score simulator).
#'
#' Parallel processing is provided by the `mirai` package (optional). Install
#' it with `install.packages("mirai")` to enable parallelisation.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' @seealso [RMiteminfit()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#' )
#' colnames(sim_data) <- paste0("Item", 1:10)
#'
#' # Run 100 iterations sequentially for a quick demo
#' cutoff_res <- RMinfitcutoff(sim_data, iterations = 100, parallel = FALSE,
#'                             seed = 42)
#' cutoff_res$item_cutoffs
#'
#' # Use the cutoffs in RMiteminfit()
#' RMiteminfit(sim_data)
#' }
RMinfitcutoff <- function(data, iterations = 250, parallel = TRUE,
                           n_cores = NULL, verbose = FALSE, seed = NULL,
                           cutoff_method = "hdi", hdi_width = 0.999) {

  cutoff_method <- match.arg(cutoff_method, c("hdi", "quantile"))

  if (cutoff_method == "hdi" && !requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required when cutoff_method = \"hdi\" but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")\n",
      "Alternatively, use cutoff_method = \"quantile\" to avoid this dependency.",
      call. = FALSE
    )
  }

  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMinfitcutoff() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  validate_response_data(data)

  # rgl workaround
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # Only complete cases (matching RMiteminfit behaviour)
  data <- stats::na.omit(data)
  if (nrow(data) == 0L) {
    stop("No complete cases in data. All rows contain at least one NA.",
         call. = FALSE)
  }

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
          "Your computer appears to have ", parallel::detectCores(), " cores available.\n",
          "Falling back to sequential (single core) processing."
        ),
        call. = FALSE
      )
      use_parallel <- FALSE
    } else {
      n_cores <- min(n_cores, iterations)
    }
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate per-iteration seeds
  sim_seeds <- sample.int(.Machine$integer.max, iterations)

  data_mat <- as.matrix(data)
  sample_n <- nrow(data_mat)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  item_names_vec <- colnames(data_mat)

  if (is_polytomous) {
    pcm_fit <- eRm::PCM(data_mat)
    pp <- eRm::person.parameter(pcm_fit)
    thetas <- as.numeric(stats::na.omit(unlist(pp$thetapar)))
    thresh_mat <- extract_item_thresholds(data_mat)
    deltaslist <- lapply(seq_len(nrow(thresh_mat)), function(i) {
      as.numeric(thresh_mat[i, !is.na(thresh_mat[i, ])])
    })
    sim_data_list <- list(
      type = "polytomous",
      thetas = thetas,
      deltaslist = deltaslist,
      n_items = ncol(data_mat),
      sample_n = sample_n,
      item_names = item_names_vec
    )
  } else {
    rm_fit <- eRm::RM(data_mat)
    pp <- eRm::person.parameter(rm_fit)
    thetas <- as.numeric(stats::na.omit(unlist(pp$thetapar)))
    item_params <- -rm_fit$betapar
    sim_data_list <- list(
      type = "dichotomous",
      thetas = thetas,
      item_params = item_params,
      n_items = ncol(data_mat),
      sample_n = sample_n,
      item_names = item_names_vec
    )
  }

  if (use_parallel) {
    results_raw <- run_infit_sim_parallel(iterations, sim_seeds, sim_data_list, n_cores, verbose)
  } else {
    results_raw <- run_infit_sim_sequential(iterations, sim_seeds, sim_data_list, verbose)
  }

  # Filter out failures (character strings indicate errors)
  ok <- vapply(results_raw, is.list, logical(1L))
  successful <- results_raw[ok]

  if (length(successful) == 0L) {
    stop("All simulation iterations failed. Check your data.", call. = FALSE)
  }

  actual_iterations <- length(successful)

  # Combine per-iteration data.frames
  iter_dfs <- lapply(seq_along(successful), function(i) {
    df <- successful[[i]]
    df$iteration <- i
    df
  })
  results_df <- do.call(rbind, iter_dfs)
  rownames(results_df) <- NULL

  # Compute per-item cutoffs
  item_names <- unique(results_df$Item)
  item_cutoffs <- do.call(rbind, lapply(item_names, function(item) {
    sub <- results_df[results_df$Item == item, ]
    if (cutoff_method == "hdi") {
      # ggdist::hdi() returns a matrix with ncol = 2: column 1 is the lower
      # bound, column 2 is the upper bound. Row 1 contains the widest/only
      # interval. For multimodal distributions `ggdist::hdi()` may return
      # multiple rows (disjoint intervals); we take only row 1 (the first
      # interval). In practice, infit/outfit MSQ distributions from parametric
      # bootstrap are unimodal, so a single interval is expected.
      infit_interval  <- ggdist::hdi(sub$InfitMSQ,  .width = hdi_width)
      outfit_interval <- ggdist::hdi(sub$OutfitMSQ, .width = hdi_width)
      data.frame(
        Item        = item,
        infit_low   = infit_interval[1L, 1L],
        infit_high  = infit_interval[1L, 2L],
        outfit_low  = outfit_interval[1L, 1L],
        outfit_high = outfit_interval[1L, 2L],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    } else {
      data.frame(
        Item        = item,
        infit_low   = stats::quantile(sub$InfitMSQ,  0.025, na.rm = TRUE),
        infit_high  = stats::quantile(sub$InfitMSQ,  0.975, na.rm = TRUE),
        outfit_low  = stats::quantile(sub$OutfitMSQ, 0.025, na.rm = TRUE),
        outfit_high = stats::quantile(sub$OutfitMSQ, 0.975, na.rm = TRUE),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }))
  rownames(item_cutoffs) <- NULL

  list(
    results           = results_df,
    item_cutoffs      = item_cutoffs,
    actual_iterations = actual_iterations,
    sample_n          = sample_n,
    sample_summary    = summary(thetas),
    item_names        = item_names_vec,
    cutoff_method     = cutoff_method,
    hdi_width         = hdi_width
  )
}

# ---------------------------------------------------------------------------
# Internal: single simulation iteration
# ---------------------------------------------------------------------------

#' Run a single infit simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMinfitcutoff()].
#' @return A data.frame with columns `Item`, `InfitMSQ`, `OutfitMSQ`, or a
#'   character string on failure.
#' @keywords internal
run_single_infit_sim <- function(seed, data_list) {
  set.seed(seed)

  thetas_res <- sample(data_list$thetas, size = data_list$sample_n, replace = TRUE)

  tryCatch({
    if (data_list$type == "dichotomous") {
      sim_mat <- psychotools::rrm(
        theta = thetas_res,
        delta = data_list$item_params
      )
      sim_df <- as.data.frame(sim_mat)
      colnames(sim_df) <- data_list$item_names

      pos_counts <- colSums(sim_df, na.rm = TRUE)
      # Fewer than 8 positive responses per item can cause numerical
      # instability in conditional MLE (eRm) fitting.
      if (any(pos_counts < 8L)) {
        return("validation_failed: fewer than 8 positive responses in at least one item")
      }

      model_fit <- eRm::RM(sim_df, se = FALSE)
    } else {
      sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
      sim_df <- as.data.frame(sim_mat)
      colnames(sim_df) <- data_list$item_names

      n_cats <- vapply(data_list$deltaslist, function(d) length(d) + 1L, integer(1L))
      for (j in seq_len(ncol(sim_df))) {
        tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
        if (any(tab == 0L)) {
          return("validation_failed: not all categories represented")
        }
      }

      model_fit <- psychotools::pcmodel(sim_df, hessian = FALSE)
    }

    cfit <- iarm::out_infit(model_fit)

    data.frame(
      Item      = data_list$item_names,
      InfitMSQ  = round(cfit$Infit,  3),
      OutfitMSQ = round(cfit$Outfit, 3),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }, error = function(e) {
    as.character(conditionMessage(e))
  })
}

# ---------------------------------------------------------------------------
# Internal: parallel runner
# ---------------------------------------------------------------------------

#' Run infit simulations in parallel using mirai
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param n_cores Number of mirai daemons.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_infit_sim_parallel <- function(iterations, sim_seeds, sim_data_list,
                                   n_cores, verbose = FALSE) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
    completed <- 0L
  }

  # Submit all tasks
  tasks <- lapply(seq_len(iterations), function(sim) {
    mirai::mirai(
      {
        run_single_infit_sim(seed, data_list)
      },
      seed = sim_seeds[sim],
      data_list = sim_data_list,
      run_single_infit_sim = run_single_infit_sim,
      sim_partial_score = sim_partial_score,
      sim_poly_item = sim_poly_item
    )
  })

  # Collect results
  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    result <- mirai::call_mirai(tasks[[sim]])$data
    if (!inherits(result, "errorValue")) {
      results[[sim]] <- result
    } else {
      results[[sim]] <- "mirai_error"
    }
    if (verbose) {
      completed <- completed + 1L
      utils::setTxtProgressBar(pb, completed)
    }
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

#' Run infit simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_infit_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                     verbose = FALSE) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_infit_sim(sim_seeds[sim], sim_data_list)
    if (verbose) {
      utils::setTxtProgressBar(pb, sim)
    }
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}
