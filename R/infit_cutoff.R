#' Simulation-Based Infit MSQ Cutoff Determination
#'
#' Uses parametric bootstrap simulation to determine appropriate cutoff values
#' for \code{\link{RMitemInfit}}. This function simulates data from a correctly fitting 
#' Rasch model that mimics your data and returns per-item empirical cutoffs.
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
#'   computed. Either `"hdci"` (default) for the Highest Density Interval via
#'   `ggdist::hdci()`, or `"quantile"` for the 2.5th/97.5th percentiles via
#'   `stats::quantile()`.
#' @param hdci_width Numeric. Width of the HDCI when `cutoff_method = "hdci"`.
#'   Default is `0.999` (99.9% HDCI). Ignored when `cutoff_method = "quantile"`.
#' @param dgp Character. Data-generating process for the parametric bootstrap.
#'   `"resample"` (default) resamples WLE person locations with replacement and
#'   simulates responses under the model (a *marginal* null). `"conditional"`
#'   simulates each respondent's pattern from the exact Rasch conditional
#'   distribution given their observed total score, item parameters fixed (a
#'   *conditional* null). Because the conditional infit/outfit statistic is
#'   itself conditional on the total score, `"conditional"` is its naturally
#'   matched null. \strong{Experimental.}
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
#'   \item{`cutoff_method`}{The method used to compute cutoffs (`"hdci"` or
#'     `"quantile"`).}
#'   \item{`hdci_width`}{The HDCI width used (only meaningful when
#'     `cutoff_method = "hdci"`).}
#'   \item{`dgp`}{The data-generating process used (`"resample"` or
#'     `"conditional"`).}
#' }
#'
#' @details
#' The generating model is CML item parameters (via `psychotools`) with WLE
#' person locations. For each iteration a dataset is simulated under the chosen
#' `dgp`, the model is refitted by CML (`psychotools::pcmodel()`, which handles
#' dichotomous and polytomous data and is accepted by `iarm`), and conditional
#' infit and outfit MSQ are computed via `iarm::out_infit()`. The distribution
#' of these statistics across iterations provides empirical critical values per
#' item. Failed iterations (e.g., degenerate simulated data) are silently
#' discarded.
#'
#' Parallel processing is provided by the `mirai` package (optional). Install
#' it with `install.packages("mirai")` to enable parallelisation.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' @seealso \code{\link{RMitemInfit}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#' )
#' colnames(sim_data) <- paste0("Item", 1:10)
#'
#' # Run 100 iterations sequentially for a quick demo
#' cutoff_res <- RMitemInfitCutoff(sim_data, iterations = 100, parallel = FALSE,
#'                             seed = 42)
#' cutoff_res$item_cutoffs
#'
#' # Use the cutoffs in RMitemInfit()
#' RMitemInfit(sim_data)
#' }
RMitemInfitCutoff <- function(data, iterations = 250, parallel = TRUE,
                           n_cores = NULL, verbose = FALSE, seed = NULL,
                           cutoff_method = "hdci", hdci_width = 0.999,
                           dgp = c("resample", "conditional")) {

  dgp           <- match.arg(dgp)
  cutoff_method <- match.arg(cutoff_method, c("hdci", "quantile"))

  if (cutoff_method == "hdci" && !requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required when cutoff_method = \"hdci\" but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")\n",
      "Alternatively, use cutoff_method = \"quantile\" to avoid this dependency.",
      call. = FALSE
    )
  }

  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMitemInfitCutoff() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  validate_response_data(data)

  # rgl workaround
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # Only complete cases (matching RMitemInfit behaviour)
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

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate per-iteration seeds
  sim_seeds <- sample.int(.Machine$integer.max, iterations)

  data_mat <- as.matrix(data)
  sample_n <- nrow(data_mat)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  item_names_vec <- colnames(data_mat)

  # Generating model: CML item thresholds (psychotools), computed once. The
  # conditional infit statistic is conditional on the total score, so the
  # "conditional" DGP (simulate each respondent's pattern given their observed
  # score) is the matched null. The "resample" DGP draws WLE person locations
  # with replacement and simulates parametrically (a marginal null).
  thr_list   <- .fit_cml_thresholds(data_mat)
  wle_thetas <- .estimate_thetas(data_mat, thr_list, method = "WLE")$theta
  wle_thetas <- wle_thetas[is.finite(wle_thetas)]

  sim_data_list <- list(
    dgp        = dgp,
    type       = if (is_polytomous) "polytomous" else "dichotomous",
    thr_list   = thr_list,
    n_items    = ncol(data_mat),
    sample_n   = sample_n,
    item_names = item_names_vec
  )
  if (dgp == "resample") {
    sim_data_list$thetas <- wle_thetas
    if (is_polytomous) {
      sim_data_list$deltaslist <- thr_list
    } else {
      sim_data_list$item_params <- unlist(thr_list, use.names = FALSE)
    }
  } else {
    sim_data_list$cond_groups <- .cond_groups(data_mat, thr_list)
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
    if (cutoff_method == "hdci") {
      # ggdist::hdci() returns a matrix with ncol = 2: column 1 is the lower
      # bound, column 2 is the upper bound. Row 1 contains the continuous
      # interval.
      infit_interval  <- ggdist::hdci(sub$InfitMSQ,  .width = hdci_width)
      outfit_interval <- ggdist::hdci(sub$OutfitMSQ, .width = hdci_width)
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
    sample_summary    = summary(wle_thetas),
    item_names        = item_names_vec,
    cutoff_method     = cutoff_method,
    hdci_width         = hdci_width,
    dgp               = dgp
  )
}

# ---------------------------------------------------------------------------
# Internal: single simulation iteration
# ---------------------------------------------------------------------------

#' Run a single infit simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMitemInfitCutoff()].
#' @return A data.frame with columns `Item`, `InfitMSQ`, `OutfitMSQ`, or a
#'   character string on failure.
#' @keywords internal
run_single_infit_sim <- function(seed, data_list) {
  set.seed(seed)

  tryCatch({
    # --- Generate one simulated dataset under the chosen DGP -----------------
    if (identical(data_list$dgp, "conditional")) {
      sim_df <- .sim_cond_dataset(data_list)
    } else if (data_list$type == "dichotomous") {
      thetas_res <- sample(data_list$thetas, size = data_list$sample_n,
                           replace = TRUE)
      sim_df <- as.data.frame(
        psychotools::rrm(theta = thetas_res, beta = data_list$item_params)$data
      )
    } else {
      thetas_res <- sample(data_list$thetas, size = data_list$sample_n,
                           replace = TRUE)
      sim_df <- as.data.frame(sim_partial_score(data_list$deltaslist, thetas_res))
    }
    colnames(sim_df) <- data_list$item_names

    # --- Validate the simulated dataset (estimable refit) --------------------
    if (data_list$type == "dichotomous") {
      if (any(colSums(sim_df, na.rm = TRUE) < 8L)) {
        return("validation_failed: fewer than 8 positive responses in at least one item")
      }
    } else {
      n_cats <- vapply(data_list$thr_list, function(d) length(d) + 1L, integer(1L))
      for (j in seq_len(ncol(sim_df))) {
        tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
        if (any(tab == 0L)) {
          return("validation_failed: not all categories represented")
        }
      }
    }

    # Conditional infit/outfit from a CML refit. psychotools::pcmodel() handles
    # both polytomous and dichotomous (a 2-category PCM is the Rasch model) and
    # is accepted by iarm::out_infit(); it matches eRm to ~1e-6 but is faster.
    model_fit <- psychotools::pcmodel(sim_df)
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
      sim_poly_item = sim_poly_item,
      # Conditional-DGP generators (shared with the Q3 cutoff).
      .sim_cond_dataset = .sim_cond_dataset,
      .sim_conditional  = .sim_conditional,
      .esf_convolve     = .esf_convolve
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
