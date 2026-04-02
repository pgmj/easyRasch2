#' Q3 Residual Correlations for Local Dependence Assessment
#'
#' Computes Yen's Q3 residual correlations between item pairs using a Rasch
#' model fitted via Marginal Maximum Likelihood (MML) in `mirt`. High
#' correlations (above the dynamic cut-off) indicate potential local dependence
#' between items.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed.
#' @param cutoff Optional single numeric value added to the mean off-diagonal Q3
#'   correlation to produce the dynamic cut-off threshold. When `NULL`
#'   (default), the raw Q3 residual correlation matrix is returned without any
#'   dynamic cut-off applied. Use [RMlocdepQ3cutoff()] to derive a
#'   simulation-based cutoff value.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying numeric data.frame.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object showing the lower triangle
#'   of the Q3 correlation matrix. When `cutoff` is provided, a footnote
#'   describing the dynamic cut-off is included and an extra `above_cutoff`
#'   column marks rows containing at least one value above the threshold.
#' * If `output = "dataframe"`: a data.frame (rounded to 2 decimal places) with
#'   the lower triangle of the Q3 correlation matrix; the upper triangle and
#'   diagonal are set to `NA`. When `cutoff` is provided, an additional logical
#'   column `above_cutoff` indicates whether the row contains any value
#'   exceeding the dynamic cut-off.
#'
#' @details
#' The Q3 statistic (Yen, 1984) is the correlation between residuals of pairs
#' of items after accounting for the latent trait. Under local independence,
#' Q3 values are expected to be around \eqn{-1/(k-1)} where \eqn{k} is the
#' number of items. When `cutoff` is supplied, the dynamic cut-off is the mean
#' of all off-diagonal Q3 values plus `cutoff`, following the approach of
#' Christensen et al. (2017). Use [RMlocdepQ3cutoff()] to obtain a
#' simulation-based cutoff recommendation.
#'
#' `mirt` is used for model fitting here because Q3 requires model-based
#' expected responses, which are most readily available from MML estimation.
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
#' # Get the underlying data.frame
#' q3_df <- RMlocdepQ3(sim_data, output = "dataframe")
#'
#' # With a simulation-based cutoff
#' cutoff_res <- RMlocdepQ3cutoff(sim_data, iterations = 100)
#' RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
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

  # Compute mean of off-diagonal correlations before zeroing upper triangle
  mean_resid <- mean(resid_mat, na.rm = TRUE)

  resid_df <- as.data.frame(resid_mat)

  # Round values
  resid_df[] <- lapply(resid_df, function(x) round(x, 2))

  # Blank out upper triangle and diagonal (set to NA)
  resid_df[upper.tri(resid_df)] <- NA
  diag(resid_df) <- NA

  # --- Add above_cutoff column when cutoff is provided ------------------------
  if (!is.null(cutoff)) {
    dyn_cutoff <- mean_resid + cutoff
    # Check each row for any value exceeding the dynamic cut-off
    resid_df$above_cutoff <- apply(
      resid_df, 1,
      function(row) any(row > dyn_cutoff, na.rm = TRUE)
    )
  }

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(resid_df)
  }

  # Build caption
  if (!is.null(cutoff)) {
    caption_text <- paste0(
      "Dynamic cut-off: ", round(dyn_cutoff, 3),
      " (mean Q3 = ", round(mean_resid, 3),
      " + ", round(cutoff, 3), ").",
      " Correlations exceeding the cut-off may indicate local dependence."
    )
  } else {
    caption_text <- "Raw Q3 residual correlations (lower triangle). Use RMlocdepQ3cutoff() to derive a cutoff."
  }

  # Replace NA with "" for display (knitr::kable doesn't render NA nicely)
  # Determine which columns are the item columns (all except above_cutoff)
  item_cols <- setdiff(names(resid_df), "above_cutoff")
  resid_display <- as.data.frame(
    lapply(resid_df[item_cols], function(x) {
      ifelse(is.na(x), "", as.character(x))
    }),
    stringsAsFactors = FALSE
  )
  rownames(resid_display) <- rownames(resid_df)

  # Add above_cutoff column for kable display
  if (!is.null(cutoff)) {
    resid_display$above_cutoff <- ifelse(resid_df$above_cutoff, "*", "")
  }

  knitr::kable(
    resid_display,
    caption = caption_text,
    format = "pipe"
  )
}

#' Simulation-Based Q3 Cutoff Determination
#'
#' Uses parametric bootstrap simulation to determine an appropriate cutoff
#' value for [RMlocdepQ3()]. Under a correctly fitting Rasch model, Q3
#' residuals have an unknown distribution; this function simulates that
#' distribution and returns empirical percentiles.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers).
#' @param iterations Integer. Number of simulation iterations (default 500).
#' @param parallel Logical. Use parallel processing via `mirai` if available
#'   (default `TRUE`).
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first. If neither is set and
#'   `parallel = TRUE`, a warning is issued and execution falls back to
#'   sequential (single core) processing.
#' @param verbose Logical. Show a progress bar (default `FALSE`).
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{`results`}{data.frame with columns `mean`, `max`, `diff` (one row
#'     per successful iteration).}
#'   \item{`actual_iterations`}{Number of successful iterations.}
#'   \item{`sample_n`}{Number of persons in the original data.}
#'   \item{`sample_summary`}{Summary statistics of estimated person parameters.}
#'   \item{`max_diff`, `sd_diff`}{Max and SD of the `diff` distribution.}
#'   \item{`p95`, `p99`, `p995`, `p999`}{Empirical percentiles of `diff`.}
#'   \item{`suggested_cutoff`}{The 99th percentile (`p99`) — recommended cutoff
#'     for [RMlocdepQ3()].}
#' }
#'
#' @details
#' For each simulation iteration, person parameters (thetas) are resampled
#' with replacement from ML estimates, response data are simulated under the
#' Rasch model, a `mirt` Rasch model is fitted to the simulated data, and Q3
#' residuals are extracted. The distribution of `max(Q3) - mean(Q3)` across
#' iterations provides empirical critical values. Failed iterations (e.g., due
#' to convergence issues) are silently discarded.
#'
#' Supports both **dichotomous** data (via `eRm::RM()` and
#' `psychotools::rrm()`) and **polytomous** data (via `eRm::PCM()` and an
#' internal partial credit score simulator).
#'
#' Parallel processing is provided by the `mirai` package (optional). Install
#' it with `install.packages("mirai")` to enable parallelisation.
#'
#' @seealso [RMlocdepQ3()]
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
#' cutoff_res <- RMlocdepQ3cutoff(sim_data, iterations = 100, parallel = FALSE,
#'                                seed = 42)
#' cutoff_res$suggested_cutoff  # 99th percentile
#'
#' # Use the cutoff in RMlocdepQ3()
#' RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
#' }
RMlocdepQ3cutoff <- function(data, iterations = 500, parallel = TRUE,
                              n_cores = NULL, verbose = FALSE, seed = NULL) {
  validate_response_data(data)

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

  if (is_polytomous) {
    # Fit PCM, extract thresholds and thetas
    pcm_fit <- eRm::PCM(data_mat)
    pp <- eRm::person.parameter(pcm_fit)
    theta_table <- pp$theta.table[["Person Parameter"]]
    raw_scores <- rowSums(data_mat, na.rm = TRUE)
    thetas <- as.numeric(stats::na.omit(theta_table[raw_scores]))
    thresh_mat <- extract_item_thresholds(data_mat)
    # Convert threshold matrix rows to a list of threshold vectors (one per item)
    deltaslist <- lapply(seq_len(nrow(thresh_mat)), function(i) {
      as.numeric(thresh_mat[i, !is.na(thresh_mat[i, ])])
    })
    sim_data_list <- list(
      type = "polytomous",
      thetas = thetas,
      deltaslist = deltaslist,
      n_items = ncol(data_mat),
      sample_n = sample_n
    )
  } else {
    # Fit RM, extract thetas
    rm_fit <- eRm::RM(data_mat)
    pp <- eRm::person.parameter(rm_fit)
    theta_table <- pp$theta.table[["Person Parameter"]]
    raw_scores <- rowSums(data_mat, na.rm = TRUE)
    thetas <- as.numeric(stats::na.omit(theta_table[raw_scores]))
    item_params <- -rm_fit$betapar
    sim_data_list <- list(
      type = "dichotomous",
      thetas = thetas,
      item_params = item_params,
      n_items = ncol(data_mat),
      sample_n = sample_n
    )
  }

  if (use_parallel) {
    results_raw <- run_q3_sim_parallel(iterations, sim_seeds, sim_data_list, n_cores, verbose)
  } else {
    results_raw <- run_q3_sim_sequential(iterations, sim_seeds, sim_data_list, verbose)
  }

  # Filter out failures (character strings indicate errors)
  ok <- vapply(results_raw, function(x) is.list(x), logical(1L))
  successful <- results_raw[ok]

  if (length(successful) == 0L) {
    stop("All simulation iterations failed. Check your data.", call. = FALSE)
  }

  actual_iterations <- length(successful)
  mean_q3 <- vapply(successful, function(x) x$mean, numeric(1L))
  max_q3  <- vapply(successful, function(x) x$max,  numeric(1L))
  diff_q3 <- max_q3 - mean_q3

  results <- data.frame(
    mean = mean_q3,
    max  = max_q3,
    diff = diff_q3
  )

  out <- list()
  out$results           <- results
  out$actual_iterations <- actual_iterations
  out$sample_n          <- sample_n
  out$sample_summary    <- summary(sim_data_list$thetas)
  out$max_diff          <- max(results$diff)
  out$sd_diff           <- stats::sd(results$diff)
  out$p95               <- stats::quantile(results$diff, 0.95)
  out$p99               <- stats::quantile(results$diff, 0.99)
  out$p995              <- stats::quantile(results$diff, 0.995)
  out$p999              <- stats::quantile(results$diff, 0.999)
  out$suggested_cutoff  <- out$p99
  out
}

# ---------------------------------------------------------------------------
# Internal: single simulation iteration
# ---------------------------------------------------------------------------

#' Run a single Q3 simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMlocdepQ3cutoff()].
#' @return A list with `mean` and `max` Q3, or a character string on failure.
#' @keywords internal
run_single_q3_sim <- function(seed, data_list) {
  set.seed(seed)

  thetas_res <- sample(data_list$thetas, size = data_list$sample_n, replace = TRUE)

  tryCatch({
    if (data_list$type == "dichotomous") {
      # Simulate dichotomous data via psychotools::rrm()
      sim_mat <- psychotools::rrm(
        theta = thetas_res,
        beta = data_list$item_params
      )
      sim_df <- as.data.frame(sim_mat$data)
      
      # Validate: every item must have at least 8 positive responses.
      # Fewer than 8 positives can cause numerical instability in mirt fitting.
      pos_counts <- colSums(sim_df, na.rm = TRUE)
      if (any(pos_counts < 8L)) {
        return("validation_failed: fewer than 8 positive responses in at least one item")
      }
    } else {
      # Simulate polytomous data
      sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
      sim_df <- as.data.frame(sim_mat)

      # Validate: all categories must be represented in each item
      n_cats <- vapply(data_list$deltaslist, function(d) length(d) + 1L, integer(1L))
      for (j in seq_len(ncol(sim_df))) {
        tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
        if (any(tab == 0L)) {
          return("validation_failed: not all categories represented")
        }
      }
    }

    # Fit Rasch model via mirt
    mirt_fit <- mirt::mirt(
      sim_df,
      model = 1,
      itemtype = "Rasch",
      verbose = FALSE,
      accelerate = "squarem"
    )

    q3_mat <- mirt::residuals(mirt_fit, type = "Q3", digits = 4, verbose = FALSE)
    diag(q3_mat) <- NA

    mean_q3 <- mean(q3_mat, na.rm = TRUE)
    max_q3  <- max(q3_mat, na.rm = TRUE)

    list(mean = mean_q3, max = max_q3)
  }, error = function(e) {
    as.character(conditionMessage(e))
  })
}

# ---------------------------------------------------------------------------
# Internal: parallel runner
# ---------------------------------------------------------------------------

#' Run Q3 simulations in parallel using mirai
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param n_cores Number of mirai daemons.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_q3_sim_parallel <- function(iterations, sim_seeds, sim_data_list,
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
        run_single_q3_sim(seed, data_list)
      },
      seed = sim_seeds[sim],
      data_list = sim_data_list,
      run_single_q3_sim = run_single_q3_sim,
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

#' Run Q3 simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_q3_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                  verbose = FALSE) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_q3_sim(sim_seeds[sim], sim_data_list)
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
