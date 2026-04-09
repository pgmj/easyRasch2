#' Partial Gamma Local Dependence Analysis
#'
#' Computes partial gamma coefficients for Local Dependence (LD) assessment
#' using \code{iarm::partgam_LD()}. Each pair of items is tested for residual
#' association, controlling for the rest score (total score minus one of the
#' items in the pair).
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed,
#'   but at least one complete case must exist.
#' @param cutoff Optional. Default `NULL` (no cutoff applied). Can be:
#'   * The return value of \code{\link{RMpgLDcutoff}} (a list with
#'     `$pair_cutoffs`): the data.frame is extracted automatically and
#'     simulation metadata is included in the kable caption.
#'   * The `$pair_cutoffs` data.frame from \code{\link{RMpgLDcutoff}}
#'     directly: must have columns `Item1`, `Item2`, `gamma_low`, `gamma_high`.
#'   When provided, adds columns `Gamma_low`, `Gamma_high`, and `Flagged`
#'   (logical; `TRUE` when the observed partial gamma falls outside the
#'   credible range) to the result.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying data.frame.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object with columns "Item 1",
#'   "Item 2", "Partial gamma", and "Adjusted p-value (BH)". When `cutoff`
#'   is provided, additional columns "Gamma low", "Gamma high", and "Flagged"
#'   are included. Two tables are returned (one per rest-score direction),
#'   combined into a single output.
#' * If `output = "dataframe"`: a list of two data.frames (one per rest-score
#'   direction) with columns `Item1`, `Item2`, `gamma`, `padj_bh`. When
#'   `cutoff` is provided, columns `gamma_low`, `gamma_high`, and `flagged`
#'   are also included.
#'
#' @details
#' Partial gamma (Christensen, Kreiner & Mesbah, 2013) measures the residual
#' association between pairs of items after controlling for the rest score
#' (total score minus one item). Because it matters which item is subtracted,
#' calculations are done for each pair in both directions, yielding two
#' data.frames.
#'
#' Values near 0 indicate no local dependence. Large positive values suggest
#' positive LD (items share variance beyond the latent trait), while large
#' negative values suggest negative LD.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' @references
#' Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.)
#' \emph{Rasch Models in Health}. Iste and Wiley (2013), pp. 133--135.
#'
#' @seealso \code{\link{RMpgLDcutoff}}, \code{\link{RMpgLDplot}}
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
#' # Default kable output
#' RMpartgamLD(sim_data)
#'
#' # Return as data.frame list
#' RMpartgamLD(sim_data, output = "dataframe")
#'
#' # With simulation-based cutoffs
#' cutoff_res <- RMpgLDcutoff(sim_data, iterations = 100, parallel = FALSE,
#'                            seed = 42)
#' RMpartgamLD(sim_data, cutoff = cutoff_res)
#' }
RMpartgamLD <- function(data, cutoff = NULL, output = "kable") {

  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMpartgamLD() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe"))

  validate_response_data(data)

  # --- Validate and normalise cutoff ------------------------------------------
  cutoff_n_iter     <- NULL
  cutoff_method     <- NULL
  cutoff_hdci_width <- NULL
  if (!is.null(cutoff)) {
    if (is.list(cutoff) && !is.data.frame(cutoff) && "pair_cutoffs" %in% names(cutoff)) {
      cutoff_n_iter     <- cutoff$actual_iterations
      cutoff_method     <- cutoff$cutoff_method
      cutoff_hdci_width <- cutoff$hdci_width
      cutoff <- cutoff$pair_cutoffs
    }
    if (!is.data.frame(cutoff)) {
      stop(
        "`cutoff` must be NULL, the return value of RMpgLDcutoff(), or its ",
        "$pair_cutoffs data.frame.",
        call. = FALSE
      )
    }
    required_cols <- c("Item1", "Item2", "gamma_low", "gamma_high")
    missing_cols <- setdiff(required_cols, names(cutoff))
    if (length(missing_cols) > 0L) {
      stop(
        "`cutoff` data.frame is missing required columns: ",
        paste(missing_cols, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # --- rgl workaround ---------------------------------------------------------
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # --- Compute partial gamma LD via iarm --------------------------------------
  sink(nullfile())
  pgam_raw <- iarm::partgam_LD(as.data.frame(data))
  sink()

  # pgam_raw is a list of two data.frames (one per rest-score direction)
  process_pgam_df <- function(raw_df) {
    df <- data.frame(
      Item1   = as.character(raw_df$Item1),
      Item2   = as.character(raw_df$Item2),
      gamma   = round(as.numeric(raw_df$gamma), 3),
      padj_bh = round(as.numeric(raw_df[[6]]), 3),
      stringsAsFactors = FALSE
    )
    df
  }

  result_list <- list(
    process_pgam_df(pgam_raw[[1]]),
    process_pgam_df(pgam_raw[[2]])
  )

  # --- Apply cutoff if provided -----------------------------------------------
  # Cutoffs are keyed by direction-1 pairs (i < j), so for direction 2 (i > j)

  # we must look up by the canonical (sorted) pair order.
  if (!is.null(cutoff)) {
    # Build a canonical key on the cutoff table
    cutoff$canonical_key <- paste(
      pmin(cutoff$Item1, cutoff$Item2),
      pmax(cutoff$Item1, cutoff$Item2),
      sep = "___"
    )
    cutoff_sub <- cutoff[, c("canonical_key", "gamma_low", "gamma_high")]

    for (idx in seq_along(result_list)) {
      result_df <- result_list[[idx]]
      # Build the same canonical key on the result (sorted pair order)
      result_df$canonical_key <- paste(
        pmin(result_df$Item1, result_df$Item2),
        pmax(result_df$Item1, result_df$Item2),
        sep = "___"
      )

      merged <- merge(result_df, cutoff_sub, by = "canonical_key", all.x = TRUE,
                      sort = FALSE)
      # Restore original row order
      merged <- merged[match(result_df$canonical_key, merged$canonical_key), ]
      rownames(merged) <- NULL
      merged$gamma_low  <- round(merged$gamma_low, 3)
      merged$gamma_high <- round(merged$gamma_high, 3)
      merged$flagged <- !is.na(merged$gamma_low) &
        (merged$gamma < merged$gamma_low | merged$gamma > merged$gamma_high)

      # Remove helper column
      merged$canonical_key <- NULL

      merged <- merged[, c("Item1", "Item2", "gamma", "padj_bh",
                           "gamma_low", "gamma_high", "flagged")]
      result_list[[idx]] <- merged
    }
    # Clean up cutoff helper column
    cutoff$canonical_key <- NULL
  }

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(result_list)
  }

  # Build caption
  n_complete <- sum(stats::complete.cases(as.data.frame(data)))
  if (is.null(cutoff)) {
    caption_text <- paste0(
      "Partial gamma LD analysis (n = ", n_complete, " complete cases). ",
      "Positive gamma indicates positive local dependence between items."
    )
  } else if (!is.null(cutoff_n_iter)) {
    method_label <- .format_gamma_cutoff_method_label(cutoff_method, cutoff_hdci_width)
    iter_part <- paste0(cutoff_n_iter, " simulation iterations")
    caption_text <- paste0(
      "Partial gamma LD analysis (n = ", n_complete, " complete cases). Cutoff values based on ",
      if (!is.null(method_label)) {
        paste0(iter_part, " (", method_label, ").")
      } else {
        paste0(iter_part, ".")
      }
    )
  } else {
    caption_text <- paste0(
      "Partial gamma LD analysis (n = ", n_complete,
      " complete cases). Simulation-based cutoff values applied."
    )
  }

  col_names_no_cutoff <- c("Item 1", "Item 2", "Partial gamma",
                           "Adjusted p-value (BH)")
  col_names_cutoff    <- c("Item 1", "Item 2", "Partial gamma",
                           "Adjusted p-value (BH)",
                           "Gamma low", "Gamma high", "Flagged")

  kable1 <- knitr::kable(
    result_list[[1]],
    format    = "pipe",
    col.names = if (is.null(cutoff)) col_names_no_cutoff else col_names_cutoff,
    caption   = paste0(caption_text,
                       "\n\nDirection 1: rest score = total - Item2")
  )

  kable2 <- knitr::kable(
    result_list[[2]],
    format    = "pipe",
    col.names = if (is.null(cutoff)) col_names_no_cutoff else col_names_cutoff,
    caption   = paste0(caption_text,
                       "\n\nDirection 2: rest score = total - Item1")
  )

  knitr::asis_output(paste(kable1, "\n\n", kable2))
}

#' Simulation-Based Partial Gamma LD Cutoff Determination
#'
#' Uses parametric bootstrap simulation to determine appropriate cutoff values
#' for partial gamma Local Dependence analysis via
#' \code{\link[iarm]{partgam_LD}}. Under a correctly fitting Rasch model where
#' items are locally independent, this function generates the expected
#' distribution of partial gamma values per item pair, providing empirical
#' critical values.
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
#'   Default is `0.99` (99\% HDCI). Ignored when
#'   `cutoff_method = "quantile"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`results`}{data.frame with columns `iteration`, `Item1`, `Item2`,
#'     and `gamma` (one row per item pair per successful iteration). Contains
#'     results from direction 1 only (rest score = total - Item2), which is
#'     the conventional direction.}
#'   \item{`pair_cutoffs`}{data.frame with per-pair cutoff summaries: `Item1`,
#'     `Item2`, `gamma_low`, `gamma_high`. Bounds are computed using the method
#'     specified by `cutoff_method`.}
#'   \item{`actual_iterations`}{Number of successful iterations.}
#'   \item{`sample_n`}{Number of complete cases used.}
#'   \item{`sample_summary`}{Summary statistics of estimated person
#'     parameters.}
#'   \item{`item_names`}{Character vector of item names from data.}
#'   \item{`cutoff_method`}{The method used to compute cutoffs (`"hdci"` or
#'     `"quantile"`).}
#'   \item{`hdci_width`}{The HDCI width used (only meaningful when
#'     `cutoff_method = "hdci"`).}
#' }
#'
#' @details
#' For each simulation iteration the function:
#' \enumerate{
#'   \item Resamples person parameters (thetas) with replacement from ML
#'     estimates.
#'   \item Simulates item response data under a Rasch model (dichotomous via
#'     `psychotools::rrm()` or polytomous via an internal partial credit
#'     simulator).
#'   \item Computes partial gamma LD statistics via
#'     `iarm::partgam_LD()`.
#' }
#'
#' Because the data are simulated under the Rasch model, items are locally
#' independent by construction. The distribution of partial gamma values
#' across iterations provides empirical critical values per item pair. Values
#' from real data that fall outside these bounds suggest local dependence that
#' exceeds what would be expected by chance. Failed iterations (e.g., due to
#' convergence issues or degenerate data) are silently discarded.
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
#' @references
#' Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.)
#' \emph{Rasch Models in Health}. Iste and Wiley (2013), pp. 133--135.
#'
#' @seealso \code{\link[iarm]{partgam_LD}}, \code{\link{RMpartgamLD}},
#'   \code{\link{RMpgLDplot}}
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
#' cutoff_res <- RMpgLDcutoff(sim_data, iterations = 100, parallel = FALSE,
#'                            seed = 42)
#' cutoff_res$pair_cutoffs
#' }
RMpgLDcutoff <- function(data, iterations = 250,
                         parallel = TRUE, n_cores = NULL,
                         verbose = FALSE, seed = NULL,
                         cutoff_method = "hdci",
                         hdci_width = 0.99) {

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
      "Package 'iarm' is required for RMpgLDcutoff() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  validate_response_data(data)

  # rgl workaround (iarm depends on vcdExtra -> rgl)
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # Only complete cases
  complete_idx <- stats::complete.cases(data)
  data <- data[complete_idx, , drop = FALSE]

  if (nrow(data) == 0L) {
    stop("No complete cases in data after removing rows with NA.",
         call. = FALSE)
  }

  # --- Parallel setup ---------------------------------------------------------
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
    theta_table <- pp$theta.table[["Person Parameter"]]
    raw_scores <- rowSums(data_mat, na.rm = TRUE)
    thetas <- as.numeric(stats::na.omit(theta_table[raw_scores]))
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
    theta_table <- pp$theta.table[["Person Parameter"]]
    raw_scores <- rowSums(data_mat, na.rm = TRUE)
    thetas <- as.numeric(stats::na.omit(theta_table[raw_scores]))
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
    results_raw <- run_partgam_LD_sim_parallel(iterations, sim_seeds,
                                               sim_data_list, n_cores, verbose)
  } else {
    results_raw <- run_partgam_LD_sim_sequential(iterations, sim_seeds,
                                                 sim_data_list, verbose)
  }

  # Filter out failures (character strings indicate errors)
  ok <- vapply(results_raw, is.data.frame, logical(1L))
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

  # Compute per-pair cutoffs
  pair_keys <- unique(paste(results_df$Item1, results_df$Item2, sep = "___"))
  pair_cutoffs <- do.call(rbind, lapply(pair_keys, function(pk) {
    parts <- strsplit(pk, "___", fixed = TRUE)[[1]]
    sub <- results_df[results_df$Item1 == parts[1] & results_df$Item2 == parts[2], ]
    if (cutoff_method == "hdci") {
      gamma_interval <- ggdist::hdci(sub$gamma, .width = hdci_width)
      data.frame(
        Item1      = parts[1],
        Item2      = parts[2],
        gamma_low  = gamma_interval[1L, 1L],
        gamma_high = gamma_interval[1L, 2L],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    } else {
      data.frame(
        Item1      = parts[1],
        Item2      = parts[2],
        gamma_low  = stats::quantile(sub$gamma, 0.025, na.rm = TRUE),
        gamma_high = stats::quantile(sub$gamma, 0.975, na.rm = TRUE),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }))
  rownames(pair_cutoffs) <- NULL

  list(
    results           = results_df,
    pair_cutoffs      = pair_cutoffs,
    actual_iterations = actual_iterations,
    sample_n          = sample_n,
    sample_summary    = summary(thetas),
    item_names        = item_names_vec,
    cutoff_method     = cutoff_method,
    hdci_width        = hdci_width
  )
}

# ---------------------------------------------------------------------------
# Internal: single simulation iteration
# ---------------------------------------------------------------------------

#' Run a single partial gamma LD simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMpgLDcutoff()].
#' @return A data.frame with columns `Item1`, `Item2`, and `gamma`, or a
#'   character string on failure.
#' @keywords internal
run_single_partgam_LD_sim <- function(seed, data_list) {
  set.seed(seed)

  thetas_res <- sample(data_list$thetas, size = data_list$sample_n,
                       replace = TRUE)

  tryCatch({
    if (data_list$type == "dichotomous") {
      sim_mat <- psychotools::rrm(
        theta = thetas_res,
        beta = data_list$item_params
      )
      sim_df <- as.data.frame(sim_mat$data)
      colnames(sim_df) <- data_list$item_names

      pos_counts <- colSums(sim_df, na.rm = TRUE)
      if (any(pos_counts < 8L)) {
        return("validation_failed: fewer than 8 positive responses in at least one item")
      }
    } else {
      sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
      sim_df <- as.data.frame(sim_mat)
      colnames(sim_df) <- data_list$item_names

      n_cats <- vapply(data_list$deltaslist, function(d) length(d) + 1L,
                       integer(1L))
      for (j in seq_len(ncol(sim_df))) {
        tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
        if (any(tab == 0L)) {
          return("validation_failed: not all categories represented")
        }
      }
    }

    # Compute partial gamma LD via iarm
    # iarm::partgam_LD returns a list of two data.frames;
    # suppress its console output with sink()
    sink(nullfile())
    pgam <- iarm::partgam_LD(sim_df)
    sink()

    # Use direction 1 (rest score = total - Item2) as the canonical direction
    pgam_df <- as.data.frame(pgam[[1]])

    data.frame(
      Item1 = as.character(pgam_df$Item1),
      Item2 = as.character(pgam_df$Item2),
      gamma = as.numeric(pgam_df$gamma),
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

#' Run partial gamma LD simulations in parallel using mirai
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param n_cores Number of mirai daemons.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_partgam_LD_sim_parallel <- function(iterations, sim_seeds, sim_data_list,
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
        run_single_partgam_LD_sim(seed, data_list)
      },
      seed = sim_seeds[sim],
      data_list = sim_data_list,
      run_single_partgam_LD_sim = run_single_partgam_LD_sim,
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

#' Run partial gamma LD simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_partgam_LD_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                          verbose = FALSE) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_partgam_LD_sim(sim_seeds[sim], sim_data_list)
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

#' Plot Distribution of Simulated Partial Gamma LD Values
#'
#' Visualises the distribution of simulation-based partial gamma LD values
#' from \code{\link{RMpgLDcutoff}}, optionally overlaying observed partial
#' gamma values computed from real data via \code{\link[iarm]{partgam_LD}}.
#'
#' Uses `ggdist::stat_dotsinterval()` (when `data` is not supplied) or
#' `ggdist::stat_dots()` (when `data` is supplied) with
#' `point_interval = "median_hdci"` and `.width = c(0.66, 0.95, 0.99)`.
#'
#' @param simfit The return value of \code{\link{RMpgLDcutoff}} (a list with
#'   components `results`, `pair_cutoffs`, `actual_iterations`, `sample_n`, and
#'   `item_names`).
#' @param data Optional. A data.frame or matrix of item responses for computing
#'   and overlaying observed partial gamma values. Items must be scored starting
#'   at 0 (non-negative integers). When provided, the plot includes orange
#'   diamond markers for the observed partial gamma alongside the simulated
#'   distribution, plus segment summaries from the cutoff intervals.
#' @param items Optional character vector of item names to include in the plot.
#'   Only item pairs where **both** items are in this vector will be shown. When
#'   `NULL` (default), all item pairs are plotted.
#'
#' @return A `ggplot` object.
#'
#' @details
#' The plot shows one row per item pair (labelled as "Item1 - Item2"). Only
#' direction 1 (rest score = total - Item2) is plotted, matching the
#' convention used in the simulation.
#'
#' When `data` is **not** supplied, the function plots the simulated partial
#' gamma distributions as dot-interval plots using
#' `ggdist::stat_dotsinterval()` with median and Highest Density Continuous
#' Interval (HDCI) summaries.
#'
#' When `data` **is** supplied, the function:
#' \enumerate{
#'   \item Computes observed partial gamma values via
#'     `iarm::partgam_LD()`.
#'   \item Overlays observed gamma values as orange diamond markers on the
#'     simulated distributions.
#'   \item Shows per-pair cutoff intervals (from `simfit$pair_cutoffs`) as
#'     black line segments, with thicker segments for the 66\% interval and
#'     black dots for the median.
#' }
#'
#' The `ggplot2`, `ggdist`, and optionally `iarm` packages must be installed
#' (they are in Suggests, not Imports).
#'
#' @seealso \code{\link{RMpgLDcutoff}}, \code{\link{RMpartgamLD}}
#'
#' @importFrom rlang .data
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
#' # Run simulation
#' cutoff_res <- RMpgLDcutoff(sim_data, iterations = 100, parallel = FALSE,
#'                            seed = 42)
#'
#' # Simulated distribution only
#' RMpgLDplot(cutoff_res)
#'
#' # With observed partial gamma overlaid
#' RMpgLDplot(cutoff_res, data = sim_data)
#'
#' # Plot only a subset of items
#' RMpgLDplot(cutoff_res, data = sim_data,
#'            items = c("Item1", "Item2", "Item3"))
#' }
RMpgLDplot <- function(simfit, data, items = NULL) {

  # --- Check required packages ------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for RMpgLDplot() but is not installed.\n",
      "Install it with: install.packages(\"ggplot2\")",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required for RMpgLDplot() but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")",
      call. = FALSE
    )
  }

  # --- Validate simfit --------------------------------------------------------
  required_names <- c("results", "pair_cutoffs", "actual_iterations",
                      "sample_n", "item_names")
  missing_names <- setdiff(required_names, names(simfit))
  if (length(missing_names) > 0L) {
    stop(
      "`simfit` is missing required components: ",
      paste(missing_names, collapse = ", "),
      ".\nExpected the return value of RMpgLDcutoff().",
      call. = FALSE
    )
  }

  results_df        <- simfit$results
  pair_cutoffs      <- simfit$pair_cutoffs
  actual_iterations <- simfit$actual_iterations
  sample_n          <- simfit$sample_n
  item_names        <- simfit$item_names

  # --- Validate items parameter -----------------------------------------------
  if (!is.null(items)) {
    unknown_items <- setdiff(items, item_names)
    if (length(unknown_items) > 0L) {
      stop(
        "Unknown item(s) in `items`: ",
        paste(unknown_items, collapse = ", "),
        ".\nAvailable items: ", paste(item_names, collapse = ", "),
        call. = FALSE
      )
    }
    if (length(items) < 2L) {
      stop("`items` must contain at least 2 item names to form a pair.",
           call. = FALSE)
    }
  }

  # Create pair labels
  results_df$Pair <- paste(results_df$Item1, "-", results_df$Item2)
  pair_cutoffs$Pair <- paste(pair_cutoffs$Item1, "-", pair_cutoffs$Item2)

  # --- Filter to selected items -----------------------------------------------
  if (!is.null(items)) {
    keep <- results_df$Item1 %in% items & results_df$Item2 %in% items
    results_df <- results_df[keep, , drop = FALSE]
    keep_cut <- pair_cutoffs$Item1 %in% items & pair_cutoffs$Item2 %in% items
    pair_cutoffs <- pair_cutoffs[keep_cut, , drop = FALSE]

    if (nrow(results_df) == 0L) {
      stop("No item pairs remain after filtering by `items`.", call. = FALSE)
    }
  }

  # Pair factor levels (reversed for plotting top-to-bottom)
  pair_levels <- rev(unique(results_df$Pair))

  # --- Compute per-pair summary intervals for segment overlays ----------------
  pair_names <- unique(results_df$Pair)
  lo_hi <- do.call(rbind, lapply(pair_names, function(pair) {
    sub <- results_df[results_df$Pair == pair, ]
    data.frame(
      Pair            = pair,
      min_gamma       = stats::quantile(sub$gamma, 0.005, na.rm = TRUE),
      max_gamma       = stats::quantile(sub$gamma, 0.995, na.rm = TRUE),
      p66lo_gamma     = stats::quantile(sub$gamma, 0.167, na.rm = TRUE),
      p66hi_gamma     = stats::quantile(sub$gamma, 0.833, na.rm = TRUE),
      median_gamma    = stats::median(sub$gamma, na.rm = TRUE),
      stringsAsFactors = FALSE,
      row.names       = NULL
    )
  }))
  rownames(lo_hi) <- NULL

  # --- Case 1: no observed data, show simulation distribution only ------------
  if (missing(data)) {

    results_plot <- data.frame(
      Pair  = results_df$Pair,
      Value = results_df$gamma,
      stringsAsFactors = FALSE
    )
    results_plot$Pair <- factor(results_plot$Pair, levels = pair_levels)

    p <- ggplot2::ggplot(
      results_plot,
      ggplot2::aes(
        x = .data$Value,
        y = .data$Pair
      )
    ) +
      ggdist::stat_dotsinterval(
        ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
        quantiles = actual_iterations,
        point_interval = "median_hdci",
        layout = "weave",
        slab_color = NA,
        .width = c(0.66, 0.95, 0.99)
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "grey50",
        linewidth = 0.4
      ) +
      ggplot2::labs(
        x = "Partial gamma",
        y = "Item pair",
        caption = paste0(
          "Note: Results from ", actual_iterations,
          " simulated datasets with ", sample_n,
          " respondents (no true local dependence)."
        )
      ) +
      ggplot2::scale_color_manual(
        values = scales::brewer_pal()(3),
        aesthetics = "slab_fill",
        guide = "none"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"))

    return(p)
  }

  # --- Case 2: observed data supplied -----------------------------------------
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required to compute observed partial gamma but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  validate_response_data(data)

  # rgl workaround
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  sink(nullfile())
  pgam_raw <- iarm::partgam_LD(as.data.frame(data))
  sink()

  # Use direction 1 (rest score = total - Item2)
  observed_df <- data.frame(
    Item1          = as.character(pgam_raw[[1]]$Item1),
    Item2          = as.character(pgam_raw[[1]]$Item2),
    observed_gamma = as.numeric(pgam_raw[[1]]$gamma),
    stringsAsFactors = FALSE
  )
  observed_df$Pair <- paste(observed_df$Item1, "-", observed_df$Item2)

  # Filter observed data to selected items
  if (!is.null(items)) {
    keep_obs <- observed_df$Item1 %in% items & observed_df$Item2 %in% items
    observed_df <- observed_df[keep_obs, , drop = FALSE]
  }

  # --- Build plot data --------------------------------------------------------
  gamma_sim <- data.frame(
    Pair  = results_df$Pair,
    Value = results_df$gamma,
    stringsAsFactors = FALSE
  )
  gamma_sim <- merge(gamma_sim, observed_df[, c("Pair", "observed_gamma")],
                     by = "Pair", sort = FALSE)
  gamma_sim$Pair <- factor(gamma_sim$Pair, levels = pair_levels)

  lo_hi$Pair_f <- factor(lo_hi$Pair, levels = pair_levels)

  caption_text <- paste0(
    "Note: Results from ", actual_iterations,
    " simulated datasets with ", sample_n, " respondents.\n",
    "Orange diamonds indicate observed partial gamma LD. ",
    "Black dots indicate median gamma from simulations."
  )

  p <- ggplot2::ggplot(
    gamma_sim,
    ggplot2::aes(
      x = .data$Value,
      y = .data$Pair
    )
  ) +
    ggdist::stat_dots(
      ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
      quantiles = actual_iterations,
      layout = "weave",
      slab_color = NA,
      .width = c(0.66, 0.95, 0.99)
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x    = .data$min_gamma,
        xend = .data$max_gamma,
        y    = .data$Pair_f,
        yend = .data$Pair_f
      ),
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x    = .data$p66lo_gamma,
        xend = .data$p66hi_gamma,
        y    = .data$Pair_f,
        yend = .data$Pair_f
      ),
      color = "black",
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = lo_hi,
      ggplot2::aes(
        x = .data$median_gamma,
        y = .data$Pair_f
      ),
      size = 3.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$observed_gamma),
      color = "sienna2",
      shape = 18,
      position = ggplot2::position_nudge(y = -0.1),
      size = 4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.4
    ) +
    ggplot2::labs(
      x = "Partial gamma",
      y = "Item pair",
      caption = caption_text
    ) +
    ggplot2::scale_color_manual(
      values = scales::brewer_pal()(3),
      aesthetics = "slab_fill",
      guide = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"))

  p
}
