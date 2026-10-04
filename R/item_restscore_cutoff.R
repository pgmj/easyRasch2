#' Simulation-Based Item-Restscore Null Distribution
#'
#' Uses parametric bootstrap simulation to build the null distribution of the
#' item-restscore statistic for \code{\link{RMitemRestscore}}. This function
#' simulates data from a correctly fitting Rasch model that mimics your data
#' and returns, per item, the simulated difference between observed and
#' expected item-restscore gamma.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Only complete cases (rows without
#'   any `NA`) are used.
#' @param iterations Integer. Number of simulation iterations (default 400).
#' @param parallel Logical. Use parallel processing via `mirai` if available
#'   (default `TRUE`).
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first. If neither is set and
#'   `parallel = TRUE`, a warning is issued and execution falls back to
#'   sequential (single core) processing.
#' @param verbose Logical. Show a progress bar (default `FALSE`).
#' @param seed Integer or `NULL`. Random seed for reproducibility. See
#'   [easyRasch2-reproducibility] for what this guarantees and how it
#'   interacts with `parallel`.
#' @param cutoff_method Character string specifying how the intervals are
#'   computed. Either `"hdci"` (default) for the Highest Density Interval via
#'   `ggdist::hdci()`, or `"quantile"` for the 2.5th/97.5th percentiles via
#'   `stats::quantile()`.
#' @param hdci_width Numeric. Width of the HDCI when `cutoff_method = "hdci"`.
#'   Default is `0.95` (95% HDCI). Ignored when `cutoff_method = "quantile"`.
#'
#'   The interval is a **description** of where a fitting item's difference is
#'   expected to fall, not a decision rule. Flagging every item outside a
#'   width-`w` interval tests all `k` items at once, so the family-wise error
#'   rate is `1 - w^k`. Decisions should come from the corrected p-value
#'   instead (\code{\link{RMitemRestscore}} with `p_value = NULL` and the full
#'   object returned here).
#' @param dgp Character. Data-generating process for the parametric bootstrap.
#'   `"conditional"` (default) simulates each respondent's pattern from the
#'   exact Rasch conditional distribution given their observed total score,
#'   item parameters fixed (a *conditional* null). The expected gamma is
#'   computed from the observed score distribution, which the conditional
#'   null holds fixed. `"resample"` resamples WLE person locations with
#'   replacement and simulates responses under the model (a *marginal* null).
#'   In simulation under a true Rasch model, the conditional null gave a
#'   family-wise rate of 5.0 percent with `correction = "fwer"` and the
#'   resample null 6.4 percent, pooled over four designs. The conditional
#'   null takes about 1.2 to 1.6 times as long.
#'   \strong{Experimental.}
#'
#' @return A list with components:
#' \describe{
#'   \item{`results`}{data.frame with columns `iteration`, `Item`, `Observed`,
#'     `Expected`, `Difference` (one row per item per successful iteration).
#'     `Difference` is `Observed - Expected`, the statistic
#'     \code{\link{RMitemRestscore}} tests.}
#'   \item{`item_cutoffs`}{data.frame with per-item interval bounds for the
#'     difference: `Item`, `diff_low`, `diff_high`. Bounds are computed using
#'     the method specified by `cutoff_method`.}
#'   \item{`actual_iterations`}{Number of successful iterations. Everything
#'     downstream rests on this rather than on `iterations`, so it is the
#'     number to report.}
#'   \item{`requested_iterations`}{The `iterations` argument, kept so callers
#'     can tell how many simulated datasets were discarded.}
#'   \item{`sample_n`}{Number of complete cases used.}
#'   \item{`sample_n_total`}{Number of respondents in the raw input data,
#'     before the complete-case filter.}
#'   \item{`sample_has_na`}{Logical. Whether the raw input data contained
#'     any missing values.}
#'   \item{`sample_summary`}{Summary statistics of estimated person parameters.}
#'   \item{`item_names`}{Character vector of item names from data.}
#'   \item{`cutoff_method`}{The method used to compute the intervals (`"hdci"`
#'     or `"quantile"`).}
#'   \item{`hdci_width`}{The HDCI width used (only meaningful when
#'     `cutoff_method = "hdci"`).}
#'   \item{`dgp`}{The data-generating process used (`"resample"` or
#'     `"conditional"`).}
#' }
#'
#' @details
#' The asymptotic test in `iarm::item_restscore()` divides the observed minus
#' expected gamma by the standard error of the observed gamma alone. The
#' expected gamma is estimated from the same data and correlates with the
#' observed one, so that standard error is too large for the difference, and
#' the difference is also biased upwards in small samples. Under a true Rasch
#' model the resulting test is liberal for dichotomous items in small or
#' mistargeted samples and conservative for polytomous items, in every case
#' flagging too few underfitting items. The bootstrap null replaces the
#' asymptotic reference distribution and absorbs both problems.
#'
#' The generating model is CML item parameters (via `psychotools`) with WLE
#' person locations. For each iteration a dataset is simulated under the chosen
#' `dgp`, the model is **refitted** by CML (`psychotools::pcmodel()`), and the
#' observed and expected item-restscore gamma are computed via
#' `iarm::item_restscore()`. The refit matters: the expected gamma varies from
#' sample to sample because the thresholds do, and holding them fixed would
#' reproduce the problem the bootstrap exists to solve. Failed iterations
#' (e.g., degenerate simulated data) are silently discarded.
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
#' Johansson, M. (2026). Simulation-based cutoffs for conditional item fit in
#' Rasch models: Iterations, multiplicity correction, and decision stability.
#' *PsyArXiv*. \doi{10.31234/osf.io/7pqz4_v2}
#'
#' @seealso \code{\link{RMitemRestscore}}, \code{\link{RMitemRestscorePlot}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("iarm", quietly = TRUE) &&
#'     requireNamespace("ggdist", quietly = TRUE)) {
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:10)
#'
#'   # Run 100 iterations sequentially for a quick demo
#'   cutoff_res <- RMitemRestscoreCutoff(sim_data, iterations = 100,
#'                                       parallel = FALSE, seed = 42)
#'   cutoff_res$item_cutoffs
#'
#'   # Flag on bootstrap p-values in RMitemRestscore()
#'   RMitemRestscore(sim_data, cutoff = cutoff_res)
#' }
#' }
RMitemRestscoreCutoff <- function(
  data,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95,
  dgp = c("conditional", "resample")
) {
  dgp <- match.arg(dgp)
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
      "Package 'iarm' is required for RMitemRestscoreCutoff() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  validate_response_data(data)

  # rgl workaround
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # Only complete cases. iarm::item_restscore() refits on complete cases when
  # the fitted object contains NA, so the observed statistic in
  # RMitemRestscore() is complete-case as well.
  n_total <- nrow(as.data.frame(data))
  has_na <- anyNA(data)
  data <- stats::na.omit(data)
  if (nrow(data) == 0L) {
    stop(
      "No complete cases in data. All rows contain at least one NA.",
      call. = FALSE
    )
  }

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

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate per-iteration seeds
  sim_seeds <- sample.int(.Machine$integer.max, iterations)

  data_mat <- as.matrix(data)
  sample_n <- nrow(data_mat)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  item_names_vec <- colnames(data_mat)

  # Generating model: CML item thresholds (psychotools), computed once, with
  # the same two DGPs as RMitemInfitCutoff().
  pool <- .wle_theta_pool(data_mat)
  thr_list <- pool$thr_list
  wle_thetas <- pool$thetas

  sim_data_list <- list(
    dgp = dgp,
    type = if (is_polytomous) "polytomous" else "dichotomous",
    thr_list = thr_list,
    n_items = ncol(data_mat),
    sample_n = sample_n,
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
    results_raw <- run_restscore_sim_parallel(
      iterations,
      sim_seeds,
      sim_data_list,
      n_cores,
      verbose
    )
  } else {
    results_raw <- run_restscore_sim_sequential(
      iterations,
      sim_seeds,
      sim_data_list,
      verbose
    )
  }

  # Filter out failures (character strings indicate errors)
  ok <- vapply(results_raw, is.list, logical(1L))
  successful <- results_raw[ok]

  if (length(successful) == 0L) {
    stop(.all_sims_failed_message(data_mat), call. = FALSE)
  }

  actual_iterations <- length(successful)

  # Combine per-iteration data.frames
  iter_dfs <- lapply(seq_along(successful), function(i) {
    df <- successful[[i]]
    df$iteration <- i
    df
  })
  results_df <- do.call(rbind, iter_dfs)
  results_df <- results_df[, c(
    "iteration",
    "Item",
    "Observed",
    "Expected",
    "Difference"
  )]
  rownames(results_df) <- NULL

  # Per-item interval for the difference
  item_cutoffs <- do.call(
    rbind,
    lapply(item_names_vec, function(item) {
      d <- results_df$Difference[results_df$Item == item]
      d <- d[is.finite(d)]
      if (cutoff_method == "hdci") {
        # ggdist::hdci() returns a matrix with ncol = 2: column 1 is the lower
        # bound, column 2 is the upper bound. Row 1 contains the continuous
        # interval.
        interval <- ggdist::hdci(d, .width = hdci_width)
        lo <- interval[1L, 1L]
        hi <- interval[1L, 2L]
      } else {
        lo <- stats::quantile(d, 0.025, names = FALSE)
        hi <- stats::quantile(d, 0.975, names = FALSE)
      }
      data.frame(
        Item = item,
        diff_low = lo,
        diff_high = hi,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    })
  )
  rownames(item_cutoffs) <- NULL

  list(
    results = results_df,
    item_cutoffs = item_cutoffs,
    actual_iterations = actual_iterations,
    requested_iterations = iterations,
    sample_n = sample_n,
    sample_n_total = n_total,
    sample_has_na = has_na,
    sample_summary = summary(wle_thetas),
    item_names = item_names_vec,
    cutoff_method = cutoff_method,
    hdci_width = hdci_width,
    dgp = dgp
  )
}

# ---------------------------------------------------------------------------
# Internal: single simulation iteration
# ---------------------------------------------------------------------------

#' Run a single item-restscore simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMitemRestscoreCutoff()].
#' @return A data.frame with columns `Item`, `Observed`, `Expected`,
#'   `Difference`, or a character string on failure.
#' @keywords internal
run_single_restscore_sim <- function(seed, data_list) {
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

  tryCatch(
    {
      # --- Generate one simulated dataset under the chosen DGP -----------------
      if (identical(data_list$dgp, "conditional")) {
        sim_df <- .sim_cond_dataset(data_list)
      } else if (data_list$type == "dichotomous") {
        thetas_res <- sample(
          data_list$thetas,
          size = data_list$sample_n,
          replace = TRUE
        )
        sim_df <- as.data.frame(
          psychotools::rrm(
            theta = thetas_res,
            beta = data_list$item_params
          )$data
        )
      } else {
        thetas_res <- sample(
          data_list$thetas,
          size = data_list$sample_n,
          replace = TRUE
        )
        sim_df <- as.data.frame(sim_partial_score(
          data_list$deltaslist,
          thetas_res
        ))
      }
      colnames(sim_df) <- data_list$item_names

      # --- Validate the simulated dataset (estimable refit) --------------------
      if (data_list$type == "dichotomous") {
        if (any(colSums(sim_df, na.rm = TRUE) < 8L)) {
          return(
            "validation_failed: fewer than 8 positive responses in at least one item"
          )
        }
      } else {
        n_cats <- vapply(
          data_list$thr_list,
          function(d) length(d) + 1L,
          integer(1L)
        )
        for (j in seq_len(ncol(sim_df))) {
          tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
          if (any(tab == 0L)) {
            return("validation_failed: not all categories represented")
          }
        }
      }

      # CML refit, then observed and expected gamma. The refit is what makes
      # the expected gamma vary across iterations as it does across samples.
      # iarm prints an empty line per call, which would otherwise repeat once
      # per iteration in the console.
      model_fit <- psychotools::pcmodel(sim_df, hessian = FALSE)
      utils::capture.output(
        res_mat <- iarm::item_restscore(model_fit, p.adj = "none")
      )
      k <- length(data_list$item_names)
      observed <- as.numeric(res_mat[seq_len(k), "observed"])
      expected <- as.numeric(res_mat[seq_len(k), "expected"])

      data.frame(
        Item = data_list$item_names,
        Observed = observed,
        Expected = expected,
        Difference = observed - expected,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    },
    error = function(e) {
      as.character(conditionMessage(e))
    }
  )
}

# ---------------------------------------------------------------------------
# Internal: parallel runner
# ---------------------------------------------------------------------------

#' Run item-restscore simulations in parallel using mirai
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param n_cores Number of mirai daemons.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_restscore_sim_parallel <- function(
  iterations,
  sim_seeds,
  sim_data_list,
  n_cores,
  verbose = FALSE
) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
    completed <- 0L
  }

  tasks <- lapply(seq_len(iterations), function(sim) {
    mirai::mirai(
      {
        options(rgl.useNULL = TRUE)
        run_single_restscore_sim(seed, data_list)
      },
      seed = sim_seeds[sim],
      data_list = sim_data_list,
      run_single_restscore_sim = run_single_restscore_sim,
      sim_partial_score = sim_partial_score,
      sim_poly_item = sim_poly_item,
      # Conditional-DGP generators (shared with the Q3 and infit cutoffs).
      .sim_cond_dataset = .sim_cond_dataset,
      .sim_conditional = .sim_conditional,
      .esf_convolve = .esf_convolve
    )
  })

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

#' Run item-restscore simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_restscore_sim_sequential <- function(
  iterations,
  sim_seeds,
  sim_data_list,
  verbose = FALSE
) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_restscore_sim(sim_seeds[sim], sim_data_list)
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
