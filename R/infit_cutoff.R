#' Simulation-based infit MSQ cutoffs
#'
#' Simulates datasets from a Rasch or partial credit model that fits and
#' resembles your data, and computes the conditional infit and outfit MSQ of
#' every item in each. Passed to [RMitemInfit()], the result gives
#' bootstrap p-values corrected across items, and an interval per item that
#' describes where a fitting item's statistic is expected to fall.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Only complete cases (rows without `NA`) are used.
#' @param iterations Integer. Number of simulated datasets. Default `400`,
#'   the recommended minimum for reproducible decisions. Use 1000 to 2000
#'   for a final analysis (Johansson, 2026).
#' @param parallel Logical. Run the simulations in parallel via \pkg{mirai}
#'   when it is installed. Default `TRUE`.
#' @param n_cores Integer or `NULL`. Number of parallel workers. When
#'   `NULL`, `getOption("mc.cores")` is used. If that is unset too, a
#'   warning is given and the simulations run sequentially.
#' @param verbose Logical. Show a progress bar. Default `FALSE`.
#' @param seed Integer or `NULL`. Random seed for reproducibility. See
#'   [easyRasch2-reproducibility] for what this guarantees and how it
#'   interacts with `parallel`.
#' @param cutoff_method Character. How the interval per item is computed:
#'   `"hdci"` (default), the highest-density continuous interval from
#'   `ggdist::hdci()` (needs \pkg{ggdist}), or `"quantile"`, the 2.5th to
#'   97.5th percentiles.
#' @param hdci_width Numeric. Width of the interval with
#'   `cutoff_method = "hdci"`. Default `0.95`. Ignored with `"quantile"`.
#'
#'   The interval describes where a fitting item's statistic is expected
#'   to fall. It is not a decision rule. Flagging every item outside a
#'   width-`w` interval tests all `k` items at once, with a family-wise
#'   error rate of `1 - w^k`, 37% for a 95% interval over nine items.
#'   Decisions come from the corrected p-value in [RMitemInfit()]. Up to
#'   version 1.1.1 the default was `0.999`, which needs roughly 5000
#'   iterations before the interval reaches its stated width. A 95%
#'   interval reaches it by about 1000 (Johansson, 2026).
#' @param dgp Character. How datasets are simulated. `"resample"`
#'   (default) draws person locations from the WLE estimates with
#'   replacement and simulates responses from the model (a marginal null).
#'   `"conditional"` simulates each respondent's response pattern from the
#'   model given their observed total score, with the item parameters fixed
#'   (a conditional null). The conditional infit is itself conditional on
#'   the total score, so `"conditional"` is its natural null.
#'   \strong{Experimental.}
#'
#' @return A list with:
#' \describe{
#'   \item{`results`}{data.frame with the columns `iteration`, `Item`,
#'     `InfitMSQ` and `OutfitMSQ`, one row per item and simulated dataset.}
#'   \item{`item_cutoffs`}{data.frame with the interval per item: `Item`,
#'     `infit_low`, `infit_high`, `outfit_low`, `outfit_high`, computed by
#'     `cutoff_method`.}
#'   \item{`actual_iterations`}{Number of simulated datasets that could be
#'     used. Everything downstream rests on this rather than on
#'     `iterations`, so it is the number to report.}
#'   \item{`requested_iterations`}{The `iterations` argument, so that
#'     callers can tell how many simulated datasets were left out.}
#'   \item{`sample_n`}{Number of complete cases used.}
#'   \item{`sample_n_total`}{Number of respondents in `data`, before the
#'     complete-case filter.}
#'   \item{`sample_has_na`}{Logical. Whether `data` contained missing
#'     values.}
#'   \item{`sample_summary`}{Summary of the WLE person locations.}
#'   \item{`item_names`}{Item names from `data`.}
#'   \item{`cutoff_method`}{`"hdci"` or `"quantile"`.}
#'   \item{`hdci_width`}{The HDCI width (used only with `"hdci"`).}
#'   \item{`dgp`}{`"resample"` or `"conditional"`.}
#' }
#'
#' @details
#' The generating model has the CML item parameters (`psychotools`) and the
#' WLE person locations of `data`. In each iteration a dataset is simulated
#' by `dgp`, the model is refitted by CML (`psychotools::pcmodel()`), and the
#' conditional infit and outfit MSQ are computed with `iarm::out_infit()`.
#' A simulated dataset that cannot be refitted, for example because a
#' response category is empty, is left out. [RMitemInfit()] reports how
#' many were left out in its caption. If none can be refitted, the function
#' stops, naming any category used only by respondents with a zero or
#' perfect total score.
#'
#' The method and the default settings come from simulation studies of
#' detection (Johansson, 2025) and of iterations and multiplicity
#' correction (Johansson, 2026).
#'
#' Parallel processing needs the optional \pkg{mirai} package. The
#' \pkg{iarm} package must be installed (it is in Suggests).
#'
#' @references
#' Johansson, M. (2025). Detecting item misfit in Rasch models.
#' *Educational Methods & Psychometrics, 3*(18).
#' \doi{10.61186/emp.2025.5}
#'
#' Johansson, M. (2026). Simulation-based cutoffs for conditional item fit in
#' Rasch models: Iterations, multiplicity correction, and decision stability.
#' *PsyArXiv*. \doi{10.31234/osf.io/7pqz4_v2}
#'
#' @seealso [RMitemInfit()] for the test, [RMitemInfitPlot()] to plot the
#'   simulated values, [RMitemInfitCutoffMI()] for multiply imputed data.
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
#'   cutoff_res <- RMitemInfitCutoff(sim_data, iterations = 100,
#'                                   parallel = FALSE, seed = 42)
#'   cutoff_res$item_cutoffs
#'
#'   # Use the simulations in RMitemInfit()
#'   RMitemInfit(sim_data, cutoff = cutoff_res)
#' }
#' }
RMitemInfitCutoff <- function(
  data,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95,
  dgp = c("resample", "conditional")
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

  # Only complete cases (matching RMitemInfit behaviour). Record the raw
  # total and whether anything was dropped so callers (e.g. RMitemInfitPlot)
  # can report the sample in the standard `n = X of Y respondents` form.
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

  # Generating model: CML item thresholds (psychotools), computed once. The
  # conditional infit statistic is conditional on the total score, so the
  # "conditional" DGP (simulate each respondent's pattern given their observed
  # score) is the matched null. The "resample" DGP draws WLE person locations
  # with replacement and simulates parametrically (a marginal null).
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
    results_raw <- run_infit_sim_parallel(
      iterations,
      sim_seeds,
      sim_data_list,
      n_cores,
      verbose
    )
  } else {
    results_raw <- run_infit_sim_sequential(
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
  rownames(results_df) <- NULL

  # Compute per-item cutoffs
  item_names <- unique(results_df$Item)
  item_cutoffs <- do.call(
    rbind,
    lapply(item_names, function(item) {
      sub <- results_df[results_df$Item == item, ]
      if (cutoff_method == "hdci") {
        # ggdist::hdci() returns a matrix with ncol = 2: column 1 is the lower
        # bound, column 2 is the upper bound. Row 1 contains the continuous
        # interval.
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

#' Run a single infit simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMitemInfitCutoff()].
#' @return A data.frame with columns `Item`, `InfitMSQ`, `OutfitMSQ`, or a
#'   character string on failure.
#' @noRd
run_single_infit_sim <- function(seed, data_list) {
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

      # Conditional infit/outfit from a CML refit. psychotools::pcmodel() handles
      # both polytomous and dichotomous (a 2-category PCM is the Rasch model) and
      # is accepted by iarm::out_infit(); it matches eRm to ~1e-6 but is faster.
      # The polytomous validation above (every category present) is stricter
      # than the category check in .cml_fit(). The dichotomous one only needs 8
      # positive responses, so an item answered 1 by everyone passes it and
      # stops in the check instead; the iteration is discarded either way.
      model_fit <- .cml_fit(sim_df)
      cfit <- iarm::out_infit(model_fit)

      data.frame(
        Item = data_list$item_names,
        InfitMSQ = round(cfit$Infit, 3),
        OutfitMSQ = round(cfit$Outfit, 3),
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

#' Run infit simulations in parallel using mirai
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param n_cores Number of mirai daemons.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @noRd
run_infit_sim_parallel <- function(
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
      .sim_conditional = .sim_conditional,
      .esf_convolve = .esf_convolve,
      .cml_fit = .cml_fit,
      .check_response_categories = .check_response_categories
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
#' @noRd
run_infit_sim_sequential <- function(
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
