#' Simulation-based item-restscore cutoffs
#'
#' Simulates datasets from a Rasch or partial credit model that fits and
#' resembles your data, and computes each item's observed minus expected
#' item-restscore gamma (Kreiner, 2011) in each. Passed to
#' [RMitemRestscore()], the result gives bootstrap p-values corrected across
#' items, and an interval per item that describes where a fitting item's
#' difference is expected to fall.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Only complete cases (rows without `NA`) are used.
#' @param iterations Integer. Number of simulated datasets. Default `400`.
#'   Use 1000 or more for a final analysis (see Reproducibility in
#'   [RMitemRestscore()]).
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
#'   The interval describes where a fitting item's difference is expected
#'   to fall. It is not a decision rule. Flagging every item outside a
#'   width-`w` interval tests all `k` items at once, with a family-wise
#'   error rate of `1 - w^k`. Decisions come from the corrected p-value in
#'   [RMitemRestscore()], as shown for infit in Johansson (2026).
#' @param dgp Character. How datasets are simulated. `"conditional"`
#'   (default) simulates each respondent's response pattern from the model
#'   given their observed total score, with the item parameters fixed (a
#'   conditional null). The expected gamma is computed from the observed
#'   score distribution, which the conditional null keeps. `"resample"`
#'   draws person locations from the WLE estimates with replacement and
#'   simulates responses from the model (a marginal null). With the Rasch
#'   model true, pooled over four designs, the family-wise rate with
#'   `correction = "fwer"` was 5.0% with the conditional null and 6.4% with
#'   the resample null. The conditional null takes about 1.2 to 1.6 times
#'   as long. \strong{Experimental.}
#'
#' @return A list with:
#' \describe{
#'   \item{`results`}{data.frame with the columns `iteration`, `Item`,
#'     `Observed`, `Expected` and `Difference`, one row per item and
#'     simulated dataset. `Difference` is `Observed - Expected`, the
#'     statistic [RMitemRestscore()] tests.}
#'   \item{`item_cutoffs`}{data.frame with the interval per item: `Item`,
#'     `diff_low`, `diff_high`, computed by `cutoff_method`.}
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
#'   \item{`dgp`}{`"conditional"` or `"resample"`.}
#' }
#'
#' @details
#' **Why simulate.** The asymptotic test in `iarm::item_restscore()` divides
#' the observed minus expected gamma by the standard error of the observed
#' gamma alone. The expected gamma is estimated from the same data and
#' correlates with the observed one, so that standard error is too large for
#' the difference, and the difference is also biased upwards in small
#' samples. With the Rasch model true, the test is liberal for dichotomous
#' items in small or mistargeted samples and conservative for polytomous
#' items, and it flags too few underfitting items in every case. The
#' simulated reference distribution replaces the asymptotic one and absorbs
#' both problems. The simulation study of the asymptotic test is reported at
#' \url{https://pgmj.github.io/rasch_book/chapters/restscore-asymptotic-null.html}.
#'
#' **Simulation.** The generating model has the CML item parameters
#' (`psychotools`) and the WLE person locations of `data`. In each
#' iteration a dataset is simulated by `dgp`, the model is refitted by CML
#' (`psychotools::pcmodel()`), and the observed and expected gamma are
#' computed as in `iarm::item_restscore()`, by a faster internal routine
#' that skips the standard errors and the rounding of the printed values.
#' The refit matters. The expected gamma varies from sample to sample
#' because the thresholds do, and holding them fixed would reproduce the
#' problem the simulation is meant to solve. A simulated dataset that cannot
#' be refitted, for example because a response category is empty, is left
#' out, and [RMitemRestscore()] reports how many in its caption.
#'
#' Parallel processing needs the optional \pkg{mirai} package. The
#' \pkg{iarm} package must be installed (it is in Suggests).
#'
#' @references
#' Johansson, M. (2026). Simulation-based cutoffs for conditional item fit in
#' Rasch models: Iterations, multiplicity correction, and decision stability.
#' *PsyArXiv*. \doi{10.31234/osf.io/7pqz4_v2}
#'
#' Kreiner, S. (2011). A Note on Item-Restscore Association in Rasch Models.
#' *Applied Psychological Measurement, 35*(7), 557-561.
#' \doi{10.1177/0146621611410227}
#'
#' @seealso [RMitemRestscore()] for the test, [RMitemRestscorePlot()] to
#'   plot the simulated values.
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
#' @noRd
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
      # The polytomous validation above (every category present) is stricter
      # than the category check in .cml_fit(). The dichotomous one only needs 8
      # positive responses, so an item answered 1 by everyone passes it and
      # stops in the check instead; the iteration is discarded either way.
      model_fit <- .cml_fit(sim_df, hessian = FALSE)
      res_mat <- .restscore_gamma(
        sim_df,
        lapply(psychotools::threshpar(model_fit), cumsum)
      )
      observed <- res_mat[, "observed"]
      expected <- res_mat[, "expected"]

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
# Internal: observed and expected item-restscore gamma
# ---------------------------------------------------------------------------

#' Observed and expected item-restscore gamma
#'
#' Computes the two gammas that `iarm::item_restscore()` reports, without its
#' standard errors and p-values, which the bootstrap does not use. iarm's
#' `pscore_poly()` recomputes two full sets of elementary symmetric functions
#' for every cell of every item's restscore-by-category table. Here they are
#' computed once for all items and once per item with that item left out,
#' which makes the call about 50 times faster. The values agree with iarm's
#' unrounded ones to machine precision. iarm itself returns them through
#' `format(digits = 3)`, so they are rounded by an amount that depends on the
#' smallest p-value in the same table.
#'
#' @param X Complete-case integer response matrix, scored from 0.
#' @param coeff List of cumulative threshold parameters per item, as
#'   `lapply(psychotools::threshpar(fit), cumsum)`.
#' @return Matrix with columns `observed` and `expected`, one row per item.
#' @keywords internal
#' @noRd
.restscore_gamma <- function(X, coeff) {
  X <- as.matrix(X)
  storage.mode(X) <- "integer"
  k <- ncol(X)
  mi <- vapply(coeff, length, integer(1L))
  # pcmodel() drops the threshold of an unobserved middle category, so an
  # item's top score can exceed its threshold count; the tables below would
  # then be silently truncated.
  if (!identical(unname(apply(X, 2L, max)), unname(mi))) {
    stop(
      "Each item's highest response must equal its number of thresholds.",
      call. = FALSE
    )
  }
  m <- sum(mi)
  score <- rowSums(X)
  score_n <- tabulate(score + 1L, nbins = m + 1L)
  g_all <- psychotools::elementary_symmetric_functions(coeff)[[1L]]

  observed <- expected <- numeric(k)
  for (i in seq_len(k)) {
    nx <- mi[i] + 1L
    nr <- m - mi[i] + 1L
    # Observed restscore (rows) by item response (columns)
    N <- matrix(
      tabulate((score - X[, i]) * nx + X[, i] + 1L, nbins = nr * nx),
      ncol = nx,
      byrow = TRUE
    )
    observed[i] <- .gk_gamma(N)

    # Expected table: score-group sizes times the conditional probability of
    # each response given the total score r + x,
    # P(x | r + x) = w_x * gamma_r(rest) / gamma_{r+x}(all).
    g_rest <- psychotools::elementary_symmetric_functions(coeff[-i])[[1L]]
    w <- exp(-c(0, coeff[[i]]))
    tot <- outer(0:(nr - 1L), 0:mi[i], `+`)
    pmat <- outer(g_rest, w) / matrix(g_all[tot + 1L], ncol = nx)
    nmat <- matrix(score_n[tot + 1L], ncol = nx)
    expected[i] <- .gk_gamma(nmat * pmat)
  }
  cbind(observed = observed, expected = expected)
}

#' Goodman-Kruskal gamma of a two-way table
#'
#' Rows and columns are both taken as ordered. Each unordered pair of
#' observations is counted once, as in `.partgam_one()`.
#'
#' @param N Numeric matrix of (possibly non-integer) counts.
#' @return Gamma, or `NA_real_` when there are no concordant and no discordant
#'   pairs.
#' @keywords internal
#' @noRd
.gk_gamma <- function(N) {
  upper <- function(d) outer(seq_len(d), seq_len(d), function(a, b) b > a) + 0
  A <- upper(nrow(N)) %*% N
  Gc <- upper(ncol(N))
  conc <- sum(N * (A %*% t(Gc)))
  disc <- sum(N * (A %*% Gc))
  if (conc + disc == 0) NA_real_ else (conc - disc) / (conc + disc)
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
#' @noRd
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
      .esf_convolve = .esf_convolve,
      .restscore_gamma = .restscore_gamma,
      .gk_gamma = .gk_gamma,
      .cml_fit = .cml_fit,
      .check_response_categories = .check_response_categories
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
#' @noRd
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
