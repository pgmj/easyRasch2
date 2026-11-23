#' Category fit test given the total score
#'
#' Tests, for each item, whether its response categories work as the
#' partial credit model predicts given the respondent's total score. Each
#' item with categories \eqn{0, \dots, m} has \eqn{m} cumulative events
#' \eqn{X_i \ge k} (a response in category \eqn{k} or higher), one per
#' step between adjacent categories. The model's prediction for an event,
#' \eqn{P(X_i \ge k \mid r)}, depends on the item parameters alone, so no
#' person locations are estimated. The test detects categories that
#' respondents do not distinguish, categories used in the wrong order, and
#' steps that separate respondents more or less than the scoring assumes,
#' also when the threshold estimates are ordered. [RMitemCatFitPlot()]
#' shows where an item departs.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Respondents with missing responses are left out (complete cases),
#'   since the total score must cover the same items for everyone.
#' @param iterations Integer, at least 20. Number of datasets simulated
#'   from the model. Default `1000L`. Below 400 a notice is given. The
#'   smallest possible p is 1 / (`iterations_used` + 1).
#' @param correction Character. Correction of the item p-values across
#'   items: `"fwer"` (default, Westfall-Young step-down, controls the
#'   family-wise error rate), `"fdr_bh"` or `"fdr_by"` (Benjamini-Hochberg
#'   or Benjamini-Yekutieli, control the false discovery rate), or `"none"`.
#'   The FDR corrections need more iterations than `"fwer"`. A warning says
#'   when there are too few for a single item to reach .05 after
#'   correction.
#' @param item_labels Optional character vector of unique item labels, one
#'   per column of `data`. Default: column names.
#' @param parallel,n_cores Logical / integer. Run the simulations in
#'   parallel via \pkg{mirai} when available. Default sequential. When
#'   `n_cores` is `NULL`, `getOption("mc.cores")` is used. If that is unset
#'   too, a warning is given and the simulations run sequentially.
#' @param seed Optional integer for reproducible simulations. See
#'   [easyRasch2-reproducibility].
#' @param output Character. `"kable"` (default), `"dataframe"`, or
#'   `"list"`. The list also holds the simulations, so that
#'   [RMitemCatFitPlot()] can draw every plot from it (`fit` argument)
#'   without simulating again.
#'
#' @details
#' **Data.** Respondents with the lowest or highest possible total score
#' are left out, since their responses are fixed by the score. Every
#' category of every item must be used by the remaining respondents. A
#' category that is empty, or used only by respondents with an extreme
#' total score, stops the function with an error that names the item,
#' since its threshold cannot be estimated.
#'
#' **Statistic.** For each item i and event \eqn{X_i \ge k}, let
#' \eqn{y_v} be 1 when respondent v's response is in category k or higher
#' and 0 otherwise, and \eqn{p_v = P(X_i \ge k \mid r_v)} the model's
#' prediction. The observed \eqn{y_v} are regressed on \eqn{p_v} by
#' pool-adjacent-violators (PAV), which gives \eqn{\hat c_v}, the best
#' non-decreasing fit. The statistic is the miscalibration component of
#' the Brier score in the CORP decomposition (Dimitriadis et al., 2021),
#' \deqn{\mathrm{MCB} = \frac{1}{n}\sum_v (p_v - y_v)^2 -
#'   \frac{1}{n}\sum_v (\hat c_v - y_v)^2,}
#' where n is the number of respondents used. MCB is zero when the PAV
#' curve lies on the diagonal of the calibration plot in
#' [RMitemCatFitPlot()] and grows with its distance from it.
#'
#' **Reference distribution.** Response patterns are simulated from the
#' model given each respondent's total score, at the item parameters
#' estimated by conditional maximum likelihood, and the parameters are
#' re-estimated in every simulated dataset. The fit of the prediction to the
#' same data is therefore part of the reference distribution. A simulated
#' dataset in which some category of some item is empty cannot be refitted.
#' It is left out of the reference distribution, and `iterations_used`
#' gives the number of simulated datasets that remain.
#'
#' **p-values.** For each event, p is the share of the observed MCB and the
#' simulated ones that lie at or above the observed value, that is
#' (1 + the number of simulated values at or above it) divided by
#' (`iterations_used` + 1). These event p-values are not corrected for the
#' several events of an item and serve to show which step departs. The test of an item is
#' its item p. To get it, each event's MCB is studentised by the mean and
#' SD of the observed and simulated values together, and the item statistic
#' is the largest of these. The item p is the share of the observed and
#' simulated datasets whose item statistic is at or above the observed one.
#' The corrected item p (`padj`) adjusts it across items. With `"fwer"`
#' this is the Westfall-Young step-down on the item statistics (Ge et al.,
#' 2003; Ferreira, 2024). With `"none"`, `padj` equals `p`.
#'
#' **Calibration and power.** In a simulation study (5 or 10 items with
#' five categories, n = 250 to 1000, 199 simulated datasets), when the
#' model held, the item p was below .05 in .039 to .051 of items, and the
#' family-wise rate over items was .030 to .074 per condition (.041 over
#' all of them). A check of the condition with the highest rate (n = 1000,
#' 5 items) with 999 simulated datasets gave .052 with ordered and .045
#' with disordered thresholds (Monte Carlo SE .007). The test detected
#' interchangeable or swapped categories and a step that discriminated
#' twice as strongly as the others, also when the thresholds were ordered.
#' It rarely detected "don't know" use in up to 10% of responses or an
#' extreme response style shared by all items. The study is described at
#' \url{https://pgmj.github.io/rasch_book/chapters/category-fit.html}.
#'
#' **What a flag means.** The test compares an item's whole table of
#' category probabilities given the total score with the model, so it is
#' not specific to the categories. An item whose categories work as scored
#' but which discriminates less or more than the other items (underfit or
#' overfit, for example from a second dimension or local dependence)
#' departs as well. In a check with 10 items, one of which had ordered,
#' distinct categories but a slope of 0.5 or 2, that item was flagged in
#' 36% and 8% of datasets at n = 300, and in 96% and 83% at n = 1000.
#' Check item fit first, for example with [RMitemInfit()] or
#' [RMitemRestscore()]. Read a flag as a problem with the categories when
#' the item otherwise fits, or when the departure is in one event and not
#' in all. In [RMitemCatFitPlot()], a misfitting item shows in every
#' category, while a category problem is local to one category or one pair
#' of adjacent categories.
#'
#' @return
#' * `output = "kable"`: a `knitr::kable()` table, one row per item, with
#'   the event p-values, the item p and the corrected item p.
#' * `output = "dataframe"`: a data.frame with columns `item`, `mcb_1`,
#'   ..., `mcb_m` (observed MCB per event, where m is the largest number of
#'   events of any item, `NA` for events an item does not have), `p_1`,
#'   ..., `p_m` (uncorrected event p), `p` (item p) and `padj` (item p
#'   corrected across items).
#' * `output = "list"`: a list of class `RMitemCatFit` with `table` (the
#'   data.frame above), `correction`, `iterations` (requested),
#'   `iterations_used` (refitted), `n_used`, `n_total`, `qualifiers`,
#'   `item_labels`, and `prep` and `sims`, the data and simulations that
#'   [RMitemCatFitPlot()] uses, and `kable` (the table from
#'   `output = "kable"`).
#'
#' @references
#' Dimitriadis, T., Gneiting, T., & Jordan, A. I. (2021). Stable reliability
#' diagrams for probabilistic classifiers. *Proceedings of the National
#' Academy of Sciences, 118*(8), e2016191118.
#' \doi{10.1073/pnas.2016191118}
#'
#' Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate' number
#' of hypotheses simultaneously. *Journal of Statistical Theory and Practice,
#' 19*(1), 6. \doi{10.1007/s42519-024-00412-4}
#'
#' Ge, Y., Dudoit, S., & Speed, T. P. (2003). Resampling-based multiple
#' testing for microarray data analysis. *Test, 12*(1), 1--77.
#' \doi{10.1007/BF02595811}
#'
#' @seealso [RMitemCatFitPlot()] for the plots behind the test.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("eRm", quietly = TRUE)) {
#'   data(pcmdat2, package = "eRm")
#'   fit <- RMitemCatFit(pcmdat2, iterations = 200, seed = 1,
#'                       output = "list")
#'   fit$table
#'   RMitemCatFitPlot(pcmdat2, fit = fit)  # calibration plot
#' }
#' }
#'
#' @export
RMitemCatFit <- function(
  data,
  iterations = 1000L,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  item_labels = NULL,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  output = c("kable", "dataframe", "list")
) {
  correction <- match.arg(correction)
  output <- match.arg(output)
  .check_iterations(iterations, min = 20L)
  iterations <- as.integer(iterations)
  prep <- .catfit_prepare(data, who = "RMitemCatFit()")
  item_labels <- .catfit_item_labels(item_labels, prep$item_names)
  if (!is.null(seed)) set.seed(seed)

  sims <- .catfit_simulate(prep, iterations, c("mcb", "cep", "rest"),
                           parallel, n_cores)
  if (sims$refit_used == 0L) {
    stop("None of the simulated datasets could be refitted, usually because ",
         "an item ended up with an empty category. Check for rarely used ",
         "categories.", call. = FALSE)
  }
  obs <- .mcb_vector(prep$X, prep$r, prep$P)
  pv <- .catfit_pvalues(obs, sims$mcb, prep$event_item, correction)

  m_max <- max(lengths(prep$thr))
  wide <- function(v) {
    out <- matrix(NA_real_, length(prep$thr), m_max)
    out[cbind(prep$event_item, prep$event_k)] <- v
    out
  }
  M <- wide(obs)
  Pe <- wide(pv$p_event)
  tab <- data.frame(item = item_labels, stringsAsFactors = FALSE)
  for (k in seq_len(m_max)) tab[[paste0("mcb_", k)]] <- M[, k]
  for (k in seq_len(m_max)) tab[[paste0("p_", k)]] <- Pe[, k]
  tab$p <- pv$p_item
  tab$padj <- pv$padj

  if (sims$refit_used < 400L) {
    .notify_low_iterations(sims$refit_used, iterations, fn = "RMitemCatFit()",
                           id = "easyRasch2_catfit_low_iterations",
                           correction = correction)
  }
  .warn_fdr_floor(sims$refit_used, length(prep$thr), correction,
                  unit = "items", fn = "RMitemCatFit()")

  if (output == "dataframe") return(tab)

  disp <- data.frame(Item = tab$item, stringsAsFactors = FALSE)
  for (k in seq_len(m_max)) {
    disp[[paste0("X \u2265 ", k)]] <- .format_p_table(tab[[paste0("p_", k)]])
  }
  disp[["Item p"]] <- .format_p_table(tab$p)
  disp[["Adjusted p"]] <- .format_p_table(tab$padj)
  caption <- paste0(
    "Calibration of the cumulative events X \u2265 k given the total score. ",
    "Columns X \u2265 k: uncorrected p for each event. Item p: the events ",
    "combined (largest studentised value). Adjusted p: item p corrected ",
    "across items with the ", .correction_label(correction), ". Reference: ",
    sims$refit_used, " datasets simulated given each respondent's total ",
    "score, with item parameters re-estimated in each. ",
    .n_caption(prep$n_used, prep$n_total, prep$quals), ".",
    .attrition_clause(sims$refit_used, iterations),
    .iteration_note(sims$refit_used, correction)
  )
  kbl <- knitr::kable(disp, format = "pipe", caption = caption)
  if (output == "kable") return(kbl)

  structure(list(
    table = tab, correction = correction, iterations = iterations,
    iterations_used = sims$refit_used, n_used = prep$n_used,
    n_total = prep$n_total, qualifiers = prep$quals,
    item_labels = item_labels, prep = prep, sims = sims, kable = kbl
  ), class = c("RMitemCatFit", "list"))
}

# ---------------------------------------------------------------------
# Internal: shared by RMitemCatFit() and RMitemCatFitPlot()
# ---------------------------------------------------------------------

#' Data, CML fit and conditional probabilities for the category functions
#'
#' Complete cases (the total score must cover the same items for everyone),
#' CML thresholds, then respondents with an extreme total score removed.
#' @param data Response data.
#' @param who Subject of the category-check message.
#' @return List with `X` (non-extreme complete cases), `r` (their totals),
#'   `thr`, `P` (`.cond_cat_probs()`), `item_names`, `event_item` and
#'   `event_k` (item and event of every cumulative event, in the order used
#'   by `.mcb_vector()`), `n_used`, `n_total`, `quals`.
#' @keywords internal
#' @noRd
.catfit_prepare <- function(data, who = "This function") {
  validate_response_data(data)
  .check_complete_case_categories(data, who = who)
  item_names <- colnames(data)
  if (is.null(item_names)) item_names <- paste0("V", seq_len(ncol(data)))
  n_total <- nrow(data)
  has_na <- anyNA(data)
  X <- as.matrix(data)
  X <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(X) == 0L) {
    stop("No complete cases in data. All rows contain at least one NA.",
         call. = FALSE)
  }
  storage.mode(X) <- "integer"
  dimnames(X) <- list(NULL, item_names)
  # A category used only by extreme scorers is empty in the data CML uses:
  # its threshold diverges and no simulated dataset can be refitted.
  ec <- .extreme_only_categories(X)
  if (nrow(ec) > 0L) {
    stop(.category_error(
      paste0(.extreme_only_message(ec),
             if (has_na) {
               paste0(" (In the ", nrow(X), " respondents with complete ",
                      "responses, which ", who, " uses.)")
             },
             " ", who, " stops here."),
      unique(ec$Item)
    ))
  }
  thr <- .fit_cml_thresholds(X)
  K <- sum(lengths(thr))
  r_all <- rowSums(X)
  keep <- r_all > 0 & r_all < K
  X <- X[keep, , drop = FALSE]
  if (nrow(X) < 10L) {
    stop("Fewer than 10 respondents with a non-extreme total score.",
         call. = FALSE)
  }
  m <- lengths(thr)
  list(
    X = X, r = r_all[keep], thr = thr, P = .cond_cat_probs(thr),
    item_names = item_names,
    event_item = rep(seq_along(m), m),
    event_k = sequence(m),
    n_used = nrow(X), n_total = n_total,
    quals = c(if (has_na) "complete cases",
              if (any(!keep)) "extreme scores excluded")
  )
}

#' Check or default the item labels
#' @keywords internal
#' @noRd
.catfit_item_labels <- function(item_labels, item_names) {
  if (is.null(item_labels)) return(item_names)
  if (length(item_labels) != length(item_names)) {
    stop("`item_labels` must have the same length as ncol(data) (",
         length(item_names), ").", call. = FALSE)
  }
  if (anyDuplicated(item_labels) > 0L) {
    stop("`item_labels` must be unique.", call. = FALSE)
  }
  as.character(item_labels)
}

#' Simulate datasets given the total scores
#'
#' One dataset per iteration, drawn by `.sim_given_totals()` at the fitted
#' thresholds. Depending on `what`:
#' * `"cep"`: PAV curves of the simulated events at the observed
#'   predictions (the calibration band, no refit needed),
#' * `"mcb"`: MCB per event after re-estimating the thresholds,
#' * `"rest"`: rest-score mean and SD differences after re-estimating.
#' Datasets with an empty category cannot be refitted and give no `"mcb"`
#' or `"rest"`.
#' @return List with matrices `cep`, `mcb`, `rest` (statistics in rows,
#'   iterations in columns, `NULL` when not requested), `requested` and
#'   `refit_used`.
#' @keywords internal
#' @noRd
.catfit_simulate <- function(prep, iterations, what, parallel, n_cores) {
  it_seeds <- sample.int(.Machine$integer.max, iterations)
  ev <- if ("cep" %in% what) .event_forecasts(prep$P, prep$r) else NULL
  refit <- any(c("mcb", "rest") %in% what)
  thr <- prep$thr
  r <- prep$r
  item_names <- prep$item_names
  one <- function(s) {
    # Generator pinned as well as the seed, so parallel and sequential runs
    # agree (see easyRasch2-reproducibility).
    set.seed(s, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    Xb <- .sim_given_totals(thr, r)
    colnames(Xb) <- item_names
    out <- list(cep = NULL, mcb = NULL, rest = NULL, refit = FALSE)
    if ("cep" %in% what) out$cep <- .pav_replicate(Xb, ev)
    if (refit) {
      full <- all(vapply(seq_along(thr), function(i) {
        all(tabulate(Xb[, i] + 1L, nbins = length(thr[[i]]) + 1L) > 0L)
      }, logical(1L)))
      tb <- if (full) {
        tryCatch(.fit_cml_thresholds(Xb), error = function(e) NULL)
      }
      if (!is.null(tb)) {
        Pb <- .cond_cat_probs(tb)
        out$refit <- TRUE
        out$thr <- tb
        if ("mcb" %in% what) out$mcb <- .mcb_vector(Xb, r, Pb)
        if ("rest" %in% what) out$rest <- .restscore_stats(Xb, r, Pb)$diff
      }
    }
    out
  }
  res <- .maybe_parallel_lapply(it_seeds, one, parallel, n_cores)
  ok <- vapply(res, function(z) isTRUE(z$refit), logical(1L))
  grab <- function(nm, rows) {
    if (!nm %in% what) return(NULL)
    v <- lapply(res[rows], `[[`, nm)
    if (length(v) == 0L) return(NULL)
    do.call(cbind, v)
  }
  list(
    cep = grab("cep", rep(TRUE, length(res))),
    mcb = grab("mcb", ok),
    rest = grab("rest", ok),
    # Seed and re-estimated thresholds of each refitted dataset, so
    # RMitemCatFitPlot(type = "residuals") can regenerate the datasets for
    # any score grouping without refitting.
    seeds = if (refit) it_seeds[ok],
    thr = if (refit) lapply(res[ok], `[[`, "thr"),
    requested = iterations,
    refit_used = if (refit) sum(ok) else NA_integer_
  )
}

#' MCB of every cumulative event
#' @return Numeric vector, events of item 1 first, in the order of
#'   `prep$event_item` and `prep$event_k`.
#' @keywords internal
#' @noRd
.mcb_vector <- function(X, r, P) {
  ev <- .event_forecasts(P, r)
  unlist(lapply(seq_along(ev), function(i) {
    vapply(seq_len(ncol(ev[[i]])), function(k) {
      p <- ev[[i]][, k]
      y <- as.numeric(X[, i] >= k)
      pv <- .pav(p, y)
      cep <- pv$cep[match(p, pv$p)]
      mean((p - y)^2) - mean((cep - y)^2)
    }, numeric(1L))
  }))
}

#' Event, item and corrected p-values from the pooled MCB values
#'
#' The observed MCB and the B simulated ones form a pool of B + 1 that is
#' exchangeable when the model holds, and every quantity is computed on the
#' pool, so the observed and simulated datasets are treated alike:
#' * event p: share of the pool at or above the observed MCB,
#' * item statistic: each event's MCB studentised by the mean and SD of the
#'   pool, maximum over the item's events (as `.bootstrap_pvalues()` with
#'   `moments = "joint"`),
#' * item p: share of the pool at or above the observed item statistic,
#' * `"fwer"`: Westfall-Young step-down on the item statistics.
#' A single-step min-p over the events was used first. Its item p cannot go
#' below m / (B + 1) and its step-down not below about m k / (B + 1), since
#' one simulated dataset holds the extreme rank of every event. On the stored
#' MCB values of the simulation study the studentised maximum held its level
#' and detected as often or more (dev/category_conditional_null,
#' R/compare_combination.R and R/run_calibration_B999.R).
#' @param obs Observed MCB per event.
#' @param sim Simulated MCB, events in rows, datasets in columns.
#' @param event_item Item of every event.
#' @param correction Correction across items.
#' @return List with `p_event`, `p_item`, `padj`.
#' @keywords internal
#' @noRd
.catfit_pvalues <- function(obs, sim, event_item, correction) {
  pool <- cbind(obs, sim)
  N <- ncol(pool)
  # Event p: share of the pool at or above the observed value.
  p_event <- vapply(seq_len(nrow(pool)), function(e) {
    mean(pool[e, ] >= pool[e, 1L])
  }, numeric(1L))
  # Studentise every event over the pool. An event with no variation in the
  # pool carries no information and gets 0.
  Z <- t(apply(pool, 1L, function(t) {
    s <- stats::sd(t)
    if (!is.finite(s) || s < 1e-12) rep(0, length(t)) else (t - mean(t)) / s
  }))
  Z <- matrix(Z, nrow = nrow(pool))
  items <- sort(unique(event_item))
  Tm <- vapply(items, function(i) {
    apply(Z[event_item == i, , drop = FALSE], 2L, max)
  }, numeric(N))
  Tm <- matrix(Tm, nrow = N)
  p_item <- colMeans(sweep(Tm, 2L, Tm[1L, ], `>=`))
  padj <- switch(
    correction,
    none = p_item,
    fdr_bh = stats::p.adjust(p_item, method = "BH"),
    fdr_by = stats::p.adjust(p_item, method = "BY"),
    fwer = pmax(as.numeric(.wy_stepdown(Tm[1L, ], Tm[-1L, , drop = FALSE])),
                p_item)
  )
  list(p_event = p_event, p_item = p_item, padj = padj)
}
