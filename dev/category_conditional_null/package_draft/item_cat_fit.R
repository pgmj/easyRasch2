#' Category Fit Test Given the Total Score
#'
#' Tests, for each item, whether its cumulative response events
#' \eqn{X_i \ge k} depend on the total score as the partial credit model
#' predicts. The prediction is \eqn{P(X_i \ge k \mid r)}, which under the
#' model depends on the item parameters alone, so no person locations are
#' estimated. The test detects categories that respondents do not
#' distinguish, categories used in the wrong order, and steps that separate
#' respondents more or less than the scoring assumes, including when the
#' threshold estimates are ordered. [RMitemCatFitPlot()] shows where an
#' item departs.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Respondents with missing responses are left out (complete cases),
#'   since the total score must cover the same items for everyone.
#' @param iterations Integer. Number of datasets simulated from the model.
#'   Default `1000L`. The p-values are counts out of `iterations + 1`.
#' @param correction Character. Correction of the item p-values across
#'   items: `"fwer"` (default, Westfall-Young step-down, controls the
#'   family-wise error rate), `"fdr_bh"`, `"fdr_by"` (false discovery rate)
#'   or `"none"`.
#' @param item_labels Optional character vector of item labels, one per
#'   column of `data`. Default: column names.
#' @param parallel,n_cores Logical / integer. Run the simulations in
#'   parallel via \pkg{mirai} when available. Default sequential. When
#'   `n_cores` is `NULL`, `getOption("mc.cores")` is used, and two workers
#'   when that is unset too.
#' @param seed Optional integer for reproducible simulations. See
#'   [easyRasch2-reproducibility].
#' @param output Character. `"kable"` (default), `"dataframe"`, or
#'   `"list"`. The list also holds the simulations, so that
#'   [RMitemCatFitPlot()] can draw every plot from it (`fit` argument)
#'   without simulating again.
#'
#' @details
#' **Statistic.** For each item i and event \eqn{X_i \ge k}, the observed
#' events are regressed on the predicted probability
#' \eqn{p = P(X_i \ge k \mid r)} by pool-adjacent-violators (PAV). The
#' statistic is the miscalibration component of the Brier score in the
#' CORP decomposition (Dimitriadis et al., 2021),
#' \deqn{\mathrm{MCB} = \frac{1}{n}\sum_v (p_v - y_v)^2 -
#'   \frac{1}{n}\sum_v (\hat c_v - y_v)^2,}
#' where \eqn{\hat c_v} is the PAV-recalibrated probability. MCB is zero
#' when the PAV curve lies on the diagonal of the calibration plot in
#' [RMitemCatFitPlot()] and grows with its distance from it. Respondents
#' with the lowest or highest possible total score are left out, since their
#' responses are fixed by the score.
#'
#' **Reference distribution.** Response patterns are simulated from the
#' model given each respondent's total score, at the item parameters
#' estimated by conditional maximum likelihood, and the parameters are
#' re-estimated in every simulated dataset. The fit of the prediction to the
#' same data is therefore part of the reference distribution. Datasets in
#' which some category of some item is empty cannot be refitted and are
#' dropped.
#'
#' **p-values.** For each event, p is the share of the observed MCB and the
#' simulated ones that lie at or above the observed value. To combine the
#' events of an item, each event's MCB is studentised by the mean and SD of
#' the observed and simulated values together, and the item statistic is the
#' largest of them. The item p is the share of the observed and simulated
#' datasets whose item statistic is at or above the observed one. With
#' `"fwer"` the item p-values are corrected across items by the
#' Westfall-Young step-down on the item statistics (Ge et al., 2003).
#'
#' In a simulation study (5 or 10 items with five categories, n = 250 to
#' 1000, 200 iterations) the item p rejected at a rate between .039 and
#' .051 at .05 when the model held, and the family-wise rate over items was
#' .030 to .074 per condition (.041 over all of them). The test detected
#' interchangeable or swapped categories and a step that discriminated twice
#' as strongly as the others, also when the thresholds were ordered. It
#' rarely detected "don't know" use in up to 10% of responses or an extreme
#' response style shared by all items.
#'
#' @return
#' * `output = "kable"`: a `knitr::kable()` table, one row per item, with
#'   the event p-values, the item p and the corrected item p.
#' * `output = "dataframe"`: a data.frame with columns `item`, `mcb_1`,
#'   ..., `mcb_m` (observed MCB per event, `NA` for events an item does not
#'   have), `p_1`, ..., `p_m` (uncorrected event p), `p` (item p) and
#'   `padj` (item p corrected across items).
#' * `output = "list"`: a list of class `RMitemCatFit` with `table` (the
#'   data.frame above), `correction`, `iterations` (requested),
#'   `iterations_used` (refitted), `n_used`, `n_total`, `qualifiers`,
#'   `item_labels`, and `prep` and `sims`, the data and simulations that
#'   [RMitemCatFitPlot()] uses.
#'
#' @references
#' Dimitriadis, T., Gneiting, T., & Jordan, A. I. (2021). Stable reliability
#' diagrams for probabilistic classifiers. *Proceedings of the National
#' Academy of Sciences, 118*(8), e2016191118.
#' \doi{10.1073/pnas.2016191118}
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
#'   RMitemCatFitPlot(pcmdat2, type = "calibration", fit = fit)
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
  if (!is.numeric(iterations) || length(iterations) != 1L ||
        !is.finite(iterations) || iterations < 20) {
    stop("`iterations` must be a single integer >= 20.", call. = FALSE)
  }
  iterations <- as.integer(iterations)
  prep <- .catfit_prepare(data)
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
  if (output == "list") {
    return(structure(list(
      table = tab, correction = correction, iterations = iterations,
      iterations_used = sims$refit_used, n_used = prep$n_used,
      n_total = prep$n_total, qualifiers = prep$quals,
      item_labels = item_labels, prep = prep, sims = sims
    ), class = c("RMitemCatFit", "list")))
  }

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
  knitr::kable(disp, format = "pipe", caption = caption)
}

# ---------------------------------------------------------------------
# Internal: shared by RMitemCatFit() and RMitemCatFitPlot()
# ---------------------------------------------------------------------

#' Data, CML fit and conditional probabilities for the category functions
#'
#' Complete cases (the total score must cover the same items for everyone),
#' CML thresholds, then respondents with an extreme total score removed.
#' @param data Response data.
#' @return List with `X` (non-extreme complete cases), `r` (their totals),
#'   `thr`, `P` (`.cond_cat_probs()`), `item_names`, `event_item` and
#'   `event_k` (item and event of every cumulative event, in the order used
#'   by `.mcb_vector()`), `n_used`, `n_total`, `quals`.
#' @keywords internal
#' @noRd
.catfit_prepare <- function(data) {
  validate_response_data(data)
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
#' R/compare_combination.R).
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
