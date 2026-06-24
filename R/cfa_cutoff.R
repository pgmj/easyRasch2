#' Simulated null distribution for one-factor CFA fit and loadings under PCM unidimensionality
#'
#' Generates a parametric-bootstrap null distribution against which observed
#' one-factor categorical-CFA results can be compared. The simulation draws
#' \code{iterations} datasets from the fitted PCM (or RM, for dichotomous
#' data) using the observed item parameters and a resampled person
#' distribution; each simulated dataset is fitted with
#' \code{lavaan::cfa(..., ordered = TRUE, estimator = "WLSMV")} and both the
#' three fit indices (CFI, RMSEA, SRMR) and the per-item standardized factor
#' loadings are recorded. Because the simulated data satisfy the PCM
#' unidimensional assumption by construction, the resulting distributions are
#' the "expected" reference for what a correctly fitting unidimensional model
#' produces at this sample size and item structure.
#'
#' This function only generates the simulated reference. To obtain the
#' observed-vs-expected tables, pass its result to \code{\link{RMdimCFA}}; for
#' the figures, pass it to \code{\link{RMdimCFAPlot}}.
#'
#' @param data A data.frame or matrix of item responses (non-negative
#'   integers, 0-based). One column per item, one row per person.
#' @param iterations Integer. Number of parametric-bootstrap iterations.
#'   Default `250`.
#' @param percentile Numeric in (50, 100). The strictness of the cutoffs.
#'   Default `99`. Fit indices use a one-sided cutoff in the unfavourable
#'   direction (CFI from below; RMSEA and SRMR from above); standardized
#'   loadings use a two-sided central interval covering `percentile`% of the
#'   simulated distribution (an item is flagged when its observed loading
#'   falls in the outer `100 - percentile`%, split across the two tails).
#' @param output Character. Only `"list"` (the default) is supported and the
#'   function always returns the simulation object. `"kable"` is retained only
#'   to raise an informative error: tables now come from
#'   \code{\link{RMdimCFA}}.
#' @param parallel Logical. If `TRUE` (default), uses parallel processing
#'   via `mirai`. Falls back to sequential if `mirai` is not installed or
#'   `n_cores` cannot be resolved.
#' @param n_cores Integer or `NULL`. Number of parallel workers. When
#'   `NULL`, `getOption("mc.cores")` is consulted; if neither is set,
#'   sequential is used.
#' @param verbose Logical. Show a progress bar (default `FALSE`).
#' @param seed Integer or `NULL`. Master seed for reproducibility.
#' @param estimator Character. The lavaan estimator passed to
#'   `lavaan::cfa()`. Default `"WLSMV"`. Other limited-information
#'   estimators that produce robust/scaled fit indices (e.g.,
#'   `"DWLS"`, `"ULSMV"`) are also accepted; full-information ML is
#'   rejected (incompatible with `ordered = TRUE`).
#'
#' @return A list (the simulation object), with components:
#' \describe{
#'   \item{`simulated`}{data.frame with one row per successful iteration
#'     and columns `iteration`, `cfi`, `rmsea`, `srmr`.}
#'   \item{`simulated_loadings`}{data.frame with one row per successful
#'     iteration: an `iteration` column followed by one column per item
#'     holding the simulated standardized loading.}
#'   \item{`percentile`}{Numeric: the strictness setting used.}
#'   \item{`cutoffs`}{Named numeric vector (`cfi`, `rmsea`, `srmr`) of
#'     one-sided fit-index cutoffs at the chosen percentile.}
#'   \item{`loading_cutoffs`}{data.frame `Item`, `low`, `high` — the
#'     two-sided expected loading interval per item.}
#'   \item{`actual_iterations`}{Number of successful MC iterations.}
#'   \item{`sample_n`}{Number of complete cases used.}
#'   \item{`n_items`}{Number of items.}
#'   \item{`item_names`}{Character vector of item names.}
#'   \item{`is_polytomous`}{Logical: was a PCM (vs RM) fitted?}
#'   \item{`estimator`}{The lavaan estimator used.}
#' }
#'
#' @details
#' \strong{Generative model.} The data-generating process for each
#' simulated dataset is the PCM (or RM) fitted to the observed data,
#' with persons drawn from the empirical theta distribution
#' (resampled with replacement). This means the simulated data perfectly
#' satisfy the PCM unidimensional assumption.
#'
#' \strong{Estimation model.} The CFA on each simulated dataset uses a
#' single-factor model with all items as ordinal indicators
#' (`F1 =~ I1 + I2 + ...`), fitted with `WLSMV` by default. Reported
#' CFI / RMSEA are the Satorra-Bentler-scaled variants (`cfi.scaled`,
#' `rmsea.scaled`) for consistency across iterations; SRMR is reported
#' unchanged. Standardized loadings are the `est.std` of the `=~` paths
#' from `lavaan::standardizedSolution()`.
#'
#' \strong{Why a null distribution.} A perfectly PCM-unidimensional
#' dataset will typically not yield CFA fit indices at their ideal
#' values (CFI = 1, RMSEA = 0), nor identical loadings across items: PCM
#' uses a logistic threshold structure while WLSMV uses a probit link via
#' the polychoric correlation matrix, and finite samples add sampling
#' variability. The simulated distributions capture both, giving a more
#' honest reference than rule-of-thumb cutoffs derived under continuous-data
#' ML.
#'
#' \strong{Iteration failures.} Some simulated datasets cause WLSMV to
#' fail (non-positive-definite polychoric matrix, boundary thresholds,
#' empty categories). Failed iterations are dropped; `actual_iterations`
#' reflects the number that succeeded.
#'
#' @references
#' Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods
#' for mean and covariance structure analysis with nonnormal missing
#' data. \emph{Sociological Methodology, 30}(1), 165-200.
#' \doi{10.1111/0081-1750.00078}
#'
#' Rosseel, Y. (2012). lavaan: An R Package for Structural Equation
#' Modeling. \emph{Journal of Statistical Software, 48}(2), 1-36.
#' \doi{10.18637/jss.v048.i02}
#'
#' @seealso \code{\link{RMdimCFA}}, \code{\link{RMdimCFAPlot}},
#'   \code{\link{RMdimResidualPCA}}, \code{\link{RMdimMartinLof}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   data("raschdat1", package = "eRm")
#'
#'   # Few iterations for a fast example; use 250+ in real analyses
#'   sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 50,
#'                         parallel = FALSE, seed = 1)
#'
#'   # Observed-vs-expected tables
#'   RMdimCFA(raschdat1[, 1:8], cutoff = sim)
#'
#'   if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plots <- RMdimCFAPlot(sim, data = raschdat1[, 1:8])
#'     plots$loadings
#'     plots$fit
#'   }
#' }
#' }
#'
#' @export
RMdimCFACutoff <- function(data,
                        iterations = 250L,
                        percentile = 99,
                        output     = c("list", "kable"),
                        parallel   = TRUE,
                        n_cores    = NULL,
                        verbose    = FALSE,
                        seed       = NULL,
                        estimator  = "WLSMV") {

  output <- match.arg(output)
  if (output == "kable") {
    stop("RMdimCFACutoff() now returns the simulation object only (a list).\n",
         "Use RMdimCFA(data, cutoff = <this object>) for the fit and loadings ",
         "tables, and RMdimCFAPlot(<this object>, data) for the figures.",
         call. = FALSE)
  }

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for RMdimCFACutoff(). ",
         "Install with: install.packages(\"lavaan\")",
         call. = FALSE)
  }

  if (!is.numeric(percentile) || length(percentile) != 1L ||
      !is.finite(percentile) || percentile <= 50 || percentile >= 100) {
    stop("`percentile` must be a single numeric in (50, 100). ",
         "Common choices: 95, 99, 99.5.", call. = FALSE)
  }

  estimator <- toupper(estimator[1L])
  if (estimator %in% c("ML", "MLR", "MLM", "MLF")) {
    stop("`estimator = \"", estimator, "\"` is not appropriate for ",
         "ordinal-CFA. Use a limited-information estimator such as ",
         "\"WLSMV\", \"DWLS\", or \"ULSMV\".", call. = FALSE)
  }

  validate_response_data(data)

  # rgl workaround
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  data <- stats::na.omit(as.data.frame(data))
  if (nrow(data) == 0L) {
    stop("No complete cases in `data`.", call. = FALSE)
  }
  if (ncol(data) < 3L) {
    stop("RMdimCFACutoff() requires at least 3 items for a one-factor CFA.",
         call. = FALSE)
  }

  # Parallel setup -- mirrors RMdimResidualPCACutoff()
  use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)
  if (parallel && !use_parallel) {
    message("Install 'mirai' for parallel processing: install.packages(\"mirai\")")
    message("Running sequentially...")
  }
  if (use_parallel) {
    if (is.null(n_cores)) n_cores <- getOption("mc.cores")
    if (is.null(n_cores)) {
      warning(
        "For parallel processing, specify n_cores or set options(mc.cores = N).\n",
        "Falling back to sequential.",
        call. = FALSE
      )
      use_parallel <- FALSE
    } else {
      n_cores <- min(n_cores, iterations)
    }
  }

  if (!is.null(seed)) set.seed(seed)
  sim_seeds <- sample.int(.Machine$integer.max, iterations)

  data_mat       <- as.matrix(data)
  sample_n       <- nrow(data_mat)
  is_polytomous  <- max(data_mat, na.rm = TRUE) > 1L
  item_names_vec <- colnames(data_mat)
  if (is.null(item_names_vec))
    item_names_vec <- paste0("V", seq_len(ncol(data_mat)))

  # Generating Rasch model: CML item thresholds (psychotools) + WLE person
  # locations, consistent with the rest of the package. The DGP is unchanged
  # (thetas resampled with replacement; data simulated under the model). The
  # CFA itself is fitted by lavaan; `estimator` is the lavaan estimator.
  thr_list   <- .fit_cml_thresholds(data_mat)
  wle_thetas <- .estimate_thetas(data_mat, thr_list, method = "WLE")$theta
  wle_thetas <- wle_thetas[is.finite(wle_thetas)]

  sim_data_list <- list(
    type       = if (is_polytomous) "polytomous" else "dichotomous",
    thetas     = wle_thetas,
    n_items    = ncol(data_mat),
    sample_n   = sample_n,
    item_names = item_names_vec,
    estimator  = estimator
  )
  if (is_polytomous) {
    sim_data_list$deltaslist <- thr_list
  } else {
    sim_data_list$item_params <- unlist(thr_list, use.names = FALSE)
  }

  # Run the simulation
  if (use_parallel) {
    results_raw <- run_cfa_sim_parallel(iterations, sim_seeds,
                                        sim_data_list, n_cores, verbose)
  } else {
    results_raw <- run_cfa_sim_sequential(iterations, sim_seeds,
                                          sim_data_list, verbose)
  }

  ok         <- vapply(results_raw, is.list, logical(1L))
  successful <- results_raw[ok]

  if (length(successful) == 0L) {
    failed_msgs <- unlist(results_raw[!ok])
    sample_msg <- if (length(failed_msgs) > 0L) {
      unique(failed_msgs)[1L]
    } else {
      "(no message captured)"
    }
    stop("All CFA simulation iterations failed. Example: ", sample_msg,
         call. = FALSE)
  }

  actual_iterations <- length(successful)

  fit_mat  <- do.call(rbind, lapply(successful, `[[`, "fit"))
  load_mat <- do.call(rbind, lapply(successful, `[[`, "loadings"))
  colnames(fit_mat)  <- c("cfi", "rmsea", "srmr")
  colnames(load_mat) <- item_names_vec

  simulated_df <- data.frame(
    iteration = seq_len(actual_iterations),
    cfi       = as.numeric(fit_mat[, "cfi"]),
    rmsea     = as.numeric(fit_mat[, "rmsea"]),
    srmr      = as.numeric(fit_mat[, "srmr"]),
    stringsAsFactors = FALSE
  )

  simulated_loadings <- data.frame(
    iteration = seq_len(actual_iterations),
    as.data.frame(load_mat, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cutoffs         <- compute_cfa_cutoffs(simulated_df, percentile)
  loading_cutoffs <- compute_cfa_loading_cutoffs(simulated_loadings,
                                                 item_names_vec, percentile)

  list(
    simulated          = simulated_df,
    simulated_loadings = simulated_loadings,
    percentile         = percentile,
    cutoffs            = cutoffs,
    loading_cutoffs    = loading_cutoffs,
    actual_iterations  = actual_iterations,
    sample_n           = sample_n,
    n_items            = ncol(data_mat),
    item_names         = item_names_vec,
    is_polytomous      = is_polytomous,
    estimator          = estimator
  )
}

# ===========================================================================
# Internal: cutoff / flag helpers
# ===========================================================================

#' Compute one-sided CFA-fit cutoffs at the given percentile
#'
#' @keywords internal
#' @noRd
compute_cfa_cutoffs <- function(simulated_df, percentile) {
  pct <- percentile / 100
  # Filter to finite values -- lavaan can return Inf for RMSEA when
  # chi-square is exactly 0 (degenerate near-perfect fit on a
  # simulated dataset), which would otherwise propagate into the cutoff.
  cfi_v   <- simulated_df$cfi[is.finite(simulated_df$cfi)]
  rmsea_v <- simulated_df$rmsea[is.finite(simulated_df$rmsea)]
  srmr_v  <- simulated_df$srmr[is.finite(simulated_df$srmr)]
  c(
    cfi   = as.numeric(stats::quantile(cfi_v,   1 - pct, na.rm = TRUE)),
    rmsea = as.numeric(stats::quantile(rmsea_v, pct,     na.rm = TRUE)),
    srmr  = as.numeric(stats::quantile(srmr_v,  pct,     na.rm = TRUE))
  )
}

#' Compute two-sided per-item expected loading intervals
#'
#' Central interval covering `percentile`% of each item's simulated
#' standardized loadings (tails of `(100 - percentile) / 2`% each).
#'
#' @keywords internal
#' @noRd
compute_cfa_loading_cutoffs <- function(simulated_loadings, item_names,
                                        percentile) {
  tail <- (1 - percentile / 100) / 2
  do.call(rbind, lapply(item_names, function(it) {
    v <- simulated_loadings[[it]]
    v <- v[is.finite(v)]
    data.frame(
      Item = it,
      low  = as.numeric(stats::quantile(v, tail,     na.rm = TRUE)),
      high = as.numeric(stats::quantile(v, 1 - tail, na.rm = TRUE)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
}

#' Tiny ordinal-suffix helper for percentile labels (1st, 2nd, 3rd, 99th)
#'
#' @keywords internal
#' @noRd
ordinal_suffix <- function(n) {
  n_int <- as.integer(round(n))
  if (n_int %% 100 %in% c(11L, 12L, 13L)) return("th")
  switch(as.character(n_int %% 10),
         "1" = "st", "2" = "nd", "3" = "rd", "th")
}

#' Apply directional fit-index flagging given observed values + cutoffs
#'
#' @keywords internal
#' @noRd
compute_cfa_flagged <- function(observed, cutoffs) {
  c(
    cfi   = !is.na(observed[["cfi"]])   && observed[["cfi"]]   < cutoffs[["cfi"]],
    rmsea = !is.na(observed[["rmsea"]]) && observed[["rmsea"]] > cutoffs[["rmsea"]],
    srmr  = !is.na(observed[["srmr"]])  && observed[["srmr"]]  > cutoffs[["srmr"]]
  )
}

#' Flag per-item loadings against the two-sided expected interval
#'
#' Returns a character vector aligned to `loading_cutoffs$Item`:
#' `"below"` (observed loading under the expected interval),
#' `"above"` (over it), or `""` (within).
#'
#' @keywords internal
#' @noRd
compute_cfa_loading_flagged <- function(observed_loadings, loading_cutoffs) {
  obs <- as.numeric(observed_loadings[loading_cutoffs$Item])
  ifelse(is.na(obs), "",
         ifelse(obs < loading_cutoffs$low, "below",
                ifelse(obs > loading_cutoffs$high, "above", "")))
}

# ===========================================================================
# Observed-vs-expected tables
# ===========================================================================

#' Observed one-factor CFA fit and loadings vs a simulated reference
#'
#' Fits the observed one-factor categorical CFA to `data` and compares its
#' fit indices and per-item standardized loadings against the simulated null
#' distribution produced by \code{\link{RMdimCFACutoff}}. Returns a list of
#' two tables: model-fit indices and per-item loadings, each with the observed
#' value, the expected reference from the simulation, and a flag.
#'
#' @param data A data.frame or matrix of item responses (non-negative
#'   integers, 0-based), the same items used for the cutoff simulation.
#' @param cutoff The list returned by \code{\link{RMdimCFACutoff}}. Required:
#'   observed CFA fit indices are not interpretable without the simulated
#'   reference, so the function errors if it is missing.
#' @param output Character. `"kable"` (default) returns each table as a
#'   `knitr::kable()`; `"dataframe"` returns plain data.frames.
#'
#' @return A named list with two elements, `fit` and `loadings`:
#' \describe{
#'   \item{`fit`}{CFI / RMSEA / SRMR with columns `Index`, `Observed`,
#'     `Cutoff`, `Direction`, `Flagged` (one-sided, in the unfavourable
#'     direction).}
#'   \item{`loadings`}{One row per item with columns `Item`, `Observed`,
#'     `Expected_low`, `Expected_high`, `Flagged` (`"below"` / `"above"` /
#'     `""`).}
#' }
#' Each element is a `knitr_kable` (when `output = "kable"`) or a data.frame
#' (when `output = "dataframe"`).
#'
#' @seealso \code{\link{RMdimCFACutoff}}, \code{\link{RMdimCFAPlot}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   data("raschdat1", package = "eRm")
#'   sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 50,
#'                         parallel = FALSE, seed = 1)
#'   tabs <- RMdimCFA(raschdat1[, 1:8], cutoff = sim)
#'   tabs$fit
#'   tabs$loadings
#' }
#' }
#'
#' @export
RMdimCFA <- function(data, cutoff, output = c("kable", "dataframe")) {

  output <- match.arg(output)

  if (missing(cutoff) || is.null(cutoff)) {
    stop("`cutoff` is required: pass the object returned by RMdimCFACutoff(). ",
         "Observed CFA fit indices are not interpretable without the ",
         "simulated reference distribution.", call. = FALSE)
  }
  req <- c("simulated", "simulated_loadings", "cutoffs", "loading_cutoffs",
           "percentile", "item_names", "estimator")
  if (!is.list(cutoff) || !all(req %in% names(cutoff))) {
    stop("`cutoff` must be the list returned by RMdimCFACutoff().",
         call. = FALSE)
  }

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for RMdimCFA(). ",
         "Install with: install.packages(\"lavaan\")", call. = FALSE)
  }
  if (output == "kable" && !requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required for output = \"kable\".", call. = FALSE)
  }

  validate_response_data(data)

  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  data <- stats::na.omit(as.data.frame(data))
  if (nrow(data) == 0L) stop("No complete cases in `data`.", call. = FALSE)

  item_names <- cutoff$item_names
  if (!setequal(names(data), item_names)) {
    stop("Item names in `data` do not match those used for the cutoff.\n",
         "  data   : ", paste(names(data),  collapse = ", "), "\n",
         "  cutoff : ", paste(item_names,   collapse = ", "), call. = FALSE)
  }
  data <- data[, item_names, drop = FALSE]

  # Observed CFA fit + standardized loadings (one lavaan fit)
  obs <- run_single_cfa_sim(
    seed      = NA,
    data_list = list(item_names = item_names, estimator = cutoff$estimator),
    obs_data  = data
  )
  if (!is.list(obs)) {
    stop("Observed CFA fit failed: ", obs, call. = FALSE)
  }
  observed_fit  <- stats::setNames(obs$fit, c("cfi", "rmsea", "srmr"))
  observed_load <- obs$loadings

  flagged_fit     <- compute_cfa_flagged(observed_fit, cutoff$cutoffs)
  loading_flagged <- compute_cfa_loading_flagged(observed_load,
                                                 cutoff$loading_cutoffs)

  fit_df  <- build_cfa_fit_df(observed_fit, cutoff$cutoffs, flagged_fit,
                              cutoff$percentile)
  load_df <- build_cfa_loadings_df(observed_load, cutoff$loading_cutoffs,
                                   loading_flagged)

  if (output == "dataframe") {
    return(list(fit = fit_df, loadings = load_df))
  }

  list(
    fit      = render_cfa_fit_kable(fit_df, cutoff, nrow(data)),
    loadings = render_cfa_loadings_kable(load_df, cutoff, nrow(data))
  )
}

#' Assemble the fit-index data.frame (Index/Observed/Cutoff/Direction/Flagged)
#'
#' @keywords internal
#' @noRd
build_cfa_fit_df <- function(observed, cutoffs, flagged, percentile) {
  cfi_pct_lbl <- 100 - percentile
  pct_suffix  <- ordinal_suffix(percentile)
  cfi_suffix  <- ordinal_suffix(cfi_pct_lbl)
  data.frame(
    Index     = c("CFI", "RMSEA", "SRMR"),
    Observed  = round(c(observed[["cfi"]], observed[["rmsea"]],
                        observed[["srmr"]]), 4),
    Cutoff    = round(c(cutoffs[["cfi"]], cutoffs[["rmsea"]],
                        cutoffs[["srmr"]]), 4),
    Direction = c(paste0("< ", cfi_pct_lbl, cfi_suffix, " pct"),
                  paste0("> ", percentile,  pct_suffix, " pct"),
                  paste0("> ", percentile,  pct_suffix, " pct")),
    Flagged   = ifelse(flagged[c("cfi", "rmsea", "srmr")], "TRUE", ""),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Assemble the loadings data.frame (Item/Observed/Expected_low/high/Flagged)
#'
#' @keywords internal
#' @noRd
build_cfa_loadings_df <- function(observed_loadings, loading_cutoffs,
                                  loading_flagged) {
  data.frame(
    Item          = loading_cutoffs$Item,
    Observed      = round(as.numeric(observed_loadings[loading_cutoffs$Item]), 3),
    Expected_low  = round(loading_cutoffs$low,  3),
    Expected_high = round(loading_cutoffs$high, 3),
    Flagged       = loading_flagged,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @keywords internal
#' @noRd
render_cfa_fit_kable <- function(fit_df, cutoff, observed_n) {
  pct_label <- cutoff$percentile
  caption <- paste0(
    if (cutoff$is_polytomous) "Partial Credit Model" else "Rasch Model",
    " posterior-predictive CFA fit-index check. ",
    "Observed CFA fit (one-factor, lavaan ", cutoff$estimator,
    ", ordered = TRUE; n = ", observed_n, ") vs simulated null under ",
    if (cutoff$is_polytomous) "PCM" else "RM",
    " unidimensionality (", cutoff$actual_iterations,
    " iterations at n = ", cutoff$sample_n, "). ",
    "Cutoffs one-sided at the ", pct_label,
    "th percentile; flagged when the observed value lies in the worst ",
    100 - pct_label, "% of the null in the unfavourable direction."
  )
  knitr::kable(fit_df, format = "pipe", caption = caption, row.names = FALSE)
}

#' @keywords internal
#' @noRd
render_cfa_loadings_kable <- function(load_df, cutoff, observed_n) {
  pct_label <- cutoff$percentile
  caption <- paste0(
    "Standardized factor loadings (one-factor, lavaan ", cutoff$estimator,
    "; n = ", observed_n, ") vs the simulated expected range under ",
    if (cutoff$is_polytomous) "PCM" else "RM",
    " unidimensionality (", cutoff$actual_iterations,
    " iterations at n = ", cutoff$sample_n, "). ",
    "Expected range is the two-sided central ", pct_label,
    "% interval of the simulated loadings; Flagged = below / above that range."
  )
  knitr::kable(
    load_df, format = "pipe", row.names = FALSE,
    col.names = c("Item", "Observed", "Expected low", "Expected high",
                  "Flagged"),
    caption = caption
  )
}

# ===========================================================================
# Plot companion
# ===========================================================================

#' Plot observed CFA fit and loadings against the simulated null
#'
#' Returns two figures (in a list) comparing the observed one-factor CFA to
#' the simulated null distribution from \code{\link{RMdimCFACutoff}}:
#' a per-item standardized-loadings plot (observed marker against each item's
#' simulated distribution and expected range, in the style of
#' \code{\link{RMitemInfitCutoffPlot}}), and a faceted plot of the CFI / RMSEA
#' / SRMR distributions with the observed value overlaid.
#'
#' @param cutoff_res The list returned by \code{\link{RMdimCFACutoff}}.
#' @param data The item-response data the CFA was run on (the same items used
#'   for the cutoff). Required: the observed values are computed from it.
#' @param percentile Numeric in (50, 100) or `NULL`. When supplied, the
#'   cutoffs and flags are recomputed at this percentile from the stored
#'   simulated distributions (no re-simulation). When `NULL` (default), the
#'   percentile from the original `RMdimCFACutoff()` call is reused.
#'
#' @return A named list of two `ggplot` objects:
#' \describe{
#'   \item{`loadings`}{Per-item standardized loadings: simulated distribution
#'     (dots), expected interval, and the observed loading as a diamond
#'     (red when flagged).}
#'   \item{`fit`}{Faceted CFI / RMSEA / SRMR simulated distributions with the
#'     observed value and cutoff overlaid.}
#' }
#'
#' @seealso \code{\link{RMdimCFACutoff}}, \code{\link{RMdimCFA}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("lavaan", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   data("raschdat1", package = "eRm")
#'   sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 50,
#'                         parallel = FALSE, seed = 1)
#'   plots <- RMdimCFAPlot(sim, data = raschdat1[, 1:8])
#'   plots$loadings
#'   plots$fit
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
RMdimCFAPlot <- function(cutoff_res, data, percentile = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for RMdimCFAPlot(). ",
         "Install with: install.packages(\"ggplot2\")", call. = FALSE)
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required for RMdimCFAPlot(). ",
         "Install with: install.packages(\"ggdist\")", call. = FALSE)
  }

  req <- c("simulated", "simulated_loadings", "cutoffs", "loading_cutoffs",
           "percentile", "item_names", "estimator")
  if (!is.list(cutoff_res) || !all(req %in% names(cutoff_res))) {
    stop("`cutoff_res` must be the result returned by RMdimCFACutoff().",
         call. = FALSE)
  }
  if (missing(data)) {
    stop("`data` is required: the observed CFA fit and loadings are computed ",
         "from it for the overlay.", call. = FALSE)
  }

  # Resolve / override percentile
  if (is.null(percentile)) {
    percentile <- cutoff_res$percentile
  } else if (!is.numeric(percentile) || length(percentile) != 1L ||
             !is.finite(percentile) || percentile <= 50 || percentile >= 100) {
    stop("`percentile` must be a single numeric in (50, 100), or NULL.",
         call. = FALSE)
  }
  cutoffs         <- compute_cfa_cutoffs(cutoff_res$simulated, percentile)
  loading_cutoffs <- compute_cfa_loading_cutoffs(cutoff_res$simulated_loadings,
                                                 cutoff_res$item_names,
                                                 percentile)

  # Observed fit + loadings from the data
  validate_response_data(data)
  data <- stats::na.omit(as.data.frame(data))
  if (!setequal(names(data), cutoff_res$item_names)) {
    stop("Item names in `data` do not match those used for the cutoff.",
         call. = FALSE)
  }
  data <- data[, cutoff_res$item_names, drop = FALSE]

  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  obs <- run_single_cfa_sim(
    seed      = NA,
    data_list = list(item_names = cutoff_res$item_names,
                     estimator  = cutoff_res$estimator),
    obs_data  = data
  )
  if (!is.list(obs)) stop("Observed CFA fit failed: ", obs, call. = FALSE)
  observed_fit  <- stats::setNames(obs$fit, c("cfi", "rmsea", "srmr"))
  observed_load <- obs$loadings

  list(
    loadings = cfa_loadings_plot(cutoff_res, observed_load, loading_cutoffs,
                                 percentile),
    fit      = cfa_fit_plot(cutoff_res, observed_fit, cutoffs, percentile)
  )
}

#' Per-item standardized-loadings plot (RMitemInfitCutoffPlot style)
#'
#' @keywords internal
#' @noRd
cfa_loadings_plot <- function(cutoff_res, observed_load, loading_cutoffs,
                              percentile) {
  item_names  <- cutoff_res$item_names
  item_levels <- rev(item_names)
  sim_load    <- cutoff_res$simulated_loadings

  # Long-format simulated loadings
  sim_long <- do.call(rbind, lapply(item_names, function(it) {
    data.frame(Item = it, Value = sim_load[[it]],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  sim_long$Item <- factor(sim_long$Item, levels = item_levels)

  # Per-item summary intervals (thin = 0.1/99.9, thick = 16.7/83.3, median)
  lo_hi <- do.call(rbind, lapply(item_names, function(it) {
    v <- sim_load[[it]]; v <- v[is.finite(v)]
    data.frame(
      Item   = it,
      lo     = stats::quantile(v, 0.001, na.rm = TRUE),
      hi     = stats::quantile(v, 0.999, na.rm = TRUE),
      lo66   = stats::quantile(v, 0.167, na.rm = TRUE),
      hi66   = stats::quantile(v, 0.833, na.rm = TRUE),
      median = stats::median(v, na.rm = TRUE),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }))
  lo_hi$Item_f <- factor(lo_hi$Item, levels = item_levels)

  flags <- compute_cfa_loading_flagged(observed_load, loading_cutoffs)
  obs_df <- data.frame(
    Item     = loading_cutoffs$Item,
    observed = as.numeric(observed_load[loading_cutoffs$Item]),
    Flagged  = flags != "",
    stringsAsFactors = FALSE
  )
  obs_df$Item  <- factor(obs_df$Item, levels = item_levels)
  obs_df$Color <- ifelse(obs_df$Flagged, "red", "sienna2")

  caption_text <- er2_caption(paste0(
    "Standardized factor loadings from ", cutoff_res$actual_iterations,
    " datasets simulated under ",
    if (cutoff_res$is_polytomous) "PCM" else "RM",
    " unidimensionality at n = ", cutoff_res$sample_n,
    ", refitted with lavaan::cfa(ordered = TRUE, estimator = \"",
    cutoff_res$estimator, "\").\n",
    "Diamonds: observed loading (red = outside the two-sided ", percentile,
    "% expected range). Black dots: simulated median."
  ))

  ggplot2::ggplot(sim_long, ggplot2::aes(x = .data$Value, y = .data$Item)) +
    ggdist::stat_dots(
      ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
      quantiles = cutoff_res$actual_iterations,
      layout = "weave", slab_color = NA, .width = c(0.666, 0.999)
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(x = .data$lo, xend = .data$hi,
                   y = .data$Item_f, yend = .data$Item_f),
      color = "black", linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(x = .data$lo66, xend = .data$hi66,
                   y = .data$Item_f, yend = .data$Item_f),
      color = "black", linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = lo_hi,
      ggplot2::aes(x = .data$median, y = .data$Item_f), size = 3.6
    ) +
    ggplot2::geom_point(
      data = obs_df,
      ggplot2::aes(x = .data$observed, y = .data$Item, colour = .data$Color),
      shape = 18, size = 4,
      position = ggplot2::position_nudge(y = -0.1)
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_color_manual(
      values = scales::brewer_pal()(3)[-1],
      aesthetics = "slab_fill", guide = "none"
    ) +
    ggplot2::labs(
      x = "Standardized factor loading", y = "Item", caption = caption_text
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
    er2_axis_margins() +
    er2_plot_caption()
}

#' Faceted CFI / RMSEA / SRMR distribution plot with observed overlay
#'
#' @keywords internal
#' @noRd
cfa_fit_plot <- function(cutoff_res, observed_fit, cutoffs, percentile) {
  flagged <- compute_cfa_flagged(observed_fit, cutoffs)

  sim_long <- data.frame(
    Index = factor(rep(c("CFI", "RMSEA", "SRMR"),
                       each = nrow(cutoff_res$simulated)),
                   levels = c("CFI", "RMSEA", "SRMR")),
    Value = c(cutoff_res$simulated$cfi,
              cutoff_res$simulated$rmsea,
              cutoff_res$simulated$srmr),
    stringsAsFactors = FALSE
  )

  obs_df <- data.frame(
    Index    = factor(c("CFI", "RMSEA", "SRMR"),
                      levels = c("CFI", "RMSEA", "SRMR")),
    Observed = c(observed_fit[["cfi"]], observed_fit[["rmsea"]],
                 observed_fit[["srmr"]]),
    Cutoff   = c(cutoffs[["cfi"]], cutoffs[["rmsea"]], cutoffs[["srmr"]]),
    Flagged  = c(flagged[["cfi"]], flagged[["rmsea"]], flagged[["srmr"]]),
    stringsAsFactors = FALSE
  )
  obs_df$Color <- ifelse(obs_df$Flagged, "red", "grey30")

  cfi_pct_lbl <- 100 - percentile
  caption <- er2_caption(paste0(
    "Histograms: ", cutoff_res$actual_iterations,
    " datasets simulated under ",
    if (cutoff_res$is_polytomous) "PCM" else "RM",
    " unidimensionality at n = ", cutoff_res$sample_n,
    ",\nrefitted with lavaan::cfa(ordered = TRUE, estimator = \"",
    cutoff_res$estimator, "\").\n",
    "Diamond: observed value (red = flagged at the ", percentile,
    "th percentile in the unfavourable direction).\n",
    "Dashed line: cutoff (CFI: ", cfi_pct_lbl,
    "th pct; RMSEA / SRMR: ", percentile, "th pct)."
  ))

  ggplot2::ggplot(sim_long, ggplot2::aes(x = .data$Value)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey80", colour = "white") +
    ggplot2::geom_vline(
      data = obs_df, ggplot2::aes(xintercept = .data$Cutoff),
      linetype = "dashed", colour = "grey40"
    ) +
    ggplot2::geom_point(
      data = obs_df,
      ggplot2::aes(x = .data$Observed, colour = .data$Color), y = 0,
      size = 4.5, shape = 18
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::facet_wrap(~ Index, scales = "free", nrow = 1) +
    ggplot2::labs(
      title   = "Observed CFA fit vs simulated null distribution",
      x       = "Fit index value",
      y       = "Count of simulated datasets",
      caption = caption
    ) +
    ggplot2::theme_bw(base_size = 13) +
    er2_axis_margins() +
    er2_plot_caption()
}

# ===========================================================================
# Internal: single iteration
# ===========================================================================

#' Run one CFA simulation iteration
#'
#' If `obs_data` is supplied, it is used directly (the observed-fit call);
#' otherwise a new dataset is simulated under the fitted PCM/RM. Returns a
#' list `list(fit = c(cfi, rmsea, srmr), loadings = <named vector>)` on
#' success, or a character message on failure.
#'
#' @keywords internal
#' @noRd
run_single_cfa_sim <- function(seed, data_list, obs_data = NULL) {

  build_sim_df <- function() {
    set.seed(seed)
    thetas_res <- sample(data_list$thetas, size = data_list$sample_n,
                         replace = TRUE)

    if (data_list$type == "dichotomous") {
      sim_mat <- psychotools::rrm(theta = thetas_res,
                                  beta  = data_list$item_params)
      sim_df <- as.data.frame(sim_mat$data)
      colnames(sim_df) <- data_list$item_names

      pos_counts <- colSums(sim_df, na.rm = TRUE)
      neg_counts <- nrow(sim_df) - pos_counts
      if (any(pos_counts < 2L) || any(neg_counts < 2L)) {
        return("validation_failed: an item has < 2 responses in one of the two categories")
      }
      return(sim_df)
    }

    sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
    sim_df  <- as.data.frame(sim_mat)
    colnames(sim_df) <- data_list$item_names

    n_cats <- vapply(data_list$deltaslist,
                     function(d) length(d) + 1L, integer(1L))
    for (j in seq_len(ncol(sim_df))) {
      tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
      if (any(tab == 0L)) {
        return("validation_failed: not all categories represented")
      }
    }
    sim_df
  }

  tryCatch({
    fit_df <- if (is.null(obs_data)) build_sim_df() else obs_data
    if (is.character(fit_df)) return(fit_df)

    fmla <- paste0("F1 =~ ",
                   paste(data_list$item_names, collapse = " + "))

    fit <- suppressWarnings(suppressMessages(
      lavaan::cfa(
        model     = fmla,
        data      = fit_df,
        ordered   = data_list$item_names,
        estimator = data_list$estimator,
        warn      = FALSE,
        verbose   = FALSE
      )
    ))

    if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
      return("convergence_failed: lavaan did not converge")
    }

    suppressWarnings(list(
      fit      = extract_cfa_fit(fit, data_list$estimator),
      loadings = extract_cfa_loadings(fit, data_list$item_names)
    ))
  }, error = function(e) as.character(conditionMessage(e)))
}

#' Pull (CFI, RMSEA, SRMR) from a fitted lavaan object
#'
#' Uses the Satorra-Bentler-SCALED CFI / RMSEA variants when available
#' (`cfi.scaled`, `rmsea.scaled`), falling back to the uncorrected
#' `cfi` / `rmsea` when not. SRMR is unaffected by the correction. The
#' `.scaled` variants are well-defined across nearly all fits (unlike the
#' `.robust` variants, which return `NA` for a non-trivial fraction at small
#' n); for a percentile comparison the binding requirement is that the same
#' metric is computed for observed and simulated iterations.
#'
#' @keywords internal
#' @noRd
extract_cfa_fit <- function(fit, estimator = NULL) {
  fm <- lavaan::fitMeasures(fit)
  pick_scaled <- function(name) {
    key_s <- paste0(name, ".scaled")
    if (key_s %in% names(fm) && is.finite(fm[[key_s]])) {
      return(as.numeric(fm[[key_s]]))
    }
    if (name %in% names(fm) && is.finite(fm[[name]])) {
      return(as.numeric(fm[[name]]))
    }
    NA_real_
  }
  c(pick_scaled("cfi"),
    pick_scaled("rmsea"),
    as.numeric(fm[["srmr"]]))
}

#' Pull standardized factor loadings (one factor) from a fitted lavaan object
#'
#' Returns the `est.std` of the `=~` paths from
#' `lavaan::standardizedSolution()`, aligned to `item_names`.
#'
#' @keywords internal
#' @noRd
extract_cfa_loadings <- function(fit, item_names) {
  ss  <- lavaan::standardizedSolution(fit)
  lam <- ss[ss$op == "=~", , drop = FALSE]
  v   <- stats::setNames(as.numeric(lam$est.std), lam$rhs)
  # Keep names so callers can index by item (observed-side lookups rely on it).
  stats::setNames(as.numeric(v[item_names]), item_names)
}

# ===========================================================================
# Internal: parallel runner
# ===========================================================================

#' @keywords internal
#' @noRd
run_cfa_sim_parallel <- function(iterations, sim_seeds, sim_data_list,
                                 n_cores, verbose = FALSE) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  tasks <- lapply(seq_len(iterations), function(i) {
    mirai::mirai(
      { run_single_cfa_sim(seed, data_list) },
      seed                  = sim_seeds[i],
      data_list             = sim_data_list,
      run_single_cfa_sim    = run_single_cfa_sim,
      extract_cfa_fit       = extract_cfa_fit,
      extract_cfa_loadings  = extract_cfa_loadings,
      sim_partial_score     = sim_partial_score,
      sim_poly_item         = sim_poly_item
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

# ===========================================================================
# Internal: sequential runner
# ===========================================================================

#' @keywords internal
#' @noRd
run_cfa_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                   verbose = FALSE) {
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)

  results <- vector("list", iterations)
  for (i in seq_len(iterations)) {
    results[[i]] <- run_single_cfa_sim(sim_seeds[i], sim_data_list)
    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}
