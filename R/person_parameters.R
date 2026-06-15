#' Person Locations for a Rasch / Partial Credit Model
#'
#' Estimates a person location (theta) and its standard error of
#' measurement (SEM) for every respondent. Estimation is performed
#' directly on each person's observed response pattern, so partial
#' missingness is handled correctly: two respondents with the same sum
#' score on different subsets of items receive different estimates. This
#' avoids the sum-score lookup used by [RMscoreSE()], which assumes
#' complete data.
#'
#' @param data A data.frame or matrix of item responses. Items must be
#'   scored starting at 0 (non-negative integers). Missing values (`NA`)
#'   are allowed.
#' @param method Character. `"WLE"` (default) for Warm's Weighted
#'   Likelihood Estimate (lower bias than MLE; Warm, 1989), or `"EAP"`
#'   for the Expected A Posteriori estimate under a normal prior.
#' @param item_params Optional pre-specified item parameters, e.g. for
#'   anchoring or equating. Either a named list of Andrich-threshold
#'   vectors (one per item) or the long-format data.frame returned by
#'   [RMitemParameters()]. When `NULL` (default) item parameters are
#'   estimated from `data`.
#' @param estimator Character. How item parameters are estimated when
#'   `item_params` is `NULL`: `"CML"` (default, via \pkg{eRm}) or
#'   `"MML"` (via \pkg{mirt}). Ignored when `item_params` is supplied.
#' @param theta_range Numeric length 2. Search range for the WLE root
#'   and bounds for the EAP quadrature grid. Default `c(-10, 10)`. Only
#'   if the WLE root falls outside this range is the boundary returned.
#' @param prior_mean Numeric. Mean of the normal prior used by
#'   `method = "EAP"`. Default `0`, the natural choice when item
#'   parameters are centred at mean difficulty zero. Ignored for WLE.
#' @param prior_sd Numeric or `NULL`. Standard deviation of the normal
#'   prior used by `method = "EAP"`. When `NULL` (default) it is
#'   estimated from the data by marginal maximum likelihood, holding the
#'   item parameters fixed. Supply a number (e.g. `1` for a standard
#'   `N(0, 1)` prior, or a value carried over from a reference sample
#'   for equating) to fix it. Ignored for WLE.
#' @param output Character. `"kable"` (default) for a formatted
#'   [knitr::kable()] table, `"dataframe"` for the underlying
#'   data.frame, or `"ggplot"` for a histogram of the estimated person
#'   locations.
#'
#' @return
#' For `output = "dataframe"`, a data.frame with one row per respondent
#' (in input order) and columns `theta` (the person location), `sem`
#' (standard error of measurement), `sum_score`, `n_answered` (number of
#' non-missing responses), and `extreme` (logical: a minimum or maximum
#' possible score given the items answered). For `output = "kable"`, the
#' same content as a `knitr_kable` object. For `output = "ggplot"`, a
#' histogram of `theta`.
#'
#' @details
#' Item parameters are obtained once (by CML or MML, or taken from
#' `item_params`) and treated as fixed. Each person's theta is then
#' estimated from the items they answered.
#'
#' \strong{WLE} solves Warm's weighted-likelihood score equation
#' `l'(theta) + J(theta) / (2 I(theta)) = 0` by bracketed root finding,
#' where the bias-correction term `J / (2 I)` keeps the equation solvable
#' at the boundaries. The SEM is `1 / sqrt(I(theta))`. Unlike the MLE,
#' Warm's estimator therefore yields *finite* locations for the minimum
#' and maximum possible scores (a sensible step beyond the next-most
#' extreme score rather than +/-Inf; cf. Warm, 1989), matching
#' established implementations such as \pkg{catR} and \pkg{TAM}. Such
#' scores are still marked in the `extreme` column because they are
#' extrapolated and carry large standard errors. Only if the root lies
#' outside `theta_range` is the boundary returned with `NA` SEM.
#'
#' \strong{EAP} integrates the pattern likelihood against a normal prior
#' over a quadrature grid; the estimate is the posterior mean and the
#' SEM is the posterior standard deviation. Unlike WLE, EAP is finite at
#' extreme scores because the prior shrinks them inward, at the cost of
#' depending on the assumed prior. The prior actually used -- including
#' the marginal-ML estimate of `prior_sd` when it is left `NULL` -- is
#' reported in the kable caption and attached to the result as
#' `attr(result, "prior")`.
#'
#' @references
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. *Psychometrika, 54*(3), 427-450.
#' \doi{10.1007/BF02294627}
#'
#' Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of
#' ability in a microcomputer environment. *Applied Psychological
#' Measurement, 6*(4), 431-444. \doi{10.1177/014662168200600405}
#'
#' Magis, D. (2015). A note on weighted likelihood and Jeffreys modal
#' estimation of proficiency levels in polytomous item response models.
#' *Psychometrika, 80*(1), 200-204. \doi{10.1007/s11336-013-9378-5}
#'
#' @seealso [RMitemParameters()], [RMscoreSE()], [RMreliability()]
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- as.data.frame(
#'   matrix(sample(0:2, 200 * 6, replace = TRUE), nrow = 200, ncol = 6)
#' )
#' colnames(dat) <- paste0("Item", 1:6)
#' # Introduce some missingness
#' dat[cbind(sample(200, 30), sample(6, 30, replace = TRUE))] <- NA
#'
#' # Default: WLE person locations
#' RMpersonParameters(dat, output = "dataframe") |> head()
#'
#' # EAP with a data-estimated prior SD
#' eap <- RMpersonParameters(dat, method = "EAP", output = "dataframe")
#' attr(eap, "prior")
#'
#' # EAP with a fixed N(0, 1) prior
#' RMpersonParameters(dat, method = "EAP", prior_sd = 1, output = "dataframe") |>
#'   head()
#' }
RMpersonParameters <- function(data,
                               method      = c("WLE", "EAP"),
                               item_params = NULL,
                               estimator   = c("CML", "MML"),
                               theta_range = c(-10, 10),
                               prior_mean  = 0,
                               prior_sd    = NULL,
                               output      = c("kable", "dataframe", "ggplot")) {

  method    <- match.arg(method)
  estimator <- match.arg(estimator)
  output    <- match.arg(output)

  validate_response_data(data)

  if (!is.numeric(theta_range) || length(theta_range) != 2L ||
      theta_range[1L] >= theta_range[2L]) {
    stop("`theta_range` must be a numeric vector of length 2 with ",
         "theta_range[1] < theta_range[2].", call. = FALSE)
  }

  data <- as.data.frame(data)

  # --- Item parameters (estimated or supplied) --------------------------------
  if (is.null(item_params)) {
    fit <- if (estimator == "CML") {
      .rasch_fit_cml(data, se = FALSE)
    } else {
      .rasch_fit_mml(data, se = FALSE)
    }
    # Standard Rasch identification (mean threshold zero), so locations
    # match RMscoreSE() and RMitemParameters().
    thr_list <- .center_thresholds(fit$thr_list)
  } else {
    thr_list <- .coerce_item_params(item_params)
  }

  data <- .align_items(data, thr_list)
  data_mat <- as.matrix(data)

  # --- Person-level bookkeeping -----------------------------------------------
  n_answered <- rowSums(!is.na(data_mat))
  sum_score  <- rowSums(data_mat, na.rm = TRUE)
  max_score  <- vapply(seq_len(nrow(data_mat)), function(p) {
    ans <- !is.na(data_mat[p, ])
    sum(vapply(thr_list[ans], length, integer(1)))
  }, numeric(1))
  extreme <- n_answered > 0L & (sum_score == 0 | sum_score == max_score)

  # --- Estimation -------------------------------------------------------------
  prior_used <- NULL
  if (method == "WLE") {
    est <- t(vapply(seq_len(nrow(data_mat)), function(p) {
      .theta_wle(data_mat[p, ], thr_list, theta_range)
    }, numeric(2)))
    theta <- est[, 1L]
    sem   <- est[, 2L]
  } else {
    grid      <- seq(theta_range[1L], theta_range[2L], length.out = 121L)
    logp_tabs <- .logp_tables(thr_list, grid)
    loglik    <- .grid_loglik(data_mat, logp_tabs, grid)

    sd_used <- if (is.null(prior_sd)) {
      .estimate_prior_sd(loglik, grid, prior_mean)
    } else {
      if (!is.numeric(prior_sd) || length(prior_sd) != 1L || prior_sd <= 0) {
        stop("`prior_sd` must be a single positive number or NULL.",
             call. = FALSE)
      }
      prior_sd
    }
    eap   <- .theta_eap(loglik, grid, prior_mean, sd_used)
    theta <- eap$theta
    sem   <- eap$sem
    prior_used <- c(mean = prior_mean, sd = sd_used,
                    estimated = as.numeric(is.null(prior_sd)))
  }

  out <- data.frame(
    theta      = theta,
    sem        = sem,
    sum_score  = sum_score,
    n_answered = n_answered,
    extreme    = extreme,
    stringsAsFactors = FALSE
  )
  rownames(out) <- rownames(data)
  if (!is.null(prior_used)) attr(out, "prior") <- prior_used

  # --- Output -----------------------------------------------------------------
  if (output == "dataframe") {
    out$theta <- round(out$theta, 4)
    out$sem   <- round(out$sem, 4)
    return(out)
  }

  if (output == "ggplot") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for output = \"ggplot\". ",
           "Install it with: install.packages(\"ggplot2\")", call. = FALSE)
    }
    plot_df <- out[is.finite(out$theta), , drop = FALSE]
    return(
      ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$theta)) +
        ggplot2::geom_histogram(bins = 30, fill = "grey70", colour = "white") +
        ggplot2::labs(x = "Person location (logits)", y = "Count") +
        ggplot2::theme_bw() +
        er2_axis_margins()
    )
  }

  display <- out
  display$theta <- round(display$theta, 3)
  display$sem   <- round(display$sem, 3)
  caption <- if (method == "WLE") {
    paste0("Person locations via Warm's WLE (", estimator,
           " item parameters). Extreme scores receive finite extrapolated ",
           "estimates (flagged in the extreme column).")
  } else {
    paste0("Person locations via EAP (", estimator,
           " item parameters; normal prior mean = ",
           round(prior_used[["mean"]], 3), ", sd = ",
           round(prior_used[["sd"]], 3),
           if (prior_used[["estimated"]] == 1) " [estimated]" else " [fixed]",
           ").")
  }
  knitr::kable(display, format = "pipe", caption = caption, row.names = TRUE)
}

# ---------------------------------------------------------------------
# Internal: coerce supplied item parameters to a threshold list
# ---------------------------------------------------------------------

#' @keywords internal
#' @noRd
.coerce_item_params <- function(item_params) {
  if (is.list(item_params) && !is.data.frame(item_params)) {
    thr_list <- lapply(item_params, function(x) as.numeric(x))
    if (is.null(names(thr_list))) {
      stop("When `item_params` is a list it must be named by item.",
           call. = FALSE)
    }
    return(thr_list)
  }
  if (is.data.frame(item_params)) {
    if (!all(c("item", "location") %in% names(item_params))) {
      stop("`item_params` data.frame must have `item` and `location` ",
           "columns (as returned by RMitemParameters(format = \"long\")).",
           call. = FALSE)
    }
    sp <- split(item_params$location, factor(item_params$item,
                                             levels = unique(item_params$item)))
    return(lapply(sp, as.numeric))
  }
  stop("`item_params` must be a named list of threshold vectors or a ",
       "data.frame from RMitemParameters().", call. = FALSE)
}

#' Align response-data columns with the item-parameter list
#'
#' @param data Response data.frame.
#' @param thr_list Named (or order-matched) list of threshold vectors.
#' @return `data` with columns ordered/subset to match `thr_list`.
#' @keywords internal
#' @noRd
.align_items <- function(data, thr_list) {
  if (!is.null(names(thr_list)) && all(names(thr_list) %in% colnames(data))) {
    return(data[, names(thr_list), drop = FALSE])
  }
  if (length(thr_list) != ncol(data)) {
    stop("Number of items in `item_params` (", length(thr_list),
         ") does not match the number of columns in `data` (",
         ncol(data), "), and names do not align.", call. = FALSE)
  }
  data
}
