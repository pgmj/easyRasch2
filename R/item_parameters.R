#' Item Parameters for a Rasch / Partial Credit Model
#'
#' Estimates item difficulty (dichotomous) or item-category threshold
#' (polytomous) parameters and returns them in long or wide format, with
#' optional standard errors and Wald confidence intervals. Item
#' parameters are estimated by conditional maximum likelihood (CML, via
#' \pkg{psychotools}) by default, with marginal maximum likelihood (MML, via
#' \pkg{mirt}) available for sparse data where CML can be unstable.
#'
#' @param data A data.frame or matrix of item responses. Items must be
#'   scored starting at 0 (non-negative integers). Missing values (`NA`)
#'   are allowed; both estimators handle them.
#' @param estimator Character. `"CML"` (default) estimates item
#'   parameters with [psychotools::pcmodel()] (a dichotomous item is a
#'   2-category PCM); `"MML"` uses [mirt::mirt()] with
#'   `itemtype = "Rasch"`. CML is preferred for
#'   Rasch measurement; MML can be more robust when data are sparse.
#'   When `estimator = "CML"` and sparse response categories are
#'   detected, a warning suggests switching to MML.
#' @param format Character. `"long"` (default) returns one row per item
#'   (dichotomous) or per item-threshold (polytomous); `"wide"` returns
#'   one row per item with threshold columns `t1`, `t2`, ... plus a mean
#'   `location` column.
#' @param se Logical. If `TRUE` (default), standard-error and
#'   confidence-interval columns are added.
#' @param ci_level Numeric in (0, 1). Confidence level for the Wald
#'   interval (`estimate +/- z * SE`). Default `0.95`.
#' @param center Logical. If `TRUE` (default), all thresholds are shifted
#'   so their grand mean is zero, the usual Rasch identification. CML
#'   estimates are already centred; the shift mainly affects the MML
#'   path, keeping the two estimators on a common scale.
#' @param output Character. `"kable"` (default) for a formatted
#'   [knitr::kable()] table, `"dataframe"` for the underlying
#'   data.frame, or `"file"` to write that data.frame to a CSV at
#'   `filename` (the data.frame is also returned invisibly).
#' @param filename Character. Path to the CSV file to write when
#'   `output = "file"`. Required in that case; ignored otherwise.
#'
#' @return
#' For `output = "dataframe"`, a data.frame. In **long** format the
#' columns are `item`, `threshold` (integer; `1` for dichotomous items),
#' `location`, and -- when `se = TRUE` -- `se`, `ci_lower`, `ci_upper`.
#' In **wide** format the columns are `item`, the threshold locations
#' (`t1`, `t2`, ... or `location` for dichotomous items), and a mean
#' `location`; when `se = TRUE`, matching `se_t1`, `se_t2`, ... columns
#' are appended. For `output = "kable"`, the same content as a
#' `knitr_kable` object.
#'
#' @details
#' Items are detected as dichotomous (maximum score 1) or polytomous
#' (maximum score > 1), and the Rasch or Partial Credit model is chosen
#' accordingly. Thresholds are reported as Andrich thresholds (the
#' person locations at which adjacent response categories are equally
#' probable) on the logit difficulty scale, matching [RMtargeting()].
#'
#' \strong{Standard errors.} For the CML path, threshold SEs are the
#' square roots of the diagonal of the threshold-parameter covariance from
#' `psychotools::threshpar(vcov = TRUE)`. For the MML path, SEs come from
#' the \pkg{mirt} parameter covariance, propagated by the delta method
#' through the linear threshold map.
#' Confidence intervals are Wald intervals and are symmetric on the
#' logit scale.
#'
#' @references
#' Andrich, D. (1978). A rating formulation for ordered response
#' categories. *Psychometrika, 43*(4), 561-573.
#' \doi{10.1007/BF02293814}
#'
#' Mair, P., & Hatzinger, R. (2007). Extended Rasch modeling: The eRm
#' package for the application of IRT models in R. *Journal of
#' Statistical Software, 20*(9), 1-20. \doi{10.18637/jss.v020.i09}
#'
#' @seealso [RMpersonParameters()], [RMscoreSE()], [RMtargeting()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' poly <- as.data.frame(
#'   matrix(sample(0:2, 250 * 5, replace = TRUE), nrow = 250, ncol = 5)
#' )
#' colnames(poly) <- paste0("Item", 1:5)
#'
#' # Default: long-format kable with SE and 95% CI
#' RMitemParameters(poly)
#'
#' # Wide format, point estimates only
#' RMitemParameters(poly, format = "wide", se = FALSE, output = "dataframe")
#'
#' # Dichotomous data
#' dich <- as.data.frame(
#'   matrix(sample(0:1, 250 * 6, replace = TRUE), nrow = 250, ncol = 6)
#' )
#' colnames(dich) <- paste0("Item", 1:6)
#' RMitemParameters(dich, output = "dataframe")
#'
#' # Write the parameter table to a CSV (also returned invisibly)
#' RMitemParameters(poly, output = "file",
#'                  filename = tempfile(fileext = ".csv"))
#' }
RMitemParameters <- function(
  data,
  estimator = c("CML", "MML"),
  format = c("long", "wide"),
  se = TRUE,
  ci_level = 0.95,
  center = TRUE,
  output = c("kable", "dataframe", "file"),
  filename = NULL
) {
  estimator <- match.arg(estimator)
  format <- match.arg(format)
  output <- match.arg(output)

  if (output == "file") {
    .validate_filename(filename)
  }

  validate_response_data(data)

  if (
    !is.numeric(ci_level) ||
      length(ci_level) != 1L ||
      ci_level <= 0 ||
      ci_level >= 1
  ) {
    stop("`ci_level` must be a single number in (0, 1).", call. = FALSE)
  }

  fit <- if (estimator == "CML") {
    .rasch_fit_cml(data, se = se)
  } else {
    .rasch_fit_mml(data, se = se)
  }

  if (center) {
    grand <- mean(unlist(fit$thr_list))
    fit$thr_list <- lapply(fit$thr_list, function(x) x - grand)
  }

  z <- stats::qnorm(1 - (1 - ci_level) / 2)

  # --- Long table -------------------------------------------------------------
  long_rows <- list()
  for (i in seq_along(fit$items)) {
    thr <- fit$thr_list[[i]]
    se_thr <- if (se) fit$se_list[[i]] else rep(NA_real_, length(thr))
    for (k in seq_along(thr)) {
      row <- data.frame(
        item = fit$items[i],
        threshold = as.integer(k),
        location = thr[k],
        stringsAsFactors = FALSE
      )
      if (se) {
        row$se <- se_thr[k]
        row$ci_lower <- thr[k] - z * se_thr[k]
        row$ci_upper <- thr[k] + z * se_thr[k]
      }
      long_rows[[length(long_rows) + 1L]] <- row
    }
  }
  long_tbl <- do.call(rbind, long_rows)
  rownames(long_tbl) <- NULL

  result <- if (format == "long") {
    long_tbl
  } else {
    .itempar_widen(fit, se)
  }

  # --- Output -----------------------------------------------------------------
  if (output %in% c("dataframe", "file")) {
    num <- vapply(result, is.numeric, logical(1))
    result[num] <- lapply(result[num], round, 4)
    if (output == "file") {
      return(.write_output_csv(result, filename, row.names = FALSE))
    }
    return(result)
  }

  caption <- paste0(
    "Item ",
    if (fit$is_polytomous) "thresholds" else "difficulties",
    " via ",
    estimator,
    if (fit$is_polytomous) {
      " (Andrich thresholds, logit scale)."
    } else {
      " (logit scale)."
    },
    if (se) {
      # Wide format shows only SE columns (se_t*); the Wald CI columns are
      # long-format only.
      if (format == "long") {
        paste0(" SE and ", round(ci_level * 100), "% Wald CI shown.")
      } else {
        " SE shown."
      }
    } else {
      ""
    }
  )
  display <- result
  num <- vapply(display, is.numeric, logical(1))
  display[num] <- lapply(display[num], round, 3)
  knitr::kable(display, format = "pipe", caption = caption, row.names = FALSE)
}

# ---------------------------------------------------------------------
# Internal: reshape long item parameters to wide
# ---------------------------------------------------------------------

#' @keywords internal
#' @noRd
.itempar_widen <- function(fit, se) {
  max_thr <- max(vapply(fit$thr_list, length, integer(1)))
  rows <- list()
  for (i in seq_along(fit$items)) {
    thr <- fit$thr_list[[i]]
    se_thr <- if (se) fit$se_list[[i]] else NULL
    row <- list(item = fit$items[i])

    if (max_thr == 1L) {
      row[["location"]] <- thr[1L]
      if (se) row[["se"]] <- se_thr[1L]
    } else {
      for (k in seq_len(max_thr)) {
        row[[paste0("t", k)]] <- if (k <= length(thr)) thr[k] else NA_real_
      }
      row[["location"]] <- mean(thr)
      if (se) {
        for (k in seq_len(max_thr)) {
          row[[paste0("se_t", k)]] <-
            if (k <= length(thr)) se_thr[k] else NA_real_
        }
      }
    }
    rows[[length(rows) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[,, drop = FALSE]
}

# ---------------------------------------------------------------------
# Internal: CML fit + threshold extraction (eRm)
# ---------------------------------------------------------------------

#' Fit a Rasch/PCM by CML and extract thresholds (+ delta-method SEs)
#'
#' @param data Validated response data.
#' @param se Logical; compute threshold standard errors.
#' @return A normalised list: `model`, `items`, `is_polytomous`,
#'   `thr_list` (Andrich thresholds per item) and, if `se`, `se_list`.
#' @keywords internal
#' @noRd
.rasch_fit_cml <- function(data, se = TRUE) {
  data <- as.data.frame(data)
  is_poly <- max(as.matrix(data), na.rm = TRUE) > 1L

  .sparsity_warning(data, is_poly)

  # CML via psychotools (a dichotomous item is a 2-category PCM). Andrich
  # thresholds and their covariance come from threshpar(); callers re-centre
  # the thresholds, so the parameterisation here is immaterial. (Previously
  # eRm::PCM/RM + eRm threshold SEs; threshold locations match to ~5e-5, the
  # reported SEs differ slightly, psychotools vs eRm vcov.)
  fit <- tryCatch(
    psychotools::pcmodel(data),
    error = function(e) {
      stop(
        "CML estimation failed (",
        conditionMessage(e),
        "). ",
        "This often indicates sparse data; try estimator = \"MML\".",
        call. = FALSE
      )
    }
  )

  items <- colnames(data)
  tp <- psychotools::threshpar(fit, vcov = se)
  thr_list <- stats::setNames(lapply(tp, as.numeric), items)

  se_list <- NULL
  if (se) {
    se_all <- sqrt(diag(attr(tp, "vcov")))
    se_list <- split(se_all, rep.int(seq_along(thr_list), lengths(thr_list)))
    se_list <- stats::setNames(lapply(se_list, as.numeric), items)
  }

  list(
    model = if (is_poly) "PCM" else "RM",
    items = items,
    is_polytomous = is_poly,
    thr_list = thr_list,
    se_list = se_list
  )
}

# ---------------------------------------------------------------------
# Internal: MML fit + threshold extraction (mirt)
# ---------------------------------------------------------------------

#' Fit a Rasch/PCM by MML and extract thresholds (+ delta-method SEs)
#'
#' @inheritParams .rasch_fit_cml
#' @return Same normalised structure as `.rasch_fit_cml()`.
#' @keywords internal
#' @noRd
.rasch_fit_mml <- function(data, se = TRUE) {
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop(
      "Package 'mirt' is required for estimator = \"MML\". ",
      "Install it with: install.packages(\"mirt\")",
      call. = FALSE
    )
  }
  data <- as.data.frame(data)
  is_poly <- max(as.matrix(data), na.rm = TRUE) > 1L
  items <- colnames(data)

  fit <- suppressMessages(
    mirt::mirt(
      data,
      model = 1,
      itemtype = "Rasch",
      SE = se,
      verbose = FALSE,
      accelerate = "squarem"
    )
  )

  cf <- mirt::coef(fit, simplify = FALSE)
  V <- if (se) {
    tryCatch(mirt::extract.mirt(fit, "vcov"), error = function(e) NULL)
  } else {
    NULL
  }

  # mirt parameterises P(X = k) proportional to exp(k * theta + d_k) with
  # d_0 = 0, so the Andrich thresholds are tau_k = d_{k-1} - d_k. For a
  # dichotomous item the single intercept "d" gives difficulty -d.
  thr_list <- vector("list", length(items))
  names(thr_list) <- items
  steps <- integer(length(items))
  for (i in seq_along(items)) {
    pars <- cf[[items[i]]]
    d_cols <- grep("^d([0-9]+)?$", colnames(pars), value = TRUE)
    d <- as.numeric(pars[1L, d_cols]) # poly: d0, d1, ...; dich: d
    tau <- if (length(d_cols) == 1L) -d else -diff(d)
    thr_list[[i]] <- tau
    steps[i] <- length(tau)
  }

  se_list <- NULL
  if (se) {
    se_list <- .mml_threshold_se(V, steps)
  }

  list(
    model = if (is_poly) "PCM" else "RM",
    items = items,
    is_polytomous = is_poly,
    thr_list = thr_list,
    se_list = se_list
  )
}

#' Delta-method threshold SEs from a mirt parameter covariance
#'
#' The free intercept parameters appear in the covariance in item order
#' (`d1.<idx>, d2.<idx>, ...`; d_0 is fixed at 0), so they are sliced
#' into per-item blocks of `steps[i]` parameters and propagated through
#' the linear threshold map `tau_k = d_{k-1} - d_k`.
#'
#' @param V Parameter covariance from `vcov(fit)`, or `NULL`.
#' @param steps Integer vector of threshold counts per item.
#' @return List of threshold-SE vectors, one per item (all `NA` if the
#'   covariance is unavailable or does not line up).
#' @keywords internal
#' @noRd
.mml_threshold_se <- function(V, steps) {
  na_out <- lapply(steps, function(m) rep(NA_real_, m))
  if (is.null(V)) {
    return(na_out)
  }

  d_rows <- grep("^d[0-9]*\\.", rownames(V))
  if (length(d_rows) != sum(steps)) {
    return(na_out)
  }

  out <- vector("list", length(steps))
  off <- 0L
  for (i in seq_along(steps)) {
    m <- steps[i]
    idx <- d_rows[(off + 1L):(off + m)]
    v_d <- V[idx, idx, drop = FALSE]

    # tau_k = d_{k-1} - d_k (d_0 fixed): linear map J on free d_1..d_m.
    J <- matrix(0, m, m)
    for (k in seq_len(m)) {
      if (k == 1L) {
        J[1L, 1L] <- -1
      } else {
        J[k, k - 1L] <- 1
        J[k, k] <- -1
      }
    }
    out[[i]] <- sqrt(pmax(diag(J %*% v_d %*% t(J)), 0))
    off <- off + m
  }
  out
}

# ---------------------------------------------------------------------
# Internal: sparsity warning for the CML path
# ---------------------------------------------------------------------

#' Warn about response-category sparsity that can destabilise CML
#'
#' @param data Response data.frame.
#' @param is_poly Logical; polytomous data.
#' @return Invisibly `NULL`; emits a warning when sparse categories or
#'   zero-variance items are found.
#' @keywords internal
#' @noRd
.sparsity_warning <- function(data, is_poly) {
  flagged <- character(0)
  for (nm in colnames(data)) {
    x <- data[[nm]]
    x <- x[!is.na(x)]
    if (length(unique(x)) <= 1L) {
      flagged <- c(flagged, nm)
      next
    }
    if (is_poly) {
      counts <- tabulate(x + 1L, nbins = max(x) + 1L)
      if (any(counts == 0L)) flagged <- c(flagged, nm)
    }
  }
  if (length(flagged) > 0L) {
    warning(
      "Sparse or zero-variance response categories detected in ",
      "item(s): ",
      paste(unique(flagged), collapse = ", "),
      ". ",
      "CML estimates may be unstable; consider estimator = \"MML\".",
      call. = FALSE
    )
  }
  invisible(NULL)
}
