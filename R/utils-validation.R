# Internal helpers for input validation across easyRasch2 functions

#' Validate response data for Rasch analysis
#'
#' Checks that the supplied object is a data.frame or matrix, that all values
#' are non-negative integers (or `NA`), and that the minimum non-missing value
#' is 0 (i.e., items are scored starting at 0).
#'
#' @param data A data.frame or matrix of item responses.
#'
#' @return Invisibly returns `TRUE` if validation passes; otherwise stops with
#'   an informative error message.
#'
#' @noRd
validate_response_data <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data.frame or matrix.", call. = FALSE)
  }

  vals <- as.vector(as.matrix(data))
  vals_nonmissing <- vals[!is.na(vals)]

  if (length(vals_nonmissing) == 0L) {
    stop("`data` contains no non-missing values.", call. = FALSE)
  }

  if (any(vals_nonmissing < 0)) {
    stop(
      "`data` contains negative values. Item responses must be non-negative integers.",
      call. = FALSE
    )
  }

  if (!all(vals_nonmissing == floor(vals_nonmissing))) {
    stop(
      "`data` contains non-integer values. Item responses must be integers.",
      call. = FALSE
    )
  }

  min_val <- min(vals_nonmissing)
  if (min_val != 0) {
    stop(
      paste0(
        "The minimum value in `data` is ",
        min_val,
        ", but Rasch analysis ",
        "requires items scored starting at 0. Please recode your data."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Response categories used only by extreme scorers
#'
#' CML conditions on the total score, so respondents with a zero or perfect
#' score on the items they answered carry no information about the item
#' parameters. A category observed only among them is empty in the data the
#' model is estimated from: its threshold diverges, the covariance matrix of
#' all thresholds degrades, and a parametric bootstrap built on those
#' thresholds never generates the category again.
#'
#' @param data Numeric response matrix or data.frame (items from 0; `NA`
#'   allowed). The top category of each item is its observed maximum.
#' @return data.frame with columns `Item` and `Category`, one row per
#'   affected category (zero rows when there are none).
#' @keywords internal
#' @noRd
.extreme_only_categories <- function(data) {
  data <- as.matrix(data)
  empty <- data.frame(Item = character(0), Category = integer(0))
  if (ncol(data) == 0L || all(is.na(data))) {
    return(empty)
  }
  item_max <- suppressWarnings(apply(data, 2L, max, na.rm = TRUE))
  item_max[!is.finite(item_max)] <- 0
  answered <- !is.na(data)
  raw <- rowSums(data, na.rm = TRUE)
  possible <- as.numeric(answered %*% item_max)
  extreme <- rowSums(answered) > 0L & (raw == 0 | raw == possible)
  items <- colnames(data)
  if (is.null(items)) {
    items <- paste0("V", seq_len(ncol(data)))
  }

  out <- lapply(seq_len(ncol(data)), function(j) {
    m <- item_max[j]
    if (m < 1) {
      return(NULL)
    }
    x <- data[, j]
    all_n <- tabulate(x[!is.na(x)] + 1L, nbins = m + 1L)
    inner_n <- tabulate(x[!is.na(x) & !extreme] + 1L, nbins = m + 1L)
    cats <- which(all_n > 0L & inner_n == 0L) - 1L
    if (length(cats) == 0L) {
      return(NULL)
    }
    data.frame(Item = items[j], Category = as.integer(cats))
  })
  out <- do.call(rbind, out)
  if (is.null(out)) empty else out
}

#' Message naming categories used only by extreme scorers
#'
#' @param ec Result of `.extreme_only_categories()`, at least one row.
#' @return A single string.
#' @keywords internal
#' @noRd
.extreme_only_message <- function(ec) {
  per_item <- vapply(
    unique(ec$Item),
    function(it) {
      cats <- ec$Category[ec$Item == it]
      paste0(
        it,
        if (length(cats) == 1L) " (category " else " (categories ",
        paste(cats, collapse = ", "),
        ")"
      )
    },
    character(1)
  )
  paste0(
    "Response categories used only by respondents with a zero or perfect ",
    "total score: ",
    paste(per_item, collapse = ", "),
    ". Conditional maximum likelihood leaves these respondents out, so the ",
    "thresholds of these categories cannot be estimated. Consider merging ",
    "each with its adjacent category."
  )
}

#' Error message when every simulation iteration failed
#'
#' Names the most common cause, a category used only by extreme scorers,
#' when it is present in the data the simulation was built from.
#'
#' @param data The response matrix the simulation was built from.
#' @param detail Optional extra text appended to the generic message.
#' @return A single string.
#' @keywords internal
#' @noRd
.all_sims_failed_message <- function(data, detail = "Check your data.") {
  ec <- .extreme_only_categories(data)
  if (nrow(ec) > 0L) {
    return(paste0("All simulation iterations failed. ", .extreme_only_message(ec)))
  }
  paste0("All simulation iterations failed. ", detail)
}

#' Build the standard estimation-sample-size caption clause
#'
#' Returns `"n = X respondents"`, appending ` of Y` only when respondents were
#' excluded (`X < Y`) and a ` (qualifiers)` parenthetical only when there is
#' something to qualify (both fully conditional, so complete data with no
#' exclusions reads simply `n = X respondents`). No trailing punctuation, so
#' callers append `.` or continue the sentence.
#'
#' @param n_used Number of respondents contributing to the estimate.
#' @param n_total Total respondents supplied (raw input rows).
#' @param qualifiers Character vector of parenthetical notes (e.g.
#'   `"incomplete responses retained"`, `"complete cases"`, `"extreme scores
#'   excluded"`); empty for none. Joined with `"; "`.
#' @param noun Respondent noun (default `"respondents"`; RMpersonFit uses
#'   `"respondents assessed"`).
#' @return A single caption-clause string.
#' @noRd
.n_caption <- function(
  n_used,
  n_total,
  qualifiers = character(),
  noun = "respondents"
) {
  of_total <- if (n_used < n_total) paste0(" of ", n_total) else ""
  paren <- if (length(qualifiers) > 0L) {
    paste0(" (", paste(qualifiers, collapse = "; "), ")")
  } else {
    ""
  }
  paste0("n = ", n_used, of_total, " ", noun, paren)
}

#' Caption clause for simulated datasets that had to be discarded
#'
#' A parametric-bootstrap iteration is dropped when its simulated dataset
#' cannot be refitted: an item with an unused response category
#' (polytomous), an item with almost no positive responses (dichotomous), or
#' a refit that fails to converge. All of these get more likely with small
#' samples, extreme item locations and rarely used categories.
#'
#' Everything downstream rests on the successful count, so that is what the
#' captions report. This clause explains the gap when the two differ, and
#' names the remedy, which is simply to ask for more iterations. Returns
#' `NULL` when nothing was lost or when the requested count is unknown
#' (cutoff objects made before the count was stored).
#'
#' @param actual Number of successful iterations.
#' @param requested Number of iterations asked for, or `NULL`.
#' @return A single string with a leading space, or `NULL`.
#' @keywords internal
#' @noRd
.attrition_clause <- function(actual, requested) {
  if (is.null(requested) || is.null(actual)) {
    return(NULL)
  }
  if (!is.finite(requested) || !is.finite(actual) || actual >= requested) {
    return(NULL)
  }
  paste0(
    " ",
    requested - actual,
    " of the ",
    requested,
    " simulated datasets could not be refitted, usually because an item ",
    "ended up with an unused response category or almost no variation, so ",
    "the results above rest on ",
    actual,
    ". Raise `iterations` to recover the intended number."
  )
}

#' Round selected columns of a data.frame for kable display
#'
#' The `output = "dataframe"` contract is *unrounded* values (rounding is a
#' presentation concern); the kable path rounds a display copy just before
#' rendering via this helper. Columns absent from `df` are silently skipped,
#' so one digits spec can serve several table variants (e.g. with and
#' without p-value columns).
#'
#' @param df A data.frame.
#' @param digits Named numeric vector: column name -> decimal places.
#' @return `df` with the named columns rounded.
#' @noRd
.round_display <- function(df, digits) {
  for (nm in intersect(names(digits), names(df))) {
    df[[nm]] <- round(df[[nm]], digits[[nm]])
  }
  df
}

#' Coerce a completed multiply-imputed dataset back to numeric responses
#'
#' `mice::complete()` returns the item columns in their original type. When the
#' items were coded as ordered factors for imputation (required by mice's
#' `method = "polr"`), the completed data are ordered factors, which the
#' downstream CML routines cannot consume. This converts any factor column back
#' to its numeric level values (via `as.character()`, so the 0-based coding is
#' preserved -- not the 1-based integer codes), leaving numeric columns
#' untouched.
#'
#' @param data A completed data.frame from `mice::complete()`.
#' @return `data` with factor columns coerced to numeric.
#' @noRd
.mi_completed_to_numeric <- function(data) {
  factor_cols <- vapply(data, is.factor, logical(1L))
  if (any(factor_cols)) {
    data[factor_cols] <- lapply(
      data[factor_cols],
      function(x) as.numeric(as.character(x))
    )
  }
  data
}

#' Drop respondents with no responses (all-NA rows)
#'
#' A row that is `NA` for every item contributes nothing to any estimate and
#' breaks CML fitting (`psychotools::pcmodel()` errors on all-NA rows). Removes
#' such rows and, matching the package's existing NA-drop messaging (e.g.
#' `RMdifTree`, `RMdifLR`), emits a one-time message when any are dropped.
#'
#' @param data A data.frame or matrix of item responses.
#' @return `data` with all-NA rows removed.
#' @noRd
.drop_empty_respondents <- function(data) {
  empty <- rowSums(!is.na(as.matrix(data))) == 0L
  if (any(empty)) {
    message(sum(empty), " respondent(s) with no responses dropped.")
    data <- data[!empty, , drop = FALSE]
  }
  data
}

#' Validate a `filename` supplied for `output = "file"`
#'
#' Called early (before any model fitting) so a bad path fails fast.
#'
#' @param filename Candidate file path.
#' @return Invisibly `TRUE`; otherwise stops with an informative message.
#' @noRd
.validate_filename <- function(filename) {
  if (
    is.null(filename) ||
      !is.character(filename) ||
      length(filename) != 1L ||
      is.na(filename) ||
      !nzchar(filename)
  ) {
    stop(
      "`filename` must be a single non-empty string when ",
      "`output = \"file\"`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Write a result data.frame to CSV for `output = "file"`
#'
#' Shared by the parameter-export functions ([RMitemParameters()],
#' [RMpersonParameters()]). Writes `df` via [utils::write.csv()], emits a
#' one-line confirmation message, and returns `df` invisibly so the value is
#' still available for assignment.
#'
#' @param df Result data.frame to write.
#' @param filename Single non-empty file path (already validated).
#' @param row.names Passed to [utils::write.csv()]: `TRUE` keeps row names
#'   (e.g. respondent identifiers), `FALSE` drops them.
#' @return Invisibly, `df`.
#' @noRd
.write_output_csv <- function(df, filename, row.names = FALSE) {
  .validate_filename(filename)
  utils::write.csv(df, file = filename, row.names = row.names)
  message("Wrote ", nrow(df), " row(s) to '", filename, "'.")
  invisible(df)
}

#' Warn when a cutoff object was simulated for a different sample
#'
#' A parametric-bootstrap null depends on the sample size it was simulated
#' at, so a cutoff object built on one dataset and applied to another gives
#' intervals and p-values for the wrong n, with no error. The consumer only
#' checks that the item names match, which misses a changed set of
#' respondents. This compares the respondents the consumer analyses with the
#' `sample_n` the cutoff object recorded, using the same missing-data policy
#' on both sides (complete cases for most functions, non-empty respondents
#' for Q3), and for DIF also the group sizes, since the DIF null draws group
#' membership with the observed proportions.
#'
#' Cutoff objects that do not record `sample_n` (made by older versions, or
#' the bare interval data.frame) are skipped silently.
#'
#' @param cutoff_n `sample_n` stored in the cutoff object, or `NULL`.
#' @param n_used Number of respondents the consumer analyses, counted under
#'   the same policy as the cutoff function.
#' @param fn Name of the cutoff function to re-run, e.g.
#'   `"RMitemInfitCutoff()"`.
#' @param policy Short description of what is counted, used in the message.
#' @param groups_cutoff,groups_used Optional integer vectors of group sizes
#'   (DIF), in the same level order.
#' @return Invisibly `TRUE` when the samples agree or cannot be compared,
#'   `FALSE` after warning.
#' @keywords internal
#' @noRd
.check_cutoff_sample <- function(
  cutoff_n,
  n_used,
  fn,
  policy = "complete cases",
  groups_cutoff = NULL,
  groups_used = NULL
) {
  if (is.null(cutoff_n) || length(cutoff_n) != 1L || !is.finite(cutoff_n)) {
    return(invisible(TRUE))
  }
  if (cutoff_n != n_used) {
    warning(
      "The cutoff was simulated for n = ", cutoff_n, " but `data` has n = ",
      n_used, " (", policy, "). The null distribution depends on sample ",
      "size, so the intervals and p-values do not apply to this data. ",
      "Re-run ", fn, " on the data being tested.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  if (
    !is.null(groups_cutoff) &&
      !is.null(groups_used) &&
      !identical(as.integer(groups_cutoff), as.integer(groups_used))
  ) {
    warning(
      "The cutoff was simulated with group sizes ",
      paste(groups_cutoff, collapse = "/"), " but `dif_var` has ",
      paste(groups_used, collapse = "/"), ". The DIF null keeps each ",
      "respondent's observed group, so the intervals and p-values do not ",
      "apply to this grouping. Re-run ", fn,
      " with the `dif_var` being tested.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  invisible(TRUE)
}
