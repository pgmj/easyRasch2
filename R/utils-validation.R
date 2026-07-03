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
