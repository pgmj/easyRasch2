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
