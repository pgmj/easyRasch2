#' Conditional Item Infit MSQ
#'
#' Computes conditional infit mean-square (MSQ) statistics for each item using
#' `iarm::out_infit()`, enriched with item locations relative to the sample
#' mean person location.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed,
#'   but at least one complete case (row with no `NA`) must be present.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying data.frame.
#' @param sort Optional character string. When `sort = "infit"`, rows are
#'   sorted by `Infit_MSQ` in descending order before output.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object (plain text table via
#'   `format = "simple"`) with columns "Item", "Infit MSQ", and "Relative
#'   location", and a caption showing the number of complete cases.
#' * If `output = "dataframe"`: a data.frame with columns `Item`,
#'   `Infit_MSQ`, and `Relative_location`.
#'
#' @details
#' Infit MSQ is a weighted fit statistic that emphasises deviations near the
#' item location. Values close to 1.0 indicate good fit. Values substantially
#' above 1.0 suggest underfit (unexpected responses), while values substantially
#' below 1.0 suggest overfit (overly predictable responses).
#'
#' Conditional infit MSQ statistics are computed via `iarm::out_infit()`, which
#' uses the conditional distribution of the sufficient statistics (Müller, 2020).
#' Only complete cases (rows without any `NA`) are used in the conditional fit
#' calculation.
#'
#' For **dichotomous** data (maximum score = 1), a Rasch model is fitted via
#' `eRm::RM()`. Item locations are the negative beta parameters. Person
#' locations are estimated via `eRm::person.parameter()`.
#'
#' For **polytomous** data (maximum score > 1), a Partial Credit Model is
#' fitted via `eRm::PCM()`. Item average locations are taken from the
#' "Location" column of the threshold parameter table returned by
#' `eRm::thresholds()`; if that column is absent, row means of the threshold
#' columns are used instead. Person locations are estimated via
#' `eRm::person.parameter()`.
#'
#' Relative item location is defined as the item's average location minus the
#' sample mean person location, providing a measure of item targeting.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' @references
#' Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust
#' them? *Journal of Statistical Distributions and Applications*, 7(5).
#' \doi{10.1186/s40488-020-00108-7}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate binary item response data (5 items, 40 persons)
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 40 * 5, replace = TRUE), nrow = 40, ncol = 5)
#' )
#' colnames(sim_data) <- paste0("Item", 1:5)
#'
#' # Default kable output
#' RMiteminfit(sim_data)
#'
#' # Sorted by infit MSQ descending
#' RMiteminfit(sim_data, sort = "infit")
#'
#' # Return as data.frame for further processing
#' df <- RMiteminfit(sim_data, output = "dataframe")
#' }
RMiteminfit <- function(data, output = "kable", sort) {

  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMiteminfit() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe"))

  validate_response_data(data)

  if (nrow(stats::na.omit(data)) == 0L) {
    stop("No complete cases in data. All rows contain at least one NA.",
         call. = FALSE)
  }

  data_mat <- as.matrix(data)

  # --- Fit Rasch model and compute item/person locations ----------------------
  if (max(data_mat, na.rm = TRUE) == 1L) {
    # Dichotomous: Rasch model
    erm_out <- eRm::RM(data)
    item_avg_locations <- coef(erm_out, "beta") * -1
    pp <- eRm::person.parameter(erm_out)
    person_avg_location <- mean(pp$theta.table[["Person Parameter"]], na.rm = TRUE)
  } else {
    # Polytomous: Partial Credit Model
    erm_out <- eRm::PCM(data)
    thresh_obj <- eRm::thresholds(erm_out)
    thresh_table <- thresh_obj$threshtable[[1]]
    if ("Location" %in% colnames(thresh_table)) {
      item_avg_locations <- thresh_table[, "Location"]
    } else {
      item_avg_locations <- rowMeans(thresh_table, na.rm = TRUE)
    }
    pp <- eRm::person.parameter(erm_out)
    person_avg_location <- mean(pp$theta.table[["Person Parameter"]], na.rm = TRUE)
  }

  relative_item_avg_locations <- item_avg_locations - person_avg_location

  # --- rgl.useNULL workaround -------------------------------------------------
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # --- Compute conditional infit MSQ via iarm ---------------------------------
  cfit <- iarm::out_infit(erm_out)

  # --- Count complete cases ---------------------------------------------------
  n_complete <- nrow(stats::na.omit(data))

  # --- Assemble result data.frame ---------------------------------------------
  item_fit_table <- data.frame(
    Item              = names(data),
    Infit_MSQ         = round(cfit$Infit, 3),
    Relative_location = round(relative_item_avg_locations, 2),
    stringsAsFactors  = FALSE,
    row.names         = NULL
  )

  # --- Sort if requested ------------------------------------------------------
  if (!missing(sort) && identical(sort, "infit")) {
    item_fit_table <- item_fit_table[order(item_fit_table$Infit_MSQ, decreasing = TRUE), ]
    rownames(item_fit_table) <- NULL
  }

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(item_fit_table)
  }

  knitr::kable(
    item_fit_table,
    format    = "simple",
    col.names = c("Item", "Infit MSQ", "Relative location"),
    caption   = paste0(
      "MSQ values based on conditional estimation (n = ",
      n_complete,
      " complete cases)."
    )
  )
}
