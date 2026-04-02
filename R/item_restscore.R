#' Item Restscore Analysis
#'
#' Computes observed and model-expected item-restscore correlations using
#' `iarm::item_restscore()`, and enriches the output with the absolute
#' difference between observed and expected values, item average locations, and
#' item locations relative to the sample mean person location.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed,
#'   but at least one complete case (row with no `NA`) must be present.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying data.frame.
#' @param sort Optional character string. When `sort = "diff"`, rows are sorted
#'   by `Absolute_difference` in descending order before output.
#' @param p.adj Character string specifying the p-value adjustment method
#'   passed to `iarm::item_restscore()`. Default `"BH"` (Benjamini-Hochberg).
#'   See [stats::p.adjust()] for available methods.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object (plain text table via
#'   `format = "pipe"`) with columns for item name, observed and expected
#'   restscore correlations, absolute difference, adjusted p-value,
#'   significance level, item location, and item location relative to the
#'   sample mean person location.
#' * If `output = "dataframe"`: a data.frame with columns `Item`, `Observed`,
#'   `Expected`, `Absolute_difference`, `p_adjusted`, `Significance`,
#'   `Location`, and `Relative_location`.
#'
#' @details
#' Item-restscore correlations using Goodman-Kruskal's gamma (Kreiner, 2011) measure
#' the association between a person's score on a single item and their total
#' score on the remaining items (the "restscore"). Under a correctly fitting
#' Rasch model, observed and model-expected correlations should agree closely.
#'
#' For **dichotomous** data (maximum score = 1), a Rasch model is fitted via
#' `eRm::RM()`. Item locations are the negative beta parameters. Person
#' locations are estimated via `eRm::person.parameter()`.
#'
#' For **polytomous** data (maximum score > 1), a Partial Credit Model is
#' fitted via `eRm::PCM()`. Item average locations are the row-means of the
#' threshold parameter table returned by `eRm::thresholds()`. Person locations
#' are estimated via `eRm::person.parameter()`.
#'
#' Relative item location is defined as the item's average location minus the
#' sample mean person location, providing a measure of item targeting.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' @references
#' Kreiner, S. (2011). A Note on Item–Restscore Association in Rasch Models. 
#' *Applied Psychological Measurement, 35*(7), 557–561. 
#' <https://doi.org/10.1177/0146621611410227>
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate binary item response data (8 items, 200 persons)
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
#' )
#' colnames(sim_data) <- paste0("Item", 1:8)
#'
#' # Default kable output
#' RMitemrestscore(sim_data)
#'
#' # Sorted by absolute difference
#' RMitemrestscore(sim_data, sort = "diff")
#'
#' # Return as data.frame for further processing
#' df <- RMitemrestscore(sim_data, output = "dataframe")
#' }
RMitemrestscore <- function(data, output = "kable", sort, p.adj = "BH") {

  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMitemrestscore() but is not installed.\n",
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
  n_items <- ncol(data)

  # --- Fit Rasch model and compute item/person locations ----------------------
  if (max(data_mat, na.rm = TRUE) == 1L) {
    # Dichotomous: Rasch model
    erm_out <- eRm::RM(data)
    item_avg_locations <- stats::coef(erm_out, "beta") * -1
    pp <- eRm::person.parameter(erm_out)
    person_avg_location <- mean(pp$theta.table[["Person Parameter"]], na.rm = TRUE)
  } else {
    # Polytomous: Partial Credit Model
    erm_out <- eRm::PCM(data)
    thresh_obj <- eRm::thresholds(erm_out)
    thresh_table <- thresh_obj$threshtable[[1]]
    item_avg_locations <- rowMeans(thresh_table, na.rm = TRUE)
    pp <- eRm::person.parameter(erm_out)
    person_avg_location <- mean(pp$theta.table[["Person Parameter"]], na.rm = TRUE)
  }

  relative_item_avg_locations <- item_avg_locations - person_avg_location

  # --- Compute item-restscore statistics via iarm ----------------------------
  # Temporarily set rgl.useNULL to avoid rgl device issues during iarm fitting
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  i1 <- iarm::item_restscore(erm_out, p.adj = p.adj)
  i1 <- as.data.frame(i1)

  # i1[[1]] is the results matrix: col 1 = observed, col 2 = expected,
  # col 3 = se, col 4 = z-statistic, col 5 = adjusted p-value, col 6 = significance
  res_mat  <- i1[[1]]
  observed  <- round(as.numeric(res_mat[seq_len(n_items), 1L]), 2)
  expected  <- round(as.numeric(res_mat[seq_len(n_items), 2L]), 2)
  p_adjusted <- round(as.numeric(res_mat[seq_len(n_items), 5L]), 3)
  significance <- as.character(res_mat[seq_len(n_items), 6L])

  # --- Assemble result data.frame --------------------------------------------
  i2 <- data.frame(
    Item                = names(data),
    Observed            = observed,
    Expected            = expected,
    Absolute_difference = round(abs(expected - observed), 3),
    p_adjusted          = p_adjusted,
    Significance        = significance,
    Location            = round(item_avg_locations, 2),
    Relative_location   = round(relative_item_avg_locations, 2),
    stringsAsFactors    = FALSE,
    row.names           = NULL
  )

  # --- Sort if requested -----------------------------------------------------
  if (!missing(sort) && identical(sort, "diff")) {
    i2 <- i2[order(i2$Absolute_difference, decreasing = TRUE), ]
    rownames(i2) <- NULL
  }

  # --- Return ----------------------------------------------------------------
  if (output == "dataframe") {
    return(i2)
  }

  knitr::kable(
    i2,
    format   = "pipe",
    col.names = c(
      "Item",
      "Observed value",
      "Model expected value",
      "Absolute difference",
      paste0("Adjusted p-value (", p.adj, ")"),
      "Statistical significance level",
      "Location",
      "Relative location"
    )
  )
}
