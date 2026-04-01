#' Conditional Item Infit MSQ
#'
#' Computes conditional infit mean-square (MSQ) statistics for each item using
#' `iarm::out_infit()`, enriched with item locations relative to the sample
#' mean person location.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed,
#'   but at least one complete case (row with no `NA`) must be present.
#' @param cutoff Optional. Default `NULL` (no cutoff applied, behaviour is
#'   identical to the current version). Can be:
#'   * The return value of [RMinfitcutoff()] (a list with `$item_cutoffs`):
#'     the data.frame is extracted automatically and the number of simulation
#'     iterations is included in the kable caption.
#'   * The `$item_cutoffs` data.frame from [RMinfitcutoff()] directly: must
#'     have columns `Item`, `infit_low`, and `infit_high`.
#'   When provided, adds columns `Infit_low`, `Infit_high`, and `Flagged`
#'   (logical; `TRUE` when `Infit_MSQ` falls outside the credible range) to
#'   the result.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying data.frame.
#' @param sort Optional character string. When `sort = "infit"`, rows are
#'   sorted by `Infit_MSQ` in descending order before output.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object (plain text table via
#'   `format = "simple"`) with columns "Item", "Infit MSQ", and "Relative
#'   location", and a caption showing the number of complete cases. When
#'   `cutoff` is provided, columns "Infit low", "Infit high", and "Flagged"
#'   are also included, and the caption notes the simulation-based cutoffs.
#' * If `output = "dataframe"`: a data.frame with columns `Item`,
#'   `Infit_MSQ`, and `Relative_location`. When `cutoff` is provided, columns
#'   `Infit_low`, `Infit_high`, and `Flagged` are also included (inserted
#'   after `Infit_MSQ`, before `Relative_location`).
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
#' @seealso [RMinfitcutoff()]
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
#'
#' # With simulation-based cutoffs from RMinfitcutoff()
#' cutoff_res <- RMinfitcutoff(sim_data, iterations = 100, parallel = FALSE,
#'                             seed = 42)
#' RMiteminfit(sim_data, cutoff = cutoff_res)
#' RMiteminfit(sim_data, cutoff = cutoff_res, output = "dataframe")
#' }
RMiteminfit <- function(data, cutoff = NULL, output = "kable", sort) {

  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMiteminfit() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe"))

  # --- Validate and normalise cutoff ------------------------------------------
  cutoff_n_iter <- NULL
  if (!is.null(cutoff)) {
    if (is.list(cutoff) && !is.data.frame(cutoff) && "item_cutoffs" %in% names(cutoff)) {
      cutoff_n_iter <- cutoff$actual_iterations
      cutoff <- cutoff$item_cutoffs
    }
    if (!is.data.frame(cutoff)) {
      stop("`cutoff` must be NULL, the return value of RMinfitcutoff(), or its ",
           "$item_cutoffs data.frame.", call. = FALSE)
    }
    required_cols <- c("Item", "infit_low", "infit_high")
    missing_cols <- setdiff(required_cols, names(cutoff))
    if (length(missing_cols) > 0L) {
      stop("`cutoff` data.frame is missing required columns: ",
           paste(missing_cols, collapse = ", "), ".", call. = FALSE)
    }
  }

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

  # --- Apply cutoff if provided ------------------------------------------------
  if (!is.null(cutoff)) {
    data_items   <- item_fit_table$Item
    cutoff_items <- cutoff$Item
    if (!setequal(data_items, cutoff_items)) {
      stop(
        "Item names in `cutoff` do not match item names in `data`.\n",
        "  data items  : ", paste(data_items,   collapse = ", "), "\n",
        "  cutoff items: ", paste(cutoff_items, collapse = ", "),
        call. = FALSE
      )
    }
    cutoff_sub <- cutoff[, c("Item", "infit_low", "infit_high")]
    item_fit_table <- merge(item_fit_table, cutoff_sub, by = "Item", sort = FALSE)
    # Restore original row order (merge may reorder)
    item_fit_table <- item_fit_table[match(data_items, item_fit_table$Item), ]
    rownames(item_fit_table) <- NULL
    item_fit_table$Infit_low  <- round(item_fit_table$infit_low,  3)
    item_fit_table$Infit_high <- round(item_fit_table$infit_high, 3)
    item_fit_table$infit_low  <- NULL
    item_fit_table$infit_high <- NULL
    item_fit_table$Flagged <- item_fit_table$Infit_MSQ < item_fit_table$Infit_low |
                              item_fit_table$Infit_MSQ > item_fit_table$Infit_high
    # Reorder columns: Item, Infit_MSQ, Infit_low, Infit_high, Flagged, Relative_location
    item_fit_table <- item_fit_table[, c("Item", "Infit_MSQ", "Infit_low",
                                         "Infit_high", "Flagged",
                                         "Relative_location")]
  }

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
    col.names = if (is.null(cutoff)) {
      c("Item", "Infit MSQ", "Relative location")
    } else {
      c("Item", "Infit MSQ", "Infit low", "Infit high", "Flagged", "Relative location")
    },
    caption   = if (is.null(cutoff)) {
      paste0(
        "MSQ values based on conditional estimation (n = ",
        n_complete,
        " complete cases)."
      )
    } else if (!is.null(cutoff_n_iter)) {
      paste0(
        "MSQ values based on conditional estimation (n = ",
        n_complete,
        " complete cases). Cutoff values based on ",
        cutoff_n_iter,
        " simulation iterations."
      )
    } else {
      paste0(
        "MSQ values based on conditional estimation (n = ",
        n_complete,
        " complete cases). Simulation-based cutoff values applied."
      )
    }
  )
}
