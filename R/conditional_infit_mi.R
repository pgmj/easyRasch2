#' Conditional Item Infit MSQ for Multiply Imputed Data
#'
#' Extends \code{\link{RMiteminfit}} to work with multiply imputed datasets
#' produced by the `mice` package. Computes conditional infit MSQ on each
#' imputed dataset and pools the results using Rubin's rules.
#'
#' @param mids_object A `mids` object (multiply imputed dataset) as returned by
#'   `mice::mice()`. Each completed dataset must contain only the item response
#'   columns to be analysed (i.e., no ID or grouping variables). Items must be
#'   scored starting at 0 (non-negative integers).
#' @param cutoff Optional. Default `NULL` (no cutoff applied). Can be:
#'   * The return value of \code{\link{RMinfitcutoff}} or
#'     \code{\link{RMinfitcutoff_mi}} (a list with `$item_cutoffs`): the
#'     data.frame is extracted automatically and metadata is included in the
#'     kable caption.
#'   * The `$item_cutoffs` data.frame directly: must have columns `Item`,
#'     `infit_low`, and `infit_high`.
#'   When provided, adds columns `Infit_low`, `Infit_high`, and `Flagged` to
#'   the result.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying data.frame.
#' @param sort Optional character string. When `sort = "infit"`, rows are
#'   sorted by `Infit_MSQ` in descending order before output.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object (plain text table via
#'   `format = "pipe"`) with columns "Item", "Infit MSQ", "Infit SE",
#'   "Relative location", and a caption noting the number of imputations and
#'   complete cases. When `cutoff` is provided, columns "Infit low",
#'   "Infit high", and "Flagged" are also included.
#' * If `output = "dataframe"`: a data.frame with columns `Item`,
#'   `Infit_MSQ`, `Infit_SE`, and `Relative_location`. When `cutoff` is
#'   provided, columns `Infit_low`, `Infit_high`, and `Flagged` are also
#'   included (inserted after `Infit_SE`, before `Relative_location`).
#'
#' @details
#' For each of the `m` imputed datasets, the function:
#' \enumerate{
#'   \item Fits a Rasch model (`eRm::RM()` for dichotomous data or
#'     `eRm::PCM()` for polytomous data).
#'   \item Computes conditional infit MSQ and its standard error via
#'     `iarm::out_infit()`.
#'   \item Computes item and person locations.
#' }
#'
#' The per-imputation estimates are then pooled using Rubin's rules:
#' \describe{
#'   \item{Pooled MSQ}{The mean of the `m` infit MSQ point estimates.}
#'   \item{Within-imputation variance}{The mean of the `m` squared standard
#'     errors.}
#'   \item{Between-imputation variance}{The sample variance of the `m` point
#'     estimates.}
#'   \item{Total variance}{Within + (1 + 1/m) * Between.}
#'   \item{Pooled SE}{The square root of the total variance.}
#' }
#'
#' Relative item location is the mean of per-imputation relative locations
#' (item location minus sample mean person location).
#'
#' Imputed datasets that cause model convergence failures are dropped with a
#' warning. If all imputations fail, the function stops with an error. At
#' least two successful imputations are required to estimate between-imputation
#' variance.
#'
#' The `mice` and `iarm` packages must be installed (they are in Suggests, not
#' Imports).
#'
#' @seealso \code{\link{RMiteminfit}}, \code{\link{RMinfitcutoff_mi}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(mice)
#'
#' # Create example data with missing values
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
#' )
#' colnames(sim_data) <- paste0("Item", 1:8)
#' sim_data[sample(length(sim_data), 0.10 * length(sim_data))] <- NA
#'
#' # Impute
#' imp <- mice(sim_data, m = 5, method = "polr", seed = 123, printFlag = FALSE)
#'
#' # Pooled infit table (no cutoffs)
#' RMiteminfit_mi(imp)
#'
#' # With simulation-based cutoffs
#' cutoff_mi <- RMinfitcutoff_mi(imp, iterations = 250, parallel = FALSE,
#'                               seed = 42)
#' RMiteminfit_mi(imp, cutoff = cutoff_mi)
#'
#' # As data.frame
#' df <- RMiteminfit_mi(imp, cutoff = cutoff_mi, output = "dataframe")
#' }
RMiteminfit_mi <- function(mids_object, cutoff = NULL, output = "kable",
                           sort) {

  # --- Check required packages ------------------------------------------------
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop(
      "Package 'mice' is required for RMiteminfit_mi() but is not installed.\n",
      "Install it with: install.packages(\"mice\")",
      call. = FALSE
    )
  }

  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMiteminfit_mi() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  if (!inherits(mids_object, "mids")) {
    stop("`mids_object` must be a 'mids' object as returned by mice::mice().",
         call. = FALSE)
  }

  output <- match.arg(output, c("kable", "dataframe"))

  # --- Validate and normalise cutoff ------------------------------------------
  cutoff_n_iter     <- NULL
  cutoff_method     <- NULL
  cutoff_hdci_width <- NULL
  cutoff_n_imp      <- NULL
  if (!is.null(cutoff)) {
    if (is.list(cutoff) && !is.data.frame(cutoff) &&
        "item_cutoffs" %in% names(cutoff)) {
      cutoff_n_iter     <- cutoff$actual_iterations
      cutoff_method     <- cutoff$cutoff_method
      cutoff_hdci_width <- cutoff$hdci_width
      cutoff_n_imp      <- cutoff$n_imputations
      cutoff <- cutoff$item_cutoffs
    }
    if (!is.data.frame(cutoff)) {
      stop(
        "`cutoff` must be NULL, the return value of RMinfitcutoff() / ",
        "RMinfitcutoff_mi(), or its $item_cutoffs data.frame.",
        call. = FALSE
      )
    }
    required_cols <- c("Item", "infit_low", "infit_high")
    missing_cols  <- setdiff(required_cols, names(cutoff))
    if (length(missing_cols) > 0L) {
      stop("`cutoff` data.frame is missing required columns: ",
           paste(missing_cols, collapse = ", "), ".", call. = FALSE)
    }
  }

  # --- rgl workaround ---------------------------------------------------------
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  m <- mids_object$m

  # --- Compute infit on each imputed dataset ----------------------------------
  per_imp <- vector("list", m)
  n_failed <- 0L
  n_complete_first <- NULL

  for (i in seq_len(m)) {
    completed_data <- mice::complete(mids_object, action = i)

    result_i <- tryCatch({
      validate_response_data(completed_data)

      data_mat <- as.matrix(completed_data)
      n_complete <- nrow(completed_data)

      # Fit model
      if (max(data_mat, na.rm = TRUE) == 1L) {
        erm_out <- eRm::RM(completed_data)
        item_avg_locations <- stats::coef(erm_out, "beta") * -1
        pp <- eRm::person.parameter(erm_out)
        person_avg_location <- mean(
          pp$theta.table[["Person Parameter"]], na.rm = TRUE
        )
      } else {
        erm_out <- eRm::PCM(completed_data)
        thresh_obj   <- eRm::thresholds(erm_out)
        thresh_table <- thresh_obj$threshtable[[1]]
        if ("Location" %in% colnames(thresh_table)) {
          item_avg_locations <- thresh_table[, "Location"]
        } else {
          item_avg_locations <- rowMeans(thresh_table, na.rm = TRUE)
        }
        pp <- eRm::person.parameter(erm_out)
        person_avg_location <- mean(
          pp$theta.table[["Person Parameter"]], na.rm = TRUE
        )
      }

      relative_locations <- item_avg_locations - person_avg_location

      # Conditional infit
      cfit <- iarm::out_infit(erm_out)

      list(
        infit_msq         = cfit$Infit,
        infit_se          = cfit$Infit.se,
        relative_location = as.numeric(relative_locations),
        n_complete        = n_complete
      )
    }, error = function(e) {
      warning(
        sprintf("Model fitting failed for imputation %d: %s", i,
                conditionMessage(e)),
        call. = FALSE
      )
      NULL
    })

    if (is.null(result_i)) {
      n_failed <- n_failed + 1L
      next
    }

    per_imp[[i]] <- result_i

    if (is.null(n_complete_first)) {
      n_complete_first <- result_i$n_complete
    }
  }

  if (n_failed == m) {
    stop("Model fitting failed for all ", m, " imputed datasets. ",
         "Check your data.", call. = FALSE)
  }

  successful <- per_imp[!vapply(per_imp, is.null, logical(1L))]
  m_ok <- length(successful)

  if (m_ok < 2L) {
    stop(
      "Only ", m_ok, " imputation(s) succeeded. At least 2 are required ",
      "to estimate between-imputation variance.",
      call. = FALSE
    )
  }

  if (n_failed > 0L) {
    warning(
      sprintf("%d of %d imputed datasets failed and were excluded.",
              n_failed, m),
      call. = FALSE
    )
  }

  # --- Pool using Rubin's rules -----------------------------------------------
  item_names <- names(mice::complete(mids_object, action = 1L))
  n_items    <- length(item_names)

  # Collect per-imputation vectors into matrices (items x imputations)
  msq_mat  <- matrix(NA_real_, nrow = n_items, ncol = m_ok)
  se_mat   <- matrix(NA_real_, nrow = n_items, ncol = m_ok)
  loc_mat  <- matrix(NA_real_, nrow = n_items, ncol = m_ok)

  for (j in seq_len(m_ok)) {
    msq_mat[, j] <- successful[[j]]$infit_msq
    se_mat[, j]  <- successful[[j]]$infit_se
    loc_mat[, j] <- successful[[j]]$relative_location
  }

  # Pooled point estimate: mean across imputations
  pooled_msq <- rowMeans(msq_mat)

  # Within-imputation variance: mean of squared SEs
  within_var <- rowMeans(se_mat^2)

  # Between-imputation variance: sample variance of point estimates
  between_var <- apply(msq_mat, 1, stats::var)

  # Total variance (Rubin's combining rule)
  total_var <- within_var + (1 + 1 / m_ok) * between_var

  # Pooled SE
  pooled_se <- sqrt(total_var)

  # Pooled relative location: simple mean

  pooled_location <- rowMeans(loc_mat)

  # --- Assemble result data.frame ---------------------------------------------
  item_fit_table <- data.frame(
    Item              = item_names,
    Infit_MSQ         = round(pooled_msq, 3),
    Infit_SE          = round(pooled_se, 3),
    Relative_location = round(pooled_location, 2),
    stringsAsFactors  = FALSE,
    row.names         = NULL
  )

  # --- Apply cutoff if provided -----------------------------------------------
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
    item_fit_table <- merge(item_fit_table, cutoff_sub, by = "Item",
                            sort = FALSE)
    # Restore original row order
    item_fit_table <- item_fit_table[match(data_items, item_fit_table$Item), ]
    rownames(item_fit_table) <- NULL
    item_fit_table$Infit_low  <- round(item_fit_table$infit_low,  3)
    item_fit_table$Infit_high <- round(item_fit_table$infit_high, 3)
    item_fit_table$infit_low  <- NULL
    item_fit_table$infit_high <- NULL
    item_fit_table$Flagged <- item_fit_table$Infit_MSQ < item_fit_table$Infit_low |
      item_fit_table$Infit_MSQ > item_fit_table$Infit_high
    # Reorder columns
    item_fit_table <- item_fit_table[, c("Item", "Infit_MSQ", "Infit_SE",
                                         "Infit_low", "Infit_high", "Flagged",
                                         "Relative_location")]
  }

  # --- Sort if requested ------------------------------------------------------
  if (!missing(sort) && identical(sort, "infit")) {
    item_fit_table <- item_fit_table[order(item_fit_table$Infit_MSQ,
                                           decreasing = TRUE), ]
    rownames(item_fit_table) <- NULL
  }

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(item_fit_table)
  }

  # Build caption
  base_caption <- paste0(
    "Pooled MSQ values from ", m_ok, " imputations (Rubin's rules). ",
    "n = ", n_complete_first, " per imputed dataset."
  )

  cutoff_caption <- NULL
  if (!is.null(cutoff)) {
    if (!is.null(cutoff_n_iter)) {
      method_label <- .format_cutoff_method_label(cutoff_method,
                                                  cutoff_hdci_width)
      iter_part <- paste0(cutoff_n_iter, " total simulation iterations")
      imp_part  <- if (!is.null(cutoff_n_imp)) {
        paste0(" across ", cutoff_n_imp, " imputations")
      } else {
        ""
      }
      cutoff_caption <- paste0(
        " Cutoff values based on ",
        if (!is.null(method_label)) {
          paste0(iter_part, imp_part, " (", method_label, ").")
        } else {
          paste0(iter_part, imp_part, ".")
        }
      )
    } else {
      cutoff_caption <- " Simulation-based cutoff values applied."
    }
  }

  caption <- paste0(base_caption, if (!is.null(cutoff_caption)) cutoff_caption)

  knitr::kable(
    item_fit_table,
    format    = "pipe",
    col.names = if (is.null(cutoff)) {
      c("Item", "Infit MSQ", "Infit SE", "Relative location")
    } else {
      c("Item", "Infit MSQ", "Infit SE", "Infit low", "Infit high",
        "Flagged", "Relative location")
    },
    caption = caption
  )
}
