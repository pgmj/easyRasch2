#' Q3 Residual Correlations for Local Dependence Assessment
#'
#' Computes Yen's Q3 residual correlations between item pairs using a Rasch
#' model fitted via Marginal Maximum Likelihood (MML) in `mirt`. High
#' correlations (above the dynamic cut-off) indicate potential local dependence
#' between items.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed.
#' @param cutoff A single numeric value added to the mean off-diagonal Q3
#'   correlation to produce the dynamic cut-off threshold. A common choice is
#'   `0.2`. Use a helper such as `RIgetResidCor()` from **easyRasch** to derive
#'   a simulation-based cutoff.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying numeric data.frame.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object showing the lower triangle
#'   of the Q3 correlation matrix together with a footnote describing the
#'   dynamic cut-off.
#' * If `output = "dataframe"`: a data.frame (rounded to 2 decimal places) with
#'   the lower triangle of the Q3 correlation matrix; the upper triangle and
#'   diagonal are set to `NA`.
#'
#' @details
#' The Q3 statistic (Yen, 1984) is the correlation between residuals of pairs
#' of items after accounting for the latent trait. Under local independence,
#' Q3 values are expected to be around \eqn{-1/(k-1)} where \eqn{k} is the
#' number of items. The dynamic cut-off used here is the mean of all
#' off-diagonal Q3 values plus `cutoff`, following the approach of Christensen
#' et al. (2017).
#'
#' `mirt` is used for model fitting here because Q3 requires model-based
#' expected responses, which are most readily available from MML estimation.
#'
#' @references
#' Yen, W. M. (1984). Effects of local item dependence on the fit and
#' equating performance of the three-parameter logistic model.
#' *Applied Psychological Measurement*, 8(2), 125–145.
#'
#' Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
#' for Yen's Q3: Identification of local dependence in the Rasch model.
#' *Applied Psychological Measurement*, 41(3), 178–194.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate binary item response data (10 items, 200 persons)
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#' )
#' colnames(sim_data) <- paste0("Item", 1:10)
#'
#' # Display Q3 table with a cutoff of 0.2 above the mean
#' RMlocdepQ3(sim_data, cutoff = 0.2)
#'
#' # Get the underlying data.frame
#' q3_df <- RMlocdepQ3(sim_data, cutoff = 0.2, output = "dataframe")
#' }
RMlocdepQ3 <- function(data, cutoff, output = "kable") {
  # --- Input validation -------------------------------------------------------
  if (missing(cutoff)) {
    stop(
      "Please provide a `cutoff` value (e.g., 0.2). ",
      "See `RIgetResidCor()` in the easyRasch package for a simulation-based approach.",
      call. = FALSE
    )
  }

  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
    stop("`cutoff` must be a single numeric value.", call. = FALSE)
  }

  output <- match.arg(output, c("kable", "dataframe"))

  validate_response_data(data)

  # --- Model fitting ----------------------------------------------------------
  mirt_model <- mirt::mirt(
    data,
    model = 1,
    itemtype = "Rasch",
    verbose = FALSE,
    accelerate = "squarem"
  )

  # --- Compute Q3 residual correlations ---------------------------------------
  resid_mat <- mirt::residuals(mirt_model, type = "Q3", digits = 2, verbose = FALSE)

  diag(resid_mat) <- NA
  resid_df <- as.data.frame(resid_mat)

  # Mean of off-diagonal correlations and dynamic cut-off
  mean_resid <- mean(as.matrix(resid_df), na.rm = TRUE)
  dyn_cutoff <- mean_resid + cutoff

  # Round values
  resid_df[] <- lapply(resid_df, function(x) round(x, 2))

  # Blank out upper triangle and diagonal (set to NA)
  resid_df[upper.tri(resid_df)] <- NA
  diag(resid_df) <- NA

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(resid_df)
  }

  # "kable" output
  footnote_text <- paste0(
    "Dynamic cut-off: ", round(dyn_cutoff, 3),
    " (mean Q3 = ", round(mean_resid, 3),
    " + ", round(cutoff, 3), ").",
    " Correlations exceeding the cut-off may indicate local dependence."
  )

  # Replace NA with "" for display (knitr::kable doesn't render NA nicely)
  resid_display <- as.data.frame(
    lapply(resid_df, function(x) {
      ifelse(is.na(x), "", as.character(x))
    }),
    stringsAsFactors = FALSE
  )
  rownames(resid_display) <- rownames(resid_df)

  knitr::kable(
    resid_display,
    caption = footnote_text,
    format = "simple"
  )
}
