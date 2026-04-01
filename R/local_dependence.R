#' Q3 Residual Correlations for Local Dependence Assessment
#'
#' Computes Yen's Q3 residual correlations between item pairs using a Rasch
#' model fitted via Marginal Maximum Likelihood (MML) in `mirt`. When a
#' `cutoff` is provided, correlations above the dynamic cut-off indicate
#' potential local dependence between items.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed.
#' @param cutoff Optional. A single numeric value added to the mean off-diagonal
#'   Q3 correlation to produce the dynamic cut-off threshold. When `NULL`
#'   (default), the raw Q3 matrix is returned without any cut-off applied.
#'   Use `RMlocdepQ3cutoff()` to derive a simulation-based cutoff value.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying numeric data.frame.
#'
#' @return
#' * If `output = "kable"` and `cutoff = NULL`: a `knitr_kable` object showing
#'   the lower triangle of the Q3 correlation matrix with no cut-off annotation.
#' * If `output = "kable"` and `cutoff` is provided: a `knitr_kable` object
#'   with a caption describing the dynamic cut-off.
#' * If `output = "dataframe"`: a data.frame (rounded to 2 decimal places) with
#'   the lower triangle of the Q3 correlation matrix; the upper triangle and
#'   diagonal are set to `NA`.
#'
#' @details
#' The Q3 statistic (Yen, 1984) is the correlation between residuals of pairs
#' of items after accounting for the latent trait. Under local independence,
#' Q3 values are expected to be around \eqn{-1/(k-1)} where \eqn{k} is the
#' number of items. When `cutoff` is supplied, the dynamic cut-off is the mean
#' of all off-diagonal Q3 values plus `cutoff`, following the approach of
#' Christensen et al. (2017).
#'
#' `mirt` is used for model fitting here because Q3 requires model-based
#' expected responses, which are most readily available from MML estimation.
#'
#' The recommended workflow is to first run `RMlocdepQ3cutoff()` to obtain a
#' simulation-based cutoff, then pass that value to `RMlocdepQ3()`.
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
#' # Display raw Q3 matrix (no cutoff)
#' RMlocdepQ3(sim_data)
#'
#' # Get simulation-based cutoff, then display with cutoff
#' cutoff_res <- RMlocdepQ3cutoff(sim_data)
#' RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
#'
#' # Get the underlying data.frame
#' q3_df <- RMlocdepQ3(sim_data, output = "dataframe")
#' }
RMlocdepQ3 <- function(data, cutoff = NULL, output = "kable") {
  # --- Input validation -------------------------------------------------------
  if (!is.null(cutoff)) {
    if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
      stop("`cutoff` must be a single numeric value or NULL.", call. = FALSE)
    }
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
  if (!is.null(cutoff)) {
    mean_resid <- mean(as.matrix(resid_df), na.rm = TRUE)
    dyn_cutoff <- mean_resid + cutoff
    caption_text <- paste0(
      "Dynamic cut-off: ", round(dyn_cutoff, 3),
      " (mean Q3 = ", round(mean_resid, 3),
      " + ", round(cutoff, 3), ").",
      " Correlations exceeding the cut-off may indicate local dependence.",
      " Use RMlocdepQ3cutoff() to derive a simulation-based cutoff."
    )
  } else {
    caption_text <- "Q3 residual correlation matrix. Consider using RMlocdepQ3cutoff() to determine a simulation-based cutoff."
  }

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
    caption = caption_text,
    format = "simple"
  )
}

#' Simulation-Based Cutoff for Q3 Residual Correlations
#'
#' Uses parametric bootstrap simulation to determine an appropriate cutoff
#' value for Yen's Q3 residual correlations. Item parameters are estimated
#' from the observed data, then many datasets are simulated under the Rasch
#' model to build a null distribution of Q3 statistics. The resulting
#' percentile values can be used as cutoffs for `RMlocdepQ3()`.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed.
#'   Both dichotomous (max = 1) and polytomous (max > 1) data are supported.
#' @param iterations Integer. Number of simulation iterations. Default `500`.
#' @param cpu Integer. Number of CPU cores for parallel processing. Default `4`.
#' @param seed Integer. Random seed for reproducibility. Default `123`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{`results`}{data.frame with columns `mean`, `max`, and `diff`
#'     (one row per successful iteration).}
#'   \item{`actual_iterations`}{Number of successful iterations (may be less
#'     than `iterations` if some simulated datasets had too few responses).}
#'   \item{`sample_n`}{Sample size (number of rows in `data`).}
#'   \item{`sample_summary`}{Summary statistics of the estimated theta
#'     distribution used for simulation.}
#'   \item{`max_diff`}{Maximum value of the `diff` column.}
#'   \item{`sd_diff`}{Standard deviation of the `diff` column.}
#'   \item{`p95`}{95th percentile of `diff` (recommended default cutoff).}
#'   \item{`p99`}{99th percentile of `diff`.}
#'   \item{`p995`}{99.5th percentile of `diff`.}
#'   \item{`p999`}{99.9th percentile of `diff`.}
#'   \item{`suggested_cutoff`}{The 95th percentile value, recommended for use
#'     as the `cutoff` argument to `RMlocdepQ3()`.}
#' }
#'
#' @details
#' For **dichotomous** data (max response = 1), item parameters are estimated
#' with `eRm::RM()` and person parameters with `eRm::person.parameter()`.
#' Simulated datasets are generated using `psychotools::rrm()`.
#'
#' For **polytomous** data (max response > 1), item threshold parameters are
#' estimated with `eRm::PCM()` and person parameters with
#' `eRm::person.parameter()`. Simulated datasets are generated using the
#' internal `sim_partial_score()` helper.
#'
#' Each simulated dataset is analysed with a Rasch model in `mirt` and Q3
#' residual correlations are computed. The `diff` column in `results` is
#' `max - mean` Q3 across items, which provides a scale-free measure of
#' the maximum local dependence expected by chance. The 95th percentile of
#' this distribution is returned as `suggested_cutoff`.
#'
#' Parallel processing uses `foreach`, `doParallel`, and `doRNG` for
#' reproducible results across platforms.
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
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#' )
#' colnames(sim_data) <- paste0("Item", 1:10)
#'
#' # Determine simulation-based cutoff
#' cutoff_res <- RMlocdepQ3cutoff(sim_data, iterations = 200, cpu = 2)
#' cutoff_res$suggested_cutoff
#'
#' # Use the suggested cutoff with RMlocdepQ3()
#' RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
#' }
RMlocdepQ3cutoff <- function(data, iterations = 500, cpu = 4, seed = 123) {
  # --- Input validation -------------------------------------------------------
  validate_response_data(data)

  if (!is.numeric(iterations) || length(iterations) != 1L || iterations < 1L) {
    stop("`iterations` must be a positive integer.", call. = FALSE)
  }
  iterations <- as.integer(iterations)

  if (!is.numeric(cpu) || length(cpu) != 1L || cpu < 1L) {
    stop("`cpu` must be a positive integer.", call. = FALSE)
  }
  cpu <- as.integer(cpu)

  if (!is.numeric(seed) || length(seed) != 1L) {
    stop("`seed` must be a single numeric value.", call. = FALSE)
  }

  sample_n <- nrow(data)
  data_mat <- as.matrix(data)
  data_max <- max(data_mat, na.rm = TRUE)
  data_min <- min(data_mat, na.rm = TRUE)

  # --- Set up parallel backend ------------------------------------------------
  doParallel::registerDoParallel(cores = cpu)
  doRNG::registerDoRNG(seed)
  on.exit(foreach::registerDoSEQ(), add = TRUE)

  message("Starting RMlocdepQ3cutoff() with ", iterations, " iterations on ",
          cpu, " cores...")

  # --- Estimate item/person parameters ----------------------------------------
  if (data_max == 1L && data_min == 0L) {
    # Dichotomous branch
    erm_out <- eRm::RM(data)
    item_locations <- erm_out$betapar * -1
    names(item_locations) <- names(as.data.frame(data))
    thetas <- eRm::person.parameter(erm_out)[["theta.table"]][["Person Parameter"]]
    item_type <- "dichotomous"
  } else {
    # Polytomous branch
    erm_out <- eRm::PCM(data)
    item_thresholds <- eRm::thresholds(erm_out)$threshpar
    thetas <- eRm::person.parameter(erm_out)[["theta.table"]][["Person Parameter"]]
    n_items <- nrow(item_thresholds)
    itemlist <- vector("list", n_items)
    for (i in seq_len(n_items)) {
      itemlist[[i]] <- list(stats::na.omit(item_thresholds[i, ]))
    }
    item_type <- "polytomous"
  }

  # --- Run simulation in parallel ---------------------------------------------
  i <- NULL  # avoid R CMD check NOTE for foreach variable
  if (item_type == "dichotomous") {
    residcor <- foreach::foreach(i = seq_len(iterations)) %dopar% {
      input_thetas <- sample(thetas, size = sample_n, replace = TRUE)
      test_data <- as.data.frame(
        psychotools::rrm(input_thetas, item_locations, return_setting = FALSE)
      )
      # Check minimum column sums (need at least 8 responses per item for
      # stable model estimation; see Hambleton & Swaminathan, 1985)
      n_resp <- colSums(as.matrix(test_data), na.rm = TRUE)
      if (min(n_resp, na.rm = TRUE) < 8L) {
        return("Missing cells in generated data.")
      }
      mirt_fit <- mirt::mirt(test_data, model = 1, itemtype = "Rasch",
                             verbose = FALSE, accelerate = "squarem")
      resid <- mirt::residuals(mirt_fit, type = "Q3", digits = 2, verbose = FALSE)
      diag(resid) <- NA
      data.frame(mean = mean(resid, na.rm = TRUE), max = max(resid, na.rm = TRUE))
    }
  } else {
    residcor <- foreach::foreach(i = seq_len(iterations)) %dopar% {
      input_thetas <- sample(thetas, size = sample_n, replace = TRUE)
      test_data <- as.data.frame(
        sim_partial_score(deltaslist = itemlist, thetavec = input_thetas)
      )
      names(test_data) <- names(as.data.frame(data))
      # Check that each item has all expected response categories
      ok <- vapply(test_data, function(col) {
        all(seq(0L, data_max) %in% col)
      }, logical(1L))
      if (!all(ok)) {
        return("Missing categories in generated data.")
      }
      mirt_fit <- mirt::mirt(test_data, model = 1, itemtype = "Rasch",
                             verbose = FALSE, accelerate = "squarem")
      resid <- mirt::residuals(mirt_fit, type = "Q3", digits = 2, verbose = FALSE)
      diag(resid) <- NA
      data.frame(mean = mean(resid, na.rm = TRUE), max = max(resid, na.rm = TRUE))
    }
  }

  # --- Post-process results ---------------------------------------------------
  nodata <- vapply(residcor, is.character, logical(1L))
  iterations_nodata <- which(nodata)
  actual_iterations <- iterations - length(iterations_nodata)

  if (actual_iterations == 0L) {
    stop("All simulation iterations failed. Check your data.", call. = FALSE)
  }

  if (actual_iterations == iterations) {
    good_results <- residcor
  } else {
    good_results <- residcor[-iterations_nodata]
  }

  results <- do.call(rbind, good_results)
  results$diff <- results$max - results$mean

  message("RMlocdepQ3cutoff() completed: ", actual_iterations, " of ",
          iterations, " iterations succeeded.")

  out <- list()
  out$results          <- results
  out$actual_iterations <- actual_iterations
  out$sample_n         <- sample_n
  out$sample_summary   <- summary(thetas)
  out$max_diff         <- max(results$diff)
  out$sd_diff          <- stats::sd(results$diff)
  out$p95              <- stats::quantile(results$diff, 0.95)
  out$p99              <- stats::quantile(results$diff, 0.99)
  out$p995             <- stats::quantile(results$diff, 0.995)
  out$p999             <- stats::quantile(results$diff, 0.999)
  out$suggested_cutoff <- unname(stats::quantile(results$diff, 0.95))
  out
}
