#' Partial Gamma DIF Analysis
#'
#' Computes partial gamma coefficients for Differential Item Functioning (DIF)
#' using \code{iarm::partgam_DIF()}. Each item is tested for association with
#' a single categorical DIF variable, controlling for the total score.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed,
#'   but at least one complete case must exist after combining `data` and
#'   `dif_var`.
#' @param dif_var A vector (factor or character) of the same length as
#'   `nrow(data)`, representing the grouping variable for DIF analysis.
#' @param cutoff Optional. Default `NULL` (no cutoff applied). Can be:
#'   * The return value of \code{\link{RMpgDIFcutoff}} (a list with
#'     `$item_cutoffs`): the data.frame is extracted automatically and
#'     simulation metadata is included in the kable caption.
#'   * The `$item_cutoffs` data.frame from \code{\link{RMpgDIFcutoff}}
#'     directly: must have columns `Item`, `gamma_low`, `gamma_high`.
#'   When provided, adds columns `Gamma_low`, `Gamma_high`, and `Flagged`
#'   (logical; `TRUE` when the observed partial gamma falls outside the
#'   credible range) to the result.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying data.frame.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object with columns "Item",
#'   "Partial gamma", "SE", "Lower CI", "Upper CI", and "Adjusted p-value
#'   (BH)". When `cutoff` is provided, additional columns "Gamma low",
#'   "Gamma high", and "Flagged" are included.
#' * If `output = "dataframe"`: a data.frame with columns `Item`, `gamma`,
#'   `se`, `lower`, `upper`, `padj_bh`. When `cutoff` is provided, columns
#'   `gamma_low`, `gamma_high`, and `flagged` are also included.
#'
#' @details
#' Partial gamma (Bjorner et al., 1998) measures the association between item
#' response and an exogenous grouping variable, controlling for the total
#' score. Values near 0 indicate no DIF. Recommended interpretive thresholds
#' (Bjorner et al., 1998):
#'
#' * **No or negligible DIF**: gamma within \eqn{[-0.21, 0.21]}, *or* gamma
#'   not significantly different from 0.
#' * **Slight to moderate DIF**: gamma within \eqn{[-0.31, 0.31]} (and outside
#'   \eqn{[-0.21, 0.21]}), *or* not significantly outside \eqn{[-0.21, 0.21]}.
#' * **Moderate to large DIF**: gamma outside \eqn{[-0.31, 0.31]}, **and**
#'   significantly outside \eqn{[-0.21, 0.21]}.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' @references
#' Bjorner, J., Kreiner, S., Ware, J., Damsgaard, M. and Bech, P.
#' Differential item functioning in the Danish translation of the SF-36.
#' *Journal of Clinical Epidemiology*, 51(11), 1998, pp. 1189-1202.
#'
#' @seealso \code{\link{RMpgDIFcutoff}}
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
#' dif_group <- factor(sample(c("A", "B"), 200, replace = TRUE))
#'
#' # Default kable output
#' RMpartgamDIF(sim_data, dif_group)
#'
#' # Return as data.frame
#' RMpartgamDIF(sim_data, dif_group, output = "dataframe")
#'
#' # With simulation-based cutoffs
#' cutoff_res <- RMpgDIFcutoff(sim_data, dif_var = dif_group,
#'                                   iterations = 100, parallel = FALSE,
#'                                   seed = 42)
#' RMpartgamDIF(sim_data, dif_group, cutoff = cutoff_res)
#' }
RMpartgamDIF <- function(data, dif_var, cutoff = NULL, output = "kable") {

  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMpartgamDIF() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe"))

  validate_response_data(data)

  # --- Validate dif_var -------------------------------------------------------
  if (length(dif_var) != nrow(data)) {
    stop(
      "`dif_var` must have the same length as `nrow(data)` (",
      nrow(data), "), but has length ", length(dif_var), ".",
      call. = FALSE
    )
  }

  # Require at least 2 levels
  dif_levels <- length(unique(stats::na.omit(dif_var)))
  if (dif_levels < 2L) {
    stop(
      "`dif_var` must have at least 2 non-missing levels, but has ", dif_levels, ".",
      call. = FALSE
    )
  }

  # --- Validate and normalise cutoff ------------------------------------------
  cutoff_n_iter     <- NULL
  cutoff_method     <- NULL
  cutoff_hdci_width <- NULL
  if (!is.null(cutoff)) {
    if (is.list(cutoff) && !is.data.frame(cutoff) && "item_cutoffs" %in% names(cutoff)) {
      cutoff_n_iter     <- cutoff$actual_iterations
      cutoff_method     <- cutoff$cutoff_method
      cutoff_hdci_width <- cutoff$hdci_width
      cutoff <- cutoff$item_cutoffs
    }
    if (!is.data.frame(cutoff)) {
      stop(
        "`cutoff` must be NULL, the return value of RMpgDIFcutoff(), or its ",
        "$item_cutoffs data.frame.",
        call. = FALSE
      )
    }
    required_cols <- c("Item", "gamma_low", "gamma_high")
    missing_cols <- setdiff(required_cols, names(cutoff))
    if (length(missing_cols) > 0L) {
      stop(
        "`cutoff` data.frame is missing required columns: ",
        paste(missing_cols, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # --- rgl workaround ---------------------------------------------------------
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # --- Compute partial gamma DIF via iarm -------------------------------------
  sink(nullfile())
  #on.exit(sink(), add = TRUE)
  pgam_raw <- iarm::partgam_DIF(as.data.frame(data), dif_var)
  sink()

  # pgam_raw is a data.frame with character columns that need numeric conversion
  pgam_df <- data.frame(
    Item    = as.character(pgam_raw$Item),
    gamma   = as.numeric(pgam_raw$gamma),
    se      = as.numeric(pgam_raw$se),
    pvalue  = as.numeric(pgam_raw$pvalue),
    padj_bh = as.numeric(pgam_raw[[6]]),
    sig     = as.character(pgam_raw$sig),
    lower   = as.numeric(pgam_raw$lower),
    upper   = as.numeric(pgam_raw$upper),
    stringsAsFactors = FALSE
  )

  # Round numeric columns
  pgam_df$gamma   <- round(pgam_df$gamma, 3)
  pgam_df$se      <- round(pgam_df$se, 3)
  pgam_df$pvalue  <- round(pgam_df$pvalue, 3)
  pgam_df$padj_bh <- round(pgam_df$padj_bh, 3)
  pgam_df$lower   <- round(pgam_df$lower, 3)
  pgam_df$upper   <- round(pgam_df$upper, 3)

  # Keep output columns
  result_df <- pgam_df[, c("Item", "gamma", "se", "lower", "upper", "padj_bh")]

  # --- Apply cutoff if provided -----------------------------------------------
  if (!is.null(cutoff)) {
    data_items   <- result_df$Item
    cutoff_items <- cutoff$Item
    if (!setequal(data_items, cutoff_items)) {
      stop(
        "Item names in `cutoff` do not match item names in `data`.\n",
        "  data items  : ", paste(data_items, collapse = ", "), "\n",
        "  cutoff items: ", paste(cutoff_items, collapse = ", "),
        call. = FALSE
      )
    }
    cutoff_sub <- cutoff[, c("Item", "gamma_low", "gamma_high")]
    result_df <- merge(result_df, cutoff_sub, by = "Item", sort = FALSE)
    # Restore original row order
    result_df <- result_df[match(data_items, result_df$Item), ]
    rownames(result_df) <- NULL
    result_df$gamma_low  <- round(result_df$gamma_low, 3)
    result_df$gamma_high <- round(result_df$gamma_high, 3)
    result_df$flagged <- result_df$gamma < result_df$gamma_low |
      result_df$gamma > result_df$gamma_high
    # Reorder columns
    result_df <- result_df[, c("Item", "gamma", "se", "lower", "upper",
                               "padj_bh", "gamma_low", "gamma_high", "flagged")]
  }

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(result_df)
  }

  # Build caption
  n_complete <- sum(stats::complete.cases(cbind(as.data.frame(data), dif_var = dif_var)))
  if (is.null(cutoff)) {
    caption_text <- paste0(
      "Partial gamma DIF analysis (n = ", n_complete, " complete cases). ",
      "Positive gamma indicates higher scores in higher DIF group levels."
    )
  } else if (!is.null(cutoff_n_iter)) {
    method_label <- .format_gamma_cutoff_method_label(cutoff_method, cutoff_hdci_width)
    iter_part <- paste0(cutoff_n_iter, " simulation iterations")
    caption_text <- paste0(
      "Partial gamma DIF analysis (n = ", n_complete, " complete cases). Cutoff values based on ",
      if (!is.null(method_label)) {
        paste0(iter_part, " (", method_label, ").")
      } else {
        paste0(iter_part, ".")
      }
    )
  } else {
    caption_text <- paste0(
      "Partial gamma DIF analysis (n = ", n_complete,
      " complete cases). Simulation-based cutoff values applied."
    )
  }

  knitr::kable(
    result_df,
    format    = "pipe",
    col.names = if (is.null(cutoff)) {
      c("Item", "Partial gamma", "SE", "Lower CI", "Upper CI", "Adjusted p-value (BH)")
    } else {
      c("Item", "Partial gamma", "SE", "Lower CI", "Upper CI",
        "Adjusted p-value (BH)", "Gamma low", "Gamma high", "Flagged")
    },
    caption   = caption_text
  )
}

# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

#' Format a human-readable label for the gamma cutoff method
#'
#' @param cutoff_method Character. `"hdci"`, `"quantile"`, or `NULL`.
#' @param hdci_width Numeric or `NULL`. HDCI width (e.g., `0.99`).
#' @return A character label, or `NULL` if the method is unknown/unset.
#' @keywords internal
.format_gamma_cutoff_method_label <- function(cutoff_method, hdci_width) {
  if (is.null(cutoff_method)) {
    return(NULL)
  }
  switch(cutoff_method,
         quantile = "2.5th/97.5th percentile",
         hdci     = if (!is.null(hdci_width)) paste0(hdci_width * 100, "% HDCI") else "HDCI",
         NULL
  )
}


#' Simulation-Based Partial Gamma DIF Cutoff Determination
#'
#' Uses parametric bootstrap simulation to determine appropriate cutoff values
#' for partial gamma DIF analysis via \code{\link[iarm]{partgam_DIF}}. Under
#' a correctly fitting Rasch model where the DIF variable is unrelated to item
#' responses (i.e., no true DIF), this function generates the expected
#' distribution of absolute partial gamma values per item, providing empirical
#' critical values.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Only complete cases (rows without
#'   any `NA`) are used.
#' @param dif_var A vector (factor, character, or integer) defining group
#'   membership for DIF analysis. Must have the same length as `nrow(data)`.
#'   The actual group labels are used to determine the number of groups and
#'   their relative sizes; during simulation, respondents are randomly assigned
#'   to groups with the same proportions, so there is no true DIF by
#'   construction.
#' @param iterations Integer. Number of simulation iterations (default 250).
#' @param parallel Logical. Use parallel processing via `mirai` if available
#'   (default `TRUE`).
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first. If neither is set and
#'   `parallel = TRUE`, a warning is issued and execution falls back to
#'   sequential (single core) processing.
#' @param verbose Logical. Show a progress bar (default `FALSE`).
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param cutoff_method Character string specifying how cutoff intervals are
#'   computed. Either `"hdci"` (default) for the Highest Density Interval via
#'   `ggdist::hdci()`, or `"quantile"` for the 2.5th/97.5th percentiles via
#'   `stats::quantile()`.
#' @param hdci_width Numeric. Width of the HDCI when `cutoff_method = "hdci"`.
#'   Default is `0.99` (99\% HDCI). Ignored when
#'   `cutoff_method = "quantile"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`results`}{data.frame with columns `iteration`, `Item`, and
#'     `gamma` (one row per item per successful iteration).}
#'   \item{`item_cutoffs`}{data.frame with per-item cutoff summaries: `Item`,
#'     `gamma_low`, `gamma_high`. Bounds are computed using the method
#'     specified by `cutoff_method`.}
#'   \item{`actual_iterations`}{Number of successful iterations.}
#'   \item{`sample_n`}{Number of complete cases used.}
#'   \item{`sample_summary`}{Summary statistics of estimated person
#'     parameters.}
#'   \item{`item_names`}{Character vector of item names from data.}
#'   \item{`dif_group_sizes`}{Named integer vector of group sizes used in the
#'     simulation (matches proportions in the observed `dif_var`).}
#'   \item{`cutoff_method`}{The method used to compute cutoffs (`"hdci"` or
#'     `"quantile"`).}
#'   \item{`hdci_width`}{The HDCI width used (only meaningful when
#'     `cutoff_method = "hdci"`).}
#' }
#'
#' @details
#' For each simulation iteration the function:
#' \enumerate{
#'   \item Resamples person parameters (thetas) with replacement from ML
#'     estimates.
#'   \item Simulates item response data under a Rasch model (dichotomous via
#'     `psychotools::rrm()` or polytomous via an internal partial credit
#'     simulator).
#'   \item Creates a random DIF variable by sampling group labels with the
#'     same proportions as the observed `dif_var`, so there is **no true DIF**
#'     by construction.
#'   \item Computes partial gamma DIF statistics via
#'     `iarm::partgam_DIF()`.
#' }
#'
#' The distribution of partial gamma values across iterations provides
#' empirical critical values per item. Values from real data that fall
#' outside these bounds suggest DIF that exceeds what would be expected by
#' chance under a correctly fitting Rasch model. Failed iterations (e.g.,
#' due to convergence issues or degenerate data) are silently discarded.
#'
#' Supports both **dichotomous** data (via `eRm::RM()` and
#' `psychotools::rrm()`) and **polytomous** data (via `eRm::PCM()` and an
#' internal partial credit score simulator).
#'
#' Parallel processing is provided by the `mirai` package (optional). Install
#' it with `install.packages("mirai")` to enable parallelisation.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' @references
#' Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., &
#' Bech, P. (1998). Differential item functioning in the Danish translation
#' of the SF-36. *Journal of Clinical Epidemiology*, 51(11), 1189--1202.
#'
#' Henninger, M., Radek, J., Sengewald, M.-A., & Strobl, C. (2024).
#' Partial credit trees meet the partial gamma coefficient for quantifying
#' DIF and DSF in polytomous items. OSF Preprints.
#' \doi{10.31234/osf.io/47sah}
#'
#' @seealso \code{\link[iarm]{partgam_DIF}}
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
#' dif_sex <- sample(c("male", "female"), 200, replace = TRUE)
#'
#' # Run 100 iterations sequentially for a quick demo
#' cutoff_res <- RMpgDIFcutoff(sim_data, dif_var = dif_sex,
#'                                   iterations = 100, parallel = FALSE,
#'                                   seed = 42)
#' cutoff_res$item_cutoffs
#' }
RMpgDIFcutoff <- function(data, dif_var, iterations = 250,
                               parallel = TRUE, n_cores = NULL,
                               verbose = FALSE, seed = NULL,
                               cutoff_method = "hdci",
                               hdci_width = 0.99) {

  cutoff_method <- match.arg(cutoff_method, c("hdci", "quantile"))

  if (cutoff_method == "hdci" && !requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required when cutoff_method = \"hdci\" but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")\n",
      "Alternatively, use cutoff_method = \"quantile\" to avoid this dependency.",
      call. = FALSE
    )
  }

  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMpgDIFcutoff() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  validate_response_data(data)

  # --- Validate dif_var -------------------------------------------------------
  if (missing(dif_var) || is.null(dif_var)) {
    stop("`dif_var` must be provided.", call. = FALSE)
  }
  if (length(dif_var) != nrow(data)) {
    stop(
      "`dif_var` must have the same length as `nrow(data)` (",
      nrow(data), "), but has length ", length(dif_var), ".",
      call. = FALSE
    )
  }

  # rgl workaround (iarm depends on vcdExtra -> rgl)
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # Only complete cases (both data and dif_var)
  complete_idx <- stats::complete.cases(data) & !is.na(dif_var)
  data <- data[complete_idx, , drop = FALSE]
  dif_var <- dif_var[complete_idx]

  if (nrow(data) == 0L) {
    stop("No complete cases in data after removing rows with NA in data or dif_var.",
         call. = FALSE)
  }

  # Determine DIF group structure (labels and proportions)
  dif_table <- table(dif_var)
  dif_levels <- names(dif_table)
  dif_proportions <- as.numeric(dif_table) / sum(dif_table)
  n_dif_groups <- length(dif_levels)

  if (n_dif_groups < 2L) {
    stop("`dif_var` must have at least 2 distinct groups, but has ",
         n_dif_groups, ".", call. = FALSE)
  }

  # --- Parallel setup ---------------------------------------------------------
  use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)

  if (parallel && !use_parallel) {
    message("Install 'mirai' package for parallel processing: install.packages(\"mirai\")")
    message("Running sequentially...")
  }

  if (use_parallel) {
    if (is.null(n_cores)) {
      n_cores <- getOption("mc.cores")
    }
    if (is.null(n_cores)) {
      warning(
        paste0(
          "For parallel processing, specify n_cores or set options(mc.cores = N).\n",
          "Your computer appears to have ", parallel::detectCores(), " cores available.\n",
          "Falling back to sequential (single core) processing."
        ),
        call. = FALSE
      )
      use_parallel <- FALSE
    } else {
      n_cores <- min(n_cores, iterations)
    }
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate per-iteration seeds
  sim_seeds <- sample.int(.Machine$integer.max, iterations)

  data_mat <- as.matrix(data)
  sample_n <- nrow(data_mat)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  item_names_vec <- colnames(data_mat)

  if (is_polytomous) {
    pcm_fit <- eRm::PCM(data_mat)
    pp <- eRm::person.parameter(pcm_fit)
    theta_table <- pp$theta.table[["Person Parameter"]]
    raw_scores <- rowSums(data_mat, na.rm = TRUE)
    thetas <- as.numeric(stats::na.omit(theta_table[raw_scores]))
    thresh_mat <- extract_item_thresholds(data_mat)
    deltaslist <- lapply(seq_len(nrow(thresh_mat)), function(i) {
      as.numeric(thresh_mat[i, !is.na(thresh_mat[i, ])])
    })
    sim_data_list <- list(
      type = "polytomous",
      thetas = thetas,
      deltaslist = deltaslist,
      n_items = ncol(data_mat),
      sample_n = sample_n,
      item_names = item_names_vec,
      dif_levels = dif_levels,
      dif_proportions = dif_proportions
    )
  } else {
    rm_fit <- eRm::RM(data_mat)
    pp <- eRm::person.parameter(rm_fit)
    theta_table <- pp$theta.table[["Person Parameter"]]
    raw_scores <- rowSums(data_mat, na.rm = TRUE)
    thetas <- as.numeric(stats::na.omit(theta_table[raw_scores]))
    item_params <- -rm_fit$betapar
    sim_data_list <- list(
      type = "dichotomous",
      thetas = thetas,
      item_params = item_params,
      n_items = ncol(data_mat),
      sample_n = sample_n,
      item_names = item_names_vec,
      dif_levels = dif_levels,
      dif_proportions = dif_proportions
    )
  }

  if (use_parallel) {
    results_raw <- run_partgam_sim_parallel(iterations, sim_seeds,
                                            sim_data_list, n_cores, verbose)
  } else {
    results_raw <- run_partgam_sim_sequential(iterations, sim_seeds,
                                              sim_data_list, verbose)
  }

  # Filter out failures (character strings indicate errors)
  ok <- vapply(results_raw, is.data.frame, logical(1L))
  successful <- results_raw[ok]

  if (length(successful) == 0L) {
    stop("All simulation iterations failed. Check your data.", call. = FALSE)
  }

  actual_iterations <- length(successful)

  # Combine per-iteration data.frames
  iter_dfs <- lapply(seq_along(successful), function(i) {
    df <- successful[[i]]
    df$iteration <- i
    df
  })
  results_df <- do.call(rbind, iter_dfs)
  rownames(results_df) <- NULL

  # Compute per-item cutoffs
  item_names <- unique(results_df$Item)
  item_cutoffs <- do.call(rbind, lapply(item_names, function(item) {
    sub <- results_df[results_df$Item == item, ]
    if (cutoff_method == "hdci") {
      gamma_interval <- ggdist::hdci(sub$gamma, .width = hdci_width)
      data.frame(
        Item       = item,
        gamma_low  = gamma_interval[1L, 1L],
        gamma_high = gamma_interval[1L, 2L],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    } else {
      data.frame(
        Item       = item,
        gamma_low  = stats::quantile(sub$gamma, 0.025, na.rm = TRUE),
        gamma_high = stats::quantile(sub$gamma, 0.975, na.rm = TRUE),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }))
  rownames(item_cutoffs) <- NULL

  list(
    results           = results_df,
    item_cutoffs      = item_cutoffs,
    actual_iterations = actual_iterations,
    sample_n          = sample_n,
    sample_summary    = summary(thetas),
    item_names        = item_names_vec,
    dif_group_sizes   = as.integer(dif_table),
    cutoff_method     = cutoff_method,
    hdci_width        = hdci_width
  )
}

# ---------------------------------------------------------------------------
# Internal: single simulation iteration
# ---------------------------------------------------------------------------

#' Run a single partial gamma DIF simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMpgDIFcutoff()].
#' @return A data.frame with columns `Item` and `gamma`, or a character string
#'   on failure.
#' @keywords internal
run_single_partgam_sim <- function(seed, data_list) {
  set.seed(seed)

  thetas_res <- sample(data_list$thetas, size = data_list$sample_n,
                       replace = TRUE)

  tryCatch({
    if (data_list$type == "dichotomous") {
      sim_mat <- psychotools::rrm(
        theta = thetas_res,
        beta = data_list$item_params
      )
      sim_df <- as.data.frame(sim_mat$data)
      colnames(sim_df) <- data_list$item_names

      pos_counts <- colSums(sim_df, na.rm = TRUE)
      if (any(pos_counts < 8L)) {
        return("validation_failed: fewer than 8 positive responses in at least one item")
      }
    } else {
      sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
      sim_df <- as.data.frame(sim_mat)
      colnames(sim_df) <- data_list$item_names

      n_cats <- vapply(data_list$deltaslist, function(d) length(d) + 1L,
                       integer(1L))
      for (j in seq_len(ncol(sim_df))) {
        tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
        if (any(tab == 0L)) {
          return("validation_failed: not all categories represented")
        }
      }
    }

    # Create a random DIF variable with the same group proportions
    # but no actual relationship to item responses (no true DIF)
    random_dif <- sample(
      data_list$dif_levels,
      size = data_list$sample_n,
      replace = TRUE,
      prob = data_list$dif_proportions
    )

    # Compute partial gamma DIF via iarm
    # iarm::partgam_DIF expects a data.frame and a DIF variable vector
    pgam <- iarm::partgam_DIF(sim_df, random_dif)

    # pgam is a data.frame with columns including "item" and "gamma"
    # (column names may vary slightly; use positional or clean names)
    pgam_df <- as.data.frame(pgam)

    data.frame(
      Item  = data_list$item_names,
      gamma = as.numeric(pgam_df[seq_along(data_list$item_names), "gamma"]),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }, error = function(e) {
    as.character(conditionMessage(e))
  })
}

# ---------------------------------------------------------------------------
# Internal: parallel runner
# ---------------------------------------------------------------------------

#' Run partial gamma DIF simulations in parallel using mirai
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param n_cores Number of mirai daemons.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_partgam_sim_parallel <- function(iterations, sim_seeds, sim_data_list,
                                     n_cores, verbose = FALSE) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
    completed <- 0L
  }

  # Submit all tasks
  tasks <- lapply(seq_len(iterations), function(sim) {
    mirai::mirai(
      {
        run_single_partgam_sim(seed, data_list)
      },
      seed = sim_seeds[sim],
      data_list = sim_data_list,
      run_single_partgam_sim = run_single_partgam_sim,
      sim_partial_score = sim_partial_score,
      sim_poly_item = sim_poly_item
    )
  })

  # Collect results
  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    result <- mirai::call_mirai(tasks[[sim]])$data
    if (!inherits(result, "errorValue")) {
      results[[sim]] <- result
    } else {
      results[[sim]] <- "mirai_error"
    }
    if (verbose) {
      completed <- completed + 1L
      utils::setTxtProgressBar(pb, completed)
    }
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}

# ---------------------------------------------------------------------------
# Internal: sequential runner
# ---------------------------------------------------------------------------

#' Run partial gamma DIF simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_partgam_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                       verbose = FALSE) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_partgam_sim(sim_seeds[sim], sim_data_list)
    if (verbose) {
      utils::setTxtProgressBar(pb, sim)
    }
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}


#' Plot Distribution of Simulated Partial Gamma DIF Values
#'
#' Visualises the distribution of simulation-based partial gamma DIF values
#' from \code{\link{RMpgDIFcutoff}}, optionally overlaying observed partial
#' gamma values computed from real data via \code{\link[iarm]{partgam_DIF}}.
#'
#' Uses `ggdist::stat_dotsinterval()` (when `data` is not supplied) or
#' `ggdist::stat_dots()` (when `data` is supplied) with
#' `point_interval = "median_hdci"` and `.width = c(0.66, 0.95, 0.99)`.
#'
#' @param simfit The return value of \code{\link{RMpgDIFcutoff}} (a list with
#'   components `results`, `item_cutoffs`, `actual_iterations`, `sample_n`, and
#'   `item_names`).
#' @param data Optional. A data.frame or matrix of item responses for computing
#'   and overlaying observed partial gamma values. Items must be scored starting
#'   at 0 (non-negative integers). When provided, the plot includes orange
#'   diamond markers for the observed partial gamma alongside the simulated
#'   distribution, plus segment summaries from the cutoff intervals.
#' @param dif_var Required when `data` is supplied. A vector (factor, character,
#'   or integer) defining group membership for the DIF analysis. Must have the
#'   same length as `nrow(data)`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' When `data` is **not** supplied, the function plots the simulated partial
#' gamma distributions as dot-interval plots using
#' `ggdist::stat_dotsinterval()` with median and Highest Density Continuous
#' Interval (HDCI) summaries.
#'
#' When `data` **is** supplied (along with `dif_var`), the function:
#' \enumerate{
#'   \item Computes observed partial gamma values via
#'     `iarm::partgam_DIF()`.
#'   \item Overlays observed gamma values as orange diamond markers on the
#'     simulated distributions.
#'   \item Shows per-item cutoff intervals (from `simfit$item_cutoffs`) as
#'     black line segments, with thicker segments for the 66\% interval and
#'     black dots for the median.
#' }
#'
#' The `ggplot2`, `ggdist`, and optionally `iarm` packages must be installed
#' (they are in Suggests, not Imports).
#'
#' @seealso \code{\link{RMpgDIFcutoff}}, \code{\link{RMpartgamDIF}}
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#' )
#' colnames(sim_data) <- paste0("Item", 1:10)
#' dif_group <- factor(sample(c("A", "B"), 200, replace = TRUE))
#'
#' # Run simulation
#' cutoff_res <- RMpgDIFcutoff(sim_data, dif_var = dif_group,
#'                              iterations = 100, parallel = FALSE, seed = 42)
#'
#' # Simulated distribution only
#' RMpgDIFplot(cutoff_res)
#'
#' # With observed partial gamma overlaid
#' RMpgDIFplot(cutoff_res, data = sim_data, dif_var = dif_group)
#' }
RMpgDIFplot <- function(simfit, data, dif_var) {

  # --- Check required packages ------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for RMpgDIFplot() but is not installed.\n",
      "Install it with: install.packages(\"ggplot2\")",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required for RMpgDIFplot() but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")",
      call. = FALSE
    )
  }

  # --- Validate simfit --------------------------------------------------------
  required_names <- c("results", "item_cutoffs", "actual_iterations",
                      "sample_n", "item_names")
  missing_names <- setdiff(required_names, names(simfit))
  if (length(missing_names) > 0L) {
    stop(
      "`simfit` is missing required components: ",
      paste(missing_names, collapse = ", "),
      ".\nExpected the return value of RMpgDIFcutoff().",
      call. = FALSE
    )
  }

  results_df        <- simfit$results
  item_cutoffs      <- simfit$item_cutoffs
  actual_iterations <- simfit$actual_iterations
  sample_n          <- simfit$sample_n
  item_names        <- simfit$item_names

  # Item factor levels (reversed for plotting top-to-bottom)
  item_levels <- rev(item_names)

  # --- Compute per-item summary intervals for segment overlays ----------------
  lo_hi <- do.call(rbind, lapply(item_names, function(item) {
    sub <- results_df[results_df$Item == item, ]
    data.frame(
      Item            = item,
      min_gamma       = stats::quantile(sub$gamma, 0.005, na.rm = TRUE),
      max_gamma       = stats::quantile(sub$gamma, 0.995, na.rm = TRUE),
      p66lo_gamma     = stats::quantile(sub$gamma, 0.167, na.rm = TRUE),
      p66hi_gamma     = stats::quantile(sub$gamma, 0.833, na.rm = TRUE),
      median_gamma    = stats::median(sub$gamma, na.rm = TRUE),
      stringsAsFactors = FALSE,
      row.names       = NULL
    )
  }))
  rownames(lo_hi) <- NULL

  # --- Case 1: no observed data, show simulation distribution only ------------
  if (missing(data)) {

    results_plot <- data.frame(
      Item  = results_df$Item,
      Value = results_df$gamma,
      stringsAsFactors = FALSE
    )
    results_plot$Item <- factor(results_plot$Item, levels = item_levels)

    p <- ggplot2::ggplot(
      results_plot,
      ggplot2::aes(
        x = .data$Value,
        y = .data$Item
      )
    ) +
      ggdist::stat_dotsinterval(
        ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
        quantiles = actual_iterations,
        point_interval = "median_hdci",
        layout = "weave",
        slab_color = NA,
        .width = c(0.66, 0.95, 0.99)
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "grey50",
        linewidth = 0.4
      ) +
      ggplot2::labs(
        x = "Partial gamma",
        y = "Item",
        caption = paste0(
          "Note: Results from ", actual_iterations,
          " simulated datasets with ", sample_n, " respondents (no true DIF)."
        )
      ) +
      ggplot2::scale_color_manual(
        values = scales::brewer_pal()(3)[-1],
        aesthetics = "slab_fill",
        guide = "none"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"))

    return(p)
  }

  # --- Case 2: observed data supplied -----------------------------------------
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required to compute observed partial gamma but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  if (missing(dif_var)) {
    stop(
      "`dif_var` must be supplied when `data` is provided.",
      call. = FALSE
    )
  }

  validate_response_data(data)

  if (length(dif_var) != nrow(data)) {
    stop(
      "`dif_var` must have the same length as `nrow(data)` (",
      nrow(data), "), but has length ", length(dif_var), ".",
      call. = FALSE
    )
  }

  # rgl workaround
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  sink(nullfile())
  pgam_raw <- iarm::partgam_DIF(as.data.frame(data), dif_var)
  sink()

  observed_df <- data.frame(
    Item           = as.character(pgam_raw$Item),
    observed_gamma = as.numeric(pgam_raw$gamma),
    stringsAsFactors = FALSE
  )

  # --- Build plot data --------------------------------------------------------
  gamma_sim <- data.frame(
    Item  = results_df$Item,
    Value = results_df$gamma,
    stringsAsFactors = FALSE
  )
  gamma_sim <- merge(gamma_sim, observed_df, by = "Item", sort = FALSE)
  gamma_sim$Item <- factor(gamma_sim$Item, levels = item_levels)

  lo_hi$Item_f <- factor(lo_hi$Item, levels = item_levels)

  caption_text <- paste0(
    "Note: Results from ", actual_iterations,
    " simulated datasets with ", sample_n, " respondents.\n",
    "Orange diamonds indicate observed partial gamma DIF. ",
    "Black dots indicate median gamma from simulations."
  )

  p <- ggplot2::ggplot(
    gamma_sim,
    ggplot2::aes(
      x = .data$Value,
      y = .data$Item
    )
  ) +
    ggdist::stat_dots(
      ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
      quantiles = actual_iterations,
      layout = "weave",
      slab_color = NA,
      .width = c(0.66, 0.95, 0.99)
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x    = .data$min_gamma,
        xend = .data$max_gamma,
        y    = .data$Item_f,
        yend = .data$Item_f
      ),
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x    = .data$p66lo_gamma,
        xend = .data$p66hi_gamma,
        y    = .data$Item_f,
        yend = .data$Item_f
      ),
      color = "black",
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = lo_hi,
      ggplot2::aes(
        x = .data$median_gamma,
        y = .data$Item_f
      ),
      size = 3.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$observed_gamma),
      color = "sienna2",
      shape = 18,
      position = ggplot2::position_nudge(y = -0.1),
      size = 4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.4
    ) +
    ggplot2::labs(
      x = "Partial gamma",
      y = "Item",
      caption = caption_text
    ) +
    ggplot2::scale_color_manual(
      values = scales::brewer_pal()(3)[-1],
      aesthetics = "slab_fill",
      guide = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"))

  p
}
