#' PCA of Standardized Rasch Residuals
#'
#' Fits a Rasch model by CML via `psychotools` (a dichotomous item is a
#' 2-category PCM), extracts the standardized residuals \eqn{(x - E)/\sqrt{Var}}
#' at WLE person locations, and runs an unrotated principal-component
#' analysis on those residuals via `stats::prcomp()`. The function reports
#' the top `n_components` eigenvalues and their proportions of unexplained
#' variance, and optionally compares the first-contrast eigenvalue against
#' a simulation-based bound from \code{\link{RMdimResidualPCACutoff}}.
#'
#' Rule-of-thumb thresholds for the first-contrast eigenvalue (e.g., the
#' "> 2" heuristic occasionally cited from Winsteps documentation) are not
#' reliable indicators of multidimensionality; the first-contrast eigenvalue
#' under a correctly fitting unidimensional model varies systematically with
#' sample size, test length, and item-parameter spread. Empirical (simulated)
#' bounds tailored to the data structure should be used instead — see
#' \code{\link{RMdimResidualPCACutoff}}, and Chou & Wang (2010) for the underlying
#' simulation argument.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Rows with any `NA` are dropped
#'   before PCA, since `prcomp()` does not accept missing values.
#' @param cutoff Optional. The list returned by \code{\link{RMdimResidualPCACutoff}} (its
#'   `suggested_cutoff` is used), or a single numeric value to use as the
#'   cutoff directly. When provided, the result includes a `Flagged` column
#'   (logical: is the eigenvalue above the simulated bound?) and the kable
#'   caption notes the cutoff.
#' @param n_components Integer. Number of eigenvalues to report. Capped at the
#'   number of items. Default `5`.
#' @param output Character. `"kable"` (default) for a formatted
#'   `knitr::kable()` table, `"dataframe"` for the underlying data.frame, or
#'   `"loadings"` for a ggplot of PC1 loadings against item locations.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object with columns Component,
#'   Eigenvalue, Proportion of variance (and `Flagged` when `cutoff` is
#'   provided). The caption gives the variance partition (% of total observed
#'   variance explained by measures vs. unexplained), the model fitted,
#'   sample size, and cutoff metadata if applicable.
#' * If `output = "dataframe"`: a data.frame with columns `Component`,
#'   `Eigenvalue`, `Proportion_of_variance` (and `Flagged` when `cutoff` is
#'   provided). The variance partition is attached as the
#'   `"variance_partition"` attribute — a list with elements `total`,
#'   `explained`, `unexplained`, `pct_explained`, `pct_unexplained`,
#'   `n_persons`. Access via
#'   `attr(result, "variance_partition")`.
#' * If `output = "loadings"`: a ggplot showing each item's PC1 loading on
#'   the x-axis and Rasch item location on the y-axis, with dashed reference
#'   lines at zero, and the variance partition in the figure caption. Item
#'   names are labelled via `ggrepel::geom_text_repel()` when `ggrepel` is
#'   installed; otherwise plain `geom_text()`.
#'
#' @details
#' The PCA is performed on the standardized residuals
#' \eqn{(x - E)/\sqrt{Var}} from the shared CML/WLE engine (CML item
#' parameters via `psychotools`, WLE person locations). The reported
#' eigenvalues are unrotated; rotation is appropriate for *interpreting* a
#' multidimensional solution but obscures the dominant first contrast that
#' dimensionality assessment is concerned with.
#'
#' Item locations on the loadings plot are the per-item mean of the CML
#' Andrich thresholds.
#'
#' The variance partition follows Linacre's convention: per-item observed
#' variance is compared to per-item *expected* variance under the fitted
#' model, summed across items. Expected scores are computed from the CML item
#' parameters and WLE person locations. WLE is finite at extreme scores, so
#' all persons are retained (the previous MLE partition dropped extreme-score
#' cases).
#'
#' @references
#' Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item
#' response models with principal component analysis on standardized
#' residuals. *Educational and Psychological Measurement, 70*(5), 717-731.
#' \doi{10.1177/0013164410379322}
#'
#' @seealso \code{\link{RMdimResidualPCACutoff}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- as.data.frame(
#'   matrix(sample(0:1, 200 * 12, replace = TRUE), nrow = 200, ncol = 12)
#' )
#' colnames(dat) <- paste0("I", 1:12)
#'
#' # Default kable output
#' RMdimResidualPCA(dat)
#'
#' # PC1 loadings vs item location plot
#' RMdimResidualPCA(dat, output = "loadings")
#'
#' # Simulation-based cutoff (use 250+ iterations in real analyses)
#' bound <- RMdimResidualPCACutoff(dat, iterations = 50, parallel = FALSE, seed = 1)
#' RMdimResidualPCA(dat, cutoff = bound)
#' }
RMdimResidualPCA <- function(data,
                          cutoff       = NULL,
                          n_components = 5L,
                          output       = "kable") {

  output <- match.arg(output, c("kable", "dataframe", "loadings"))

  validate_response_data(data)

  # Drop incomplete rows (prcomp can't handle NA in standardized residuals)
  data <- stats::na.omit(as.data.frame(data))
  if (nrow(data) == 0L) {
    stop("No complete cases in `data` after na.omit.", call. = FALSE)
  }

  # --- Normalise the cutoff argument -----------------------------------------
  cutoff_value           <- NULL
  cutoff_n_iter          <- NULL
  cutoff_percentile_pct  <- NULL
  if (!is.null(cutoff)) {
    if (is.list(cutoff) && !is.data.frame(cutoff) &&
        "suggested_cutoff" %in% names(cutoff)) {
      cutoff_value           <- cutoff$suggested_cutoff
      cutoff_n_iter          <- cutoff$actual_iterations
      cutoff_percentile_pct  <- cutoff$suggested_cutoff_percentile
    } else if (is.numeric(cutoff) && length(cutoff) == 1L) {
      cutoff_value <- as.numeric(cutoff)
    } else {
      stop("`cutoff` must be NULL, the return value of RMdimResidualPCACutoff(), or a ",
           "single numeric value.", call. = FALSE)
    }
  }

  # rgl workaround
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # --- Fit Rasch model and get residuals -------------------------------------
  data_mat      <- as.matrix(data)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  # CML item thresholds (psychotools) and WLE-based standardized residuals,
  # consistent with the rest of the package. Item locations are per-item mean
  # thresholds; the standardized residuals (x - E)/sqrt(Var) come from the
  # shared CML/WLE engine. (Previously eRm CML + eRm::itemfit() MLE residuals.)
  thr_list       <- .fit_cml_thresholds(data)
  item_locations <- vapply(thr_list, mean, numeric(1L))
  names(item_locations) <- names(data)

  st_resids <- .rasch_std_residuals(data, method = "WLE")
  if (anyNA(st_resids)) {
    keep_rows <- stats::complete.cases(st_resids)
    st_resids <- st_resids[keep_rows, , drop = FALSE]
  }

  # --- Variance partition (Linacre-style; CML items, WLE thetas) -------------
  # Compares Var(observed) to Var(model-expected) per item, summed across
  # items. WLE locations are finite at extreme scores, so all persons are
  # retained (unlike the MLE partition, which dropped extreme-score cases).
  thetas_all    <- .estimate_thetas(data_mat, thr_list, method = "WLE")$theta
  finite_thetas <- is.finite(thetas_all)

  if (sum(finite_thetas) >= 2L) {
    th <- thetas_all[finite_thetas]
    expected_mat <- vapply(seq_along(thr_list), function(i) {
      cats <- 0:length(thr_list[[i]])
      P    <- vapply(th, function(z) .pcm_cat_probs(z, thr_list[[i]]),
                     numeric(length(cats)))
      as.numeric(crossprod(cats, P))                 # E[score] per person
    }, numeric(length(th)))
    data_finite     <- data_mat[finite_thetas, , drop = FALSE]
    var_total       <- sum(apply(data_finite,  2L, stats::var, na.rm = TRUE))
    var_explained   <- sum(apply(expected_mat, 2L, stats::var, na.rm = TRUE))
    var_unexplained <- max(var_total - var_explained, 0)
    pct_explained   <- if (var_total > 0) var_explained / var_total else NA_real_
    pct_unexplained <- if (var_total > 0) var_unexplained / var_total else NA_real_
    n_partition     <- sum(finite_thetas)
    partition_avail <- TRUE
  } else {
    var_total <- var_explained <- var_unexplained <- NA_real_
    pct_explained <- pct_unexplained <- NA_real_
    n_partition <- 0L
    partition_avail <- FALSE
  }

  partition_text <- if (partition_avail) {
    paste0(
      "Total observed variance: ",
      round(pct_explained * 100, 1), "% explained by measures, ",
      round(pct_unexplained * 100, 1),
      "% unexplained\n(basis for PCA; n = ", n_partition,
      " non-extreme cases)."
    )
  } else {
    "Variance partition unavailable (too few persons with finite theta MLEs)."
  }

  # --- Run unrotated PCA -----------------------------------------------------
  pca_fit  <- stats::prcomp(st_resids)
  eigvals  <- pca_fit$sdev^2
  total    <- sum(eigvals)
  prop_var <- eigvals / total

  k_show <- min(as.integer(n_components), length(eigvals))

  result_df <- data.frame(
    Component              = paste0("PC", seq_len(k_show)),
    Eigenvalue             = round(eigvals[seq_len(k_show)], 3),
    Proportion_of_variance = round(prop_var[seq_len(k_show)], 3),
    stringsAsFactors       = FALSE,
    row.names              = NULL
  )

  if (!is.null(cutoff_value)) {
    result_df$Flagged <- result_df$Eigenvalue > cutoff_value
  }

  # Variance partition is metadata on the dataframe rather than a separate
  # list element. Access via attr(result, "variance_partition").
  attr(result_df, "variance_partition") <- list(
    total           = var_total,
    explained       = var_explained,
    unexplained     = var_unexplained,
    pct_explained   = pct_explained,
    pct_unexplained = pct_unexplained,
    n_persons       = n_partition
  )

  # --- Loadings plot ---------------------------------------------------------
  if (output == "loadings") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for output = \"loadings\".",
           call. = FALSE)
    }
    loadings <- as.data.frame(pca_fit$rotation)
    loadings$Item     <- rownames(loadings)
    loadings$Location <- as.numeric(item_locations[loadings$Item])

    p <- ggplot2::ggplot(
      loadings,
      ggplot2::aes(x = .data$PC1, y = .data$Location, label = .data$Item)
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::scale_x_continuous(limits = c(-1, 1))

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(size = 3.5, max.overlaps = Inf)
    } else {
      p <- p + ggplot2::geom_text(nudge_y = 0.05, size = 3.5)
    }

    p <- p +
      ggplot2::labs(
        x       = "Loading on first residual contrast (PC1)",
        y       = "Item location (logit scale)",
        caption = er2_caption(partition_text)
      ) +
      ggplot2::theme_bw(base_size = 13) +
      er2_axis_margins() +
      er2_plot_caption()

    return(p)
  }

  if (output == "dataframe") {
    return(result_df)
  }

  # --- kable -----------------------------------------------------------------
  caption_parts <- c(
    paste0(
      if (is_polytomous) "Partial Credit Model" else "Rasch model",
      " (", nrow(data), " complete cases, ", ncol(data), " items)."
    ),
    partition_text
  )
  if (!is.null(cutoff_value)) {
    iter_part <- if (!is.null(cutoff_n_iter)) {
      paste0(" based on ", cutoff_n_iter, " simulation iterations")
    } else {
      ""
    }
    pct_part <- if (!is.null(cutoff_percentile_pct)) {
      paste0(" (", cutoff_percentile_pct, "th percentile)")
    } else {
      ""
    }
    caption_parts <- c(
      caption_parts,
      paste0(
        "First-contrast cutoff = ", round(cutoff_value, 3),
        iter_part, pct_part, "."
      )
    )
  }

  knitr::kable(
    result_df,
    format = "pipe",
    col.names = c("Component", "Eigenvalue", "Proportion of variance",
                  if (!is.null(cutoff_value)) "Flagged"),
    caption = paste(caption_parts, collapse = " ")
  )
}

# ===========================================================================

#' Simulation-based Cutoff for First-Contrast Eigenvalue
#'
#' Parametric bootstrap producing an empirical upper percentile for the
#' largest eigenvalue from a PCA of standardized residuals under a correctly
#' fitting Rasch model. For each iteration the function resamples theta
#' values from the data-based estimates, simulates response data, refits the
#' appropriate Rasch model, and extracts the first-contrast eigenvalue.
#' Several upper-tail percentiles of the resulting distribution are returned;
#' the 99th percentile is reported as the suggested cutoff (matching the
#' convention used by \code{\link{RMlocdepQ3Cutoff}}).
#'
#' Rule-of-thumb cutoffs for the first-contrast eigenvalue depend strongly on
#' sample size, test length, and item-parameter spread; simulation-based
#' cutoffs tailored to the data are more defensible (Chou & Wang, 2010).
#'
#' @param data A data.frame or matrix of item responses (0-based,
#'   non-negative integers).
#' @param iterations Integer. Number of simulation iterations. Default `250`.
#' @param parallel Logical. Use parallel processing via `mirai`. Default
#'   `TRUE`.
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first; if neither is set,
#'   `parallel = TRUE` falls back to sequential with a warning.
#' @param verbose Logical. Show a progress bar. Default `FALSE`.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{`results`}{data.frame: `iteration`, `eigenvalue`.}
#'   \item{`p95`, `p99`, `p995`, `p999`}{Empirical percentiles of
#'     `eigenvalue`.}
#'   \item{`max`}{The largest simulated eigenvalue.}
#'   \item{`suggested_cutoff`}{The 99th percentile (`p99`) — pass this list
#'     back into \code{\link{RMdimResidualPCA}} via `cutoff = `, or use
#'     `suggested_cutoff` directly.}
#'   \item{`suggested_cutoff_percentile`}{The percentile used for
#'     `suggested_cutoff`, currently always `99`.}
#'   \item{`actual_iterations`}{Number of successful iterations.}
#'   \item{`sample_n`}{Number of complete cases used.}
#'   \item{`item_names`}{Character vector of item names.}
#' }
#'
#' @details
#' The generating model uses CML item parameters (`psychotools`) and WLE person
#' locations. Per iteration: theta values are sampled with replacement from the
#' WLE estimates; response data are simulated under the model
#' (`psychotools::rrm` for dichotomous, partial-credit simulator for
#' polytomous); CML/WLE standardized residuals are computed (the same engine as
#' the observed analysis); `prcomp()` is run; the first-contrast eigenvalue is
#' recorded.
#'
#' Iterations that fail (e.g., due to a degenerate simulated dataset where
#' some category isn't represented) are silently dropped. The `iarm` package
#' is **not** required.
#'
#' @references
#' Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item
#' response models with principal component analysis on standardized
#' residuals. *Educational and Psychological Measurement, 70*(5), 717-731.
#' \doi{10.1177/0013164410379322}
#'
#' @seealso \code{\link{RMdimResidualPCA}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- as.data.frame(
#'   matrix(sample(0:1, 200 * 12, replace = TRUE), nrow = 200, ncol = 12)
#' )
#' colnames(dat) <- paste0("I", 1:12)
#'
#' # Few iterations for a fast example; use 250+ in real analyses
#' bound <- RMdimResidualPCACutoff(dat, iterations = 50, parallel = FALSE, seed = 1)
#' bound$suggested_cutoff
#'
#' RMdimResidualPCA(dat, cutoff = bound)
#' }
RMdimResidualPCACutoff <- function(data,
                        iterations = 250,
                        parallel   = TRUE,
                        n_cores    = NULL,
                        verbose    = FALSE,
                        seed       = NULL) {

  validate_response_data(data)

  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  data <- stats::na.omit(as.data.frame(data))
  if (nrow(data) == 0L) {
    stop("No complete cases in `data`.", call. = FALSE)
  }

  use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)
  if (parallel && !use_parallel) {
    message("Install 'mirai' for parallel processing: install.packages(\"mirai\")")
    message("Running sequentially...")
  }
  if (use_parallel) {
    if (is.null(n_cores)) n_cores <- getOption("mc.cores")
    if (is.null(n_cores)) {
      warning(
        "For parallel processing, specify n_cores or set options(mc.cores = N).\n",
        "Falling back to sequential.",
        call. = FALSE
      )
      use_parallel <- FALSE
    } else {
      n_cores <- min(n_cores, iterations)
    }
  }

  if (!is.null(seed)) set.seed(seed)
  sim_seeds <- sample.int(.Machine$integer.max, iterations)

  data_mat       <- as.matrix(data)
  sample_n       <- nrow(data_mat)
  is_polytomous  <- max(data_mat, na.rm = TRUE) > 1L
  item_names_vec <- colnames(data_mat)

  # Generating model: CML item thresholds (psychotools) + WLE person locations,
  # consistent with the rest of the package. The DGP is unchanged: thetas are
  # resampled with replacement and data simulated parametrically.
  thr_list   <- .fit_cml_thresholds(data_mat)
  wle_thetas <- .estimate_thetas(data_mat, thr_list, method = "WLE")$theta
  wle_thetas <- wle_thetas[is.finite(wle_thetas)]

  sim_data_list <- list(
    type       = if (is_polytomous) "polytomous" else "dichotomous",
    thetas     = wle_thetas,
    n_items    = ncol(data_mat),
    sample_n   = sample_n,
    item_names = item_names_vec
  )
  if (is_polytomous) {
    sim_data_list$deltaslist <- thr_list
  } else {
    sim_data_list$item_params <- unlist(thr_list, use.names = FALSE)
  }

  if (use_parallel) {
    results_raw <- run_pca_sim_parallel(iterations, sim_seeds, sim_data_list,
                                        n_cores, verbose)
  } else {
    results_raw <- run_pca_sim_sequential(iterations, sim_seeds, sim_data_list,
                                          verbose)
  }

  ok <- vapply(results_raw, is.numeric, logical(1L))
  successful <- results_raw[ok]

  if (length(successful) == 0L) {
    failed_msgs <- unlist(results_raw[!ok])
    sample_msg <- if (length(failed_msgs) > 0L) {
      unique(failed_msgs)[1L]
    } else {
      "(no message captured)"
    }
    stop("All simulation iterations failed. Example: ", sample_msg,
         call. = FALSE)
  }

  actual_iterations <- length(successful)
  eig_vec <- as.numeric(unlist(successful))

  results_df <- data.frame(
    iteration  = seq_len(actual_iterations),
    eigenvalue = eig_vec,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  p95  <- as.numeric(stats::quantile(eig_vec, 0.95,  na.rm = TRUE))
  p99  <- as.numeric(stats::quantile(eig_vec, 0.99,  na.rm = TRUE))
  p995 <- as.numeric(stats::quantile(eig_vec, 0.995, na.rm = TRUE))
  p999 <- as.numeric(stats::quantile(eig_vec, 0.999, na.rm = TRUE))

  list(
    results                     = results_df,
    p95                         = p95,
    p99                         = p99,
    p995                        = p995,
    p999                        = p999,
    max                         = max(eig_vec, na.rm = TRUE),
    suggested_cutoff            = p99,
    suggested_cutoff_percentile = 99,
    actual_iterations           = actual_iterations,
    sample_n                    = sample_n,
    item_names                  = item_names_vec
  )
}

# ===========================================================================
# Internal: single iteration
# ===========================================================================

#' Run a single PCA-eigenvalue simulation iteration
#'
#' @keywords internal
#' @noRd
run_single_pca_sim <- function(seed, data_list) {
  set.seed(seed)
  thetas_res <- sample(data_list$thetas, size = data_list$sample_n,
                       replace = TRUE)

  tryCatch({
    if (data_list$type == "dichotomous") {
      sim_mat <- psychotools::rrm(theta = thetas_res,
                                  beta  = data_list$item_params)
      sim_df <- as.data.frame(sim_mat$data)
      colnames(sim_df) <- data_list$item_names

      pos_counts <- colSums(sim_df, na.rm = TRUE)
      if (any(pos_counts < 8L)) {
        return("validation_failed: fewer than 8 positive responses in at least one item")
      }
      neg_counts <- nrow(sim_df) - pos_counts
      if (any(neg_counts < 8L)) {
        return("validation_failed: fewer than 8 negative responses in at least one item")
      }
    } else {
      sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
      sim_df  <- as.data.frame(sim_mat)
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

    # CML/WLE standardized residuals (same engine as the observed analysis).
    st_resids <- .rasch_std_residuals(sim_df, method = "WLE")

    if (anyNA(st_resids)) {
      keep <- stats::complete.cases(st_resids)
      st_resids <- st_resids[keep, , drop = FALSE]
    }
    if (nrow(st_resids) < ncol(st_resids)) {
      return("validation_failed: too few rows in residual matrix for PCA")
    }

    pca_fit <- stats::prcomp(st_resids)
    as.numeric(pca_fit$sdev[1L]^2)
  }, error = function(e) as.character(conditionMessage(e)))
}

# ===========================================================================
# Internal: parallel runner
# ===========================================================================

#' @keywords internal
#' @noRd
run_pca_sim_parallel <- function(iterations, sim_seeds, sim_data_list,
                                 n_cores, verbose = FALSE) {
  mirai::daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)

  if (verbose) {
    message(sprintf("Starting %d daemons...", n_cores))
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  tasks <- lapply(seq_len(iterations), function(i) {
    mirai::mirai(
      { run_single_pca_sim(seed, data_list) },
      seed                = sim_seeds[i],
      data_list           = sim_data_list,
      run_single_pca_sim  = run_single_pca_sim,
      sim_partial_score   = sim_partial_score,
      sim_poly_item       = sim_poly_item,
      # CML/WLE residual engine used by .rasch_std_residuals() in the daemon.
      .rasch_std_residuals = .rasch_std_residuals,
      .fit_cml_thresholds  = .fit_cml_thresholds,
      .estimate_thetas     = .estimate_thetas,
      .theta_wle           = .theta_wle,
      .pcm_cat_probs       = .pcm_cat_probs,
      .center_thresholds   = .center_thresholds
    )
  })

  results <- vector("list", iterations)
  for (i in seq_len(iterations)) {
    res <- mirai::call_mirai(tasks[[i]])$data
    results[[i]] <- if (inherits(res, "errorValue")) "mirai_error" else res
    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}

# ===========================================================================
# Internal: sequential runner
# ===========================================================================

# ===========================================================================
# Internal: PCM expected-score matrix (for the variance partition)
# ===========================================================================

#' Expected scores per (person, item) under a Partial Credit Model
#'
#' Computes \deqn{E[X_ij | theta_i, tau_j]} under the PCM:
#' \deqn{P(X = x | theta, tau) = exp(x*theta - sum_{l=1}^x tau_l) /
#'                              sum_{m=0}^M exp(m*theta - sum_{l=1}^m tau_l)}
#' Loops over items but is vectorised across persons within an item.
#'
#' @param thetas Numeric vector of person locations (length n).
#' @param thresh_mat Numeric matrix of thresholds (n_items rows; columns are
#'   threshold positions, NA-padded for items with fewer thresholds).
#'
#' @return Numeric matrix (n persons x n_items) of expected scores.
#'
#' @keywords internal
#' @noRd
.pcm_expected_scores <- function(thetas, thresh_mat) {
  n       <- length(thetas)
  n_items <- nrow(thresh_mat)
  expected <- matrix(NA_real_, nrow = n, ncol = n_items)

  for (j in seq_len(n_items)) {
    taus <- thresh_mat[j, !is.na(thresh_mat[j, ])]
    K_j  <- length(taus)
    cumsum_taus <- c(0, cumsum(taus))      # length K_j + 1
    log_num <- outer(thetas, 0:K_j) -
      matrix(cumsum_taus, nrow = n, ncol = K_j + 1L, byrow = TRUE)
    # Numerical-stability shift before exponentiating
    log_num <- log_num - apply(log_num, 1L, max)
    probs <- exp(log_num)
    probs <- probs / rowSums(probs)
    expected[, j] <- as.numeric(probs %*% (0:K_j))
  }
  expected
}

# ===========================================================================

#' @keywords internal
#' @noRd
run_pca_sim_sequential <- function(iterations, sim_seeds, sim_data_list,
                                   verbose = FALSE) {
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)

  results <- vector("list", iterations)
  for (i in seq_len(iterations)) {
    results[[i]] <- run_single_pca_sim(sim_seeds[i], sim_data_list)
    if (verbose) utils::setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("")
  }

  results
}

