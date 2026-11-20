#' Yen's Q3 residual correlations
#'
#' Computes Yen's \eqn{Q_3}, the correlations between the model residuals of
#' each pair of items, to assess local dependence. By default the item
#' parameters are estimated by conditional maximum likelihood and the person
#' locations by WLE (`estimator = "CML"`). With the result of
#' [RMlocdepQ3Cutoff()], it also tests each pair against simulations, by
#' default with a bootstrap p-value corrected across pairs.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Missing values (`NA`) are allowed. \eqn{Q_3} is computed from the
#'   respondents who answered both items of a pair, and respondents with no
#'   responses are dropped.
#' @param cutoff Optional. Either of:
#'   * `NULL` (default): the \eqn{Q_3} matrix.
#'   * A single number: the matrix with a global dynamic cutoff, the mean
#'     off-diagonal \eqn{Q_3} plus this number (for example
#'     `suggested_cutoff` from [RMlocdepQ3Cutoff()]).
#'   * The full result of [RMlocdepQ3Cutoff()]: a list of two tables, the
#'     matrix with the global cutoff and a table per pair (see Value).
#' @param output Character. `"kable"` (default) for `knitr::kable()`
#'   tables, or `"dataframe"`.
#' @param n_pairs Integer or `NULL` (default). With the full cutoff object,
#'   shows only the `n_pairs` pairs that depart most from their expected
#'   range in the pair table. `NULL` shows all pairs. The correction always
#'   covers all pairs.
#' @param p_value Logical or `NULL`. Whether to flag pairs on one-sided
#'   bootstrap p-values (`p_q3`, `padj_q3`) for excess positive dependence.
#'   `NULL` (default) uses them whenever `cutoff` is the full
#'   [RMlocdepQ3Cutoff()] result, so a numeric cutoff or none keeps the
#'   interval. `FALSE` gives the behaviour before version 1.2.0. The
#'   interval describes where a pair's \eqn{Q_3} is expected to fall. Used
#'   as a decision rule over m pairs, a width `w` gives a family-wise error
#'   rate of `1 - w^m`, and the number of pairs grows with the square of
#'   the number of items (Johansson, 2026).
#' @param correction Character. Correction of the p-values across all pairs:
#'   `"fwer"` (default, Westfall-Young step-down), `"fdr_bh"` or `"fdr_by"`
#'   (Benjamini-Hochberg or Benjamini-Yekutieli), or `"none"`.
#' @param alpha Numeric in (0, 1). A pair is flagged when its corrected
#'   p-value is below `alpha`. Default `0.05`. Used only with p-values.
#' @param estimator Character. How the residuals are computed. `"CML"`
#'   (default) uses CML item parameters (`psychotools`) and WLE person
#'   locations, which are finite at extreme scores. `"MML"` uses marginal
#'   maximum likelihood and EAP (`mirt`), kept for older analyses. It must
#'   match the estimator of [RMlocdepQ3Cutoff()]. With the full cutoff
#'   object, its stored estimator is used, with a warning if `estimator`
#'   differs.
#'
#' @return
#' With `cutoff = NULL` or a number, one table (`knitr_kable` or
#' data.frame) with the lower triangle of the \eqn{Q_3} matrix. A number
#' adds a flag for pairs above the dynamic cutoff and a caption describing
#' it.
#'
#' With the full [RMlocdepQ3Cutoff()] result, a named list:
#' \describe{
#'   \item{`matrix`}{the \eqn{Q_3} matrix with the global dynamic cutoff
#'     (mean off-diagonal \eqn{Q_3} plus `suggested_cutoff`), as above.}
#'   \item{`pairs`}{one row per item pair with `Item1`, `Item2`, `Observed`
#'     (\eqn{Q_3}), `Low` and `High` (the pair's simulated interval) and
#'     `Flagged`: `"above"` (above the interval, local dependence),
#'     `"below"` or `""`. Sorted by the departure from the pair's simulated
#'     median and cut to `n_pairs`. With p-values, also `p_q3` and
#'     `padj_q3`, and `Flagged` marks `padj_q3 < alpha`, only as
#'     `"above"`.}
#' }
#' The \eqn{Q_3} heatmap is drawn by [RMlocdepQ3Plot()] (its `matrix`
#' element), so the tables and plots share the `matrix` and `pairs`
#' structure.
#'
#' @details
#' \strong{Statistic.} \eqn{Q_3} (Yen, 1984) is the correlation between the
#' residuals of two items after the latent variable is accounted for. Under
#' local independence it is around \eqn{-1/(k-1)} for k items. With a
#' `cutoff`, the dynamic cutoff is the mean of all off-diagonal \eqn{Q_3}
#' plus `cutoff` (Christensen et al., 2017).
#'
#' \strong{Estimation.} \eqn{Q_3} is the correlation matrix of the
#' standardized residuals \eqn{(x - E)/\sqrt{Var}}. With `"MML"` the residuals
#' come from `mirt`. The two estimators give very similar \eqn{Q_3} (the
#' off-diagonal values typically correlate above .95). What matters for
#' inference is that the observed \eqn{Q_3} and the simulations in
#' [RMlocdepQ3Cutoff()] use the same estimator, which the functions enforce.
#'
#' \strong{Two views.} With the full [RMlocdepQ3Cutoff()] result, two tables
#' are returned. `matrix` applies one global cutoff, the 99th percentile of
#' the simulated maximum minus mean \eqn{Q_3} (Christensen et al., 2017),
#' an overview of whether there is any local dependence. `pairs` compares
#' each observed \eqn{Q_3} with its own simulated interval, so that single
#' dependent pairs can be read off and ranked.
#'
#' \strong{Bootstrap p-values.} Each observed \eqn{Q_3} is tested against
#' its simulated values (`cutoff$pair_results`), one-sided for excess local
#' dependence. The statistic is the deviation from the mean of the observed
#' and simulated values together, divided by their SD. The uncorrected p is
#' `(1 + #{Q3* >= Q3}) / (B + 1)` for B simulated datasets, and
#' `correction` adjusts across the \eqn{k(k-1)/2} pairs (Westfall-Young
#' step-down, Ferreira, 2024, or FDR). As for item fit, the family-wise
#' error rate is controlled at any number of iterations, but decisions
#' depend on the seed. Below 400 iterations a message says so once per
#' session, and between 400 and 1000 the caption notes it. A false
#' discovery rate correction needs far more iterations than the
#' family-wise one, since a bootstrap p-value cannot fall below
#' `1/(B + 1)` and Benjamini-Hochberg compares the smallest with
#' `alpha/m`. Over the 36 pairs of a nine-item scale that takes 720
#' iterations, and over the 190 pairs of a twenty-item scale 3799, against
#' 19 for Westfall-Young whatever the number of pairs.
#'
#' @references
#' Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
#' for Yen's \eqn{Q_3}: Identification of local dependence in the Rasch model.
#' *Applied Psychological Measurement, 41*(3), 178--194.
#' \doi{10.1177/0146621616677520}
#'
#' Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate' number
#' of hypotheses simultaneously. *Journal of Statistical Theory and Practice,
#' 19*(1), 6. \doi{10.1007/s42519-024-00412-4}
#'
#' Johansson, M. (2026). Simulation-based cutoffs for conditional item fit in
#' Rasch models: Iterations, multiplicity correction, and decision stability.
#' *PsyArXiv*. \doi{10.31234/osf.io/7pqz4_v2}
#'
#' Yen, W. M. (1984). Effects of local item dependence on the fit and
#' equating performance of the three-parameter logistic model.
#' *Applied Psychological Measurement, 8*(2), 125--145.
#' \doi{10.1177/014662168400800201}
#'
#' @inheritSection RMitemInfit Multiple comparisons
#'
#' @seealso [RMlocdepQ3Cutoff()] for the simulations, [RMlocdepQ3Plot()]
#'   for the plots, [RMlocdepGamma()] for partial gamma.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Simulate binary item response data (10 items, 200 persons)
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#' )
#' colnames(sim_data) <- paste0("Item", 1:10)
#'
#' # Raw Q3 matrix (no cutoff)
#' RMlocdepQ3(sim_data)
#'
#' # Get the underlying data.frame
#' q3_df <- RMlocdepQ3(sim_data, output = "dataframe")
#'
#' # Simulation-based cutoff (use 1000 or more iterations in real analyses)
#' if (requireNamespace("ggdist", quietly = TRUE)) {
#'   cutoff_res <- RMlocdepQ3Cutoff(sim_data, iterations = 50, parallel = FALSE)
#'
#'   # A single number gives the matrix only
#'   RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
#'
#'   # The full object gives two tables, $matrix and $pairs
#'   res <- RMlocdepQ3(sim_data, cutoff = cutoff_res, output = "dataframe")
#'   res$pairs
#'
#'   # Top 5 pairs, with bootstrap p-values (use iterations >= 1000 in practice)
#'   RMlocdepQ3(sim_data, cutoff = cutoff_res, n_pairs = 5, p_value = TRUE,
#'              output = "dataframe")$pairs
#' }
#' }
RMlocdepQ3 <- function(
  data,
  cutoff = NULL,
  output = "kable",
  n_pairs = NULL,
  p_value = NULL,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  estimator = c("CML", "MML")
) {
  # --- Input validation -------------------------------------------------------
  cutoff_full <- NULL # full object (carries simulated $pair_results)
  if (!is.null(cutoff)) {
    # Accept the full RMlocdepQ3Cutoff() return list as well as a bare numeric.
    if (
      is.list(cutoff) &&
        !is.data.frame(cutoff) &&
        "suggested_cutoff" %in% names(cutoff)
    ) {
      cutoff_full <- cutoff
      cutoff <- as.numeric(cutoff$suggested_cutoff)
    }
    if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
      stop(
        "`cutoff` must be a single numeric value, NULL, or the list ",
        "returned by RMlocdepQ3Cutoff().",
        call. = FALSE
      )
    }
  }

  output <- match.arg(output, c("kable", "dataframe"))
  correction <- match.arg(correction)
  estimator <- match.arg(estimator)
  # The observed Q3 must use the same estimator as the simulated cut-off, so
  # the estimator stored in the cutoff object takes precedence over the arg.
  if (
    !is.null(cutoff_full) &&
      !is.null(cutoff_full$estimator) &&
      !identical(cutoff_full$estimator, estimator)
  ) {
    warning(
      "Using estimator \"",
      cutoff_full$estimator,
      "\" from the ",
      "RMlocdepQ3Cutoff() object (overriding estimator = \"",
      estimator,
      "\") so the observed Q3 matches the simulated cut-off.",
      call. = FALSE
    )
    estimator <- cutoff_full$estimator
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  if (
    !is.null(n_pairs) &&
      (!is.numeric(n_pairs) || length(n_pairs) != 1L || n_pairs < 1)
  ) {
    stop("`n_pairs` must be NULL or a single positive integer.", call. = FALSE)
  }

  # --- Resolve p_value --------------------------------------------------------
  # NULL means "use the corrected p-value when the simulations are available".
  # The interval describes where a pair's Q3 is expected to fall and makes a
  # poor decision rule, because its width sets a family-wise error rate over
  # every pair at once, and pairs grow quadratically in items (Johansson,
  # 2026). A numeric cutoff, the bare $pair_cutoffs, or no cutoff at all
  # leaves nothing to compute a p-value from, so those resolve to FALSE.
  if (!is.null(p_value) && (!is.logical(p_value) || length(p_value) != 1L)) {
    stop("`p_value` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }
  have_sims <- !is.null(cutoff_full) && !is.null(cutoff_full$pair_results)
  if (is.null(p_value)) {
    p_value <- have_sims
  }
  n_pairs_total <- if (!is.null(cutoff_full$pair_cutoffs)) {
    nrow(cutoff_full$pair_cutoffs)
  } else {
    NULL
  }

  if (p_value) {
    if (!have_sims) {
      stop(
        "`p_value = TRUE` requires the full RMlocdepQ3Cutoff() object (it ",
        "carries the simulated per-pair distributions in $pair_results); a ",
        "numeric cutoff or NULL is not sufficient.",
        call. = FALSE
      )
    }
    # Below 400 decisions depend on the seed enough to say so on the
    # console. Between 400 and 1000 the table caption reports it instead.
    if (
      !is.null(cutoff_full$actual_iterations) &&
        cutoff_full$actual_iterations < 400L
    ) {
      .notify_low_iterations(
        cutoff_full$actual_iterations,
        cutoff_full$requested_iterations,
        fn = "RMlocdepQ3Cutoff()",
        id = "easyRasch2_low_iterations_locdep",
        correction = correction
      )
    }
    if (!is.null(n_pairs_total)) {
      .warn_fdr_floor(
        cutoff_full$actual_iterations,
        n_pairs_total,
        correction,
        alpha,
        unit = "item pairs",
        fn = "RMlocdepQ3Cutoff()"
      )
    }
  } else if (!is.null(cutoff_full) && !is.null(n_pairs_total)) {
    .notify_band_flagging(
      if (identical(cutoff_full$cutoff_method, "quantile")) {
        0.95
      } else {
        cutoff_full$hdci_width
      },
      n_pairs_total,
      unit = "item pairs",
      fn = "RMlocdepQ3Cutoff()",
      id = "easyRasch2_band_flagging_locdep"
    )
  }

  validate_response_data(data)

  data <- as.data.frame(data)
  n_total <- nrow(data)
  has_na <- anyNA(data)
  data <- .drop_empty_respondents(data)
  n_used <- nrow(data)
  if (!is.null(cutoff_full)) {
    .check_cutoff_sample(
      cutoff_full$sample_n,
      n_used,
      "RMlocdepQ3Cutoff()",
      policy = "respondents with at least one response"
    )
  }

  # --- Compute Q3 residual correlations ---------------------------------------
  # .q3_residual_matrix() returns a symmetric Q3 matrix with an NA diagonal,
  # computed under the chosen estimator (CML/WLE residual correlations, or
  # mirt's MML Q3). The same routine is used for the simulated Q3 in
  # RMlocdepQ3Cutoff(), so observed and simulated values are comparable.
  resid_mat <- .q3_residual_matrix(data, estimator = estimator)
  # Warn about categories used only by extreme scorers (they carry no
  # information in CML), after the fit so that data stopped by the category
  # check do not warn first.
  if (estimator == "CML") {
    .sparsity_warning(data)
  }

  # Mean of the off-diagonal correlations (diagonal already NA).
  mean_resid <- mean(resid_mat, na.rm = TRUE)

  # Matrix view (global dynamic cutoff = mean + cutoff, or raw when NULL).
  matrix_out <- .q3_matrix_output(
    resid_mat,
    mean_resid,
    cutoff,
    output,
    cutoff_full,
    n_used = n_used,
    n_total = n_total,
    has_na = has_na
  )

  # NULL or bare-numeric cutoff: just the matrix (single object).
  if (is.null(cutoff_full)) {
    return(matrix_out)
  }

  # Full RMlocdepQ3Cutoff() object: the matrix view and the per-pair table.
  # The Q3 tile heatmap now lives in RMlocdepQ3Plot()$matrix, so this table
  # function no longer emits a plot (consistent with the other *Cutoff/table
  # functions, which return data and leave rendering to the *Plot functions).
  pairs_out <- .q3_pairs_table(
    resid_mat,
    cutoff_full,
    n_pairs,
    p_value,
    correction,
    alpha,
    output
  )
  list(matrix = matrix_out, pairs = pairs_out)
}

#' Lower-triangle heatmap of the observed Q3 matrix
#'
#' Viridis tile plot of the Q3 residual correlations (lower triangle),
#' styled to match the package's other ggplots. Pairs whose Q3 exceeds the
#' global dynamic cut-off are outlined in red. Returns `NULL` (with a
#' message) when `ggplot2` is not installed, so the rest of the list output
#' is unaffected.
#'
#' @param resid_mat Observed Q3 matrix (symmetric, `NA` diagonal, item names).
#' @param dyn_cutoff Numeric global dynamic cut-off, or `NULL`.
#' @param estimator Character estimator label for the caption.
#' @param actual_iterations Number of simulation iterations (for the caption),
#'   or `NULL`.
#' @return A `ggplot` object, or `NULL`.
#' @keywords internal
#' @noRd
.q3_tile_plot <- function(
  resid_mat,
  dyn_cutoff = NULL,
  estimator = "CML",
  actual_iterations = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message(
      "Package 'ggplot2' is required for the Q3 heatmap ($plot); ",
      "returning NULL. Install it with: install.packages(\"ggplot2\")"
    )
    return(NULL)
  }

  items <- colnames(resid_mat)
  if (is.null(items)) {
    items <- as.character(seq_len(ncol(resid_mat)))
  }
  lt <- which(lower.tri(resid_mat), arr.ind = TRUE)
  df <- data.frame(
    row = factor(items[lt[, "row"]], levels = items),
    col = factor(items[lt[, "col"]], levels = items),
    Q3 = resid_mat[lt],
    stringsAsFactors = FALSE
  )
  df$flagged <- if (!is.null(dyn_cutoff)) df$Q3 > dyn_cutoff else FALSE
  df$label <- formatC(df$Q3, format = "f", digits = 2)

  # Diverging fill centred on the mean off-diagonal Q3 (the value expected
  # under local independence), following RASCHplot's ggQ3star(): the neutral
  # colour marks the baseline, warm (red) is above it, cool (blue) below.
  # Symmetric limits about the centre keep the two sides comparable in
  # intensity.
  mean_resid <- mean(resid_mat, na.rm = TRUE)
  # Round the half-range up to the next 0.01 (as RASCHplot does) so the most
  # extreme cell sits strictly inside the limits rather than on the boundary,
  # where floating-point rounding could clip it to NA (grey).
  max_dev <- ceiling(max(abs(df$Q3 - mean_resid), na.rm = TRUE) * 100) / 100

  est_label <- switch(
    estimator,
    CML = "CML item / WLE person",
    MML = "MML item / EAP person",
    estimator
  )
  cap <- paste0(
    "Yen's Q3 residual correlations (lower triangle; ",
    est_label,
    "). Diverging fill centred on the mean off-diagonal Q3 (",
    formatC(mean_resid, format = "f", digits = 3),
    "), the value expected under local independence.",
    if (!is.null(dyn_cutoff)) {
      paste0(
        " Dynamic cut-off ",
        formatC(dyn_cutoff, format = "f", digits = 3),
        "; ",
        sum(df$flagged, na.rm = TRUE),
        " pair(s) above it outlined in black."
      )
    } else {
      ""
    },
    if (!is.null(actual_iterations)) {
      paste0(" Cut-off from ", actual_iterations, " simulation iterations.")
    } else {
      ""
    }
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$col, y = .data$row, fill = .data$Q3)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5)
  if (any(df$flagged, na.rm = TRUE)) {
    # Black outline (the RdYlBu high end is red, so a red outline would clash).
    p <- p +
      ggplot2::geom_tile(
        data = df[df$flagged %in% TRUE, , drop = FALSE],
        color = "black",
        linewidth = 1.1,
        fill = NA,
        width = 0.92,
        height = 0.92
      )
  }
  p +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      color = "grey15",
      size = 3
    ) +
    # RdYlBu anchor colours, diverging about the local-independence baseline.
    # high Q3 -> red (potential local dependence), low -> blue.
    ggplot2::scale_fill_gradient2(
      name = "Q3",
      low = "#4575b4",
      mid = "#ffffbf",
      high = "#d73027",
      midpoint = mean_resid,
      limits = mean_resid + c(-1, 1) * max_dev
    ) +
    ggplot2::scale_x_discrete(limits = utils::head(items, -1L)) +
    ggplot2::scale_y_discrete(limits = rev(items[-1L])) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, caption = er2_caption(cap)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    er2_axis_margins() +
    er2_plot_caption()
}

#' Compute the Yen's Q3 residual-correlation matrix under a chosen estimator
#'
#' Single entry point used by `RMlocdepQ3()`, the per-iteration simulation in
#' `RMlocdepQ3Cutoff()`, and the observed overlay in `RMlocdepQ3Plot()`, so
#' observed and simulated Q3 are always computed the same way.
#'
#' `"CML"` (default) returns the correlation matrix of the CML/WLE standardized
#' residuals from `.rasch_std_residuals()`; `"MML"` returns `mirt`'s Q3
#' residuals from a marginal-ML Rasch fit. The matrix is symmetric with an `NA`
#' diagonal and carries the item names.
#'
#' @param data Numeric response matrix or data.frame (items from 0; `NA` ok).
#' @param estimator `"CML"` or `"MML"`.
#' @param fast Logical. For `"MML"` only, use the faster (lower-precision)
#'   `mirt` settings (`quadpts = 29`, `TOL = 0.005`) appropriate for the
#'   simulation loop. Ignored for `"CML"`.
#' @return Symmetric Q3 matrix with `NA` diagonal.
#' @keywords internal
#' @noRd
.q3_residual_matrix <- function(
  data,
  estimator = c("CML", "MML"),
  fast = FALSE
) {
  estimator <- match.arg(estimator)
  if (estimator == "MML") {
    # .mml_fit() runs the category check (mirt renumbers items with an empty
    # lowest or interior category); the CML path gets it through
    # .fit_cml_thresholds().
    mirt_fit <- if (fast) {
      .mml_fit(data, accelerate = "squarem", quadpts = 29, TOL = 0.005)
    } else {
      .mml_fit(data, accelerate = "squarem")
    }
    m <- mirt::residuals(mirt_fit, type = "Q3", digits = 4, verbose = FALSE)
    diag(m) <- NA
    return(m)
  }
  # CML/WLE: Q3 is the column-wise correlation of standardized residuals.
  m <- stats::cor(
    .rasch_std_residuals(data, method = "WLE"),
    use = "pairwise.complete.obs"
  )
  diag(m) <- NA
  m
}

#' Render the Q3 lower-triangle matrix with an optional global cutoff
#'
#' @param resid_mat Observed Q3 matrix (symmetric, `NA` diagonal).
#' @param mean_resid Mean of the off-diagonal Q3 values.
#' @param cutoff Numeric global cutoff (added to `mean_resid`), or `NULL`.
#' @param output `"kable"` or `"dataframe"`.
#' @param cutoff_full Optional full cutoff object, used to enrich the caption.
#' @keywords internal
#' @noRd
.q3_matrix_output <- function(
  resid_mat,
  mean_resid,
  cutoff,
  output,
  cutoff_full = NULL,
  n_used = NULL,
  n_total = NULL,
  has_na = FALSE
) {
  resid_df <- as.data.frame(resid_mat)
  resid_df[upper.tri(resid_df)] <- NA
  diag(resid_df) <- NA

  if (!is.null(cutoff)) {
    dyn_cutoff <- mean_resid + cutoff
    resid_df$above_cutoff <- apply(
      resid_df,
      1L,
      function(row) any(row > dyn_cutoff, na.rm = TRUE)
    )
  }

  if (output == "dataframe") {
    return(resid_df)
  }

  n_clause <- if (!is.null(n_used) && !is.null(n_total)) {
    paste0(
      " ",
      .n_caption(
        n_used,
        n_total,
        if (has_na) "incomplete responses retained" else character()
      ),
      "."
    )
  } else {
    ""
  }
  if (!is.null(cutoff)) {
    caption_text <- paste0(
      "Dynamic cut-off: ",
      round(dyn_cutoff, 3),
      " (mean Q3 ",
      round(mean_resid, 3),
      " + ",
      round(cutoff, 3),
      ").",
      if (!is.null(cutoff_full)) {
        paste0(
          " Global simulation cutoff (99th pctl of max-mean Q3) from ",
          cutoff_full$actual_iterations,
          " iterations."
        )
      } else {
        ""
      },
      " Correlations exceeding the cut-off may indicate local dependence; ",
      "see the per-pair table for detail.",
      n_clause
    )
  } else {
    caption_text <- paste0(
      "Raw Q3 residual correlations (lower triangle). Use RMlocdepQ3Cutoff() ",
      "to derive a cutoff.",
      n_clause
    )
  }

  item_cols <- setdiff(names(resid_df), "above_cutoff")
  resid_display <- as.data.frame(
    lapply(resid_df[item_cols], function(x) {
      ifelse(is.na(x), "", as.character(round(x, 2)))
    }),
    stringsAsFactors = FALSE
  )
  rownames(resid_display) <- rownames(resid_df)
  if (!is.null(cutoff)) {
    resid_display$above_cutoff <- ifelse(resid_df$above_cutoff, "*", "")
  }

  knitr::kable(resid_display, caption = caption_text, format = "pipe")
}

#' Build the per-pair Q3 table (observed vs expected range, flagged)
#'
#' One row per item pair: observed Q3, the per-pair expected range (the
#' simulated `Low`/`High` bounds), and a `Flagged` label -- `"above"` (Q3
#' above the upper bound, local dependence), `"below"` (below the lower
#' bound), or `""`. Rows are sorted by absolute departure from the per-pair
#' simulated median; `n_pairs` keeps the top rows. With `p_value = TRUE`,
#' one-sided bootstrap p-values (`p_q3`, `padj_q3`) are added and `Flagged`
#' reflects `padj_q3 < alpha`.
#'
#' @param resid_mat Observed Q3 matrix (symmetric, `NA` diagonal).
#' @param cutoff_full Full [RMlocdepQ3Cutoff()] object.
#' @param n_pairs Integer or `NULL`; keep the top-`n_pairs` pairs by departure.
#' @param p_value Logical; add bootstrap p-values and flag on them.
#' @param correction,alpha Passed to `.bootstrap_pvalues()` / flagging.
#' @param output `"kable"` or `"dataframe"`.
#' @return A per-pair table (kable or data.frame).
#' @keywords internal
#' @noRd
.q3_pairs_table <- function(
  resid_mat,
  cutoff_full,
  n_pairs,
  p_value,
  correction,
  alpha,
  output
) {
  pc <- cutoff_full$pair_cutoffs # Item1, Item2, Q3_low, Q3_high
  pr <- cutoff_full$pair_results # Item1, Item2, Q3, iteration
  # The family is every pair, not the `n_pairs` shown, so the error rate the
  # interval implies is computed before the table is truncated.
  n_pairs_total <- nrow(pc)
  pr$key <- paste(pr$Item1, pr$Item2, sep = "___")
  keys <- paste(pc$Item1, pc$Item2, sep = "___")

  observed <- mapply(function(a, b) resid_mat[a, b], pc$Item1, pc$Item2)
  median_q3 <- tapply(pr$Q3, pr$key, stats::median)[keys] # expected (sort/dir)
  low <- pc$Q3_low
  high <- pc$Q3_high

  # Unrounded values; the kable path rounds a display copy before rendering.
  tbl <- data.frame(
    Item1 = pc$Item1,
    Item2 = pc$Item2,
    Observed = as.numeric(observed),
    Low = as.numeric(low),
    High = as.numeric(high),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (p_value) {
    sim_mat <- tapply(pr$Q3, list(pr$iteration, pr$key), function(x) x[1L])
    obs_named <- stats::setNames(as.numeric(observed), keys)
    pv <- .bootstrap_pvalues(
      obs_named,
      sim_mat,
      correction = correction,
      tail = "upper"
    )
    idx <- match(keys, pv$name)
    tbl$p_q3 <- pv$p[idx]
    tbl$padj_q3 <- pv$padj[idx]
    is_flagged <- !is.na(tbl$padj_q3) & tbl$padj_q3 < alpha
  } else {
    is_flagged <- observed > high | observed < low
  }

  # Direction relative to the per-pair simulated median.
  tbl$Flagged <- ifelse(
    !is_flagged,
    "",
    ifelse(observed > median_q3, "above", "below")
  )

  # Sort by absolute departure from expected, then keep the top n_pairs.
  ord <- order(abs(observed - median_q3), decreasing = TRUE)
  tbl <- tbl[ord, , drop = FALSE]
  rownames(tbl) <- NULL
  if (!is.null(n_pairs)) {
    tbl <- tbl[seq_len(min(as.integer(n_pairs), nrow(tbl))), , drop = FALSE]
  }

  if (output == "dataframe") {
    return(tbl)
  }

  # Kable display rounding (the dataframe output above stays unrounded)
  tbl <- .round_display(tbl, c(
    Observed = 3, Low = 3, High = 3, p_q3 = 4, padj_q3 = 4
  ))

  width_pct <- round(100 * cutoff_full$hdci_width, 1)
  n_iter <- cutoff_full$actual_iterations
  caption <- paste0(
    "Q3 by item pair, sorted by departure from the expected range. ",
    "Expected range = ",
    width_pct,
    "% interval of the simulated Q3 per pair (",
    n_iter,
    " iterations)."
  )
  if (p_value) {
    corr_label <- .correction_label(correction)
    caption <- paste0(
      caption,
      " p_q3/padj_q3: one-sided bootstrap p-values (upper tail, local ",
      "dependence), ",
      corr_label,
      ". Flagged at padj < ",
      alpha,
      ". The interval is shown as description and is not the decision rule, ",
      "so a pair below the lower bound is not flagged.",
      .iteration_note(n_iter, correction),
      .attrition_clause(n_iter, cutoff_full$requested_iterations)
    )
    col.names <- c(
      "Item 1",
      "Item 2",
      "Observed Q3",
      "Exp. low",
      "Exp. high",
      "p",
      "p (adj)",
      "Flagged"
    )
  } else {
    caption <- paste0(
      caption,
      " Flagged: above = Q3 above the upper bound (local dependence), ",
      "below = below the lower bound.",
      .band_error_clause(
        if (identical(cutoff_full$cutoff_method, "quantile")) {
          0.95
        } else {
          cutoff_full$hdci_width
        },
        n_pairs_total,
        unit = "item pairs"
      ),
      .attrition_clause(n_iter, cutoff_full$requested_iterations)
    )
    col.names <- c(
      "Item 1",
      "Item 2",
      "Observed Q3",
      "Exp. low",
      "Exp. high",
      "Flagged"
    )
  }
  knitr::kable(
    tbl,
    format = "pipe",
    col.names = col.names,
    caption = caption,
    row.names = FALSE
  )
}

#' Simulation-based Q3 cutoffs
#'
#' Simulates datasets from a Rasch or partial credit model that fits and
#' resembles your data, and computes \eqn{Q_3} for every item pair in each.
#' Returns percentiles of the maximum minus mean \eqn{Q_3} for a global
#' cutoff (Christensen et al., 2017) and an interval per pair. Passed to
#' [RMlocdepQ3()], the result gives bootstrap p-values corrected across
#' pairs.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Respondents with no responses are dropped, incomplete ones are kept.
#' @param iterations Integer. Number of simulated datasets. Default `400`,
#'   the recommended minimum for reproducible Westfall-Young decisions and
#'   the number a 95% interval needs to converge (Johansson, 2026). Use 1000
#'   to 2000 for a final analysis.
#' @param parallel Logical. Run the simulations in parallel via \pkg{mirai}
#'   when it is installed. Default `TRUE`.
#' @param n_cores Integer or `NULL`. Number of parallel workers. When
#'   `NULL`, `getOption("mc.cores")` is used. If that is unset too, a
#'   warning is given and the simulations run sequentially.
#' @param verbose Logical. Show a progress bar. Default `FALSE`.
#' @param seed Integer or `NULL`. Random seed for reproducibility. See
#'   [easyRasch2-reproducibility] for what this guarantees and how it
#'   interacts with `parallel`.
#' @param cutoff_method Character. How the interval per pair in
#'   `pair_cutoffs` is computed: `"hdci"` (default), the highest-density
#'   continuous interval from `ggdist::hdci()` (needs \pkg{ggdist}), or
#'   `"quantile"`, the 2.5th to 97.5th percentiles. The global
#'   `suggested_cutoff` is not affected.
#' @param hdci_width Numeric in (0, 1). Width of the interval with
#'   `cutoff_method = "hdci"`. Default `0.95` (`0.99` before 1.2.0). The
#'   interval describes where a fitting pair's \eqn{Q_3} is expected to
#'   fall and is no longer the default decision rule, so the width is
#'   chosen to converge at the default number of iterations. Ignored with
#'   `"quantile"`.
#' @param estimator Character. `"CML"` (default, CML item parameters and
#'   WLE person locations) or `"MML"` (`mirt`). It is stored in the result
#'   and used by [RMlocdepQ3()] and [RMlocdepQ3Plot()].
#' @param dgp Character. How datasets are simulated. `"resample"` (default)
#'   draws person locations from the WLE estimates with replacement and
#'   simulates responses from the model (a marginal null). Each drawn
#'   respondent keeps their pattern of missing responses, so the simulated
#'   \eqn{Q_3} of a pair rests on as many respondents as the observed one.
#'   `"conditional"` simulates each respondent's pattern from the model
#'   given their observed total score and answered items, with the item
#'   parameters fixed (a conditional null). It needs no latent
#'   distribution and avoids the over-dispersion of resampled point
#'   estimates. The two give different cutoffs. \strong{Experimental.}
#'
#' @return A list with:
#' \describe{
#'   \item{`results`}{data.frame with `mean`, `max` and `diff` (maximum
#'     minus mean \eqn{Q_3}), one row per simulated dataset.}
#'   \item{`pair_results`}{data.frame with `Item1`, `Item2`, `Q3` and
#'     `iteration`, one row per item pair and simulated dataset.}
#'   \item{`pair_cutoffs`}{data.frame with the interval per pair: `Item1`,
#'     `Item2`, `Q3_low`, `Q3_high`, computed by `cutoff_method`.}
#'   \item{`actual_iterations`}{Number of simulated datasets that could be
#'     used.}
#'   \item{`requested_iterations`}{The `iterations` argument.}
#'   \item{`sample_n`}{Number of respondents used (those with at least one
#'     response).}
#'   \item{`sample_n_total`}{Number of respondents in `data`.}
#'   \item{`sample_has_na`}{Logical. Whether `data` contained missing
#'     values.}
#'   \item{`sample_summary`}{Summary of the WLE person locations.}
#'   \item{`item_names`}{Item names.}
#'   \item{`max_diff`, `sd_diff`}{Maximum and SD of `diff`.}
#'   \item{`p95`, `p99`, `p995`, `p999`}{Percentiles of `diff`.}
#'   \item{`suggested_cutoff`}{The 99th percentile (`p99`), the global
#'     cutoff for [RMlocdepQ3()].}
#'   \item{`cutoff_method`}{`"hdci"` or `"quantile"`.}
#'   \item{`hdci_width`}{The HDCI width (used only with `"hdci"`).}
#'   \item{`estimator`}{`"CML"` or `"MML"`.}
#'   \item{`dgp`}{`"resample"` or `"conditional"`.}
#' }
#'
#' @details
#' The generating model is fitted once, with CML item parameters
#' (`psychotools`) and WLE person locations. In each iteration a dataset is
#' simulated by `dgp`, the model is refitted, and \eqn{Q_3} is computed
#' with `estimator`. A simulated dataset that cannot be refitted is left
#' out, and [RMlocdepQ3()] reports how many in its caption.
#'
#' Parallel processing needs the optional \pkg{mirai} package.
#'
#' @references
#' Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
#' for Yen's \eqn{Q_3}: Identification of local dependence in the Rasch model.
#' *Applied Psychological Measurement, 41*(3), 178--194.
#' \doi{10.1177/0146621616677520}
#'
#' Johansson, M. (2026). Simulation-based cutoffs for conditional item fit in
#' Rasch models: Iterations, multiplicity correction, and decision stability.
#' *PsyArXiv*. \doi{10.31234/osf.io/7pqz4_v2}
#'
#' @seealso [RMlocdepQ3()] for the tables, [RMlocdepQ3Plot()] for the plots.
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggdist", quietly = TRUE)) {
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:10)
#'
#'   # Few iterations for a fast example. Use 1000 or more in real analyses.
#'   cutoff_res <- RMlocdepQ3Cutoff(sim_data, iterations = 50, parallel = FALSE,
#'                                  seed = 42)
#'   cutoff_res$suggested_cutoff  # 99th percentile
#'
#'   # Use the cutoff in RMlocdepQ3()
#'   RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
#' }
#' }
RMlocdepQ3Cutoff <- function(
  data,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95,
  estimator = c("CML", "MML"),
  dgp = c("resample", "conditional")
) {
  validate_response_data(data)

  estimator <- match.arg(estimator)
  dgp <- match.arg(dgp)
  cutoff_method <- match.arg(cutoff_method, c("hdci", "quantile"))
  if (cutoff_method == "hdci" && !requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required when cutoff_method = \"hdci\" but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")\n",
      "Alternatively, use cutoff_method = \"quantile\" to avoid this dependency.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(hdci_width) ||
      length(hdci_width) != 1L ||
      !is.finite(hdci_width) ||
      hdci_width <= 0 ||
      hdci_width >= 1
  ) {
    stop("`hdci_width` must be a single number in (0, 1).", call. = FALSE)
  }

  use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)

  if (parallel && !use_parallel) {
    message(
      "Install 'mirai' package for parallel processing: install.packages(\"mirai\")"
    )
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
          "(Use `parallel::detectCores()` to see how many cores are available.)\n",
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

  # Incomplete responses are retained (they feed the per-pattern WLE pool and
  # the conditional DGP), but respondents with no responses at all must go:
  # they contribute nothing and break the CML fit (psychotools errors on
  # all-NA rows). Record the raw total and the missingness flag first so
  # callers (e.g. RMlocdepQ3Plot) can report the sample in the standard
  # `n = X of Y respondents (policy)` form, matching RMlocdepQ3().
  sample_n_total <- nrow(as.data.frame(data))
  sample_has_na <- anyNA(data)
  data <- .drop_empty_respondents(data)

  data_mat <- as.matrix(data)
  sample_n <- nrow(data_mat)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  # Preserve item names so the per-iteration pair_q3 frames use the user's
  # labels (e.g., "q1", "q2") rather than psychotools / mirt's auto-generated
  # "V1", "V2", ... .
  item_names_sim <- colnames(data_mat)
  if (is.null(item_names_sim)) {
    item_names_sim <- paste0("V", seq_len(ncol(data_mat)))
  }

  # Generating model: CML item thresholds (psychotools), computed once. WLE
  # person locations are computed for the sample summary and, under the
  # "resample" DGP, as the generating theta pool.
  pool <- .wle_theta_pool(data_mat)
  thr_list <- pool$thr_list
  wle_thetas <- pool$thetas

  sim_data_list <- list(
    dgp = dgp,
    type = if (is_polytomous) "polytomous" else "dichotomous",
    thr_list = thr_list,
    n_items = ncol(data_mat),
    sample_n = sample_n,
    item_names = item_names_sim,
    estimator = estimator
  )
  if (dgp == "resample") {
    # Marginal DGP: resample respondents with replacement, then simulate data
    # parametrically under the model (items and thetas on a common scale).
    # A respondent carries both their WLE theta and their missingness pattern,
    # so that incomplete respondents stay incomplete in the simulated data.
    # The observed Q3 is a pairwise-complete correlation, and a complete
    # simulated dataset would give every pair more respondents than the
    # observed one has, making the null too narrow. The WLE pool keeps one
    # finite theta per row (the mask is subset the same way), and for complete
    # data the draws are identical to resampling the thetas directly.
    sim_data_list$thetas <- wle_thetas
    if (anyNA(data_mat)) {
      theta_rows <- which(is.finite(
        .estimate_thetas(data_mat, thr_list, method = "WLE")$theta
      ))
      sim_data_list$na_mask <- is.na(data_mat)[theta_rows, , drop = FALSE]
    }
    if (is_polytomous) {
      sim_data_list$deltaslist <- thr_list
    } else {
      sim_data_list$item_params <- unlist(thr_list, use.names = FALSE)
    }
  } else {
    # Conditional DGP: simulate each respondent's pattern from the exact Rasch
    # conditional distribution given their observed total score (and answered
    # items), item parameters fixed. No latent distribution is estimated, so
    # the sufficient statistic -- not a resampled point estimate -- defines the
    # null. Respondents are grouped by (answered-set, score) so each group is
    # drawn in one batch.
    sim_data_list$cond_groups <- .cond_groups(data_mat, thr_list)
  }

  if (use_parallel) {
    results_raw <- run_q3_sim_parallel(
      iterations,
      sim_seeds,
      sim_data_list,
      n_cores,
      verbose
    )
  } else {
    results_raw <- run_q3_sim_sequential(
      iterations,
      sim_seeds,
      sim_data_list,
      verbose
    )
  }

  # Filter out failures (character strings indicate errors)
  ok <- vapply(results_raw, function(x) is.list(x), logical(1L))
  successful <- results_raw[ok]

  if (length(successful) == 0L) {
    stop(.all_sims_failed_message(data_mat), call. = FALSE)
  }

  actual_iterations <- length(successful)
  mean_q3 <- vapply(successful, function(x) x$mean, numeric(1L))
  max_q3 <- vapply(successful, function(x) x$max, numeric(1L))
  diff_q3 <- max_q3 - mean_q3

  results <- data.frame(
    mean = mean_q3,
    max = max_q3,
    diff = diff_q3
  )

  # --- Per-pair aggregation -------------------------------------------------
  # Stack per-iteration pair_q3 frames into one long data.frame and add an
  # `iteration` column (matching the partial-gamma family).
  item_names_vec <- colnames(data)
  if (is.null(item_names_vec)) {
    item_names_vec <- as.character(seq_len(ncol(data)))
  }
  pair_iter_dfs <- lapply(seq_along(successful), function(i) {
    df <- successful[[i]]$pair_q3
    df$iteration <- i
    df
  })
  pair_results <- do.call(rbind, pair_iter_dfs)
  rownames(pair_results) <- NULL

  # Per-pair cutoff intervals
  pair_keys <- unique(paste(
    pair_results$Item1,
    pair_results$Item2,
    sep = "___"
  ))
  pair_cutoffs <- do.call(
    rbind,
    lapply(pair_keys, function(pk) {
      parts <- strsplit(pk, "___", fixed = TRUE)[[1L]]
      sub <- pair_results[
        pair_results$Item1 == parts[1L] &
          pair_results$Item2 == parts[2L],
      ]
      if (cutoff_method == "hdci") {
        q3_interval <- ggdist::hdci(sub$Q3, .width = hdci_width)
        data.frame(
          Item1 = parts[1L],
          Item2 = parts[2L],
          Q3_low = q3_interval[1L, 1L],
          Q3_high = q3_interval[1L, 2L],
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      } else {
        data.frame(
          Item1 = parts[1L],
          Item2 = parts[2L],
          Q3_low = stats::quantile(sub$Q3, 0.025, na.rm = TRUE),
          Q3_high = stats::quantile(sub$Q3, 0.975, na.rm = TRUE),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
    })
  )
  rownames(pair_cutoffs) <- NULL

  out <- list()
  out$results <- results
  out$pair_results <- pair_results
  out$pair_cutoffs <- pair_cutoffs
  out$actual_iterations <- actual_iterations
  out$requested_iterations <- iterations
  out$sample_n <- sample_n
  out$sample_n_total <- sample_n_total
  out$sample_has_na <- sample_has_na
  out$sample_summary <- summary(wle_thetas)
  out$item_names <- item_names_vec
  out$max_diff <- max(results$diff)
  out$sd_diff <- stats::sd(results$diff)
  out$p95 <- stats::quantile(results$diff, 0.95)
  out$p99 <- stats::quantile(results$diff, 0.99)
  out$p995 <- stats::quantile(results$diff, 0.995)
  out$p999 <- stats::quantile(results$diff, 0.999)
  out$suggested_cutoff <- out$p99
  out$cutoff_method <- cutoff_method
  out$hdci_width <- hdci_width
  out$estimator <- estimator
  out$dgp <- dgp
  out
}

# ---------------------------------------------------------------------------
# Internal: single simulation iteration
# ---------------------------------------------------------------------------

#' Group respondents by sufficient statistic for the conditional DGP
#'
#' Shared by the conditional-DGP cutoff simulators (`RMlocdepQ3Cutoff()`,
#' `RMitemInfitCutoff()`, ...). Returns one element per distinct (answered-item
#' set, total score) so each group can be drawn in a single
#' `.sim_conditional()` call.
#'
#' @param data_mat Observed response matrix.
#' @param thr_list Item thresholds (unused here but kept for symmetry).
#' @return List of `list(rows, ans, score)` groups.
#' @keywords internal
#' @noRd
.cond_groups <- function(data_mat, thr_list) {
  n <- nrow(data_mat)
  ans_key <- apply(!is.na(data_mat), 1L, function(z) {
    paste0(which(z), collapse = ",")
  })
  score <- rowSums(data_mat, na.rm = TRUE)
  key <- paste(ans_key, score, sep = "|")
  unname(lapply(split(seq_len(n), key), function(rows) {
    ans <- which(!is.na(data_mat[rows[1L], ]))
    list(rows = rows, ans = ans, score = sum(data_mat[rows[1L], ans]))
  }))
}

#' Simulate one conditional-DGP dataset
#'
#' For each (answered-set, score) group, draw the required number of response
#' patterns from the Rasch conditional distribution given the total score
#' (`.sim_conditional()`, in person_fit.R); extreme scores (0 or maximum) yield
#' their unique degenerate pattern. The observed missingness pattern is
#' preserved (unanswered cells stay `NA`). Shared by the conditional-DGP cutoff
#' simulators.
#'
#' @param data_list A simulation list with `dgp = "conditional"` (carries
#'   `thr_list`, `cond_groups`, `sample_n`, `n_items`).
#' @return A data.frame of simulated responses.
#' @keywords internal
#' @noRd
.sim_cond_dataset <- function(data_list) {
  N <- data_list$sample_n
  K <- data_list$n_items
  out <- matrix(NA_integer_, N, K)
  for (g in data_list$cond_groups) {
    thr_sub <- data_list$thr_list[g$ans]
    m <- length(g$rows)
    pat <- .sim_conditional(thr_sub, g$score, m)
    if (is.null(pat)) {
      # Extreme score: deterministic pattern (all-0 or all-maximum-category).
      steps <- vapply(thr_sub, length, integer(1L))
      pat <- if (g$score <= 0L) {
        matrix(0L, m, length(g$ans))
      } else {
        matrix(steps, nrow = m, ncol = length(g$ans), byrow = TRUE)
      }
    }
    out[g$rows, g$ans] <- pat
  }
  as.data.frame(out)
}

#' Run a single Q3 simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMlocdepQ3Cutoff()].
#' @return A list with `mean` and `max` Q3, or a character string on failure.
#' @noRd
run_single_q3_sim <- function(seed, data_list) {
  # The RNG kind is pinned, not just the seed: mirai daemons start under
  # L'Ecuyer-CMRG while the calling session uses the Mersenne-Twister
  # default, so seeding alone would make the parallel and sequential paths
  # draw different streams from the same `seed`.
  set.seed(
    seed,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )

  tryCatch(
    {
      # --- Generate one simulated dataset under the chosen DGP -----------------
      if (identical(data_list$dgp, "conditional")) {
        sim_df <- .sim_cond_dataset(data_list)
      } else {
        # Resample respondents: row indices, so each draw carries its theta
        # and (when the data are incomplete) its missingness pattern.
        # sample(x, size, replace = TRUE) draws x[sample.int(length(x), ...)],
        # so complete data reproduce the earlier theta-only draws exactly.
        idx <- sample.int(
          length(data_list$thetas),
          size = data_list$sample_n,
          replace = TRUE
        )
        thetas_res <- data_list$thetas[idx]
        sim_df <- if (data_list$type == "dichotomous") {
          as.data.frame(
            psychotools::rrm(
              theta = thetas_res,
              beta = data_list$item_params
            )$data
          )
        } else {
          as.data.frame(sim_partial_score(
            data_list$deltaslist,
            thetas_res
          ))
        }
        if (!is.null(data_list$na_mask)) {
          sim_df[data_list$na_mask[idx, , drop = FALSE]] <- NA
        }
      }

      # --- Validate the simulated dataset (estimable refit) --------------------
      if (data_list$type == "dichotomous") {
        # Every item must have at least 8 positive responses (numerical stability).
        if (any(colSums(sim_df, na.rm = TRUE) < 8L)) {
          return(
            "validation_failed: fewer than 8 positive responses in at least one item"
          )
        }
      } else {
        # Every item must show all categories.
        n_cats <- vapply(
          data_list$thr_list,
          function(d) length(d) + 1L,
          integer(1L)
        )
        for (j in seq_len(ncol(sim_df))) {
          tab <- tabulate(sim_df[[j]] + 1L, nbins = n_cats[j])
          if (any(tab == 0L)) {
            return("validation_failed: not all categories represented")
          }
        }
      }

      # Preserve the user's item labels so the Q3 matrix below uses them
      if (
        !is.null(data_list$item_names) &&
          length(data_list$item_names) == ncol(sim_df)
      ) {
        colnames(sim_df) <- data_list$item_names
      }

      # Q3 matrix under the same estimator as the observed analysis. For "MML"
      # the faster (lower-precision) mirt settings are used in the loop.
      q3_mat <- .q3_residual_matrix(
        sim_df,
        estimator = data_list$estimator,
        fast = TRUE
      )

      mean_q3 <- mean(q3_mat, na.rm = TRUE)
      max_q3 <- max(q3_mat, na.rm = TRUE)

      # Extract upper triangle as a long vector for per-pair retention.
      # mirt::residuals returns a symmetric matrix with item names on rows/cols.
      item_names_q3 <- colnames(q3_mat)
      if (is.null(item_names_q3)) {
        item_names_q3 <- as.character(seq_len(ncol(q3_mat)))
      }
      upper_idx <- which(upper.tri(q3_mat), arr.ind = TRUE)
      pair_q3 <- data.frame(
        Item1 = item_names_q3[upper_idx[, "row"]],
        Item2 = item_names_q3[upper_idx[, "col"]],
        Q3 = q3_mat[upper_idx],
        stringsAsFactors = FALSE,
        row.names = NULL
      )

      list(mean = mean_q3, max = max_q3, pair_q3 = pair_q3)
    },
    error = function(e) {
      as.character(conditionMessage(e))
    }
  )
}

# ---------------------------------------------------------------------------
# Internal: parallel runner
# ---------------------------------------------------------------------------

#' Run Q3 simulations in parallel using mirai
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param n_cores Number of mirai daemons.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @noRd
run_q3_sim_parallel <- function(
  iterations,
  sim_seeds,
  sim_data_list,
  n_cores,
  verbose = FALSE
) {
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
        run_single_q3_sim(seed, data_list)
      },
      seed = sim_seeds[sim],
      data_list = sim_data_list,
      run_single_q3_sim = run_single_q3_sim,
      sim_partial_score = sim_partial_score,
      sim_poly_item = sim_poly_item,
      # Conditional-DGP generators.
      .sim_cond_dataset = .sim_cond_dataset,
      .sim_conditional = .sim_conditional,
      .esf_convolve = .esf_convolve,
      # Engine helpers needed by .q3_residual_matrix() inside the daemon.
      .q3_residual_matrix = .q3_residual_matrix,
      .mml_fit = .mml_fit,
      .rasch_std_residuals = .rasch_std_residuals,
      .fit_cml_thresholds = .fit_cml_thresholds,
      .cml_threshold_summary = .cml_threshold_summary,
      .cml_fit = .cml_fit,
      .check_response_categories = .check_response_categories,
      .estimate_thetas = .estimate_thetas,
      .theta_wle = .theta_wle,
      .pcm_cat_probs = .pcm_cat_probs,
      .center_thresholds = .center_thresholds
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

#' Run Q3 simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @noRd
run_q3_sim_sequential <- function(
  iterations,
  sim_seeds,
  sim_data_list,
  verbose = FALSE
) {
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  results <- vector("list", iterations)
  for (sim in seq_len(iterations)) {
    results[[sim]] <- run_single_q3_sim(sim_seeds[sim], sim_data_list)
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
