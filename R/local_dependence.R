#' \eqn{Q_3} Residual Correlations for Local Dependence Assessment
#'
#' Computes Yen's \eqn{Q_3} residual correlations between item pairs. By default
#' the Rasch model is fitted by conditional maximum likelihood with WLE person
#' locations (`estimator = "CML"`); marginal ML via `mirt` is available with
#' `estimator = "MML"`. High correlations (above the dynamic cut-off) indicate
#' potential local dependence between items. See \code{\link{RMlocdepQ3Cutoff}}
#' for how to determine the appropriate dynamic cut-off for your data.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed.
#' @param cutoff Optional. `NULL` (default) returns the raw \eqn{Q_3}matrix. A single
#'   numeric value returns the \eqn{Q_3} matrix with a global dynamic cut-off (the
#'   value added to the mean off-diagonal \eqn{Q_3}). The **full list** returned by
#'   \code{\link{RMlocdepQ3Cutoff}} returns a *list of two tables* (see Value):
#'   the cut-off matrix plus a per-pair table.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for formatted `knitr::kable()` table(s), or
#'   `"dataframe"` for the underlying data.frame(s).
#' @param n_pairs Integer or `NULL` (default). When the full cutoff object is
#'   supplied, limits the per-pair table to the `n_pairs` pairs with the
#'   largest departure from their expected range. `NULL` shows all pairs.
#' @param p_value Logical. If `TRUE` (requires the full
#'   \code{\link{RMlocdepQ3Cutoff}} object), the per-pair table also reports
#'   one-sided bootstrap p-values (`p_q3`, `padj_q3`) and flags `above` pairs
#'   (only) on `padj_q3 < alpha` instead of the expected range. Default `FALSE`.
#' @param correction Character. Multiple-comparison correction across item
#'   pairs when `p_value = TRUE`: `"fwer"` (default) for the Westfall-Young
#'   studentised-max step-down, `"fdr_bh"` / `"fdr_by"` for Benjamini-Hochberg
#'   / Benjamini-Yekutieli, or `"none"`.
#' @param alpha Numeric in (0, 1). Significance level for `Flagged` when
#'   `p_value = TRUE`. Default `0.05`.
#' @param estimator Character. Estimation engine for the \eqn{Q_3} residual
#'   correlations. `"CML"` (default) fits item parameters by conditional
#'   maximum likelihood (`psychotools`) and person locations by Warm's
#'   weighted likelihood (WLE) -- true to the Rasch tradition and finite
#'   at extreme scores. `"MML"` uses the marginal-ML / EAP engine
#'   (`mirt`), retained for backward compatibility. The estimator must
#'   match the one used by \code{\link{RMlocdepQ3Cutoff}}; when the full
#'   cutoff object is supplied, its stored estimator takes precedence (a
#'   mismatching `estimator` argument is overridden with a warning).
#'
#' @return
#' With `cutoff = NULL` or a bare numeric cut-off, a single object (`kable` or
#' data.frame) holding the lower triangle of the \eqn{Q_3} matrix; a numeric cut-off
#' adds a `Flagged` row flag and a caption describing the dynamic
#' cut-off.
#'
#' With the **full `RMlocdepQ3Cutoff()` object**, a named list of two:
#' \describe{
#'   \item{`$matrix`}{the \eqn{Q_3} lower-triangle matrix with the *global* dynamic
#'     cut-off (mean off-diagonal \eqn{Q_3} + suggested cut-off), as above.}
#'   \item{`$pairs`}{one row per item pair: `Item1`, `Item2`, `Observed` (\eqn{Q_3}),
#'     `Low`/`High` (the per-pair expected range, i.e. the simulated bounds),
#'     and `Flagged` -- `"above"` (\eqn{Q_3} above the upper bound, indicating local
#'     dependence), `"below"` (below the lower bound), or `""`. Sorted by
#'     absolute departure from the per-pair simulated median and truncated to
#'     `n_pairs`. With `p_value = TRUE`, columns `p_q3` and `padj_q3` are added
#'     and `Flagged` reflects `padj_q3 < alpha` and **only flags `"above"`**.}
#' }
#' The \eqn{Q_3} tile heatmap that earlier versions returned as `$plot` is now produced
#' by \code{\link{RMlocdepQ3Plot}} (as its `$matrix` element), so the table and
#' plot outputs share the same `$matrix`/`$pairs` structure.
#'
#' @details
#' The \eqn{Q_3} statistic (Yen, 1984) is the correlation between residuals of pairs
#' of items after accounting for the latent trait. Under local independence,
#' \eqn{Q_3} values are expected to be around \eqn{-1/(k-1)} where \eqn{k} is the
#' number of items. When `cutoff` is supplied, the dynamic cut-off is the mean
#' of all off-diagonal \eqn{Q_3} values plus `cutoff`, following the approach of
#' Christensen et al. (2017). Use \code{\link{RMlocdepQ3Cutoff}} to obtain a
#' simulation-based cutoff recommendation.
#'
#' \eqn{Q_3} is the column-wise correlation matrix of the model standardized
#' residuals \eqn{(x - E)/\sqrt{Var}}. By default (`estimator = "CML"`) item
#' parameters are estimated by conditional maximum likelihood (`psychotools`)
#' and person locations by Warm's weighted likelihood; `estimator = "MML"`
#' instead uses `mirt`'s marginal-ML model and its built-in \eqn{Q_3} residuals.
#' The two estimators give very similar \eqn{Q_3} values (off-diagonal correlation
#' typically > 0.95); what matters for inference is that the observed \eqn{Q_3} and
#' the simulated cut-off in \code{\link{RMlocdepQ3Cutoff}} use the *same*
#' estimator, which the functions enforce.
#'
#' \strong{Two views of local dependence.} Given the full
#' \code{\link{RMlocdepQ3Cutoff}} object, two complementary tables are returned.
#' The `$matrix` applies a single *global* cut-off (the Christensen et al.
#' approach: the 99th percentile of the simulated max-minus-mean \eqn{Q_3}) -- a
#' family-wise "is there any local dependence" overview. The `$pairs` table is
#' the per-comparison view: each observed \eqn{Q_3} against its own simulated expected
#' range (the `Low`/`High` bounds), so individual dependent pairs can be read
#' off and ranked.
#'
#' \strong{Bootstrap p-values.} When `p_value = TRUE`, the `$pairs` table also
#' tests each observed \eqn{Q_3} against its simulated null (from `cutoff$pair_results`)
#' with a one-sided (upper-tail) test for excess local dependence. The pair
#' statistic is studentised by the bootstrap mean and SD; the marginal p-value
#' is `(1 + #{Q3* >= Q3}) / (B + 1)`, and `correction` applies the family-wise
#' (Westfall-Young step-down) or FDR adjustment across the \eqn{k(k-1)/2} pairs.
#' As for item fit, the family-wise correction is liberal when the simulation is
#' small, so >= 1000 `iterations` in [RMlocdepQ3Cutoff()] are recommended (a
#' warning is issued otherwise).
#'
#' @references
#' Yen, W. M. (1984). Effects of local item dependence on the fit and
#' equating performance of the three-parameter logistic model.
#' *Applied Psychological Measurement, 8*(2), 125--145.
#' \doi{10.1177/014662168400800201}
#'
#' Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
#' for Yen's \eqn{Q_3}: Identification of local dependence in the Rasch model.
#' *Applied Psychological Measurement, 41*(3), 178--194.
#' \doi{10.1177/0146621616677520}
#'
#' Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate' number
#' of hypotheses simultaneously. *Journal of Statistical Theory and Practice,
#' 19*(6). \doi{10.1007/s42519-024-00412-4}
#'
#' @inheritSection RMitemInfit Multiple comparisons
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
#' # Simulation-based cutoff (use 500+ iterations in real analyses)
#' if (requireNamespace("ggdist", quietly = TRUE)) {
#'   cutoff_res <- RMlocdepQ3Cutoff(sim_data, iterations = 50, parallel = FALSE)
#'
#'   # Bare numeric cutoff -> just the matrix
#'   RMlocdepQ3(sim_data, cutoff = cutoff_res$suggested_cutoff)
#'
#'   # Full object -> list of two tables: $matrix and $pairs
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
  p_value = FALSE,
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

  if (p_value) {
    if (is.null(cutoff_full) || is.null(cutoff_full$pair_results)) {
      stop(
        "`p_value = TRUE` requires the full RMlocdepQ3Cutoff() object (it ",
        "carries the simulated per-pair distributions in $pair_results); a ",
        "numeric cutoff or NULL is not sufficient.",
        call. = FALSE
      )
    }
    if (
      !is.null(cutoff_full$actual_iterations) &&
        cutoff_full$actual_iterations < 1000L
    ) {
      warning(
        "Bootstrap p-values are based on only ",
        cutoff_full$actual_iterations,
        " simulation iterations. With few ",
        "iterations the studentised-max (FWER) correction is liberal and ",
        "small p-values are imprecise; use iterations >= 1000 in ",
        "RMlocdepQ3Cutoff() for reliable p-values.",
        call. = FALSE
      )
    }
  }

  validate_response_data(data)

  data <- as.data.frame(data)
  n_total <- nrow(data)
  has_na <- anyNA(data)
  data <- .drop_empty_respondents(data)
  n_used <- nrow(data)

  # The CML engine (psychotools) can choke or destabilise on sparse / zero-
  # variance response categories; warn (and point to estimator = "MML") before
  # fitting, consistent with RMitemParameters(). Not needed for the MML path.
  if (estimator == "CML") {
    .sparsity_warning(data, is_poly = max(as.matrix(data), na.rm = TRUE) > 1L)
  }

  # --- Compute Q3 residual correlations ---------------------------------------
  # .q3_residual_matrix() returns a symmetric Q3 matrix with an NA diagonal,
  # computed under the chosen estimator (CML/WLE residual correlations, or
  # mirt's MML Q3). The same routine is used for the simulated Q3 in
  # RMlocdepQ3Cutoff(), so observed and simulated values are comparable.
  resid_mat <- .q3_residual_matrix(data, estimator = estimator)

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
    mirt_fit <- if (fast) {
      mirt::mirt(
        data,
        model = 1,
        itemtype = "Rasch",
        verbose = FALSE,
        accelerate = "squarem",
        quadpts = 29,
        TOL = 0.005
      )
    } else {
      mirt::mirt(
        data,
        model = 1,
        itemtype = "Rasch",
        verbose = FALSE,
        accelerate = "squarem"
      )
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
  pr$key <- paste(pr$Item1, pr$Item2, sep = "___")
  keys <- paste(pc$Item1, pc$Item2, sep = "___")

  observed <- mapply(function(a, b) resid_mat[a, b], pc$Item1, pc$Item2)
  median_q3 <- tapply(pr$Q3, pr$key, stats::median)[keys] # expected (sort/dir)
  low <- pc$Q3_low
  high <- pc$Q3_high

  tbl <- data.frame(
    Item1 = pc$Item1,
    Item2 = pc$Item2,
    Observed = round(as.numeric(observed), 3),
    Low = round(as.numeric(low), 3),
    High = round(as.numeric(high), 3),
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
    tbl$p_q3 <- round(pv$p[idx], 4)
    tbl$padj_q3 <- round(pv$padj[idx], 4)
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

  width_pct <- round(100 * cutoff_full$hdci_width, 1)
  caption <- paste0(
    "Q3 by item pair, sorted by departure from the expected range. ",
    "Expected range = ",
    width_pct,
    "% interval of the simulated Q3 per pair (",
    cutoff_full$actual_iterations,
    " iterations). Flagged: above = Q3 above ",
    "the upper bound (local dependence); below = below the lower bound."
  )
  if (p_value) {
    corr_label <- .correction_label(correction)
    caption <- paste0(
      caption,
      " p_q3/padj_q3: one-sided bootstrap p-values, ",
      corr_label,
      "; flagged at padj < ",
      alpha,
      "."
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

#' Simulation-Based \eqn{Q_3} Cutoff Determination
#'
#' Uses parametric bootstrap simulation to determine an appropriate cutoff
#' value for \code{\link{RMlocdepQ3}}. Under a correctly fitting Rasch model,
#' \eqn{Q_3} residuals have an unknown distribution; this function simulates
#' that distribution and returns empirical percentiles.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers).
#' @param iterations Integer. Number of simulation iterations (default 500).
#' @param parallel Logical. Use parallel processing via `mirai` if available
#'   (default `TRUE`).
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first. If neither is set and
#'   `parallel = TRUE`, a warning is issued and execution falls back to
#'   sequential (single core) processing.
#' @param verbose Logical. Show a progress bar (default `FALSE`).
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param cutoff_method Character. Method used to compute per-pair \eqn{Q_3}
#'   credible intervals in `pair_cutoffs`. One of `"hdci"` (the default,
#'   Highest Density Continuous Interval via `ggdist::hdci()`) or
#'   `"quantile"` (symmetric 2.5th / 97.5th percentiles). Only affects
#'   `pair_cutoffs`; the global `$suggested_cutoff` (99th percentile of
#'   `max(Q3) - mean(Q3)`) is unaffected.
#' @param hdci_width Numeric in (0, 1). Width of the HDCI when
#'   `cutoff_method = "hdci"`. Default `0.99`. Ignored when
#'   `cutoff_method = "quantile"`.
#' @param estimator Character. Estimation engine for the simulated \eqn{Q_3} values,
#'   passed through to the per-iteration computation. `"CML"` (default) uses
#'   CML item parameters and WLE person locations; `"MML"` uses `mirt`. This
#'   must match the `estimator` later given to \code{\link{RMlocdepQ3}}; the
#'   value is stored in the returned object's `$estimator` and reused
#'   automatically.
#' @param dgp Character. Data-generating process for the parametric bootstrap.
#'   `"resample"` (default) draws person locations by resampling the WLE
#'   estimates with replacement and simulates responses under the model -- a
#'   *marginal* null. `"conditional"` instead simulates each respondent's
#'   pattern from the exact Rasch conditional distribution given their observed
#'   total score (and answered items), with item parameters fixed -- a
#'   *conditional* null that fixes the score margin and needs no latent
#'   distribution, avoiding the over-dispersion of resampled point estimates.
#'   The two give different cut-offs; see the package's comparison study.
#'   \strong{Experimental.}
#'
#' @return A list with components:
#' \describe{
#'   \item{`results`}{data.frame with columns `mean`, `max`, `diff` (one row
#'     per successful iteration).}
#'   \item{`pair_results`}{Long data.frame with columns `Item1`, `Item2`,
#'     `Q3`, `iteration` --- one row per item pair per successful iteration.
#'     Used by \code{\link{RMlocdepQ3Plot}}.}
#'   \item{`pair_cutoffs`}{data.frame with per-pair cutoff summaries:
#'     `Item1`, `Item2`, `Q3_low`, `Q3_high`. Boundaries are computed via
#'     the method specified by `cutoff_method`.}
#'   \item{`actual_iterations`}{Number of successful iterations.}
#'   \item{`sample_n`}{Number of persons in the original data.}
#'   \item{`sample_n_total`}{Equal to `sample_n`: no respondents are dropped
#'     (incomplete responses are retained). Stored for consistency with the
#'     other `*Cutoff()` objects.}
#'   \item{`sample_has_na`}{Logical. Whether the data contained any missing
#'     values.}
#'   \item{`sample_summary`}{Summary statistics of estimated person parameters.}
#'   \item{`item_names`}{Character vector of item names from `data`.}
#'   \item{`max_diff`, `sd_diff`}{Max and SD of the `diff` distribution.}
#'   \item{`p95`, `p99`, `p995`, `p999`}{Empirical percentiles of `diff`.}
#'   \item{`suggested_cutoff`}{The 99th percentile (`p99`) --- recommended
#'     scalar cutoff for \code{\link{RMlocdepQ3}}.}
#'   \item{`cutoff_method`}{The method used for `pair_cutoffs`
#'     (`"hdci"` or `"quantile"`).}
#'   \item{`hdci_width`}{The HDCI width used (only meaningful when
#'     `cutoff_method = "hdci"`).}
#'   \item{`estimator`}{The estimator used for the simulated \eqn{Q_3} (`"CML"` or
#'     `"MML"`); reused by \code{\link{RMlocdepQ3}} and
#'     \code{\link{RMlocdepQ3Plot}}.}
#'   \item{`dgp`}{The data-generating process used (`"resample"` or
#'     `"conditional"`).}
#' }
#'
#' @details
#' The generating model is fitted once: CML item parameters (via
#' `psychotools`) and WLE person locations. For each simulation iteration,
#' those WLE thetas are resampled with replacement, response data are simulated
#' under the Rasch / Partial Credit model, the model is refitted, and \eqn{Q_3}
#' residuals are computed under `estimator`. The distribution of
#' `max(Q3) - mean(Q3)` across iterations provides empirical critical values.
#' Failed iterations (e.g., due to convergence issues) are silently discarded.
#'
#' Supports both **dichotomous** data (simulated via `psychotools::rrm()`) and
#' **polytomous** data (via an internal partial credit score simulator).
#'
#' Parallel processing is provided by the `mirai` package (optional). Install
#' it with `install.packages("mirai")` to enable parallelisation.
#'
#' @seealso \code{\link{RMlocdepQ3}}
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
#'   # Few iterations for a fast example; use 500+ in real analyses
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
  iterations = 500,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.99,
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

  data_mat <- as.matrix(data)
  sample_n <- nrow(data_mat)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  # No rows are dropped here (incomplete responses feed the per-pattern WLE
  # pool and the conditional DGP), so the raw total equals `sample_n`; record
  # the missingness flag so callers (e.g. RMlocdepQ3Plot) can report the
  # sample in the standard `n = X respondents (policy)` form.
  sample_has_na <- anyNA(data_mat)

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
    # Marginal DGP: resample WLE thetas with replacement, then simulate data
    # parametrically under the model (items and thetas on a common scale).
    sim_data_list$thetas <- wle_thetas
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
    stop("All simulation iterations failed. Check your data.", call. = FALSE)
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
  out$sample_n <- sample_n
  out$sample_n_total <- sample_n
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
#' @keywords internal
run_single_q3_sim <- function(seed, data_list) {
  set.seed(seed)

  tryCatch(
    {
      # --- Generate one simulated dataset under the chosen DGP -----------------
      if (identical(data_list$dgp, "conditional")) {
        sim_df <- .sim_cond_dataset(data_list)
      } else if (data_list$type == "dichotomous") {
        thetas_res <- sample(
          data_list$thetas,
          size = data_list$sample_n,
          replace = TRUE
        )
        sim_df <- as.data.frame(
          psychotools::rrm(
            theta = thetas_res,
            beta = data_list$item_params
          )$data
        )
      } else {
        thetas_res <- sample(
          data_list$thetas,
          size = data_list$sample_n,
          replace = TRUE
        )
        sim_df <- as.data.frame(sim_partial_score(
          data_list$deltaslist,
          thetas_res
        ))
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
#' @keywords internal
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
      # Engine helpers needed by .q3_residual_matrix() inside the daemon
      # (the CML/WLE path; the MML path uses namespaced mirt::).
      .q3_residual_matrix = .q3_residual_matrix,
      .rasch_std_residuals = .rasch_std_residuals,
      .fit_cml_thresholds = .fit_cml_thresholds,
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
#' @keywords internal
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
