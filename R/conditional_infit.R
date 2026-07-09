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
#'   * The return value of \code{\link{RMitemInfitCutoff}} (a list with `$item_cutoffs`):
#'     the data.frame is extracted automatically and the number of simulation
#'     iterations, `cutoff_method`, and `hdci_width` are included in the kable
#'     caption.
#'   * The `$item_cutoffs` data.frame from \code{\link{RMitemInfitCutoff}} directly: must
#'     have columns `Item`, `infit_low`, and `infit_high`.
#'   When provided, adds columns `Infit_low`, `Infit_high`, and `Flagged`
#'   to the result. `Flagged` labels the misfit direction: `"overfit"`
#'   (infit below the range -- more predictable than the model expects),
#'   `"underfit"` (above the range -- noisier than expected), or `""` (within
#'   range, no misfit).
#' @param p_value Logical. If `TRUE`, bootstrap p-values are computed from the
#'   simulated null distribution and added to the output. This requires
#'   `cutoff` to be the **full** \code{\link{RMitemInfitCutoff}} object (which
#'   carries the per-item simulated values in its `$results` element); the
#'   summarised `$item_cutoffs` data.frame is not sufficient. Default `FALSE`,
#'   in which case behaviour is unchanged.
#' @param correction Character. Multiple-comparison correction applied across
#'   items when `p_value = TRUE`: `"fwer"` (default) for the Westfall-Young
#'   studentised-max step-down (family-wise error rate), `"fdr_bh"` /
#'   `"fdr_by"` for Benjamini-Hochberg / Benjamini-Yekutieli false discovery
#'   rate control, or `"none"` for uncorrected per-item p-values.
#' @param alpha Numeric in (0, 1). Significance level for the `Flagged` column
#'   when `p_value = TRUE` (an item is flagged when its corrected p-value is
#'   below `alpha`). Default `0.05`.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying data.frame.
#' @param sort Optional character string. When `sort = "infit"`, rows are
#'   sorted by `Infit_MSQ` in descending order before output.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object (plain text table via
#'   `format = "pipe"`) with columns "Item", "Infit MSQ", and "Relative
#'   location", and a caption showing the number of complete cases. When
#'   `cutoff` is provided, columns "Infit low", "Infit high", and "Flagged"
#'   are also included, and the caption notes the simulation-based cutoffs.
#' * If `output = "dataframe"`: a data.frame with columns `Item`,
#'   `Infit_MSQ`, and `Relative_location`. When `cutoff` is provided, columns
#'   `Infit_low`, `Infit_high`, and `Flagged` are also included (inserted
#'   after `Infit_MSQ`, before `Relative_location`). `Flagged` is a character
#'   column with values `"overfit"`, `"underfit"`, or `""` (not the previous
#'   logical). When `p_value = TRUE`, columns `p_infit` (marginal two-sided
#'   p-value) and `padj_infit` (corrected p-value) are added and `Flagged`
#'   reflects items with `padj_infit < alpha` (direction from `Infit_MSQ`
#'   relative to 1).
#'
#' @details
#' Infit MSQ is a weighted fit statistic that emphasises deviations near the
#' item location. Values close to 1.0 indicate good fit. Values substantially
#' above 1.0 suggest underfit (unexpected responses), while values substantially
#' below 1.0 suggest overfit (overly predictable responses). The definition of
#' "substantially" depends on several factors such as sample size, and needs to
#' be determined by simulation using \code{\link{RMitemInfitCutoff}}. There is no general
#' rule-of-thumb value that is correct.
#'
#' Conditional infit MSQ statistics are computed via `iarm::out_infit()`, which
#' uses the conditional distribution of the sufficient statistics (Müller, 2020).
#' Only complete cases (rows without any `NA`) are used in the conditional fit
#' calculation.
#'
#' Item parameters are estimated by conditional maximum likelihood via
#' `psychotools::pcmodel()` (a dichotomous item is a 2-category PCM); the
#' conditional infit/outfit MSQ comes from `iarm::out_infit()` and is invariant
#' to the estimation engine. Per-item average locations are the means of the
#' CML thresholds, and the person-location reference is the mean of the Warm
#' WLE estimates.
#'
#' Relative item location is defined as the item's average location minus the
#' sample mean person location, providing a measure of item targeting.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' \strong{Bootstrap p-values.} When `p_value = TRUE`, each item's observed
#' infit is compared against its simulated null distribution (from
#' `cutoff$results`). The per-item statistic is the residual studentised by
#' the bootstrap mean and SD -- deliberately the empirical SD rather than the
#' Wilson-Hilferty / ZSTD transform, which is uninformative for conditional
#' MSQ (Müller, 2020). The marginal p-value is the two-sided Monte-Carlo
#' p-value `(1 + #{|t*| >= |t|}) / (B + 1)`. For `correction = "fwer"` the
#' family-wise adjustment uses the Westfall-Young studentised-max step-down,
#' which exploits the bootstrap dependence among items and is more powerful
#' than Bonferroni/Holm (Ferreira, 2024); its validity rests on subset
#' pivotality. `"fdr_bh"`/`"fdr_by"` apply Benjamini-Hochberg / Benjamini-
#' Yekutieli instead. These are model-conditional, sample-size-sensitive
#' p-values and are reported alongside the simulated effect-size band, not in
#' place of it. p-values can be no smaller than `1 / (B + 1)`, and the
#' studentised-max (FWER) correction is *liberal* when the simulation is small
#' (the bootstrap mean/SD used for studentisation are then too noisy). At
#' least 1000 `iterations` in [RMitemInfitCutoff()] are recommended -- in
#' simulations the family-wise error rate is then controlled at the nominal
#' level -- and a warning is issued when the simulation is smaller. The
#' marginal p-values are well calibrated even at a few hundred iterations.
#'
#' @section Multiple comparisons:
#' The marginal p-value controls the error rate of a *single* comparison: for
#' one item (or item pair) decided on in advance it is the relevant value. But
#' scanning all *k* comparisons and flagging whichever fall below `alpha` tests
#' *k* hypotheses at once, so the chance of at least one false flag inflates to
#' roughly \eqn{1 - (1 - \alpha)^k} (e.g. about 34% for *k* = 8 at
#' `alpha = 0.05`) -- even when every marginal p-value is correctly calibrated.
#' The corrected (adjusted) p-value controls this: `correction = "fwer"` bounds
#' the probability of *any* false flag (strict, lower power), while `"fdr_bh"` /
#' `"fdr_by"` bound the expected *proportion* of false flags among those raised
#' (a more lenient middle ground). Rule of thumb: use the marginal p-value for a
#' single pre-specified comparison, and a corrected p-value when screening the
#' whole table -- the usual workflow.
#'
#' @references
#' Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust
#' them? *Journal of Statistical Distributions and Applications*, 7(5).
#' \doi{10.1186/s40488-020-00108-7}
#'
#' Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate' number
#' of hypotheses simultaneously. *Journal of Statistical Theory and Practice,
#' 19*(6). \doi{10.1007/s42519-024-00412-4}
#'
#' Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple Testing*.
#' Wiley.
#'
#' @seealso \code{\link{RMitemInfitCutoff}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("iarm", quietly = TRUE)) {
#'   # Simulate binary item response data (5 items, 40 persons)
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:1, 40 * 5, replace = TRUE), nrow = 40, ncol = 5)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:5)
#'
#'   # Default kable output
#'   RMitemInfit(sim_data)
#'
#'   # Sorted by infit MSQ descending
#'   RMitemInfit(sim_data, sort = "infit")
#'
#'   # Return as data.frame for further processing
#'   df <- RMitemInfit(sim_data, output = "dataframe")
#'
#'   # Simulation-based cutoffs (100 Monte-Carlo iterations)
#'   if (requireNamespace("ggdist", quietly = TRUE)) {
#'     cutoff_res <- RMitemInfitCutoff(sim_data, iterations = 100,
#'                                     parallel = FALSE, seed = 42)
#'     RMitemInfit(sim_data, cutoff = cutoff_res)
#'     RMitemInfit(sim_data, cutoff = cutoff_res, output = "dataframe")
#'
#'     # Bootstrap p-values with family-wise (Westfall-Young) correction
#'     # (use iterations >= 1000 in real analyses for stable p-values)
#'     RMitemInfit(sim_data, cutoff = cutoff_res, p_value = TRUE,
#'                 output = "dataframe")
#'   }
#' }
#' }
RMitemInfit <- function(
  data,
  cutoff = NULL,
  p_value = FALSE,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable",
  sort
) {
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMitemInfit() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe"))
  correction <- match.arg(correction)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }

  # --- Validate and normalise cutoff ------------------------------------------
  cutoff_n_iter <- NULL
  cutoff_method <- NULL
  cutoff_hdci_width <- NULL
  cutoff_full <- NULL # full object (carries simulated $results for p-values)
  if (!is.null(cutoff)) {
    if (
      is.list(cutoff) &&
        !is.data.frame(cutoff) &&
        "item_cutoffs" %in% names(cutoff)
    ) {
      cutoff_full <- cutoff
      cutoff_n_iter <- cutoff$actual_iterations
      cutoff_method <- cutoff$cutoff_method
      cutoff_hdci_width <- cutoff$hdci_width
      cutoff <- cutoff$item_cutoffs
    }
    if (!is.data.frame(cutoff)) {
      stop(
        "`cutoff` must be NULL, the return value of RMitemInfitCutoff(), or its ",
        "$item_cutoffs data.frame.",
        call. = FALSE
      )
    }
    required_cols <- c("Item", "infit_low", "infit_high")
    missing_cols <- setdiff(required_cols, names(cutoff))
    if (length(missing_cols) > 0L) {
      stop(
        "`cutoff` data.frame is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  # --- p-value prerequisites --------------------------------------------------
  if (p_value) {
    if (is.null(cutoff_full) || is.null(cutoff_full$results)) {
      stop(
        "`p_value = TRUE` requires the full RMitemInfitCutoff() object (it ",
        "carries the simulated distributions in $results); a NULL cutoff or ",
        "the bare $item_cutoffs data.frame is not sufficient.",
        call. = FALSE
      )
    }
    if (!is.null(cutoff_n_iter) && cutoff_n_iter < 1000L) {
      warning(
        "Bootstrap p-values are based on only ",
        cutoff_n_iter,
        " simulation iterations. With few iterations the studentised-max ",
        "(FWER) correction is liberal and small p-values are imprecise; ",
        "use iterations >= 1000 in RMitemInfitCutoff() for reliable p-values.",
        call. = FALSE
      )
    }
  }

  validate_response_data(data)

  if (nrow(stats::na.omit(data)) == 0L) {
    stop(
      "No complete cases in data. All rows contain at least one NA.",
      call. = FALSE
    )
  }

  # Respondents with no responses at all contribute nothing and break the CML
  # fit (psychotools errors on all-NA rows); drop them, keeping the raw totals
  # for the caption.
  n_total <- nrow(as.data.frame(data))
  has_na <- anyNA(data)
  data <- .drop_empty_respondents(data)

  data_mat <- as.matrix(data)

  # --- Fit Rasch model and compute item/person locations ----------------------
  # CML item parameters (psychotools; a dichotomous item is a 2-category PCM)
  # and WLE person locations, consistent with the rest of the package. The
  # conditional infit/outfit statistic from iarm is engine-invariant; only the
  # relative-location reference shifts slightly (WLE vs eRm MLE person mean).
  fit <- psychotools::pcmodel(data)
  thr_list <- .center_thresholds(lapply(
    psychotools::threshpar(fit),
    as.numeric
  ))
  item_avg_locations <- vapply(thr_list, mean, numeric(1L))
  names(item_avg_locations) <- names(data)
  person_avg_location <- mean(
    .estimate_thetas(data_mat, thr_list, method = "WLE")$theta,
    na.rm = TRUE
  )

  relative_item_avg_locations <- item_avg_locations - person_avg_location

  # --- rgl.useNULL workaround -------------------------------------------------
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # --- Compute conditional infit MSQ via iarm ---------------------------------
  cfit <- iarm::out_infit(fit)

  # --- Count complete cases ---------------------------------------------------
  n_complete <- nrow(stats::na.omit(data))
  n_clause <- .n_caption(
    n_complete,
    n_total,
    if (has_na) "complete cases" else character()
  )

  # --- Assemble result data.frame (unrounded; kable rounds for display) -------
  item_fit_table <- data.frame(
    Item = names(data),
    Infit_MSQ = as.numeric(cfit$Infit),
    Relative_location = as.numeric(relative_item_avg_locations),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # --- Apply cutoff if provided ------------------------------------------------
  if (!is.null(cutoff)) {
    data_items <- item_fit_table$Item
    cutoff_items <- cutoff$Item
    if (!setequal(data_items, cutoff_items)) {
      stop(
        "Item names in `cutoff` do not match item names in `data`.\n",
        "  data items  : ",
        paste(data_items, collapse = ", "),
        "\n",
        "  cutoff items: ",
        paste(cutoff_items, collapse = ", "),
        call. = FALSE
      )
    }
    cutoff_sub <- cutoff[, c("Item", "infit_low", "infit_high")]
    item_fit_table <- merge(
      item_fit_table,
      cutoff_sub,
      by = "Item",
      sort = FALSE
    )
    # Restore original row order (merge may reorder)
    item_fit_table <- item_fit_table[match(data_items, item_fit_table$Item), ]
    rownames(item_fit_table) <- NULL
    item_fit_table$Infit_low <- item_fit_table$infit_low
    item_fit_table$Infit_high <- item_fit_table$infit_high
    item_fit_table$infit_low <- NULL
    item_fit_table$infit_high <- NULL
    # Flagged labels the misfit direction: infit below the expected range =
    # overfit (more predictable than the model expects), above = underfit
    # (noisier than expected); "" when within range.
    item_fit_table$Flagged <- ifelse(
      item_fit_table$Infit_MSQ < item_fit_table$Infit_low,
      "overfit",
      ifelse(
        item_fit_table$Infit_MSQ > item_fit_table$Infit_high,
        "underfit",
        ""
      )
    )

    if (p_value) {
      # Compare observed infit to its simulated null (cutoff_full$results).
      sim_items <- unique(cutoff_full$results$Item)
      if (!setequal(data_items, sim_items)) {
        stop(
          "Item names in the cutoff simulations ($results) do not match ",
          "`data`.",
          call. = FALSE
        )
      }
      sim_mat <- tapply(
        cutoff_full$results$InfitMSQ,
        list(cutoff_full$results$iteration, cutoff_full$results$Item),
        function(x) x[1L]
      )
      observed <- stats::setNames(as.numeric(cfit$Infit), names(data))
      pv <- .bootstrap_pvalues(observed, sim_mat, correction = correction)
      idx <- match(item_fit_table$Item, pv$name)
      item_fit_table$p_infit <- pv$p[idx]
      item_fit_table$padj_infit <- pv$padj[idx]
      sig <- item_fit_table$padj_infit < alpha
      item_fit_table$Flagged <- ifelse(
        is.na(sig) | !sig,
        "",
        ifelse(item_fit_table$Infit_MSQ > 1, "underfit", "overfit")
      )
      item_fit_table <- item_fit_table[, c(
        "Item",
        "Infit_MSQ",
        "Infit_low",
        "Infit_high",
        "p_infit",
        "padj_infit",
        "Flagged",
        "Relative_location"
      )]
    } else {
      # Reorder: Item, Infit_MSQ, Infit_low, Infit_high, Flagged, Relative_location
      item_fit_table <- item_fit_table[, c(
        "Item",
        "Infit_MSQ",
        "Infit_low",
        "Infit_high",
        "Flagged",
        "Relative_location"
      )]
    }
  }

  # --- Sort if requested ------------------------------------------------------
  if (!missing(sort) && identical(sort, "infit")) {
    item_fit_table <- item_fit_table[
      order(item_fit_table$Infit_MSQ, decreasing = TRUE),
    ]
    rownames(item_fit_table) <- NULL
  }

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(item_fit_table)
  }

  # Kable display rounding (the dataframe output above stays unrounded)
  item_fit_table <- .round_display(item_fit_table, c(
    Infit_MSQ = 3, Infit_low = 3, Infit_high = 3,
    p_infit = 4, padj_infit = 4, Relative_location = 2
  ))

  if (p_value) {
    kbl_colnames <- c(
      "Item",
      "Infit MSQ",
      "Infit low",
      "Infit high",
      "p",
      "p (adj)",
      "Flagged",
      "Relative location"
    )
    corr_label <- .correction_label(correction)
    kbl_caption <- paste0(
      "MSQ values based on conditional estimation. ",
      n_clause,
      ". Two-sided bootstrap p-values from ",
      cutoff_n_iter,
      " iterations; multiplicity correction: ",
      corr_label,
      "; flagged at alpha = ",
      alpha,
      ". p-values cannot be smaller than ",
      "1/(",
      cutoff_n_iter,
      "+1) = ",
      round(1 / (cutoff_n_iter + 1), 4),
      ".",
      " Flagged: underfit (infit > 1, noisier) / overfit (infit < 1, more predictable)."
    )
  } else if (is.null(cutoff)) {
    kbl_colnames <- c("Item", "Infit MSQ", "Relative location")
    kbl_caption <- paste0(
      "MSQ values based on conditional estimation. ",
      n_clause,
      "."
    )
  } else {
    kbl_colnames <- c(
      "Item",
      "Infit MSQ",
      "Infit low",
      "Infit high",
      "Flagged",
      "Relative location"
    )
    if (!is.null(cutoff_n_iter)) {
      method_label <- .format_cutoff_method_label(
        cutoff_method,
        cutoff_hdci_width
      )
      iter_part <- paste0(cutoff_n_iter, " simulation iterations")
      kbl_caption <- paste0(
        "MSQ values based on conditional estimation. ",
        n_clause,
        ". Cutoff values based on ",
        if (!is.null(method_label)) {
          paste0(iter_part, " (", method_label, ").")
        } else {
          paste0(iter_part, ".")
        }
      )
    } else {
      kbl_caption <- paste0(
        "MSQ values based on conditional estimation. ",
        n_clause,
        ". Simulation-based cutoff values applied."
      )
    }
    kbl_caption <- paste0(
      kbl_caption,
      " Flagged: overfit = infit below range (more predictable); ",
      "underfit = above range (noisier)."
    )
  }

  knitr::kable(
    item_fit_table,
    format = "pipe",
    col.names = kbl_colnames,
    caption = kbl_caption
  )
}

# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

#' Format a human-readable label for the cutoff method
#'
#' @param cutoff_method Character. `"hdci"`, `"quantile"`, or `NULL`.
#' @param hdci_width Numeric or `NULL`. HDCI width (e.g., `0.999`).
#' @return A character label, or `NULL` if the method is unknown/unset.
#' @keywords internal
.format_cutoff_method_label <- function(cutoff_method, hdci_width) {
  if (is.null(cutoff_method)) {
    return(NULL)
  }
  switch(
    cutoff_method,
    quantile = "2.5th/97.5th percentile",
    hdci = if (!is.null(hdci_width)) {
      paste0(hdci_width * 100, "% HDCI")
    } else {
      "HDCI"
    },
    NULL
  )
}
