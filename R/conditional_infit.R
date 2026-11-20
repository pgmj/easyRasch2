#' Conditional item infit MSQ
#'
#' Computes the conditional infit mean-square (MSQ) of each item with
#' `iarm::out_infit()`, and each item's location relative to the mean
#' person location. Given the result of [RMitemInfitCutoff()], it also
#' flags misfitting items, by default on a bootstrap p-value corrected
#' across items.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Missing values (`NA`) are allowed, but the infit statistic uses
#'   complete cases only (see Details), so at least one is needed.
#' @param cutoff Optional. Default `NULL`, no flagging. Either of:
#'   * The full result of [RMitemInfitCutoff()]. This enables the p-values
#'     (see `p_value`), and the caption reports the number of iterations
#'     and the interval method. A result of [RMitemInfitCutoffMI()] gives
#'     an error, since its simulations come from imputed data. Use it with
#'     [RMitemInfitMI()].
#'   * Its `$item_cutoffs` data.frame alone, with columns `Item`,
#'     `infit_low` and `infit_high` (or `outfit_low` and `outfit_high` with
#'     `statistic = "outfit"`). Items are then flagged against the interval.
#'
#'   Either adds the columns `Infit_low`, `Infit_high` and `Flagged`.
#'   `Flagged` gives the direction of misfit: `"overfit"` (below 1, or
#'   below the interval), `"underfit"` (above) or `""` (no misfit).
#' @param p_value Logical or `NULL`. Whether to flag on bootstrap p-values
#'   computed from the simulated values in `cutoff`. `NULL` (default) uses
#'   them whenever `cutoff` is the full [RMitemInfitCutoff()] result, and
#'   flags against the interval otherwise. `TRUE` without the full result is
#'   an error.
#'
#'   With `FALSE`, items are flagged against the interval and a one-time
#'   message gives the family-wise error rate that implies. The interval
#'   describes where a fitting item's statistic is expected to fall. Used
#'   as a decision rule over `k` items, a width `w` gives a family-wise
#'   error rate of `1 - w^k`, 37% for the default 95% interval over nine
#'   items. The corrected p-value keeps the error rate at `alpha` and needs
#'   far fewer iterations to do so (Johansson, 2026). Versions up to 1.1.1
#'   flagged against the interval by default, so older scripts that pass
#'   the full cutoff object may now flag different items.
#' @param correction Character. Correction of the p-values across items
#'   when they are used: `"fwer"` (default, Westfall-Young step-down,
#'   controls the family-wise error rate), `"fdr_bh"` or `"fdr_by"`
#'   (Benjamini-Hochberg or Benjamini-Yekutieli, control the false
#'   discovery rate), or `"none"`.
#' @param alpha Numeric in (0, 1). An item is flagged when its corrected
#'   p-value is below `alpha`. Default `0.05`. Used only with p-values.
#' @param output Character. `"kable"` (default) for a `knitr::kable()`
#'   table, or `"dataframe"`.
#' @param sort Optional. `sort = "infit"` (or `"outfit"`, matching
#'   `statistic`) sorts the rows by MSQ, largest first. Other values leave
#'   the order unchanged.
#' @param statistic Character. Which conditional fit statistic to report
#'   and test: `"infit"` (default) or `"outfit"`. Both come from the same
#'   `iarm::out_infit()` call and the same simulations, and the column
#'   names, interval, p-values and `Flagged` follow the choice. The
#'   correction is applied across the k items for the chosen statistic.
#'   Testing both statistics and reporting whichever flags is a family of
#'   2k tests, larger than either call corrects for.
#'
#' @return
#' * `output = "kable"`: a `knitr_kable` table with the columns "Item",
#'   "Infit MSQ" and "Relative location", and a caption with the sample
#'   size. With `cutoff`, also "Infit low", "Infit high" and "Flagged",
#'   and with p-values "p" and "p (adj)".
#' * `output = "dataframe"`: a data.frame with the columns `Item`,
#'   `Infit_MSQ` and `Relative_location`. With `cutoff`, also
#'   `Infit_low`, `Infit_high` and `Flagged` (after `Infit_MSQ`). With
#'   p-values, also `p_infit` (uncorrected two-sided p) and `padj_infit`
#'   (corrected p). `Flagged` then marks items with `padj_infit < alpha`,
#'   with the direction from `Infit_MSQ` relative to 1.
#'
#' With `statistic = "outfit"` the columns are named after it:
#' `Outfit_MSQ`, `Outfit_low`, `Outfit_high`, `p_outfit`, `padj_outfit`.
#'
#' @details
#' **Statistic.** Infit MSQ is a weighted fit statistic that emphasises
#' unexpected responses near the item's location. Values near 1 indicate
#' fit. Values above 1 indicate underfit (responses less predictable than
#' the model expects), and values below 1 overfit (responses more
#' predictable). How far from 1 is too far depends on the sample size and
#' the items, so it is judged against simulations from
#' [RMitemInfitCutoff()]. No fixed value is right in general. The
#' conditional statistic from `iarm::out_infit()` uses the distribution of
#' the responses given the total score (Müller, 2020).
#'
#' **Data and model.** The item parameters are estimated by conditional
#' maximum likelihood with `psychotools::pcmodel()` (a dichotomous item is
#' a 2-category partial credit model). The conditional infit uses complete
#' cases only, and every response category must be used in them, or the
#' function stops and names the item. The relative item location is the
#' mean of the item's CML thresholds minus the mean WLE person location,
#' a description of targeting. The `iarm` package must be installed (it is
#' in Suggests).
#'
#' **Bootstrap p-values.** Each item's observed statistic is compared with
#' its simulated values from `cutoff$results`. The test statistic is the
#' deviation from the mean of the observed and simulated values together,
#' divided by their SD. The ZSTD (Wilson-Hilferty) transform is not used,
#' since it does not work for the conditional MSQ (Müller, 2020). The
#' uncorrected p is two-sided, `(1 + #{|t*| >= |t|}) / (B + 1)` for B
#' simulated datasets, so it can be no smaller than `1 / (B + 1)`. With
#' `"fwer"` the correction is the Westfall-Young step-down on the
#' studentised maximum, which uses the dependence between the items'
#' statistics and has more power than Bonferroni or Holm (Westfall &
#' Young, 1993; Ferreira, 2024). Its validity rests on subset pivotality.
#' `"fdr_bh"` and `"fdr_by"` apply Benjamini-Hochberg or
#' Benjamini-Yekutieli. The p-values depend on the model and the sample
#' size, and are reported with the simulated interval, not instead of it.
#'
#' **Iterations.** The observed value is part of the mean and SD it is
#' studentised by, so the observed and simulated statistics are
#' exchangeable under the model, and the family-wise error rate is at the
#' nominal level at any number of iterations in [RMitemInfitCutoff()].
#' More iterations make decisions reproducible. Two analysts using
#' different seeds disagree about at least one item roughly 10% of the
#' time at 400 iterations and 4% at 2000 (Johansson, 2026). A message says
#' so below 400 iterations, and the caption below 1000. Use 1000 to 2000
#' for a final analysis. The uncorrected p-values are well calibrated even
#' at a few hundred iterations.
#'
#' @section Multiple comparisons:
#' The uncorrected p-value controls the error rate of a single comparison.
#' For one item (or item pair) chosen in advance it is the relevant value.
#' Flagging whichever of k comparisons fall below `alpha` tests k
#' hypotheses at once, and the chance of at least one false flag rises to
#' about \eqn{1 - (1 - \alpha)^k}, for example 34% for k = 8 at
#' `alpha = 0.05`, even when every uncorrected p-value is calibrated. The
#' corrected p-value controls this. `correction = "fwer"` bounds the
#' probability of any false flag (strict, less power), and `"fdr_bh"` and
#' `"fdr_by"` bound the expected share of false flags among those raised.
#' Use the uncorrected p-value for a single comparison chosen in advance,
#' and a corrected one when screening the whole table, which is the usual
#' case.
#'
#' @section Flags depend on the other items:
#' Item infit and item-restscore both judge each item against expectations
#' from a Rasch model fitted to all items, so misfitting items change what
#' is expected of the others. A flag in one direction can therefore
#' produce flags in the other. Overfitting items can push fitting items
#' into the underfit range, and underfitting items, for example from a
#' second dimension, can make fitting items look overfit. How large the
#' effect is depends on how strong the misfit is and how many items share
#' it. In simulation with seven polytomous items, two overfitting items led
#' to about 6% of the fitting items being flagged as underfit, and one or
#' two mildly multidimensional items to 1% to 2% being flagged as overfit.
#' When items are flagged in both directions, consider removing the
#' clearest misfit and testing again before interpreting the rest.
#'
#' @references
#' Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate' number
#' of hypotheses simultaneously. *Journal of Statistical Theory and Practice,
#' 19*(1), 6. \doi{10.1007/s42519-024-00412-4}
#'
#' Johansson, M. (2026). Simulation-based cutoffs for conditional item fit in
#' Rasch models: Iterations, multiplicity correction, and decision stability.
#' *PsyArXiv*. \doi{10.31234/osf.io/7pqz4_v2}
#'
#' Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust
#' them? *Journal of Statistical Distributions and Applications*, 7(5).
#' \doi{10.1186/s40488-020-00108-7}
#'
#' Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple Testing*.
#' Wiley.
#'
#' @seealso [RMitemInfitCutoff()] for the simulations, [RMitemInfitPlot()]
#'   to plot them with the observed values, [RMitemInfitMI()] for multiply
#'   imputed data, and [RMitemRestscore()] for a second item fit test.
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
#'
#'     # The same for conditional outfit MSQ
#'     RMitemInfit(sim_data, cutoff = cutoff_res, p_value = TRUE,
#'                 statistic = "outfit", output = "dataframe")
#'   }
#' }
#' }
RMitemInfit <- function(
  data,
  cutoff = NULL,
  p_value = NULL,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable",
  sort,
  statistic = "infit"
) {
  if (
    !is.null(p_value) &&
      (!is.logical(p_value) || length(p_value) != 1L || is.na(p_value))
  ) {
    stop("`p_value` must be NULL, TRUE or FALSE.", call. = FALSE)
  }
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMitemInfit() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe"))
  correction <- match.arg(correction)
  statistic <- match.arg(statistic, c("infit", "outfit"))
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }

  # Column names follow the reported statistic. Infit and outfit come from the
  # same iarm::out_infit() call and the same bootstrap, so only the naming and
  # the source column differ.
  stat_title <- if (statistic == "outfit") "Outfit" else "Infit"
  msq_col <- paste0(stat_title, "_MSQ") # Infit_MSQ / Outfit_MSQ
  low_col <- paste0(stat_title, "_low")
  high_col <- paste0(stat_title, "_high")
  cut_low <- paste0(statistic, "_low") # infit_low / outfit_low (in $item_cutoffs)
  cut_high <- paste0(statistic, "_high")
  p_col <- paste0("p_", statistic)
  padj_col <- paste0("padj_", statistic)
  sim_col <- paste0(stat_title, "MSQ") # InfitMSQ / OutfitMSQ (in $results)

  # --- Validate and normalise cutoff ------------------------------------------
  cutoff_n_iter <- NULL
  cutoff_req_iter <- NULL
  cutoff_method <- NULL
  cutoff_hdci_width <- NULL
  cutoff_full <- NULL # full object (carries simulated $results for p-values)
  if (!is.null(cutoff)) {
    if (
      is.list(cutoff) &&
        !is.data.frame(cutoff) &&
        "item_cutoffs" %in% names(cutoff)
    ) {
      # An RMitemInfitCutoffMI() result stacks simulations from imputed
      # datasets. Comparing complete-case infit with them is not calibrated
      # (Johansson, 2026, covers complete data only), so refuse it.
      if (!is.null(cutoff$n_imputations)) {
        stop(
          "`cutoff` is a result of RMitemInfitCutoffMI(), made from imputed ",
          "data. Use it with RMitemInfitMI(), or use RMitemInfitCutoff() on ",
          "the data passed to RMitemInfit().",
          call. = FALSE
        )
      }
      cutoff_full <- cutoff
      cutoff_n_iter <- cutoff$actual_iterations
      # NULL for cutoff objects made before 1.2.0 stored the requested count
      cutoff_req_iter <- cutoff$requested_iterations
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
    required_cols <- c("Item", cut_low, cut_high)
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

  # --- Resolve p_value --------------------------------------------------------
  # NULL means "use the corrected p-value when the simulations are available".
  # The interval is a description of where a fitting item's statistic is
  # expected to fall; it makes a poor decision rule, because its width sets a
  # family-wise error rate implicitly (Johansson, 2026). Passing the bare
  # $item_cutoffs data.frame, or no cutoff at all, leaves nothing to compute a
  # p-value from, so those resolve to FALSE.
  have_sims <- !is.null(cutoff_full) && !is.null(cutoff_full$results)
  if (is.null(p_value)) {
    p_value <- have_sims
  }

  # --- p-value prerequisites --------------------------------------------------
  if (p_value) {
    if (!have_sims) {
      stop(
        "`p_value = TRUE` requires the full RMitemInfitCutoff() object (it ",
        "carries the simulated distributions in $results); a NULL cutoff or ",
        "the bare $item_cutoffs data.frame is not sufficient.",
        call. = FALSE
      )
    }
    # Below 400 decisions depend on the seed enough to say so on the
    # console. Between 400 and 1000 the table caption reports it instead.
    if (!is.null(cutoff_n_iter) && cutoff_n_iter < 400L) {
      .notify_low_iterations(cutoff_n_iter, cutoff_req_iter,
                             correction = correction)
    }
  } else if (!is.null(cutoff)) {
    .notify_band_flagging(
      if (identical(cutoff_method, "quantile")) 0.95 else cutoff_hdci_width,
      nrow(cutoff)
    )
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
  # iarm's conditional infit refits on complete cases, so check them too.
  .check_complete_case_categories(data)

  data_mat <- as.matrix(data)

  # --- Fit Rasch model and compute item/person locations ----------------------
  # CML item parameters (psychotools; a dichotomous item is a 2-category PCM)
  # and WLE person locations, consistent with the rest of the package. The
  # conditional infit/outfit statistic from iarm is engine-invariant; only the
  # relative-location reference shifts slightly (WLE vs eRm MLE person mean).
  fit <- .cml_fit(data) # shared category check, then pcmodel()
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
  if (!is.null(cutoff_full)) {
    .check_cutoff_sample(cutoff_full$sample_n, n_complete, "RMitemInfitCutoff()")
  }
  n_clause <- .n_caption(
    n_complete,
    n_total,
    if (has_na) "complete cases" else character()
  )

  # --- Assemble result data.frame (unrounded; kable rounds for display) -------
  observed_msq <- as.numeric(if (statistic == "outfit") cfit$Outfit else cfit$Infit)
  item_fit_table <- data.frame(
    Item = names(data),
    MSQ = observed_msq,
    Relative_location = as.numeric(relative_item_avg_locations),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  names(item_fit_table)[names(item_fit_table) == "MSQ"] <- msq_col

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
    cutoff_sub <- cutoff[, c("Item", cut_low, cut_high)]
    item_fit_table <- merge(
      item_fit_table,
      cutoff_sub,
      by = "Item",
      sort = FALSE
    )
    # Restore original row order (merge may reorder)
    item_fit_table <- item_fit_table[match(data_items, item_fit_table$Item), ]
    rownames(item_fit_table) <- NULL
    item_fit_table[[low_col]] <- item_fit_table[[cut_low]]
    item_fit_table[[high_col]] <- item_fit_table[[cut_high]]
    item_fit_table[[cut_low]] <- NULL
    item_fit_table[[cut_high]] <- NULL
    # Flagged labels the misfit direction: MSQ below the expected range =
    # overfit, above = underfit
    # ; "" when within range.
    item_fit_table$Flagged <- ifelse(
      item_fit_table[[msq_col]] < item_fit_table[[low_col]],
      "overfit",
      ifelse(
        item_fit_table[[msq_col]] > item_fit_table[[high_col]],
        "underfit",
        ""
      )
    )

    if (p_value) {
      # Compare the observed statistic to its simulated null
      # (cutoff_full$results).
      sim_items <- unique(cutoff_full$results$Item)
      if (!setequal(data_items, sim_items)) {
        stop(
          "Item names in the cutoff simulations ($results) do not match ",
          "`data`.",
          call. = FALSE
        )
      }
      if (!sim_col %in% names(cutoff_full$results)) {
        stop(
          "The cutoff simulations ($results) have no `",
          sim_col,
          "` column, which `statistic = \"",
          statistic,
          "\"` requires.",
          call. = FALSE
        )
      }
      sim_mat <- tapply(
        cutoff_full$results[[sim_col]],
        list(cutoff_full$results$iteration, cutoff_full$results$Item),
        function(x) x[1L]
      )
      observed <- stats::setNames(observed_msq, names(data))
      pv <- .bootstrap_pvalues(observed, sim_mat, correction = correction)
      idx <- match(item_fit_table$Item, pv$name)
      item_fit_table[[p_col]] <- pv$p[idx]
      item_fit_table[[padj_col]] <- pv$padj[idx]
      sig <- item_fit_table[[padj_col]] < alpha
      item_fit_table$Flagged <- ifelse(
        is.na(sig) | !sig,
        "",
        ifelse(item_fit_table[[msq_col]] > 1, "underfit", "overfit")
      )
      item_fit_table <- item_fit_table[, c(
        "Item",
        msq_col,
        low_col,
        high_col,
        p_col,
        padj_col,
        "Flagged",
        "Relative_location"
      )]
    } else {
      # Reorder: Item, <stat>_MSQ, <stat>_low, <stat>_high, Flagged,
      # Relative_location
      item_fit_table <- item_fit_table[, c(
        "Item",
        msq_col,
        low_col,
        high_col,
        "Flagged",
        "Relative_location"
      )]
    }
  }

  # --- Sort if requested ------------------------------------------------------
  if (!missing(sort) && identical(sort, statistic)) {
    item_fit_table <- item_fit_table[
      order(item_fit_table[[msq_col]], decreasing = TRUE),
    ]
    rownames(item_fit_table) <- NULL
  }

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(item_fit_table)
  }

  # Kable display rounding (the dataframe output above stays unrounded)
  item_fit_table <- .round_display(
    item_fit_table,
    stats::setNames(
      c(3, 3, 3, 4, 4, 2),
      c(msq_col, low_col, high_col, p_col, padj_col, "Relative_location")
    )
  )

  if (p_value) {
    kbl_colnames <- c(
      "Item",
      paste(stat_title, "MSQ"),
      paste(stat_title, "low"),
      paste(stat_title, "high"),
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
      " iterations, multiplicity correction: ",
      corr_label,
      ". Flagged on the corrected p-value at alpha = ",
      alpha,
      ". p-values cannot be smaller than ",
      "1/(",
      cutoff_n_iter,
      "+1) = ",
      .fmt_small(1 / (cutoff_n_iter + 1)),
      ".",
      # Two tiers about reproducibility. What the error rate claim says
      # follows `correction` (.error_rate_sentence()). Below 400 the console
      # also reports it. Between 400 and 1000 it is a caption matter rather
      # than something to interrupt over.
      if (is.null(cutoff_n_iter)) {
        ""
      } else if (cutoff_n_iter < 400L) {
        paste0(
          " This is below the recommended 400. ",
          .error_rate_sentence(correction, "any"),
          " Decisions depend on the seed, so use 1000 to 2000 for a final",
          " analysis (Johansson, 2026)."
        )
      } else if (cutoff_n_iter < 1000L) {
        paste0(
          " ", .error_rate_sentence(correction, "this"),
          " Decisions are still somewhat seed-dependent. Two analysts using",
          " different seeds disagree about at least one item roughly 10% of the time at",
          " 400 iterations against 4% at 2000, so use 1000 to 2000 for a final",
          " analysis (Johansson, 2026)."
        )
      } else {
        ""
      },
      .attrition_clause(cutoff_n_iter, cutoff_req_iter),
      " Flagged: underfit (",
      statistic,
      " > 1) / overfit (",
      statistic,
      " < 1)."
    )
  } else if (is.null(cutoff)) {
    kbl_colnames <- c("Item", paste(stat_title, "MSQ"), "Relative location")
    kbl_caption <- paste0(
      "MSQ values based on conditional estimation. ",
      n_clause,
      "."
    )
  } else {
    kbl_colnames <- c(
      "Item",
      paste(stat_title, "MSQ"),
      paste(stat_title, "low"),
      paste(stat_title, "high"),
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
      " Flagged: overfit = ",
      statistic,
      " below range, ",
      "underfit = above range.",
      .band_error_clause(
        if (identical(cutoff_method, "quantile")) 0.95 else cutoff_hdci_width,
        nrow(item_fit_table)
      ),
      .attrition_clause(cutoff_n_iter, cutoff_req_iter)
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
#' @noRd
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
