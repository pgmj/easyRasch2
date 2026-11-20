#' Item-restscore association
#'
#' Computes, for each item, the observed and model-expected association
#' between the item and the rest score (the total score on the other
#' items), measured by Goodman-Kruskal's gamma with `iarm::item_restscore()`.
#' Reports their signed difference, the item's location relative to the mean
#' person location, and a misfit flag. With the result of
#' [RMitemRestscoreCutoff()], items are flagged on a bootstrap p-value
#' corrected across items, which is recommended over the default
#' asymptotic test (see Details).
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Missing values (`NA`) are allowed, but the statistic uses complete
#'   cases only, so at least one is needed.
#' @param cutoff Optional. Default `NULL`, which tests each item with the
#'   asymptotic normal test of `iarm::item_restscore()`, adjusted by
#'   `p_adj`. Either of:
#'   * The full result of [RMitemRestscoreCutoff()]. This enables the
#'     bootstrap p-values (see `p_value`) and adds the columns `Diff_low`
#'     and `Diff_high`, the interval of the simulated differences.
#'   * Its `$item_cutoffs` data.frame alone (columns `Item`, `diff_low`,
#'     `diff_high`). Items are then flagged against the interval.
#'
#'   With a `cutoff`, the asymptotic `p_adjusted` column is dropped and
#'   `p_adj` is ignored.
#' @param p_value Logical or `NULL`. Whether to flag on bootstrap p-values
#'   from the simulated values in `cutoff`. `NULL` (default) uses them
#'   whenever `cutoff` is the full [RMitemRestscoreCutoff()] result.
#'   `TRUE` without it is an error. With `FALSE` and a `cutoff`, items are
#'   flagged against the interval, and a one-time message gives the
#'   family-wise error rate that implies.
#' @param correction Character. Correction of the bootstrap p-values across
#'   items: `"fwer"` (default, Westfall-Young step-down, controls the
#'   family-wise error rate), `"fdr_bh"` or `"fdr_by"` (Benjamini-Hochberg
#'   or Benjamini-Yekutieli, control the false discovery rate), or
#'   `"none"`.
#' @param alpha Numeric in (0, 1). An item is flagged when its corrected
#'   bootstrap p-value is below `alpha`. Default `0.05`. The asymptotic
#'   test without a `cutoff` always flags at .05.
#' @param output Character. `"kable"` (default) for a `knitr::kable()`
#'   table, or `"dataframe"`.
#' @param sort Optional. `sort = "diff"` sorts the rows by the absolute
#'   `Difference`, largest first, so that over- and underfitting items both
#'   come first. Other values leave the order unchanged.
#' @param p_adj Character. Adjustment of the asymptotic p-values across
#'   items, passed to `iarm::item_restscore()`. Default `"BH"`
#'   (Benjamini-Hochberg). `"none"` gives unadjusted p-values, and any
#'   method of [stats::p.adjust()] can be used. Only used when
#'   `cutoff = NULL`.
#'
#' @return
#' * `output = "kable"`: a `knitr_kable` table with the item name, the
#'   observed and expected gamma, their difference, the adjusted p-value,
#'   the `Flagged` label and the relative item location.
#' * `output = "dataframe"`: a data.frame with the columns `Item`,
#'   `Observed`, `Expected`, `Difference`, `p_adjusted`, `Flagged` and
#'   `Relative_location`. `Flagged` is `"overfit"` (observed above
#'   expected, adjusted p < .05), `"underfit"` (below, adjusted p < .05) or
#'   `""`.
#'
#' With a `cutoff`, `p_adjusted` is replaced by `Diff_low` and `Diff_high`,
#' followed by `p_restscore` (uncorrected equal-tailed bootstrap p) and
#' `padj_restscore` (corrected p) when bootstrap p-values are used.
#' `Flagged` then marks `padj_restscore < alpha`, or values outside the
#' interval with `p_value = FALSE`.
#'
#' `Difference` is observed minus expected gamma. Positive values mean the
#' item is more strongly associated with the rest score than the model
#' predicts (overfit, often from local dependence). Negative values mean a
#' weaker association (underfit, often from multidimensionality or noise).
#'
#' @details
#' **Statistic.** The item-restscore gamma (Kreiner, 2011) measures the
#' association between the responses to one item and the rest score.
#' Under the Rasch model the observed and expected gamma agree, apart from
#' sampling error. The statistic is conditional on the total score.
#'
#' **Data and model.** The item parameters are estimated by conditional
#' maximum likelihood with `psychotools::pcmodel()` (a dichotomous item is
#' a 2-category partial credit model). The statistic uses complete cases,
#' and every response category must be used in them, or the function stops
#' and names the item. The relative item location is the mean of the
#' item's CML thresholds minus the mean WLE person location, a description
#' of targeting. The \pkg{iarm} package must be installed (it is in
#' Suggests).
#'
#' **The asymptotic p-value is miscalibrated.** Without a `cutoff`, each
#' item is tested with `(observed - expected) / SE` against the standard
#' normal, where the SE is that of the observed gamma alone. The expected
#' gamma is estimated from the same data and correlates with the observed
#' one, so the SE is too large for the difference. The difference is also
#' biased upwards in small samples. When the Rasch model holds, the test
#' therefore flags too many items as overfit in small samples and too few
#' as underfit at every sample size, and the adjustment in `p_adj` corrects
#' neither. In simulation, 20 dichotomous items at n = 150 gave at least one
#' BH-flagged item in 13% of datasets (30% when 1.5 logits off target),
#' almost all of them overfit, while 9 polytomous items at n = 1000 gave
#' about 1%. Pass the result of [RMitemRestscoreCutoff()] as `cutoff` for
#' p-values from a parametric bootstrap instead. The simulation study is
#' reported at
#' \url{https://pgmj.github.io/rasch_book/chapters/restscore-asymptotic-null.html}.
#'
#' **Bootstrap p-values.** Each item's observed `Difference` is compared
#' with its simulated values from `cutoff$results`, studentised by the mean
#' and SD of the observed and simulated values together. Subtracting that
#' mean removes the small-sample bias of the difference, and dividing by
#' that SD accounts for the correlation between the observed and expected
#' gamma.
#'
#' The test is equal-tailed: overfit and underfit each get `alpha / 2`.
#' Gamma is bounded at 1, so the simulated `Difference` is skewed to the
#' left, and a two-sided test on `|t|` would reject more often in the long
#' (underfit) tail. The uncorrected p is twice the smaller of the one-sided
#' Monte Carlo p-values, `(1 + #{t* >= t}) / (B + 1)` and
#' `(1 + #{t* <= t}) / (B + 1)` for B simulated datasets, so it cannot be
#' smaller than `2 / (B + 1)`. With `correction = "fwer"`, the
#' Westfall-Young step-down on the studentised maximum (Ferreira, 2024) is
#' run once in each direction, and `padj_restscore` is twice the smaller
#' of the two corrected p-values. `"fdr_bh"` and `"fdr_by"` adjust the
#' equal-tailed uncorrected p-values.
#'
#' In simulation with the Rasch model true, over four conditions (20
#' dichotomous or 9 polytomous items, n = 150 to 1000, on or off target),
#' `correction = "fwer"` flagged at least one item in 5.0% of datasets,
#' 2.4% underfit and 2.75% overfit. A test on `|t|` also held 5.0%, but
#' split it 3.7% and 1.6%. The equal-tailed test also gave fewer false
#' underfit flags when other items overfit (4.7% against 6.2% with two
#' overfitting items), at a cost of 2 to 3 percentage points of power to
#' detect underfitting items. Before easyRasch2 1.4.0.9004 the test was on
#' `|t|`.
#'
#' **Reproducibility.** The bootstrap p-values depend on the simulations, so
#' two analyses of the same data with different seeds can disagree about
#' items near the decision boundary. In simulation, in conditions chosen to
#' include such items, two seeds disagreed about the flag of at least one
#' item in about 10% of analyses at 400 iterations and 6% at 1000 with
#' `correction = "fwer"`, and in 13% and 11% with `"fdr_bh"`. Set `seed` in
#' [RMitemRestscoreCutoff()] for a reproducible analysis, and use 1000 or
#' more iterations for a final one.
#'
#' The direction in `Flagged` follows the studentised value, not the sign
#' of `Difference`. At small sample sizes the simulated `Difference` has a
#' slightly positive mean, so an item can, rarely, be flagged `"underfit"`
#' while its `Difference` is just above zero.
#'
#' @inheritSection RMitemInfit Multiple comparisons
#' @inheritSection RMitemInfit Flags depend on the other items
#'
#' @references
#' Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate' number
#' of hypotheses simultaneously. *Journal of Statistical Theory and Practice,
#' 19*(1), 6. \doi{10.1007/s42519-024-00412-4}
#'
#' Kreiner, S. (2011). A Note on Item–Restscore Association in Rasch Models.
#' *Applied Psychological Measurement, 35*(7), 557–561.
#' \doi{10.1177/0146621611410227}
#'
#' @seealso [RMitemRestscoreCutoff()] for the simulations,
#'   [RMitemRestscorePlot()] to plot them, [RMitemRestscoreBoot()] for
#'   large samples, and [RMitemInfit()] for a second item fit test.
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("iarm", quietly = TRUE)) {
#'   # Simulate binary item response data (8 items, 200 persons)
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:1, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:8)
#'
#'   # Default kable output
#'   RMitemRestscore(sim_data)
#'
#'   # Sorted by absolute difference
#'   RMitemRestscore(sim_data, sort = "diff")
#'
#'   # Return as data.frame for further processing
#'   df <- RMitemRestscore(sim_data, output = "dataframe")
#'
#'   # Bootstrap null distribution, flagging on Westfall-Young corrected
#'   # p-values (use 1000 or more iterations in a final analysis)
#'   if (requireNamespace("ggdist", quietly = TRUE)) {
#'     cutoff_res <- RMitemRestscoreCutoff(sim_data, iterations = 100,
#'                                         parallel = FALSE, seed = 42)
#'     RMitemRestscore(sim_data, cutoff = cutoff_res)
#'   }
#' }
#' }
RMitemRestscore <- function(
  data,
  cutoff = NULL,
  p_value = NULL,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable",
  sort,
  p_adj = "BH"
) {
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMitemRestscore() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe"))
  correction <- match.arg(correction)
  if (
    !is.null(p_value) &&
      (!is.logical(p_value) || length(p_value) != 1L || is.na(p_value))
  ) {
    stop("`p_value` must be NULL, TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }

  # --- Validate and normalise cutoff ------------------------------------------
  cutoff_n_iter <- NULL
  cutoff_req_iter <- NULL
  cutoff_method <- NULL
  cutoff_hdci_width <- NULL
  cutoff_dgp <- NULL
  cutoff_full <- NULL # full object (carries simulated $results for p-values)
  if (!is.null(cutoff)) {
    if (
      is.list(cutoff) &&
        !is.data.frame(cutoff) &&
        "item_cutoffs" %in% names(cutoff)
    ) {
      cutoff_full <- cutoff
      cutoff_n_iter <- cutoff$actual_iterations
      cutoff_req_iter <- cutoff$requested_iterations
      cutoff_method <- cutoff$cutoff_method
      cutoff_hdci_width <- cutoff$hdci_width
      cutoff_dgp <- cutoff$dgp
      cutoff <- cutoff$item_cutoffs
    }
    if (!is.data.frame(cutoff)) {
      stop(
        "`cutoff` must be NULL, the return value of RMitemRestscoreCutoff(), ",
        "or its $item_cutoffs data.frame.",
        call. = FALSE
      )
    }
    missing_cols <- setdiff(c("Item", "diff_low", "diff_high"), names(cutoff))
    if (length(missing_cols) > 0L) {
      stop(
        "`cutoff` data.frame is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    if (!missing(p_adj)) {
      warning(
        "`p_adj` is ignored when `cutoff` is supplied; the asymptotic ",
        "p-value it adjusts is replaced by the bootstrap null.",
        call. = FALSE
      )
    }
  }

  # --- Resolve p_value --------------------------------------------------------
  # NULL means "use the corrected p-value when the simulations are available",
  # as in RMitemInfit().
  have_sims <- !is.null(cutoff_full) && !is.null(cutoff_full$results)
  if (is.null(p_value)) {
    p_value <- have_sims
  }
  if (p_value) {
    if (!have_sims) {
      stop(
        "`p_value = TRUE` requires the full RMitemRestscoreCutoff() object ",
        "(it carries the simulated distributions in $results); a NULL cutoff ",
        "or the bare $item_cutoffs data.frame is not sufficient.",
        call. = FALSE
      )
    }
    if (!is.null(cutoff_n_iter) && cutoff_n_iter < 400L) {
      .notify_low_iterations(
        cutoff_n_iter,
        cutoff_req_iter,
        fn = "RMitemRestscoreCutoff()",
        id = "easyRasch2_low_iterations_restscore",
        correction = correction
      )
    }
    # The equal-tailed test's floor is 2/(B + 1), so an FDR correction needs
    # about twice the iterations it would under |t|.
    .warn_fdr_floor(
      cutoff_n_iter,
      nrow(cutoff),
      correction,
      alpha,
      unit = "items",
      fn = "RMitemRestscoreCutoff()",
      floor_mult = 2L
    )
  } else if (!is.null(cutoff)) {
    .notify_band_flagging(
      if (identical(cutoff_method, "quantile")) 0.95 else cutoff_hdci_width,
      nrow(cutoff),
      fn = "RMitemRestscoreCutoff()",
      id = "easyRasch2_band_flagging_restscore"
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

  data_mat <- as.matrix(data)
  n_items <- ncol(data)
  n_complete <- nrow(stats::na.omit(data))
  # iarm's item-restscore test refits on complete cases, so check them too.
  .check_complete_case_categories(data)
  if (!is.null(cutoff_full)) {
    .check_cutoff_sample(
      cutoff_full$sample_n,
      n_complete,
      "RMitemRestscoreCutoff()"
    )
  }

  # --- Fit Rasch model and compute item/person locations ----------------------
  # CML item parameters (psychotools; a dichotomous item is a 2-category PCM)
  # and WLE person locations, consistent with the rest of the package. The
  # item-restscore statistic from iarm is conditional and engine-invariant; the
  # relative-location reference shifts only slightly (WLE vs eRm MLE person
  # mean), since item and person locations move together with the scale.
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

  # --- Compute item-restscore statistics via iarm ----------------------------
  # Temporarily set rgl.useNULL to avoid rgl device issues during iarm fitting
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  i1 <- iarm::item_restscore(fit, p.adj = p_adj)
  i1 <- as.data.frame(i1)

  # i1[[1]] is the results matrix. iarm::item_restscore() appends an
  # adjusted-p column named "padj.<method>" (e.g. "padj.BH") only when
  # p.adj != "none"; with p.adj = "none" that column is absent and the
  # fixed position 5 is instead the significance-stars column ("sig"),
  # whose "***"/"." strings coerce to NA. Select the p-value column by
  # name: the adjusted column when present, otherwise the raw "pvalue".
  res_mat <- i1[[1]]
  cn <- colnames(res_mat)
  padj_idx <- grep("^padj", cn)
  p_col <- if (length(padj_idx) == 1L) padj_idx else match("pvalue", cn)
  observed <- as.numeric(res_mat[seq_len(n_items), 1L])
  expected <- as.numeric(res_mat[seq_len(n_items), 2L])
  p_adjusted <- as.numeric(res_mat[seq_len(n_items), p_col])

  # Flagged labels the misfit direction (only when adj. p < .05): observed
  # above expected = over-discrimination ("overfit", often local dependence);
  # below = under-discrimination ("underfit", often multidimensionality/noise);
  # "" otherwise. Note the value direction is opposite to infit (where a high
  # statistic is underfit).
  difference <- observed - expected
  flagged <- ifelse(
    !is.na(p_adjusted) & p_adjusted < 0.05 & difference > 0,
    "overfit",
    ifelse(
      !is.na(p_adjusted) & p_adjusted < 0.05 & difference < 0,
      "underfit",
      ""
    )
  )

  # --- Assemble result data.frame --------------------------------------------
  i2 <- data.frame(
    Item = names(data),
    Observed = observed,
    Expected = expected,
    Difference = difference,
    p_adjusted = p_adjusted,
    Flagged = flagged,
    Relative_location = as.numeric(relative_item_avg_locations),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # --- Apply cutoff if provided ----------------------------------------------
  # The asymptotic p-value is replaced, never shown alongside: the interval
  # and (when available) the bootstrap p-value take its place.
  if (!is.null(cutoff)) {
    # Recompute the observed statistic the way RMitemRestscoreCutoff() computes
    # its null, so both sides are unrounded (iarm returns them through
    # format(digits = 3)). Complete cases, refitted when needed, as iarm does.
    # .check_complete_case_categories() above has already stopped on complete
    # cases that reshape an item, so the refit below is aligned. The tryCatch
    # keeps the iarm values if .restscore_gamma() still refuses the data.
    cc <- stats::na.omit(data_mat)
    # check = FALSE: .check_complete_case_categories() has checked `cc`.
    fit_cc <- if (has_na) .cml_fit(cc, hessian = FALSE, check = FALSE) else fit
    rs <- tryCatch(
      .restscore_gamma(cc, lapply(psychotools::threshpar(fit_cc), cumsum)),
      error = function(e) NULL
    )
    if (!is.null(rs)) {
      i2$Observed <- rs[, "observed"]
      i2$Expected <- rs[, "expected"]
      i2$Difference <- i2$Observed - i2$Expected
    }

    data_items <- i2$Item
    if (!setequal(data_items, cutoff$Item)) {
      stop(
        "Item names in `cutoff` do not match item names in `data`.\n",
        "  data items  : ",
        paste(data_items, collapse = ", "),
        "\n",
        "  cutoff items: ",
        paste(cutoff$Item, collapse = ", "),
        call. = FALSE
      )
    }
    idx <- match(data_items, cutoff$Item)
    i2$Diff_low <- cutoff$diff_low[idx]
    i2$Diff_high <- cutoff$diff_high[idx]
    # Interval flagging, used when p_value = FALSE
    i2$Flagged <- ifelse(
      i2$Difference > i2$Diff_high,
      "overfit",
      ifelse(i2$Difference < i2$Diff_low, "underfit", "")
    )

    if (p_value) {
      sim_items <- unique(cutoff_full$results$Item)
      if (!setequal(data_items, sim_items)) {
        stop(
          "Item names in the cutoff simulations ($results) do not match ",
          "`data`.",
          call. = FALSE
        )
      }
      sim_mat <- tapply(
        cutoff_full$results$Difference,
        list(cutoff_full$results$iteration, cutoff_full$results$Item),
        function(x) x[1L]
      )
      obs_diff <- stats::setNames(i2$Difference, data_items)
      # Equal-tailed: the null of Difference is left-skewed (gamma is bounded
      # at 1), so a test on |t| would give the underfit tail more than
      # alpha/2. See dev/HANDOFF-restscore-equal-tails.md.
      pv <- .bootstrap_pvalues(
        obs_diff, sim_mat,
        correction = correction,
        tail = "split"
      )
      pidx <- match(data_items, pv$name)
      i2$p_restscore <- pv$p[pidx]
      i2$padj_restscore <- pv$padj[pidx]
      # Direction from the tail that rejected, not the sign of Difference: the
      # null mean of the difference is positive in small samples, so the sign
      # alone can disagree with the side of the null an item falls on.
      sig <- i2$padj_restscore < alpha
      i2$Flagged <- ifelse(
        is.na(sig) | !sig,
        "",
        ifelse(pv$direction[pidx] == "upper", "overfit", "underfit")
      )
      i2 <- i2[, c(
        "Item", "Observed", "Expected", "Difference", "Diff_low", "Diff_high",
        "p_restscore", "padj_restscore", "Flagged", "Relative_location"
      )]
    } else {
      i2 <- i2[, c(
        "Item", "Observed", "Expected", "Difference", "Diff_low", "Diff_high",
        "Flagged", "Relative_location"
      )]
    }
  }

  # --- Sort if requested -----------------------------------------------------
  if (!missing(sort) && identical(sort, "diff")) {
    # Sort by absolute magnitude so both over- and underfit items rise
    # to the top; the signed value is still what the user sees.
    i2 <- i2[order(abs(i2$Difference), decreasing = TRUE), ]
    rownames(i2) <- NULL
  }

  # --- Return ----------------------------------------------------------------
  if (output == "dataframe") {
    return(i2)
  }

  # Kable display rounding (the dataframe output above stays unrounded)
  i2 <- .round_display(i2, c(
    Observed = 2, Expected = 2, Difference = 3, p_adjusted = 3,
    Diff_low = 3, Diff_high = 3, p_restscore = 4, padj_restscore = 4,
    Relative_location = 2
  ))

  n_clause <- .n_caption(
    n_complete,
    n_total,
    if (has_na) "complete cases" else character()
  )
  direction_clause <- paste0(
    "overfit = observed above expected ",
    "(over-discrimination, often local dependence); underfit = below ",
    "(under-discrimination, often multidimensionality/noise)."
  )

  if (!is.null(cutoff)) {
    dgp_part <- if (is.null(cutoff_dgp)) "" else paste0(", ", cutoff_dgp, " DGP")
    if (p_value) {
      kbl_colnames <- c(
        "Item", "Observed", "Expected", "Difference", "Diff low",
        "Diff high", "p", "p (adj)", "Flagged", "Rel. location"
      )
      kbl_caption <- paste0(
        "Item-restscore associations. ",
        n_clause,
        ". Equal-tailed parametric bootstrap p-values for the difference from ",
        cutoff_n_iter,
        " iterations",
        dgp_part,
        ", multiplicity correction: ",
        .correction_label(correction),
        ". p-values cannot be smaller than 2/(",
        cutoff_n_iter,
        "+1) = ",
        .fmt_small(2 / (cutoff_n_iter + 1)),
        ".",
        if (2 / (cutoff_n_iter + 1) >= alpha) {
          paste0(" No item can be flagged at alpha = ", alpha, ".")
        },
        if (cutoff_n_iter < 400L) {
          paste0(
            " This is below the recommended 400, where decisions depend on",
            " the seed (Johansson, 2026)."
          )
        } else if (cutoff_n_iter < 1000L) {
          paste0(
            " Decisions are somewhat seed-dependent below 1000 iterations,",
            " so use 1000 to 2000 for a final analysis (Johansson, 2026)."
          )
        } else {
          ""
        },
        .attrition_clause(cutoff_n_iter, cutoff_req_iter),
        " Flagged (adj. p < ",
        alpha,
        "): ",
        direction_clause
      )
    } else {
      kbl_colnames <- c(
        "Item", "Observed", "Expected", "Difference", "Diff low",
        "Diff high", "Flagged", "Rel. location"
      )
      method_label <- .format_cutoff_method_label(
        cutoff_method,
        cutoff_hdci_width
      )
      kbl_caption <- paste0(
        "Item-restscore associations. ",
        n_clause,
        ". Interval for the difference ",
        if (is.null(cutoff_n_iter)) {
          "from simulation."
        } else {
          paste0(
            "from ",
            cutoff_n_iter,
            " simulation iterations",
            if (!is.null(method_label)) paste0(" (", method_label, ")") else "",
            "."
          )
        },
        " Flagged against the interval: ",
        direction_clause,
        .band_error_clause(
          if (identical(cutoff_method, "quantile")) 0.95 else cutoff_hdci_width,
          nrow(i2)
        ),
        .attrition_clause(cutoff_n_iter, cutoff_req_iter)
      )
    }
    return(knitr::kable(
      i2,
      format = "pipe",
      col.names = kbl_colnames,
      caption = kbl_caption
    ))
  }

  p_header <- if (identical(p_adj, "none")) {
    "p-value"
  } else {
    paste0("Adj. p-value (", p_adj, ")")
  }

  knitr::kable(
    i2,
    format = "pipe",
    col.names = c(
      "Item",
      "Observed",
      "Expected",
      "Difference",
      p_header,
      "Flagged",
      "Rel. location"
    ),
    caption = paste0(
      "Item-restscore associations. ",
      n_clause,
      ". Flagged (",
      if (identical(p_adj, "none")) "p" else "adj. p",
      " < .05): ",
      direction_clause,
      " The asymptotic p-values are miscalibrated under the Rasch model,",
      " flagging too many items as overfit in small samples and too few as",
      " underfit at any sample size. Use RMitemRestscoreCutoff() for",
      " simulation-based p-values."
    )
  )
}
