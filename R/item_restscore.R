#' Item Restscore Analysis
#'
#' Computes observed and model-expected item-restscore correlations using
#' `iarm::item_restscore()`, and enriches the output with the absolute
#' difference between observed and expected values, item average locations, and
#' item locations relative to the sample mean person location.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed,
#'   but at least one complete case (row with no `NA`) must be present.
#' @param cutoff Optional. Default `NULL`, which tests each item against the
#'   asymptotic normal reference distribution from `iarm::item_restscore()`,
#'   adjusted by `p_adj`. Can be:
#'   * The return value of \code{\link{RMitemRestscoreCutoff}}: the parametric
#'     bootstrap null distribution of the observed minus expected gamma.
#'     Adds the columns `Diff_low` and `Diff_high` (the descriptive interval),
#'     and, with `p_value` resolving to `TRUE`, `p_restscore` and
#'     `padj_restscore`.
#'   * The `$item_cutoffs` data.frame from \code{\link{RMitemRestscoreCutoff}}
#'     directly (columns `Item`, `diff_low`, `diff_high`). Items are then
#'     flagged against the interval.
#'
#'   With a `cutoff`, the asymptotic `p_adjusted` column is dropped and
#'   `p_adj` is ignored.
#' @param p_value Logical or `NULL`. Whether to compute bootstrap p-values
#'   from the simulated null distribution and flag on them.
#'
#'   `NULL` (the default) means **use them when they are available**: `TRUE`
#'   when `cutoff` is the full \code{\link{RMitemRestscoreCutoff}} object,
#'   which carries the per-item simulated values in its `$results` element,
#'   and `FALSE` otherwise. Explicit `TRUE` without the full object is an
#'   error rather than a silent downgrade. With `FALSE` and a `cutoff`, items
#'   are flagged against the interval instead, and a one-time message reports
#'   the family-wise error rate that implies.
#' @param correction Character. Multiple-comparison correction applied across
#'   items when `p_value = TRUE`: `"fwer"` (default) for the Westfall-Young
#'   studentised-max step-down (family-wise error rate), `"fdr_bh"` /
#'   `"fdr_by"` for Benjamini-Hochberg / Benjamini-Yekutieli false discovery
#'   rate control, or `"none"` for uncorrected per-item p-values.
#' @param alpha Numeric in (0, 1). Significance level for the `Flagged` column
#'   when `p_value = TRUE`. Default `0.05`.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying data.frame.
#' @param sort Optional character string. When `sort = "diff"`, rows are sorted
#'   by the absolute magnitude of `Difference` in descending order, so that
#'   both over- and underfitting items appear near the top.
#' @param p_adj Character string specifying the p-value adjustment method
#'   passed to `iarm::item_restscore()`. Default `"BH"` (Benjamini-Hochberg);
#'   use `"none"` for unadjusted p-values. Run `?stats::p.adjust` for the list
#'   of available methods. Only used when `cutoff = NULL`.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object (plain text table via
#'   `format = "pipe"`) with columns for item name, observed and expected
#'   restscore correlations, the signed difference (observed minus
#'   expected), adjusted p-value, the `Flagged` misfit label, and item
#'   location relative to the sample mean person location.
#' * If `output = "dataframe"`: a data.frame with columns `Item`, `Observed`,
#'   `Expected`, `Difference`, `p_adjusted`, `Flagged`, and
#'   `Relative_location`. `Flagged` is `"overfit"` (observed above expected,
#'   adj. p < .05), `"underfit"` (below, adj. p < .05), or `""` (not flagged).
#'
#' With a `cutoff`, `p_adjusted` is replaced by `Diff_low` and `Diff_high`,
#' followed by `p_restscore` (marginal two-sided bootstrap p-value) and
#' `padj_restscore` (corrected p-value) when `p_value` resolves to `TRUE`.
#' `Flagged` then reflects `padj_restscore < alpha`, or the interval when
#' `p_value = FALSE`.
#'
#' The `Difference` column is signed (observed minus expected):
#' *positive* values indicate that the item correlates more strongly with
#' the rest-score than the Rasch model predicts (over-discrimination /
#' *overfit*, often associated with local dependence), and *negative*
#' values indicate weaker-than-expected association (under-discrimination
#' / *underfit*, often associated with multidimensionality or noise).
#'
#' @details
#' Item-restscore correlations using Goodman-Kruskal's gamma (Kreiner, 2011) measure
#' the association between a person's score on a single item and their total
#' score on the remaining items (the "restscore"). Under a correctly fitting
#' Rasch model, observed and model-expected correlations should agree closely.
#'
#' Item parameters are estimated by conditional maximum likelihood via
#' `psychotools::pcmodel()` (a dichotomous item is a 2-category PCM); the
#' item-restscore statistic itself comes from `iarm::item_restscore()` and is
#' conditional on the total score, so it is invariant to the estimation engine.
#' Per-item average locations are the means of the CML thresholds, and the
#' person-location reference is the mean of the Warm WLE estimates.
#'
#' Relative item location is defined as the item's average location minus the
#' sample mean person location, providing a measure of item targeting.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' \strong{The asymptotic p-value is miscalibrated.} Without a `cutoff`, each
#' item is tested with `(observed - expected) / SE` against the standard
#' normal, where the SE is that of the observed gamma alone. The expected
#' gamma is estimated from the same data and correlates with the observed
#' one, so the SE is too large for the difference. The difference is also
#' biased upwards in small samples. Under a true Rasch model the test
#' therefore flags too many items as overfit in small samples and too few as
#' underfit at every sample size, and the adjustment chosen in `p_adj` does
#' not correct either. In simulation, 20 dichotomous items at n = 150 gave at
#' least one BH-flagged item in 13 percent of datasets (30 percent when 1.5
#' logits off target), almost all of them overfit, while 9 polytomous items
#' at n = 1000 gave about 1 percent. Pass the object from
#' [RMitemRestscoreCutoff()] as `cutoff` for p-values from a parametric
#' bootstrap null instead.
#'
#' \strong{Bootstrap p-values.} When `p_value = TRUE`, each item's observed
#' `Difference` is compared against its simulated null distribution (from
#' `cutoff$results`), studentised by the bootstrap mean and SD. Because the
#' bootstrap mean is subtracted, the small-sample upward bias of the
#' difference is removed, and because the bootstrap SD is used, the
#' correlation between observed and expected gamma is accounted for. The
#' marginal p-value is the two-sided Monte-Carlo p-value
#' `(1 + #{|t*| >= |t|}) / (B + 1)`, and `correction = "fwer"` applies the
#' Westfall-Young studentised-max step-down across items (Ferreira, 2024).
#'
#' The two directions do not get equal shares of alpha. Gamma is bounded at
#' 1, so the null of `Difference` is left-skewed and a two-sided test on
#' `|t|` rejects more often in the long (underfit) tail. In simulation under
#' a true Rasch model the total rate was nominal, but the per-item marginal
#' rate was about 3 percent for underfit and 2 percent for overfit with 20
#' dichotomous items at n = 150 and 1.5 logits off target, narrowing toward
#' equal shares with polytomous items and larger samples. With
#' `correction = "fwer"`, over four such conditions, at least one item was
#' flagged as underfit in 3.7 percent of datasets and as overfit in 1.6
#' percent, 5.0 percent in total. An equal-tailed test, with alpha/2 for each
#' direction, was evaluated and not adopted: it detected 2 to 3 percentage
#' points fewer underfitting items and no more overfitting ones.
#'
#' \strong{Reproducibility.} The bootstrap p-values depend on the simulated
#' null, so two analyses of the same data with different seeds can disagree
#' about items near the decision boundary. In simulation, in conditions
#' chosen to include such items, two seeds disagreed about the flag of at
#' least one item in about 10 percent of analyses at 400 iterations and about
#' 6 percent at 1000 with `correction = "fwer"`, and in 13 and 11 percent with
#' `correction = "fdr_bh"`. Set `seed` in [RMitemRestscoreCutoff()] for a
#' reproducible analysis, and use 1000 or more iterations for a final one.
#'
#' The direction in `Flagged` follows the studentised value, not the sign of
#' `Difference`. At small sample sizes the null mean of `Difference` is
#' slightly positive, so an item can be flagged `"underfit"` while its
#' `Difference` is still just above zero. This is rare.
#'
#' @inheritSection RMitemInfit Multiple comparisons
#' @inheritSection RMitemInfit Flags depend on the other items
#'
#' @references
#' Kreiner, S. (2011). A Note on Item–Restscore Association in Rasch Models.
#' *Applied Psychological Measurement, 35*(7), 557–561.
#' \doi{10.1177/0146621611410227}
#'
#' Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate' number
#' of hypotheses simultaneously. *Journal of Statistical Theory and Practice,
#' 19*(6). \doi{10.1007/s42519-024-00412-4}
#'
#' @seealso \code{\link{RMitemRestscoreCutoff}},
#'   \code{\link{RMitemRestscorePlot}}
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
        id = "easyRasch2_low_iterations_restscore"
      )
    }
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
    # With an unobserved middle category .restscore_gamma() refuses the data,
    # and the iarm values (expected NA for that item) are kept.
    cc <- stats::na.omit(data_mat)
    fit_cc <- if (has_na) psychotools::pcmodel(cc, hessian = FALSE) else fit
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
      pv <- .bootstrap_pvalues(obs_diff, sim_mat, correction = correction)
      pidx <- match(data_items, pv$name)
      i2$p_restscore <- pv$p[pidx]
      i2$padj_restscore <- pv$padj[pidx]
      # Direction from the studentised value: the null mean of the difference
      # is positive in small samples, so the sign of Difference alone can
      # disagree with the side of the null distribution an item falls on.
      null_mean <- colMeans(sim_mat, na.rm = TRUE)[data_items]
      sig <- i2$padj_restscore < alpha
      i2$Flagged <- ifelse(
        is.na(sig) | !sig,
        "",
        ifelse(i2$Difference > null_mean, "overfit", "underfit")
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
        ". Two-sided parametric bootstrap p-values for the difference from ",
        cutoff_n_iter,
        " iterations",
        dgp_part,
        ", multiplicity correction: ",
        .correction_label(correction),
        ". p-values cannot be smaller than 1/(",
        cutoff_n_iter,
        "+1) = ",
        round(1 / (cutoff_n_iter + 1), 4),
        ".",
        if (1 / (cutoff_n_iter + 1) >= alpha) {
          paste0(" No item can be flagged at alpha = ", alpha, ".")
        },
        if (cutoff_n_iter < 400L) {
          paste0(
            " This is below the calibrated floor of 400, where the",
            " correction is mildly liberal (Johansson, 2026)."
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
