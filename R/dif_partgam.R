#' Partial gamma DIF
#'
#' Computes the partial gamma coefficient of each item, the association
#' between the item and a grouping variable given the total score, to
#' assess differential item functioning (DIF). With the result of
#' [RMdifGammaCutoff()], it tests each item against simulations without
#' DIF, by default with a bootstrap p-value corrected across items.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Only complete cases (in `data` and `dif_var`) are used, so at least
#'   one is needed.
#' @param dif_var A vector (factor or character) of length `nrow(data)`
#'   with the groups, two or more.
#' @param cutoff Optional. Default `NULL`, the asymptotic test of
#'   `iarm::partgam_DIF()`. Either of:
#'   * The full result of [RMdifGammaCutoff()]. This enables the bootstrap
#'     p-values (see `p_value`), and the caption reports the simulation
#'     settings.
#'   * Its `$item_cutoffs` data.frame alone, with columns `Item`,
#'     `gamma_low` and `gamma_high`. Items are flagged against the
#'     interval.
#'
#'   Either adds the columns `Gamma_low`, `Gamma_high` and a logical
#'   `Flagged` (`TRUE` when the partial gamma is outside the simulated
#'   interval, or with p-values when the item is flagged).
#' @param p_value Logical or `NULL`. Whether to flag on two-sided bootstrap
#'   p-values (`p_gamma`, `padj_gamma`) from the simulations. The
#'   asymptotic p-value and star columns are then dropped, so that the
#'   table has one family of p-values, and the simulated interval is kept
#'   as a reference. `TRUE` needs the full [RMdifGammaCutoff()] result.
#'   `NULL` (default) uses the p-values whenever `cutoff` is that result,
#'   so a call without a cutoff keeps the asymptotic table. `FALSE` flags
#'   against the interval, a poor decision rule, since over k items a width
#'   `w` gives a family-wise error rate of `1 - w^k` (Johansson, 2026).
#' @param correction Character. Correction of the bootstrap p-values across
#'   items: `"fwer"` (default, Westfall-Young step-down), `"fdr_bh"`,
#'   `"fdr_by"` or `"none"`. Used only with p-values.
#' @param alpha Numeric in (0, 1). An item is flagged when its corrected
#'   p-value is below `alpha`. Default `0.05`. Used only with p-values.
#' @param output Character. `"kable"` (default) for a `knitr::kable()`
#'   table, or `"dataframe"`.
#'
#' @return
#' * `output = "kable"`: a `knitr_kable` table with the columns "Item",
#'   "Partial gamma", "SE", "Lower CI", "Upper CI", "Adj. p-value (BH)" and
#'   "p-value sign." (stars for the adjusted p-value, as iarm shows them),
#'   plus "Gamma low", "Gamma high" and "Flagged" with a `cutoff`. With
#'   p-values, the asymptotic columns are replaced by "p" and "p (adj)".
#' * `output = "dataframe"`: a data.frame with the columns `Item`, `gamma`,
#'   `se`, `lower`, `upper`, `padj_bh` and `Significance`, plus
#'   `gamma_low`, `gamma_high` and `flagged` with a `cutoff`. With p-values,
#'   `padj_bh` and `Significance` are replaced by `p_gamma` and
#'   `padj_gamma`.
#'
#' @details
#' \strong{Statistic.} The partial gamma (Bjorner et al., 1998) measures the
#' association between the item response and an external grouping variable
#' given the total score. Values near 0 indicate no DIF. Bjorner et al.
#' (1998) also proposed classes of DIF size, which [RMdifTree()] uses.
#' They are conventions, not values calibrated to the sample and items at
#' hand, and this function does not use them. The \pkg{iarm} package must
#' be installed (it is in Suggests).
#'
#' \strong{Asymptotic p-values.} The adjusted p-value is a
#' Benjamini-Hochberg adjustment, made by easyRasch2, of the p-values from
#' `iarm::partgam_DIF()` over the items. iarm's own adjusted column is a
#' Bonferroni correction whatever method is named, since iarm adjusts one
#' p-value at a time, so it is not used.
#'
#' \strong{Bootstrap p-values.} Each item's observed partial gamma is
#' compared with its simulated values (`cutoff$results`, simulated without
#' DIF). The statistic is the deviation from the mean of the observed and
#' simulated values together, divided by their SD. The uncorrected p is
#' two-sided, `(1 + #\{|t*| >= |t|\}) / (B + 1)` for B simulated datasets,
#' so it can be no smaller than `1 / (B + 1)`. It is two-sided since DIF
#' can favour either group, unlike [RMlocdepGamma()], which tests
#' one-sided for excess positive local dependence. `correction = "fwer"`
#' uses the Westfall-Young step-down, which uses the dependence between
#' the items (Westfall & Young, 1993, Ferreira, 2024). It controls the
#' family-wise error rate at any number of iterations in
#' [RMdifGammaCutoff()], but below 400 decisions depend on the seed and a
#' message says so. Unlike the asymptotic p-values, these are calibrated
#' against simulations from the Rasch model. They depend on the model and
#' the sample size, and are reported with the simulated interval, not
#' instead of it. The observed partial gamma they test is computed by the
#' same internal routine as the simulations, which reproduces
#' `iarm::partgam_DIF()` exactly.
#'
#' \strong{Real and artificial DIF.} Partial gamma conditions on the total
#' score, which includes every item. Real DIF in one item then induces
#' artificial DIF in the other items, in the opposite direction, since at a
#' given total score the item scores must add up to that score (Hagquist &
#' Andrich, 2017). An item can therefore be flagged only because another
#' item has DIF. The simulations in [RMdifGammaCutoff()] do not guard
#' against this, since they have no DIF in any item. Resolve DIF items one
#' at a time, starting with the largest, and repeat the analysis before
#' deciding on the next (Hagquist & Andrich, 2017). [RMdifAnchor()] avoids
#' the problem by linking the groups through anchor items.
#'
#' @inheritSection RMitemInfit Multiple comparisons
#'
#' @references
#' Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., &
#' Bech, P. (1998). Differential item functioning in the Danish translation
#' of the SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189--1202.
#' \doi{10.1016/S0895-4356(98)00111-5}
#'
#' Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate' number
#' of hypotheses simultaneously. *Journal of Statistical Theory and Practice,
#' 19*(1), 6. \doi{10.1007/s42519-024-00412-4}
#'
#' Hagquist, C., & Andrich, D. (2017). Recent advances in analysis of
#' differential item functioning in health research using the Rasch model.
#' *Health and Quality of Life Outcomes, 15*, 181.
#' \doi{10.1186/s12955-017-0755-0}
#'
#' Johansson, M. (2026). Simulation-based cutoffs for conditional item fit in
#' Rasch models: Iterations, multiplicity correction, and decision stability.
#' *PsyArXiv*. \doi{10.31234/osf.io/7pqz4_v2}
#'
#' Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple Testing*.
#' Wiley.
#'
#' @seealso [RMdifGammaCutoff()] for the simulations, [RMdifGammaPlot()]
#'   for the plot, [RMdifAnchor()] for anchored tests, [RMdifImpact()] for
#'   the size of the impact.
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("iarm", quietly = TRUE)) {
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:10)
#'   dif_group <- factor(sample(c("A", "B"), 200, replace = TRUE))
#'
#'   # Default kable output
#'   RMdifGamma(sim_data, dif_group)
#'
#'   # Return as data.frame
#'   RMdifGamma(sim_data, dif_group, output = "dataframe")
#'
#'   # Simulation-based cutoffs (100 iterations to keep the example quick)
#'   if (requireNamespace("ggdist", quietly = TRUE)) {
#'     cutoff_res <- RMdifGammaCutoff(sim_data, dif_var = dif_group,
#'                                    iterations = 100, parallel = FALSE,
#'                                    seed = 42)
#'     RMdifGamma(sim_data, dif_group, cutoff = cutoff_res)
#'
#'     # The full cutoff object switches on bootstrap p-values with the
#'     # family-wise (Westfall-Young) correction. Use 1000 or more
#'     # iterations in real analyses.
#'     RMdifGamma(sim_data, dif_group, cutoff = cutoff_res,
#'                output = "dataframe")
#'   }
#' }
#' }
RMdifGamma <- function(
  data,
  dif_var,
  cutoff = NULL,
  p_value = NULL,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable"
) {
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMdifGamma() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe"))
  correction <- match.arg(correction)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }

  validate_response_data(data)

  # --- Validate dif_var -------------------------------------------------------
  if (length(dif_var) != nrow(data)) {
    stop(
      "`dif_var` must have the same length as `nrow(data)` (",
      nrow(data),
      "), but has length ",
      length(dif_var),
      ".",
      call. = FALSE
    )
  }

  # Require at least 2 levels
  dif_levels <- length(unique(stats::na.omit(dif_var)))
  if (dif_levels < 2L) {
    stop(
      "`dif_var` must have at least 2 non-missing levels, but has ",
      dif_levels,
      ".",
      call. = FALSE
    )
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
        "`cutoff` must be NULL, the return value of RMdifGammaCutoff(), or its ",
        "$item_cutoffs data.frame.",
        call. = FALSE
      )
    }
    required_cols <- c("Item", "gamma_low", "gamma_high")
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
  # NULL means "use the corrected p-value when the simulations are available",
  # as in RMlocdepGamma() and RMitemInfit(). The bare $item_cutoffs, or no
  # cutoff at all, leaves nothing to compute a p-value from and resolves to
  # FALSE, which keeps the asymptotic table as it was.
  if (!is.null(p_value) && (!is.logical(p_value) || length(p_value) != 1L ||
    is.na(p_value))) {
    stop("`p_value` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }
  have_sims <- !is.null(cutoff_full) && !is.null(cutoff_full$results)
  if (is.null(p_value)) {
    p_value <- have_sims
  }
  n_items_total <- if (!is.null(cutoff)) nrow(cutoff) else NULL

  if (p_value) {
    if (!have_sims) {
      stop(
        "`p_value = TRUE` requires the full RMdifGammaCutoff() object (it ",
        "carries the simulated distributions in $results); a NULL cutoff or ",
        "the bare $item_cutoffs data.frame is not sufficient.",
        call. = FALSE
      )
    }
    # Below 400 decisions depend on the seed enough to say so on the
    # console. Between 400 and 1000 the table caption reports it instead.
    if (!is.null(cutoff_n_iter) && cutoff_n_iter < 400L) {
      .notify_low_iterations(
        cutoff_n_iter,
        cutoff_full$requested_iterations,
        fn = "RMdifGammaCutoff()",
        id = "easyRasch2_low_iterations_dif",
        correction = correction
      )
    }
    .warn_fdr_floor(
      cutoff_n_iter,
      n_items_total,
      correction,
      alpha,
      unit = "items",
      fn = "RMdifGammaCutoff()"
    )
  } else if (!is.null(cutoff_full)) {
    .notify_band_flagging(
      if (identical(cutoff_method, "quantile")) 0.95 else cutoff_hdci_width,
      n_items_total,
      unit = "items",
      fn = "RMdifGammaCutoff()",
      id = "easyRasch2_band_flagging_dif"
    )
  }

  # --- rgl workaround ---------------------------------------------------------
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # --- Same sample as the cutoff? ---------------------------------------------
  # The DIF null keeps each respondent's observed group, so the group sizes
  # have to match as well as the total.
  if (!is.null(cutoff_full)) {
    cc_dif <- stats::complete.cases(as.data.frame(data)) & !is.na(dif_var)
    .check_cutoff_sample(
      cutoff_full$sample_n,
      sum(cc_dif),
      "RMdifGammaCutoff()",
      policy = "complete cases in data and dif_var",
      groups_cutoff = cutoff_full$dif_group_sizes,
      groups_used = as.integer(table(droplevels(factor(dif_var[cc_dif]))))
    )
  }

  # --- Compute partial gamma DIF via iarm -------------------------------------
  sink(nullfile())
  #on.exit(sink(), add = TRUE)
  pgam_raw <- iarm::partgam_DIF(as.data.frame(data), dif_var)
  sink()

  # pgam_raw is a data.frame with character columns that need numeric conversion
  pgam_df <- data.frame(
    Item = as.character(pgam_raw$Item),
    gamma = as.numeric(pgam_raw$gamma),
    se = as.numeric(pgam_raw$se),
    pvalue = as.numeric(pgam_raw$pvalue),
    # iarm adjusts one p-value at a time (p.adjust() with n = k on a single
    # value), which makes its "BH" column a Bonferroni correction whatever
    # method is named. The Benjamini-Hochberg adjustment is applied here
    # across all items instead, and the stars are recomputed from it.
    padj_bh = stats::p.adjust(as.numeric(pgam_raw$pvalue), method = "BH"),
    Significance = "",
    lower = as.numeric(pgam_raw$lower),
    upper = as.numeric(pgam_raw$upper),
    stringsAsFactors = FALSE
  )
  pgam_df$Significance <- .p_stars(pgam_df$padj_bh)

  # Keep output columns (unrounded; the kable path rounds a display copy)
  result_df <- pgam_df[, c(
    "Item",
    "gamma",
    "se",
    "lower",
    "upper",
    "padj_bh",
    "Significance"
  )]

  # --- Apply cutoff if provided -----------------------------------------------
  if (!is.null(cutoff)) {
    data_items <- result_df$Item
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
    cutoff_sub <- cutoff[, c("Item", "gamma_low", "gamma_high")]
    result_df <- merge(result_df, cutoff_sub, by = "Item", sort = FALSE)
    # Restore original row order
    result_df <- result_df[match(data_items, result_df$Item), ]
    rownames(result_df) <- NULL
    result_df$flagged <- result_df$gamma < result_df$gamma_low |
      result_df$gamma > result_df$gamma_high

    if (p_value) {
      # Compare the observed gamma to its simulated null (cutoff_full$results).
      # Two-sided: DIF in either direction matters. The observed value comes
      # from the same vectorised routine as the simulated ones, so the two
      # cannot diverge. It equals iarm's gamma exactly
      # (dev/dif_partgam_fast_validation.qmd).
      sim_items <- unique(cutoff_full$results$Item)
      if (!setequal(data_items, sim_items)) {
        stop(
          "Item names in the cutoff simulations ($results) do not match ",
          "`data`.",
          call. = FALSE
        )
      }
      sim_mat <- tapply(
        cutoff_full$results$gamma,
        list(cutoff_full$results$iteration, cutoff_full$results$Item),
        function(x) x[1L]
      )
      obs_fast <- .partgam_dif_gamma(as.data.frame(data), dif_var)
      observed <- stats::setNames(obs_fast$gamma, obs_fast$Item)
      pv <- .bootstrap_pvalues(
        observed,
        sim_mat,
        correction = correction,
        tail = "two.sided"
      )
      idx <- match(result_df$Item, pv$name)
      result_df$p_gamma <- pv$p[idx]
      result_df$padj_gamma <- pv$padj[idx]
      result_df$flagged <- !is.na(result_df$padj_gamma) &
        result_df$padj_gamma < alpha
      # Drop the asymptotic p-value pair; the bootstrap p-values replace it.
      result_df <- result_df[, c(
        "Item",
        "gamma",
        "se",
        "lower",
        "upper",
        "gamma_low",
        "gamma_high",
        "p_gamma",
        "padj_gamma",
        "flagged"
      )]
    } else {
      # Reorder columns
      result_df <- result_df[, c(
        "Item",
        "gamma",
        "se",
        "lower",
        "upper",
        "padj_bh",
        "Significance",
        "gamma_low",
        "gamma_high",
        "flagged"
      )]
    }
  }

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(result_df)
  }

  # Kable display rounding (the dataframe output above stays unrounded)
  result_df <- .round_display(result_df, c(
    gamma = 3, se = 3, lower = 3, upper = 3, padj_bh = 3,
    gamma_low = 3, gamma_high = 3, p_gamma = 4, padj_gamma = 4
  ))

  # Build caption
  n_complete <- sum(stats::complete.cases(cbind(
    as.data.frame(data),
    dif_var = dif_var
  )))
  has_na <- anyNA(as.data.frame(data)) || anyNA(dif_var)
  n_clause <- .n_caption(
    n_complete,
    nrow(as.data.frame(data)),
    if (has_na) "complete cases" else character()
  )
  if (p_value) {
    caption_text <- paste0(
      "Partial gamma DIF analysis. ",
      n_clause,
      ". Two-sided bootstrap p-values from ",
      cutoff_n_iter,
      " iterations (replacing the asymptotic BH p-values); multiplicity ",
      "correction: ",
      .correction_label(correction),
      "; flagged at padj < ",
      alpha,
      ". p-values cannot be smaller than 1/(",
      cutoff_n_iter,
      "+1) = ",
      .fmt_small(1 / (cutoff_n_iter + 1)),
      ". The interval is shown as description and is not the decision rule.",
      " Positive gamma indicates higher scores in higher DIF group levels.",
      .iteration_note(cutoff_n_iter, correction),
      .attrition_clause(cutoff_n_iter, cutoff_full$requested_iterations)
    )
  } else if (is.null(cutoff)) {
    caption_text <- paste0(
      "Partial gamma DIF analysis. ",
      n_clause,
      ". Positive gamma indicates higher scores in higher DIF group levels."
    )
  } else if (!is.null(cutoff_n_iter)) {
    method_label <- .format_gamma_cutoff_method_label(
      cutoff_method,
      cutoff_hdci_width
    )
    iter_part <- paste0(cutoff_n_iter, " simulation iterations")
    caption_text <- paste0(
      "Partial gamma DIF analysis. ",
      n_clause,
      ". Cutoff values based on ",
      if (!is.null(method_label)) {
        paste0(iter_part, " (", method_label, ").")
      } else {
        paste0(iter_part, ".")
      },
      .band_error_clause(
        if (identical(cutoff_method, "quantile")) 0.95 else cutoff_hdci_width,
        n_items_total,
        unit = "items"
      ),
      .attrition_clause(cutoff_n_iter, cutoff_full$requested_iterations)
    )
  } else {
    caption_text <- paste0(
      "Partial gamma DIF analysis. ",
      n_clause,
      ". Simulation-based cutoff values applied."
    )
  }

  knitr::kable(
    result_df,
    format = "pipe",
    col.names = if (p_value) {
      c(
        "Item",
        "Partial gamma",
        "SE",
        "Lower CI",
        "Upper CI",
        "Gamma low",
        "Gamma high",
        "p",
        "p (adj)",
        "Flagged"
      )
    } else if (is.null(cutoff)) {
      c(
        "Item",
        "Partial gamma",
        "SE",
        "Lower CI",
        "Upper CI",
        "Adj. p-value (BH)",
        "p-value sign."
      )
    } else {
      c(
        "Item",
        "Partial gamma",
        "SE",
        "Lower CI",
        "Upper CI",
        "Adj. p-value (BH)",
        "p-value sign.",
        "Gamma low",
        "Gamma high",
        "Flagged"
      )
    },
    caption = caption_text
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
#' @noRd
.format_gamma_cutoff_method_label <- function(cutoff_method, hdci_width) {
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


#' Simulation-based partial gamma DIF cutoffs
#'
#' Simulates the partial gamma DIF coefficient (Bjorner et al., 1998,
#' Henninger et al., 2025) of each item when there is no DIF. Passed to [RMdifGamma()], the result
#' gives bootstrap p-values corrected across items, and an interval per
#' item.
#' Every simulated dataset keeps each respondent's observed group membership
#' and total score, so the null reflects the observed group sizes and any
#' difference between the groups' latent distributions.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Only complete cases (in `data` and `dif_var`) are used.
#' @param dif_var A vector (factor, character or integer) of length
#'   `nrow(data)` with the groups.
#'   Group membership is held fixed in every simulated dataset, so the group
#'   sizes and each group's distribution of total scores are those observed.
#' @param iterations Integer. Number of simulated datasets. Default `400`.
#'   Use 1000 to 2000 for a final analysis.
#' @param parallel Logical. Run the simulations in parallel via \pkg{mirai}
#'   when it is installed. Default `TRUE`.
#' @param n_cores Integer or `NULL`. Number of parallel workers. When
#'   `NULL`, `getOption("mc.cores")` is used. If that is unset too, a
#'   warning is given and the simulations run sequentially.
#' @param verbose Logical. Show a progress bar. Default `FALSE`.
#' @param seed Integer or `NULL`. Random seed for reproducibility. See
#'   [easyRasch2-reproducibility] for what this guarantees and how it
#'   interacts with `parallel`.
#' @param cutoff_method Character. How the interval per item is computed:
#'   `"hdci"` (default), the highest-density continuous interval from
#'   `ggdist::hdci()` (needs \pkg{ggdist}), or `"quantile"`, the 2.5th to
#'   97.5th percentiles.
#' @param hdci_width Numeric. Width of the interval with
#'   `cutoff_method = "hdci"`. Default `0.95`, as in [RMitemInfitCutoff()].
#'   Ignored with `"quantile"`.
#' @param dgp Character. How the datasets without DIF are simulated.
#'   `"conditional"` (default) draws each respondent's response pattern from
#'   the model given their observed total score, with the CML item thresholds
#'   fixed, as the conditional option of [RMitemRestscoreCutoff()] does. `"permutation"` keeps the
#'   observed responses and permutes group membership among respondents with
#'   the same total score, the Monte Carlo form of an exact conditional test of
#'   item-group independence given the score (Kreiner, 1987). It needs no item
#'   parameters, keeps any misfit elsewhere in the observed responses, and is
#'   25 to 60 times faster per iteration. Both held the nominal error rate in
#'   simulation, including when the groups differed by 1 SD.
#'
#' @return A list with:
#' \describe{
#'   \item{`results`}{data.frame with `iteration`, `Item` and `gamma`, one
#'     row per item and simulated dataset.}
#'   \item{`item_cutoffs`}{data.frame with the interval per item: `Item`,
#'     `gamma_low`, `gamma_high`, computed by `cutoff_method`.}
#'   \item{`actual_iterations`}{Number of simulated datasets that could be
#'     used.}
#'   \item{`requested_iterations`}{Number of iterations asked for.}
#'   \item{`sample_n`}{Number of complete cases used.}
#'   \item{`sample_n_total`}{Number of respondents in the raw input data,
#'     before removing rows with `NA` in `data` or `dif_var`.}
#'   \item{`sample_has_na`}{Logical. Whether `data` or `dif_var` contained
#'     any missing values.}
#'   \item{`sample_summary`}{Summary statistics of the WLE person locations,
#'     or `NULL` when the model could not be fitted with
#'     `dgp = "permutation"`.}
#'   \item{`item_names`}{Character vector of item names from data.}
#'   \item{`dif_group_sizes`}{Integer vector of group sizes, held fixed in
#'     every simulated dataset.}
#'   \item{`cutoff_method`}{The method used to compute cutoffs (`"hdci"` or
#'     `"quantile"`).}
#'   \item{`hdci_width`}{The HDCI width used (only meaningful when
#'     `cutoff_method = "hdci"`).}
#'   \item{`dgp`}{The data-generating process used.}
#' }
#'
#' @details
#' Partial gamma conditions on the total score, and under the Rasch model the
#' response to an item is independent of group membership given the total
#' score. Both ways of simulating keep each respondent's group and total
#' score, so the simulations have the observed distribution of groups over
#' the total scores. Assigning simulated respondents to groups at random, as
#' an earlier version did, gives both groups the same latent distribution.
#' When a minority group differed from the majority by 1 SD, that was too
#' narrow and flagged at least one item in about 14% of datasets without
#' DIF, against a nominal 5%.
#'
#' In each iteration one dataset is simulated (see `dgp`) and the partial
#' gamma of every item computed. The computation reproduces
#' `iarm::partgam_DIF()` exactly but is vectorised, so \pkg{iarm} is not
#' needed here. Datasets that cannot be used are left out, and
#' `actual_iterations` gives the number used.
#'
#' The conditional data-generating process uses CML item thresholds from
#' `psychotools::pcmodel()` (a dichotomous item is a 2-category PCM). The
#' group order follows the levels of `dif_var` (alphabetical for character
#' vectors), which sets the sign of gamma as in `iarm::partgam_DIF()`.
#'
#' The simulations have no DIF in any item, so the cutoffs do not guard
#' against artificial DIF, which real DIF in one item induces in the others
#' (Hagquist & Andrich, 2017). See "Real and artificial DIF" in
#' [RMdifGamma()].
#'
#' Parallel processing needs the optional \pkg{mirai} package.
#'
#' @references
#' Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., &
#' Bech, P. (1998). Differential item functioning in the Danish translation
#' of the SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189--1202.
#' \doi{10.1016/S0895-4356(98)00111-5}
#'
#' Hagquist, C., & Andrich, D. (2017). Recent advances in analysis of
#' differential item functioning in health research using the Rasch model.
#' *Health and Quality of Life Outcomes, 15*, 181.
#' \doi{10.1186/s12955-017-0755-0}
#'
#' Henninger, M., Radek, J., Sengewald, M.-A., & Strobl, C. (2025).
#' Partial credit trees meet the partial gamma coefficient for quantifying
#' DIF and DSF in polytomous items. *Behaviormetrika, 52*, 221--257.
#' \doi{10.1007/s41237-024-00252-3}
#'
#' Kreiner, S. (1987). Analysis of multidimensional contingency tables by
#' exact conditional tests: Techniques and strategies. *Scandinavian Journal
#' of Statistics, 14*(2), 97--112.
#'
#' @seealso [RMdifGamma()] for the test, [RMdifGammaPlot()] for the plot,
#'   [iarm::partgam_DIF()].
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
#'   dif_sex <- sample(c("male", "female"), 200, replace = TRUE)
#'
#'   # 100 iterations, sequentially, for a quick example
#'   cutoff_res <- RMdifGammaCutoff(sim_data, dif_var = dif_sex,
#'                                  iterations = 100, parallel = FALSE,
#'                                  seed = 42)
#'   cutoff_res$item_cutoffs
#'
#'   # The permutation null needs no item parameters and is much faster
#'   perm_res <- RMdifGammaCutoff(sim_data, dif_var = dif_sex,
#'                                iterations = 100, parallel = FALSE,
#'                                seed = 42, dgp = "permutation")
#'   perm_res$item_cutoffs
#' }
#' }
RMdifGammaCutoff <- function(
  data,
  dif_var,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95,
  dgp = c("conditional", "permutation")
) {
  .check_iterations(iterations)
  cutoff_method <- match.arg(cutoff_method, c("hdci", "quantile"))
  dgp <- match.arg(dgp)

  if (cutoff_method == "hdci" && !requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required when cutoff_method = \"hdci\" but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")\n",
      "Alternatively, use cutoff_method = \"quantile\" to avoid this dependency.",
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
      nrow(data),
      "), but has length ",
      length(dif_var),
      ".",
      call. = FALSE
    )
  }

  # Only complete cases (both data and dif_var). Record the raw total and
  # whether anything was dropped so callers (e.g. RMdifGammaPlot) can report
  # the sample in the standard `n = X of Y respondents` form.
  n_total <- nrow(data)
  has_na <- anyNA(data) || anyNA(dif_var)
  complete_idx <- stats::complete.cases(data) & !is.na(dif_var)
  data <- data[complete_idx, , drop = FALSE]
  dif_var <- dif_var[complete_idx]

  if (nrow(data) == 0L) {
    stop(
      "No complete cases in data after removing rows with NA in data or dif_var.",
      call. = FALSE
    )
  }

  # Group structure. factor() keeps the level order of a factor and sorts a
  # character or numeric vector, which is the coding iarm::partgam_DIF() uses
  # and so sets the sign of gamma.
  dif_factor <- factor(dif_var)
  dif_factor <- droplevels(dif_factor)
  dif_table <- table(dif_factor)
  n_dif_groups <- length(dif_table)

  if (n_dif_groups < 2L) {
    stop(
      "`dif_var` must have at least 2 distinct groups, but has ",
      n_dif_groups,
      ".",
      call. = FALSE
    )
  }

  # --- Parallel setup ---------------------------------------------------------
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
  storage.mode(data_mat) <- "integer"
  sample_n <- nrow(data_mat)
  item_names_vec <- colnames(data_mat)
  if (is.null(item_names_vec)) {
    item_names_vec <- paste0("V", seq_len(ncol(data_mat)))
    colnames(data_mat) <- item_names_vec
  }

  # The conditional DGP needs the CML thresholds. The permutation DGP needs no
  # model, so a failed fit there only loses the person-location summary.
  pool <- if (dgp == "conditional") {
    .wle_theta_pool(data_mat)
  } else {
    tryCatch(.wle_theta_pool(data_mat), error = function(e) NULL)
  }

  sim_data_list <- list(
    dgp = dgp,
    sample_n = sample_n,
    n_items = ncol(data_mat),
    item_names = item_names_vec,
    dif_factor = dif_factor
  )
  if (dgp == "conditional") {
    # Each respondent keeps their group and total score; their response
    # pattern is redrawn from the Rasch conditional distribution given that
    # score. Respondents are grouped by score so each group is drawn in one
    # batch, and the rows (and so the group labels) stay where they are.
    sim_data_list$thr_list <- pool$thr_list
    sim_data_list$cond_groups <- .cond_groups(data_mat, pool$thr_list)
  } else {
    # Responses stay as observed; group labels are permuted within strata of
    # the total score. Strata with a single respondent cannot change.
    sim_data_list$data_mat <- data_mat
    strata <- split(seq_len(sample_n), rowSums(data_mat))
    sim_data_list$strata <- strata[lengths(strata) > 1L]
  }

  if (use_parallel) {
    results_raw <- run_partgam_sim_parallel(
      iterations,
      sim_seeds,
      sim_data_list,
      n_cores,
      verbose
    )
  } else {
    results_raw <- run_partgam_sim_sequential(
      iterations,
      sim_seeds,
      sim_data_list,
      verbose
    )
  }

  # Filter out failures (character strings indicate errors)
  ok <- vapply(results_raw, is.data.frame, logical(1L))
  successful <- results_raw[ok]

  if (length(successful) == 0L) {
    stop(.all_sims_failed_message(data_mat), call. = FALSE)
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
  item_cutoffs <- do.call(
    rbind,
    lapply(item_names, function(item) {
      sub <- results_df[results_df$Item == item, ]
      g <- sub$gamma[is.finite(sub$gamma)]
      if (length(g) == 0L) {
        return(data.frame(
          Item = item,
          gamma_low = NA_real_,
          gamma_high = NA_real_,
          stringsAsFactors = FALSE,
          row.names = NULL
        ))
      }
      if (cutoff_method == "hdci") {
        gamma_interval <- ggdist::hdci(g, .width = hdci_width)
        data.frame(
          Item = item,
          gamma_low = gamma_interval[1L, 1L],
          gamma_high = gamma_interval[1L, 2L],
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      } else {
        data.frame(
          Item = item,
          gamma_low = stats::quantile(g, 0.025),
          gamma_high = stats::quantile(g, 0.975),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
    })
  )
  rownames(item_cutoffs) <- NULL

  list(
    results = results_df,
    item_cutoffs = item_cutoffs,
    actual_iterations = actual_iterations,
    requested_iterations = iterations,
    sample_n = sample_n,
    sample_n_total = n_total,
    sample_has_na = has_na,
    sample_summary = if (is.null(pool)) NULL else summary(pool$thetas),
    item_names = item_names_vec,
    dif_group_sizes = as.integer(dif_table),
    cutoff_method = cutoff_method,
    hdci_width = hdci_width,
    dgp = dgp
  )
}

# ---------------------------------------------------------------------------
# Internal: vectorised partial gamma for DIF
# ---------------------------------------------------------------------------

#' Partial gamma DIF coefficient for every item, coefficient only
#'
#' Vectorised implementation of Davis's (1967) partial gamma between each item
#' and the grouping variable, conditioned on the total score. It returns the
#' coefficient and nothing else, which is all the bootstrap needs, and is used
#' for both the simulated null in [RMdifGammaCutoff()] and the observed
#' statistic that [RMdifGamma()] tests against it, so the two cannot diverge.
#' The `gamma`, `se`, `lower` and `upper` columns shown to users still come
#' from `iarm::partgam_DIF()`.
#'
#' The coding follows `iarm::partgam_DIF()`, and agreement is exact
#' (`dev/dif_partgam_fast_validation.qmd`, asserted in
#' `tests/testthat/test-partgam_dif.R`): complete cases on the items and the
#' grouping variable, the total score over all items (the item included), and
#' groups ordered by the levels of `factor(dif_var)`, so three or more groups
#' are treated as ordinal. Positive gamma means that respondents in
#' higher-coded groups score higher on the item at the same total score.
#'
#' @details
#' For one stratum with count matrix \eqn{N} (rows: item categories, columns:
#' groups), \eqn{G_x} and \eqn{G_y} the strictly-upper indicator matrices of
#' matching size and \eqn{A = G_x N}, the concordant count is
#' \eqn{\sum N \circ (A G_y^T)} and the discordant count
#' \eqn{\sum N \circ (A G_y)}. Each unordered pair of respondents is counted
#' once, and the counts are pooled over strata.
#'
#' @param data data.frame or matrix of item responses scored from 0.
#' @param dif_var Grouping variable, same length as `nrow(data)`.
#' @return data.frame with `Item` and `gamma`, one row per item in column
#'   order. `gamma` is `NA_real_` for an item with no concordant and no
#'   discordant pair (a constant item), where `iarm::partgam_DIF()` errors for
#'   the whole dataset.
#' @keywords internal
#' @noRd
.partgam_dif_gamma <- function(data, dif_var) {
  X <- as.matrix(data)
  storage.mode(X) <- "integer"
  items <- colnames(X)
  if (is.null(items)) {
    items <- paste0("I", seq_len(ncol(X)))
  }
  f <- if (is.factor(dif_var)) dif_var else factor(dif_var)
  ok <- stats::complete.cases(X) & !is.na(f)
  X <- X[ok, , drop = FALSE]
  y <- as.integer(f[ok]) - 1L
  my <- nlevels(f)
  mx <- max(X) + 1L
  score <- rowSums(X)
  upper <- function(m) {
    outer(seq_len(m), seq_len(m), function(a, b) as.numeric(b > a))
  }
  Gx <- upper(mx)
  Gy <- upper(my)
  tGy <- t(Gy)
  data.frame(
    Item = items,
    gamma = vapply(
      seq_len(ncol(X)),
      function(i) .partgam_dif_one(X[, i], y, score, mx, my, Gx, Gy, tGy),
      numeric(1L)
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Partial gamma between one item and the grouping variable
#'
#' @param x Integer item responses scored from 0.
#' @param y Integer group codes from 0.
#' @param z Integer conditioning variable (the total score).
#' @param mx,my Number of item categories and of groups.
#' @param Gx,Gy,tGy Strictly-upper indicator matrices and the transpose of
#'   `Gy`.
#' @return The partial gamma, or `NA_real_` when no pair is concordant or
#'   discordant.
#' @keywords internal
#' @noRd
.partgam_dif_one <- function(x, y, z, mx, my, Gx, Gy, tGy) {
  zc <- z - min(z) + 1L
  nz <- max(zc)
  counts <- tabulate(
    (zc - 1L) * mx * my + y * mx + x + 1L,
    nbins = mx * my * nz
  )
  dim(counts) <- c(mx, my, nz)
  conc <- 0
  disc <- 0
  for (s in seq_len(nz)) {
    N <- counts[, , s]
    # A stratum holding fewer than two observations contributes no pairs.
    if (sum(N) < 2L) next
    A <- Gx %*% N
    conc <- conc + sum(N * (A %*% tGy))
    disc <- disc + sum(N * (A %*% Gy))
  }
  total <- conc + disc
  if (total == 0) NA_real_ else (conc - disc) / total
}

# ---------------------------------------------------------------------------
# Internal: single simulation iteration
# ---------------------------------------------------------------------------

#' Run a single partial gamma DIF simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMdifGammaCutoff()].
#' @return A data.frame with columns `Item` and `gamma`, or a character string
#'   on failure.
#' @noRd
run_single_partgam_sim <- function(seed, data_list) {
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
      if (identical(data_list$dgp, "conditional")) {
        sim_mat <- as.matrix(.sim_cond_dataset(data_list))
        groups <- data_list$dif_factor
      } else {
        sim_mat <- data_list$data_mat
        groups <- data_list$dif_factor
        for (s in data_list$strata) {
          groups[s] <- groups[s[sample.int(length(s))]]
        }
      }
      colnames(sim_mat) <- data_list$item_names
      .partgam_dif_gamma(sim_mat, groups)
    },
    error = function(e) {
      as.character(conditionMessage(e))
    }
  )
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
#' @noRd
run_partgam_sim_parallel <- function(
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
        run_single_partgam_sim(seed, data_list)
      },
      seed = sim_seeds[sim],
      data_list = sim_data_list,
      run_single_partgam_sim = run_single_partgam_sim,
      .partgam_dif_gamma = .partgam_dif_gamma,
      .partgam_dif_one = .partgam_dif_one,
      # Conditional-DGP generators.
      .sim_cond_dataset = .sim_cond_dataset,
      .sim_conditional = .sim_conditional,
      .esf_convolve = .esf_convolve
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
#' @noRd
run_partgam_sim_sequential <- function(
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


#' Plot the simulated partial gamma DIF values
#'
#' Plots the simulated partial gamma of each item from [RMdifGammaCutoff()],
#' and optionally the observed values in your data, computed with the same
#' routine as the simulations, which reproduces `iarm::partgam_DIF()`.
#'
#' @param simfit The result of [RMdifGammaCutoff()].
#' @param data Optional. The item responses the simulations were made from,
#'   items scored from 0. The observed partial gamma is then drawn as orange
#'   diamonds.
#' @param dif_var Needed with `data`. The vector of groups used in
#'   [RMdifGammaCutoff()], of length `nrow(data)`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' Each item's simulated values are drawn as dots, with the median, an
#' inner 66% interval and an outer interval. The outer interval is the one
#' in `simfit$item_cutoffs`, which [RMdifGamma()] tabulates: the HDCI of
#' width `hdci_width`, or the 2.5th to 97.5th percentiles with
#' `cutoff_method = "quantile"`. The inner interval uses the same method.
#' Without `data` the intervals are drawn with
#' `ggdist::stat_dotsinterval()`. With `data` they are black segments with
#' a black dot for the median, and the observed value is an orange diamond.
#'
#' The intervals are per item and not corrected for testing several items.
#' With the full cutoff object, [RMdifGamma()] flags items by adjusted
#' p-values instead, so an item just outside its interval need not be
#' flagged there. The plot shows no p-values.
#'
#' The \pkg{ggplot2} and \pkg{ggdist} packages must be installed (they are
#' in Suggests).
#'
#' @seealso [RMdifGammaCutoff()], [RMdifGamma()]
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggdist", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:10)
#'   dif_group <- factor(sample(c("A", "B"), 200, replace = TRUE))
#'
#'   # Run the simulation
#'   cutoff_res <- RMdifGammaCutoff(sim_data, dif_var = dif_group,
#'                                  iterations = 100, parallel = FALSE,
#'                                  seed = 42)
#'
#'   # Simulated distribution only
#'   RMdifGammaPlot(cutoff_res)
#'
#'   # With observed partial gamma overlaid
#'   RMdifGammaPlot(cutoff_res, data = sim_data, dif_var = dif_group)
#' }
#' }
RMdifGammaPlot <- function(simfit, data, dif_var) {
  # --- Check required packages ------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for RMdifGammaPlot() but is not installed.\n",
      "Install it with: install.packages(\"ggplot2\")",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required for RMdifGammaPlot() but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")",
      call. = FALSE
    )
  }

  # --- Validate simfit --------------------------------------------------------
  required_names <- c(
    "results",
    "item_cutoffs",
    "actual_iterations",
    "sample_n",
    "item_names"
  )
  missing_names <- setdiff(required_names, names(simfit))
  if (length(missing_names) > 0L) {
    stop(
      "`simfit` is missing required components: ",
      paste(missing_names, collapse = ", "),
      ".\nExpected the return value of RMdifGammaCutoff().",
      call. = FALSE
    )
  }

  results_df <- simfit$results
  item_cutoffs <- simfit$item_cutoffs
  actual_iterations <- simfit$actual_iterations
  sample_n <- simfit$sample_n
  item_names <- simfit$item_names

  # Standard sample-size clause for the (complete-case) simulation sample.
  # `sample_n_total` / `sample_has_na` are absent in cutoff objects made by
  # older versions, so fall back to the plain count.
  sample_clause <- .n_caption(
    sample_n,
    if (is.null(simfit$sample_n_total)) sample_n else simfit$sample_n_total,
    if (isTRUE(simfit$sample_has_na)) "complete cases" else character()
  )

  # Item factor levels (reversed for plotting top-to-bottom)
  item_levels <- rev(item_names)

  # --- Compute per-item summary intervals for segment overlays ----------------
  # The outer interval is the one stored in `simfit$item_cutoffs` (the one
  # RMdifGamma() tabulates) and the inner 66% interval uses the same method.
  spec <- .sim_interval_spec(simfit)
  interval_note <- .sim_interval_caption(spec, "item", "RMdifGamma()")
  key_names <- item_names
  cut_idx <- match(key_names, item_cutoffs$Item)
  iv <- .sim_plot_intervals(
    results_df$gamma, results_df$Item, key_names,
    item_cutoffs$gamma_low[cut_idx], item_cutoffs$gamma_high[cut_idx],
    spec$method
  )
  lo_hi <- data.frame(
    Item = key_names,
    min_gamma = iv$lo,
    max_gamma = iv$hi,
    p66lo_gamma = iv$lo66,
    p66hi_gamma = iv$hi66,
    median_gamma = iv$median,
    stringsAsFactors = FALSE
  )

  # --- Case 1: no observed data, show simulation distribution only ------------
  if (missing(data)) {
    results_plot <- data.frame(
      Item = results_df$Item,
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
        point_interval = .sim_point_interval(spec),
        layout = "weave",
        slab_color = NA,
        .width = c(0.66, spec$width)
      ) +
      ggplot2::labs(
        x = "Partial gamma",
        y = "Item",
        caption = er2_caption(paste0(
          "Results from ",
          actual_iterations,
          " simulated datasets. ",
          sample_clause,
          " per dataset.\n",
          interval_note
        ))
      ) +
      ggplot2::scale_color_manual(
        values = scales::brewer_pal()(3)[-1],
        aesthetics = "slab_fill",
        guide = "none"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
      er2_axis_margins() +
      er2_plot_caption()

    return(p)
  }

  # --- Case 2: observed data supplied -----------------------------------------
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
      nrow(data),
      "), but has length ",
      length(dif_var),
      ".",
      call. = FALSE
    )
  }

  cc_dif <- stats::complete.cases(as.data.frame(data)) & !is.na(dif_var)
  .check_cutoff_sample(
    simfit$sample_n,
    sum(cc_dif),
    "RMdifGammaCutoff()",
    policy = "complete cases in data and dif_var",
    groups_cutoff = simfit$dif_group_sizes,
    groups_used = as.integer(table(droplevels(factor(dif_var[cc_dif]))))
  )

  obs <- .partgam_dif_gamma(as.data.frame(data), dif_var)
  observed_df <- data.frame(
    Item = obs$Item,
    observed_gamma = obs$gamma,
    stringsAsFactors = FALSE
  )

  # --- Build plot data --------------------------------------------------------
  gamma_sim <- data.frame(
    Item = results_df$Item,
    Value = results_df$gamma,
    stringsAsFactors = FALSE
  )
  gamma_sim <- merge(gamma_sim, observed_df, by = "Item", sort = FALSE)
  gamma_sim$Item <- factor(gamma_sim$Item, levels = item_levels)

  lo_hi$Item_f <- factor(lo_hi$Item, levels = item_levels)

  caption_text <- er2_caption(paste0(
    "Results from ",
    actual_iterations,
    " simulated datasets. ",
    sample_clause,
    " per dataset.\n",
    "Orange diamonds indicate observed partial gamma DIF. ",
    "Black dots indicate median gamma from simulations.\n",
    interval_note
  ))

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
      point_interval = .sim_point_interval(spec),
      layout = "weave",
      slab_color = NA,
      .width = c(0.66, spec$width)
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x = .data$min_gamma,
        xend = .data$max_gamma,
        y = .data$Item_f,
        yend = .data$Item_f
      ),
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x = .data$p66lo_gamma,
        xend = .data$p66hi_gamma,
        y = .data$Item_f,
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
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
    er2_axis_margins() +
    er2_plot_caption()

  p
}
