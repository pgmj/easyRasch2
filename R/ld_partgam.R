#' Partial Gamma Local Dependence Analysis
#'
#' Computes partial gamma coefficients for Local Dependence (LD) assessment
#' using \code{iarm::partgam_LD()}. Each pair of items is tested for residual
#' association, controlling for the rest score (total score minus one of the
#' items in the pair).
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed,
#'   but at least one complete case must exist.
#' @param cutoff Optional. Default `NULL` (no cutoff applied). Can be:
#'   * The return value of \code{\link{RMlocdepGammaCutoff}} (a list with
#'     `$pair_cutoffs`): the data.frame is extracted automatically and
#'     simulation metadata is included in the kable caption.
#'   * The `$pair_cutoffs` data.frame from \code{\link{RMlocdepGammaCutoff}}
#'     directly: must have columns `Item1`, `Item2`, `gamma_low`, `gamma_high`.
#'   When provided, adds columns `Gamma_low`, `Gamma_high`, and `Flagged`
#'   (logical; `TRUE` when the observed partial gamma falls outside the
#'   credible range) to the result.
#' @param p_value Logical or `NULL`. When `TRUE`, adds one-sided bootstrap
#'   p-values for *excess positive* local dependence (`p_gamma`,
#'   `padj_gamma`), matching the `p_value` semantics of
#'   \code{\link{RMlocdepQ3}}, and `flagged` reflects `padj_gamma < alpha`
#'   (positive deviations only) instead of the credible range. One test per
#'   item pair: the p-value is computed on `gamma_pair`, the larger of the two
#'   conditioning directions, and repeated in the direction-2 table for the
#'   same pair. The asymptotic adjusted p-value and star columns from
#'   `iarm::partgam_LD()` are **dropped** in this mode, and the simulated
#'   `gamma_low` / `gamma_high` band is kept as the effect-size reference.
#'   `NULL`, the default, means `TRUE` when `cutoff` is the **full**
#'   \code{\link{RMlocdepGammaCutoff}} object (it carries the simulated
#'   distributions in `$results`) and `FALSE` otherwise, so calling with no
#'   cutoff keeps the asymptotic table unchanged. Pass `FALSE` for the
#'   pre-1.2.0 behaviour. The interval makes a poor decision rule, since its
#'   width sets a family-wise error rate of `1 - width^m` over all \eqn{m}
#'   pairs at once (Johansson, 2026).
#' @param correction Character. Multiplicity correction for the bootstrap
#'   p-values, applied over the family of all item pairs (before any
#'   `n_pairs` display filter): `"fwer"` (default; Westfall-Young
#'   studentised-max step-down), `"fdr_bh"`, `"fdr_by"`, or `"none"`. Ignored
#'   when `p_value = FALSE`.
#' @param alpha Numeric in (0, 1). Significance level used to flag pairs on
#'   the corrected p-value. Default `0.05`. Ignored when `p_value = FALSE`.
#' @param output Character string controlling the return value. Either
#'   `"kable"` (default) for a formatted `knitr::kable()` table, or
#'   `"dataframe"` for the underlying data.frame.
#' @param n_pairs Optional positive integer. When supplied, only the
#'   `n_pairs` item pairs with the largest absolute partial-gamma values
#'   (i.e., strongest residual dependence in either direction) are
#'   retained per rest-score direction, sorted by `|gamma|` descending.
#'   When `NULL` (default), all pairs are returned in `iarm`'s native
#'   ordering. Values larger than the total number of pairs are silently
#'   capped at that total.
#'
#' @return
#' * If `output = "kable"`: an object of class `"RMlocdepGamma"`. Internally
#'   a list with two `knitr_kable` elements, `$direction1` and
#'   `$direction2`. In both, the rest score is the total score minus
#'   Item 2 (the second column); the two elements list each item pair in
#'   the two possible orders, so together they cover both rest-score
#'   directions for every pair. Each has columns "Item 1", "Item 2",
#'   "Partial gamma",
#'   "Adj. p-value (BH)", and "p-value sign." (a star-string indicator
#'   from `iarm::partgam_LD()`). When `cutoff` is provided, additional
#'   columns "Gamma low", "Gamma high", and "Flagged" are included.
#'
#'   The object has custom `print()` and `knitr::knit_print()` methods:
#'   in the R console it prints the two tables stacked vertically; in a
#'   Quarto / R Markdown chunk it renders as two distinct pipe tables.
#'   Access the individual tables explicitly as `result$direction1` and
#'   `result$direction2` if needed.
#' * If `output = "dataframe"`: a named list of two data.frames
#'   (`$direction1`, `$direction2`) with columns `Item1`, `Item2`,
#'   `gamma`, `se`, `lower`, `upper` (95% Wald CI), `padj_bh`,
#'   `Significance`. When `cutoff` is provided, columns `gamma_low`,
#'   `gamma_high`, and `flagged` are also included.
#'   With `p_value = TRUE`, `padj_bh` and `Significance` are replaced by
#'   `p_gamma` and `padj_gamma` (identical for a pair in both directions).
#'
#' @details
#' Partial gamma (Christensen, Kreiner & Mesbah, 2013) measures the residual
#' association between pairs of items after controlling for the rest score
#' (total score minus one item). Because it matters which item is subtracted,
#' calculations are done for each pair in both directions, yielding two
#' data.frames.
#'
#' Values near 0 indicate no local dependence. Large positive values suggest
#' positive LD (items share variance beyond the latent trait), while large
#' negative values suggest negative LD.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' \strong{Bootstrap p-values.} When `p_value = TRUE`, each pair's observed
#' partial gamma (canonical direction) is compared against its simulated null
#' distribution (from `cutoff$results`, simulated under local independence).
#' The per-pair statistic is the residual studentised by the bootstrap mean
#' and SD; the marginal p-value is the one-sided Monte-Carlo p-value
#' `(1 + #\{t* >= t\}) / (B + 1)` for excess *positive* LD (redundancy, the
#' diagnostic target — matching \code{\link{RMlocdepQ3}}), so it can be no
#' smaller than `1 / (B + 1)`. The band still shows both bounds for
#' reference. `correction = "fwer"` uses the Westfall-Young studentised-max
#' step-down over the family of all pairs, which exploits the bootstrap
#' dependence among them (Ferreira, 2024); it is liberal when the simulation
#' is small, so at least 1000 `iterations` in [RMlocdepGammaCutoff()] are
#' recommended (a warning is issued below that). Unlike the asymptotic
#' p-values from `iarm::partgam_LD()`, these are calibrated against the
#' *simulated Rasch null* rather than the asymptotic SE; they are
#' model-conditional and sample-size-sensitive, and are reported alongside
#' the simulated effect-size band, not in place of it.
#'
#' @inheritSection RMitemInfit Multiple comparisons
#'
#' @references
#' Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013).
#' \emph{Rasch Models in Health}, pp. 133--135. ISTE & Wiley.
#' \doi{10.1002/9781118574454}
#'
#' Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate' number
#' of hypotheses simultaneously. *Journal of Statistical Theory and Practice,
#' 19*(6). \doi{10.1007/s42519-024-00412-4}
#'
#' Westfall, P. H., & Young, S. S. (1993). *Resampling-Based Multiple Testing*.
#' Wiley.
#'
#' @seealso \code{\link{RMlocdepGammaCutoff}}, \code{\link{RMlocdepGammaPlot}}
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
#'
#'   # Default kable output
#'   RMlocdepGamma(sim_data)
#'
#'   # Return as data.frame list
#'   RMlocdepGamma(sim_data, output = "dataframe")
#'
#'   # Simulation-based cutoffs (slow): 100+ Monte-Carlo iterations
#'   if (requireNamespace("ggdist", quietly = TRUE)) {
#'     cutoff_res <- RMlocdepGammaCutoff(sim_data, iterations = 100,
#'                                       parallel = FALSE, seed = 42)
#'     RMlocdepGamma(sim_data, cutoff = cutoff_res)
#'
#'     # Bootstrap p-values with family-wise (Westfall-Young) correction
#'     # (use iterations >= 1000 in real analyses for stable p-values)
#'     RMlocdepGamma(sim_data, cutoff = cutoff_res, p_value = TRUE,
#'                   output = "dataframe")
#'   }
#' }
#' }
RMlocdepGamma <- function(
  data,
  cutoff = NULL,
  p_value = NULL,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  alpha = 0.05,
  output = "kable",
  n_pairs = NULL
) {
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required for RMlocdepGamma() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  output <- match.arg(output, c("kable", "dataframe"))
  correction <- match.arg(correction)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }

  # --- Validate n_pairs -------------------------------------------------------
  if (!is.null(n_pairs)) {
    if (
      !is.numeric(n_pairs) ||
        length(n_pairs) != 1L ||
        !is.finite(n_pairs) ||
        n_pairs < 1 ||
        n_pairs != as.integer(n_pairs)
    ) {
      stop(
        "`n_pairs` must be a single positive integer or NULL.",
        call. = FALSE
      )
    }
    n_pairs <- as.integer(n_pairs)
  }

  validate_response_data(data)

  # --- Validate and normalise cutoff ------------------------------------------
  cutoff_n_iter <- NULL
  cutoff_method <- NULL
  cutoff_hdci_width <- NULL
  cutoff_full <- NULL # full object (carries simulated $results for p-values)
  if (!is.null(cutoff)) {
    if (
      is.list(cutoff) &&
        !is.data.frame(cutoff) &&
        "pair_cutoffs" %in% names(cutoff)
    ) {
      cutoff_full <- cutoff
      cutoff_n_iter <- cutoff$actual_iterations
      cutoff_method <- cutoff$cutoff_method
      cutoff_hdci_width <- cutoff$hdci_width
      cutoff <- cutoff$pair_cutoffs
    }
    if (!is.data.frame(cutoff)) {
      stop(
        "`cutoff` must be NULL, the return value of RMlocdepGammaCutoff(), or its ",
        "$pair_cutoffs data.frame.",
        call. = FALSE
      )
    }
    required_cols <- c("Item1", "Item2", "gamma_low", "gamma_high")
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
  # The interval describes where a pair's coefficient is expected to fall and
  # makes a poor decision rule, because its width sets a family-wise error rate
  # over every pair at once (Johansson, 2026). The bare $pair_cutoffs, or no
  # cutoff at all, leaves nothing to compute a p-value from and resolves to
  # FALSE, which keeps the asymptotic path this function shows by default
  # exactly as it was.
  if (!is.null(p_value) && (!is.logical(p_value) || length(p_value) != 1L)) {
    stop("`p_value` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }
  have_sims <- !is.null(cutoff_full) && !is.null(cutoff_full$results)
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
        "`p_value = TRUE` requires the full RMlocdepGammaCutoff() object (it ",
        "carries the simulated per-pair distributions in $results); a NULL ",
        "cutoff or the bare $pair_cutoffs data.frame is not sufficient.",
        call. = FALSE
      )
    }
    # Below 400 the correction itself is off. Between 400 and 1000 only
    # reproducibility improves, which the table caption reports instead.
    if (!is.null(cutoff_n_iter) && cutoff_n_iter < 400L) {
      .notify_low_iterations(
        cutoff_n_iter,
        cutoff_full$requested_iterations,
        fn = "RMlocdepGammaCutoff()",
        id = "easyRasch2_low_iterations_locdep"
      )
    }
    if (!is.null(n_pairs_total)) {
      .warn_fdr_floor(
        cutoff_n_iter,
        n_pairs_total,
        correction,
        alpha,
        unit = "item pairs",
        fn = "RMlocdepGammaCutoff()"
      )
    }
  } else if (!is.null(cutoff_full) && !is.null(n_pairs_total)) {
    .notify_band_flagging(
      if (identical(cutoff_method, "quantile")) 0.95 else cutoff_hdci_width,
      n_pairs_total,
      unit = "item pairs",
      fn = "RMlocdepGammaCutoff()",
      id = "easyRasch2_band_flagging_locdep"
    )
  }

  # --- rgl workaround ---------------------------------------------------------
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # --- Compute partial gamma LD via iarm --------------------------------------
  sink(nullfile())
  pgam_raw <- iarm::partgam_LD(as.data.frame(data))
  sink()

  # pgam_raw is a list of two data.frames (one per rest-score direction).
  # Columns from iarm::partgam_LD():
  #   1 Item1, 2 Item2, 3 gamma, 4 se, 5 pvalue, 6 padj.BH, 7 sig,
  #   8 lower, 9 upper. The `sig` column is a star string like " ***" /
  #   " **" / " *" / "  ." / "    "; we trim whitespace for display.
  process_pgam_df <- function(raw_df) {
    df <- data.frame(
      Item1 = as.character(raw_df$Item1),
      Item2 = as.character(raw_df$Item2),
      gamma = as.numeric(raw_df$gamma),
      se = as.numeric(raw_df$se),
      lower = as.numeric(raw_df$lower),
      upper = as.numeric(raw_df$upper),
      padj_bh = as.numeric(raw_df[[6]]),
      Significance = trimws(as.character(raw_df[[7]])),
      stringsAsFactors = FALSE
    )
    df
  }

  result_list <- list(
    direction1 = process_pgam_df(pgam_raw[[1]]),
    direction2 = process_pgam_df(pgam_raw[[2]])
  )

  # --- The tested pair statistic ----------------------------------------------
  # A pair is tested once, on the larger of its two conditioning directions,
  # because local dependence violates both of the conditional independence
  # hypotheses the Rasch model implies for it. The `gamma` column of each table
  # remains that direction's own coefficient, which is what a reader wants to
  # see, but every decision below is taken on `gamma_pair`, so no p-value or
  # flag is ever attributed to a coefficient it was not computed from.
  data_complete <- data[stats::complete.cases(data), , drop = FALSE]
  .obs1 <- .partgam_ld_gamma(data_complete, direction = 1L)
  .obs2 <- .partgam_ld_gamma(data_complete, direction = 2L)
  .pkey <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "___")
  gamma_pair <- stats::setNames(
    pmax(.obs1$gamma,
         .obs2$gamma[match(.pkey(.obs1$Item1, .obs1$Item2),
                           .pkey(.obs2$Item1, .obs2$Item2))]),
    .pkey(.obs1$Item1, .obs1$Item2)
  )
  for (idx in seq_along(result_list)) {
    df <- result_list[[idx]]
    df$gamma_pair <- as.numeric(gamma_pair[.pkey(df$Item1, df$Item2)])
    result_list[[idx]] <- df
  }

  # --- Apply cutoff if provided -----------------------------------------------
  # Cutoffs are keyed by direction-1 pairs (i < j), so for direction 2 (i > j)

  # we must look up by the canonical (sorted) pair order.
  if (!is.null(cutoff)) {
    # Build a canonical key on the cutoff table
    cutoff$canonical_key <- paste(
      pmin(cutoff$Item1, cutoff$Item2),
      pmax(cutoff$Item1, cutoff$Item2),
      sep = "___"
    )
    cutoff_sub <- cutoff[, c("canonical_key", "gamma_low", "gamma_high")]

    for (idx in seq_along(result_list)) {
      result_df <- result_list[[idx]]
      # Build the same canonical key on the result (sorted pair order)
      result_df$canonical_key <- paste(
        pmin(result_df$Item1, result_df$Item2),
        pmax(result_df$Item1, result_df$Item2),
        sep = "___"
      )

      merged <- merge(
        result_df,
        cutoff_sub,
        by = "canonical_key",
        all.x = TRUE,
        sort = FALSE
      )
      # Restore original row order
      merged <- merged[match(result_df$canonical_key, merged$canonical_key), ]
      rownames(merged) <- NULL
      # Flagged on the pair statistic, which is what the cutoff band describes.
      merged$flagged <- !is.na(merged$gamma_low) &
        (merged$gamma_pair < merged$gamma_low |
           merged$gamma_pair > merged$gamma_high)

      # Remove helper column
      merged$canonical_key <- NULL

      merged <- merged[, c(
        "Item1",
        "Item2",
        "gamma",
        "gamma_pair",
        "se",
        "lower",
        "upper",
        "padj_bh",
        "Significance",
        "gamma_low",
        "gamma_high",
        "flagged"
      )]
      result_list[[idx]] <- merged
    }
    # Clean up cutoff helper column
    cutoff$canonical_key <- NULL
  }

  # --- Bootstrap p-values (one test per pair, canonical direction) ------------
  # Computed BEFORE the n_pairs display filter so the multiplicity correction
  # always runs over the full family of pairs. The simulated null holds the
  # direction-1 gammas only, so the observed statistic is the direction-1
  # (canonical) gamma; the resulting p-value is keyed by the unordered pair
  # and repeated in the direction-2 table.
  if (p_value) {
    canon_key <- function(a, b) {
      paste(pmin(a, b), pmax(a, b), sep = "___")
    }
    sim_res <- cutoff_full$results
    sim_key <- canon_key(sim_res$Item1, sim_res$Item2)
    obs_key <- canon_key(
      as.character(pgam_raw[[1L]]$Item1),
      as.character(pgam_raw[[1L]]$Item2)
    )
    if (!setequal(obs_key, unique(sim_key))) {
      stop(
        "Item pairs in the cutoff simulations ($results) do not match the ",
        "pairs in `data`.",
        call. = FALSE
      )
    }
    sim_mat <- tapply(
      sim_res$gamma,
      list(sim_res$iteration, sim_key),
      function(x) x[1L]
    )
    # The tested statistic is the pair maximum computed above, from the same
    # code path that produced the simulated null, so the two cannot diverge.
    observed <- gamma_pair[obs_key]
    # One-sided: excess positive LD (redundancy), matching RMlocdepQ3.
    pv <- .bootstrap_pvalues(
      observed,
      sim_mat,
      correction = correction,
      tail = "upper"
    )
    p_lookup <- stats::setNames(pv$p, pv$name)
    padj_lookup <- stats::setNames(pv$padj, pv$name)

    for (idx in seq_along(result_list)) {
      df <- result_list[[idx]]
      key <- canon_key(df$Item1, df$Item2)
      df$p_gamma <- as.numeric(p_lookup[key])
      df$padj_gamma <- as.numeric(padj_lookup[key])
      df$flagged <- !is.na(df$padj_gamma) & df$padj_gamma < alpha
      # Drop the asymptotic p-value pair; the bootstrap p-values replace it.
      df <- df[, c(
        "Item1",
        "Item2",
        "gamma",
        "gamma_pair",
        "se",
        "lower",
        "upper",
        "gamma_low",
        "gamma_high",
        "p_gamma",
        "padj_gamma",
        "flagged"
      )]
      result_list[[idx]] <- df
    }
  }

  # --- Top-N filter by |gamma| per direction ----------------------------------
  total_pairs <- nrow(result_list[[1L]])
  filter_applied <- FALSE
  if (!is.null(n_pairs)) {
    keep_n <- min(n_pairs, total_pairs)
    if (keep_n < total_pairs) {
      filter_applied <- TRUE
      for (idx in seq_along(result_list)) {
        df <- result_list[[idx]]
        ord <- order(abs(df$gamma), decreasing = TRUE)
        df <- df[ord[seq_len(keep_n)], , drop = FALSE]
        rownames(df) <- NULL
        result_list[[idx]] <- df
      }
    }
  }

  # --- Return -----------------------------------------------------------------
  if (output == "dataframe") {
    return(result_list)
  }

  # Kable display rounding (the dataframe output above stays unrounded).
  # The se / lower / upper columns are dataframe-only (added for
  # programmatic use, e.g. the jamovi module); the kable keeps its
  # previous column set.
  ld_digits <- c(
    gamma = 3, padj_bh = 3, gamma_low = 3, gamma_high = 3,
    p_gamma = 4, padj_gamma = 4
  )
  result_list <- lapply(result_list, function(d) {
    .round_display(d[, setdiff(names(d), c("se", "lower", "upper")),
                     drop = FALSE], digits = ld_digits)
  })

  # Build caption
  n_complete <- sum(stats::complete.cases(as.data.frame(data)))
  n_clause <- .n_caption(
    n_complete,
    nrow(as.data.frame(data)),
    if (anyNA(as.data.frame(data))) "complete cases" else character()
  )
  filter_suffix <- if (filter_applied) {
    paste0(" Showing top ", n_pairs, " of ", total_pairs, " pairs by |gamma|.")
  } else {
    ""
  }
  if (p_value) {
    caption_text <- paste0(
      "Partial gamma LD analysis. ",
      n_clause,
      ". One-sided bootstrap p-values for excess positive LD from ",
      cutoff_n_iter,
      " iterations, computed on gamma_pair, the larger of the two ",
      "conditioning directions, and repeated across both tables (replacing ",
      "the asymptotic p-values). Multiplicity correction: ",
      .correction_label(correction),
      ". Flagged at padj < ",
      alpha,
      ". p-values cannot be smaller than 1/(",
      cutoff_n_iter,
      "+1) = ",
      round(1 / (cutoff_n_iter + 1), 4),
      ". The interval is shown as description and is not the decision rule, ",
      "so a pair below the lower bound is not flagged.",
      .iteration_note(cutoff_n_iter),
      .attrition_clause(cutoff_n_iter, cutoff_full$requested_iterations),
      filter_suffix
    )
  } else if (is.null(cutoff)) {
    caption_text <- paste0(
      "Partial gamma LD analysis. ",
      n_clause,
      ". Positive gamma indicates positive local dependence between items.",
      filter_suffix
    )
  } else if (!is.null(cutoff_n_iter)) {
    method_label <- .format_gamma_cutoff_method_label(
      cutoff_method,
      cutoff_hdci_width
    )
    iter_part <- paste0(cutoff_n_iter, " simulation iterations")
    caption_text <- paste0(
      "Partial gamma LD analysis. ",
      n_clause,
      ". Cutoff values based on ",
      if (!is.null(method_label)) {
        paste0(iter_part, " (", method_label, ").")
      } else {
        paste0(iter_part, ".")
      },
      .band_error_clause(
        if (identical(cutoff_method, "quantile")) 0.95 else cutoff_hdci_width,
        total_pairs,
        unit = "item pairs"
      ),
      .attrition_clause(cutoff_n_iter, cutoff_full$requested_iterations),
      filter_suffix
    )
  } else {
    caption_text <- paste0(
      "Partial gamma LD analysis. ",
      n_clause,
      ". Simulation-based cutoff values applied.",
      filter_suffix
    )
  }

  # One header per displayed column, in the order the tables carry them. The
  # `gamma_pair` column (the larger of the two conditioning directions, and the
  # statistic that is tested) sits fourth once a cutoff is supplied and last on
  # the asymptotic path, so the three vectors are not interchangeable.
  col_names_no_cutoff <- c(
    "Item 1",
    "Item 2",
    "Partial gamma",
    "Adj. p-value (BH)",
    "p-value sign.",
    "Gamma pair"
  )
  col_names_cutoff <- c(
    "Item 1",
    "Item 2",
    "Partial gamma",
    "Gamma pair",
    "Adj. p-value (BH)",
    "p-value sign.",
    "Gamma low",
    "Gamma high",
    "Flagged"
  )
  col_names_pvalue <- c(
    "Item 1",
    "Item 2",
    "Partial gamma",
    "Gamma pair",
    "Gamma low",
    "Gamma high",
    "p",
    "p (adj)",
    "Flagged"
  )
  col_names_used <- if (p_value) {
    col_names_pvalue
  } else if (is.null(cutoff)) {
    col_names_no_cutoff
  } else {
    col_names_cutoff
  }

  kable1 <- knitr::kable(
    result_list$direction1,
    format = "pipe",
    col.names = col_names_used,
    caption = paste0(caption_text, " Direction 1: rest score = total - Item2.")
  )

  kable2 <- knitr::kable(
    result_list$direction2,
    format = "pipe",
    col.names = col_names_used,
    caption = paste0(
      caption_text,
      " Direction 2: rest score = total - Item2 ",
      "(item pairs shown in the reverse order to ",
      "direction 1)."
    )
  )

  # `knitr::kable(format = "pipe")` returns a multi-line *character vector*
  # (one element per line), so `paste(kable1, "\n\n", kable2)` would do
  # element-wise concatenation, interleaving the two tables row-by-row.
  # Collapse each table to a single string first, then join with a
  # blank line between to render as two distinct kable tables.
  combined <- paste(
    paste(kable1, collapse = "\n"),
    paste(kable2, collapse = "\n"),
    sep = "\n\n"
  )

  # Return a custom-class list so the two display contexts (R console and
  # knitr chunk) can be handled separately. See print.RMlocdepGamma() and
  # knit_print.RMlocdepGamma() below.
  out <- list(
    direction1 = kable1,
    direction2 = kable2,
    .combined = combined
  )
  class(out) <- c("RMlocdepGamma", "list")
  out
}

#' Print method for RMlocdepGamma kable output
#'
#' Prints the two rest-score direction tables stacked vertically with a
#' blank line between them. Each table renders via `knitr_kable`'s own
#' print method as a clean pipe-markdown table.
#'
#' @param x An object of class `"RMlocdepGamma"` returned by
#'   \code{\link{RMlocdepGamma}} with `output = "kable"`.
#' @param ... Further arguments (currently unused).
#' @return Invisibly returns `x`.
#' @keywords internal
#' @exportS3Method base::print
print.RMlocdepGamma <- function(x, ...) {
  print(x$direction1)
  cat("\n")
  print(x$direction2)
  invisible(x)
}

#' knitr knit_print method for RMlocdepGamma kable output
#'
#' Inside a knitr / Quarto / R Markdown chunk, returns the pre-combined
#' two-table asis string so pandoc renders them as two distinct pipe
#' tables. Outside knitr, R's normal dispatch falls back to
#' `print.RMlocdepGamma()`.
#'
#' @param x An object of class `"RMlocdepGamma"`.
#' @param ... Further arguments passed to `knitr::asis_output()`.
#' @return A `knit_asis` character object.
#' @keywords internal
#' @exportS3Method knitr::knit_print
knit_print.RMlocdepGamma <- function(x, ...) {
  knitr::asis_output(x$.combined)
}

#' Simulation-Based Partial Gamma LD Cutoff Determination
#'
#' Uses parametric bootstrap simulation to determine appropriate cutoff values
#' for partial gamma Local Dependence analysis via
#' \code{\link[iarm]{partgam_LD}}. Under a correctly fitting Rasch model where
#' items are locally independent, this function generates the expected
#' distribution of partial gamma values per item pair, providing empirical
#' critical values.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Only complete cases (rows without
#'   any `NA`) are used.
#' @param iterations Integer. Number of simulation iterations (default 400,
#'   was 250 before 1.2.0). 400 is the calibrated floor for the Westfall-Young
#'   correction (Johansson, 2026) and the count a 95\% interval needs to
#'   converge. Use 1000 to 2000 for a final analysis.
#' @param parallel Logical. Use parallel processing via `mirai` if available
#'   (default `TRUE`).
#' @param n_cores Integer or `NULL`. Number of parallel workers. When `NULL`,
#'   `getOption("mc.cores")` is checked first. If neither is set and
#'   `parallel = TRUE`, a warning is issued and execution falls back to
#'   sequential (single core) processing.
#' @param verbose Logical. Show a progress bar (default `FALSE`).
#' @param seed Integer or `NULL`. Random seed for reproducibility. See
#'   [easyRasch2-reproducibility] for what this guarantees and how it
#'   interacts with `parallel`.
#' @param cutoff_method Character string specifying how cutoff intervals are
#'   computed. Either `"hdci"` (default) for the Highest Density Interval via
#'   `ggdist::hdci()`, or `"quantile"` for the 2.5th/97.5th percentiles via
#'   `stats::quantile()`.
#' @param hdci_width Numeric. Width of the HDCI when `cutoff_method = "hdci"`.
#'   Default is `0.95` (95\% HDCI), was `0.99` before 1.2.0. The interval
#'   describes where a fitting pair's coefficient is expected to fall and is
#'   no longer the default decision rule, so the width is chosen to converge
#'   at the default iteration count rather than to imply an error rate.
#'   Ignored when `cutoff_method = "quantile"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`results`}{data.frame with columns `iteration`, `Item1`, `Item2`,
#'     and `gamma` (one row per item pair per successful iteration). Contains
#'     results from direction 1 only (rest score = total - Item2), which is
#'     the conventional direction.}
#'   \item{`pair_cutoffs`}{data.frame with per-pair cutoff summaries: `Item1`,
#'     `Item2`, `gamma_low`, `gamma_high`. Bounds are computed using the method
#'     specified by `cutoff_method`.}
#'   \item{`actual_iterations`}{Number of successful iterations.}
#'   \item{`sample_n`}{Number of complete cases used.}
#'   \item{`sample_n_total`}{Number of respondents in the raw input data,
#'     before the complete-case filter.}
#'   \item{`sample_has_na`}{Logical. Whether the raw input data contained
#'     any missing values.}
#'   \item{`sample_summary`}{Summary statistics of estimated person
#'     parameters.}
#'   \item{`item_names`}{Character vector of item names from data.}
#'   \item{`cutoff_method`}{The method used to compute cutoffs (`"hdci"` or
#'     `"quantile"`).}
#'   \item{`hdci_width`}{The HDCI width used (only meaningful when
#'     `cutoff_method = "hdci"`).}
#' }
#'
#' @details
#' For each simulation iteration the function:
#' \enumerate{
#'   \item Resamples person parameters (thetas) with replacement from the
#'     WLE person locations.
#'   \item Simulates item response data under a Rasch model (dichotomous via
#'     `psychotools::rrm()` or polytomous via an internal partial credit
#'     simulator).
#'   \item Computes partial gamma for every item pair in the canonical
#'     rest-score direction. The coefficients are identical to those of
#'     `iarm::partgam_LD()`, but are computed by a vectorised internal, since
#'     `iarm` also derives the asymptotic standard error and confidence
#'     interval that a simulated null does not need and costs roughly two
#'     orders of magnitude more per iteration.
#' }
#'
#' Because the data are simulated under the Rasch model, items are locally
#' independent by construction. The distribution of partial gamma values
#' across iterations provides empirical critical values per item pair. Values
#' from real data that fall outside these bounds suggest local dependence that
#' exceeds what would be expected by chance. Failed iterations (e.g., due to
#' convergence issues or degenerate data) are silently discarded.
#'
#' The generating model uses CML item thresholds via `psychotools::pcmodel()`
#' (a dichotomous item is a 2-category PCM) and WLE person locations,
#' consistent with the rest of the package; responses are simulated with
#' `psychotools::rrm()` (dichotomous) or an internal partial credit score
#' simulator (polytomous).
#'
#' Parallel processing is provided by the `mirai` package (optional). Install
#' it with `install.packages("mirai")` to enable parallelisation.
#'
#' The `iarm` package must be installed (it is in Suggests, not Imports).
#'
#' @references
#' Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013).
#' \emph{Rasch Models in Health}, pp. 133--135. ISTE & Wiley.
#' \doi{10.1002/9781118574454}
#'
#' @seealso \code{\link[iarm]{partgam_LD}}, \code{\link{RMlocdepGamma}},
#'   \code{\link{RMlocdepGammaPlot}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("iarm", quietly = TRUE) &&
#'     requireNamespace("ggdist", quietly = TRUE)) {
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:10)
#'
#'   # Run 100 iterations sequentially for a quick demo
#'   cutoff_res <- RMlocdepGammaCutoff(sim_data, iterations = 100,
#'                                     parallel = FALSE, seed = 42)
#'   cutoff_res$pair_cutoffs
#' }
#' }
RMlocdepGammaCutoff <- function(
  data,
  iterations = 400,
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL,
  cutoff_method = "hdci",
  hdci_width = 0.95
) {
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
      "Package 'iarm' is required for RMlocdepGammaCutoff() but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }

  validate_response_data(data)

  # rgl workaround (iarm depends on vcdExtra -> rgl)
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # Only complete cases. Record the raw total and whether anything was
  # dropped so callers (e.g. RMlocdepGammaPlot) can report the sample in the
  # standard `n = X of Y respondents` form.
  n_total <- nrow(data)
  has_na <- anyNA(data)
  complete_idx <- stats::complete.cases(data)
  data <- data[complete_idx, , drop = FALSE]

  if (nrow(data) == 0L) {
    stop(
      "No complete cases in data after removing rows with NA.",
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
  sample_n <- nrow(data_mat)
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  item_names_vec <- colnames(data_mat)

  # Generating model: CML item thresholds (psychotools) + WLE person locations,
  # consistent with the rest of the package, replacing eRm CML +
  # eRm::person.parameter() (MLE). Thetas form the pool resampled with
  # replacement to build each simulated dataset; a dichotomous item is a
  # 2-category PCM, so its centred threshold is the item difficulty for rrm().
  pool <- .wle_theta_pool(data_mat)
  thr_list <- pool$thr_list
  thetas <- pool$thetas

  if (is_polytomous) {
    sim_data_list <- list(
      type = "polytomous",
      thetas = thetas,
      deltaslist = thr_list,
      n_items = ncol(data_mat),
      sample_n = sample_n,
      item_names = item_names_vec
    )
  } else {
    sim_data_list <- list(
      type = "dichotomous",
      thetas = thetas,
      item_params = unlist(thr_list, use.names = FALSE),
      n_items = ncol(data_mat),
      sample_n = sample_n,
      item_names = item_names_vec
    )
  }

  if (use_parallel) {
    results_raw <- run_partgam_LD_sim_parallel(
      iterations,
      sim_seeds,
      sim_data_list,
      n_cores,
      verbose
    )
  } else {
    results_raw <- run_partgam_LD_sim_sequential(
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

  # Compute per-pair cutoffs
  pair_keys <- unique(paste(results_df$Item1, results_df$Item2, sep = "___"))
  pair_cutoffs <- do.call(
    rbind,
    lapply(pair_keys, function(pk) {
      parts <- strsplit(pk, "___", fixed = TRUE)[[1]]
      sub <- results_df[
        results_df$Item1 == parts[1] & results_df$Item2 == parts[2],
      ]
      if (cutoff_method == "hdci") {
        gamma_interval <- ggdist::hdci(sub$gamma, .width = hdci_width)
        data.frame(
          Item1 = parts[1],
          Item2 = parts[2],
          gamma_low = gamma_interval[1L, 1L],
          gamma_high = gamma_interval[1L, 2L],
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      } else {
        data.frame(
          Item1 = parts[1],
          Item2 = parts[2],
          gamma_low = stats::quantile(sub$gamma, 0.025, na.rm = TRUE),
          gamma_high = stats::quantile(sub$gamma, 0.975, na.rm = TRUE),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
    })
  )
  rownames(pair_cutoffs) <- NULL

  list(
    results = results_df,
    pair_cutoffs = pair_cutoffs,
    actual_iterations = actual_iterations,
    requested_iterations = iterations,
    sample_n = sample_n,
    sample_n_total = n_total,
    sample_has_na = has_na,
    sample_summary = summary(thetas),
    item_names = item_names_vec,
    cutoff_method = cutoff_method,
    hdci_width = hdci_width
  )
}

# ---------------------------------------------------------------------------
# Internal: vectorised partial gamma
# ---------------------------------------------------------------------------

#' Partial gamma for every item pair, coefficient only
#'
#' Vectorised implementation of Davis's (1967) partial gamma. It returns the
#' coefficient and nothing else, which is all a parametric bootstrap needs.
#' `iarm::partgam_LD()` additionally computes the Goodman-Kruskal delta-method
#' variance and loops over pairs, cells and strata in R, costing roughly 300 ms
#' per call irrespective of sample size, so calling it once per bootstrap
#' iteration dominates everything else by an order of magnitude.
#'
#' Used for the simulated null in [RMlocdepGammaCutoff()] and for the observed
#' statistic that [RMlocdepGamma()] tests against it, so that the two cannot
#' diverge. The `gamma`, `se`, `lower` and `upper` columns shown to users still
#' come from `iarm::partgam_LD()`. Agreement with `iarm` is exact and is
#' asserted in `tests/testthat/test-ld_partgam_gamma.R`.
#'
#' @details
#' Partial gamma pools concordant and discordant pair counts over strata of the
#' conditioning variable, here the rest score:
#' \deqn{\gamma = \frac{\sum_k C_k - \sum_k D_k}{\sum_k C_k + \sum_k D_k}.}
#' For a stratum with an \eqn{m \times m} count matrix \eqn{N} and \eqn{G} the
#' strictly-upper indicator (\eqn{G_{ab} = 1} iff \eqn{b > a}), writing
#' \eqn{A = GN} gives \eqn{A_{ij'} = \sum_{i' > i} N_{i'j'}}, and then
#' \eqn{C = A G^{T}} sums over \eqn{j' > j} while \eqn{D = A G} sums over
#' \eqn{j' < j}. Each unordered observation pair is counted once, so the
#' halving `iarm` applies is not needed here.
#'
#' @param data data.frame or matrix of item responses scored from 0, with no
#'   missing values. `iarm::partgam_LD()` applies `complete.cases()` internally;
#'   this function does not, so the caller must filter first.
#' @param direction `1` enumerates pairs with `Item1` before `Item2` in column
#'   order, `2` the reverse. The rest score always excludes `Item2`, so the two
#'   directions give the two conditional independence hypotheses of Kreiner and
#'   Christensen (2004). Direction 1 is the canonical one stored by
#'   [RMlocdepGammaCutoff()].
#' @param strata When `TRUE`, adds the stratum sign-homogeneity columns
#'   described in `.partgam_strata_summary()`. Off by default, since the
#'   bootstrap needs the coefficient alone and calls this once per iteration.
#' @return data.frame with `Item1`, `Item2` and `gamma`, one row per pair, in
#'   the same order as `iarm::partgam_LD()[[direction]]`. `gamma` is `NA_real_`
#'   for a pair with no concordant and no discordant observations, which
#'   `iarm::partgam_LD()` cannot return at all because it errors on the whole
#'   data set when an item is constant. With `strata = TRUE` the columns
#'   `n_strata`, `n_pos`, `n_neg`, `n_zero`, `homogeneous` and `w_opposing` are
#'   appended, and the per-stratum detail behind them is attached as the
#'   `"strata"` attribute, a list with one element per pair.
#' @keywords internal
#' @noRd
.partgam_ld_gamma <- function(data, direction = 1L, strata = FALSE) {
  X <- as.matrix(data)
  storage.mode(X) <- "integer"

  if (ncol(X) < 2L) {
    stop("`data` must have at least two items.", call. = FALSE)
  }
  if (anyNA(X)) {
    stop(
      "`data` must not contain missing values; filter to complete cases first.",
      call. = FALSE
    )
  }

  items <- colnames(X)
  if (is.null(items)) {
    items <- paste0("V", seq_len(ncol(X)))
  }
  k <- ncol(X)
  m <- max(X) + 1L
  score <- rowSums(X)

  # G[a, b] = 1 iff b > a. Its transpose is the strictly-lower counterpart.
  G <- outer(seq_len(m), seq_len(m), function(a, b) as.numeric(b > a))
  tG <- t(G)

  # iarm enumerates with i in the outer loop and j in the inner one, sending
  # i < j to the first table and i > j to the second. Reproduced here so the
  # row order matches.
  pairs <- do.call(
    rbind,
    lapply(seq_len(k), function(i) {
      js <- if (direction == 1L) {
        seq_len(k)[seq_len(k) > i]
      } else {
        seq_len(k)[seq_len(k) < i]
      }
      if (length(js) == 0L) NULL else cbind(i = i, j = js)
    })
  )

  if (!strata) {
    gammas <- vapply(
      seq_len(nrow(pairs)),
      function(p) {
        i <- pairs[p, "i"]
        j <- pairs[p, "j"]
        .partgam_one(X[, i], X[, j], score - X[, j], m, G, tG)
      },
      numeric(1L)
    )
    return(data.frame(
      Item1 = items[pairs[, "i"]],
      Item2 = items[pairs[, "j"]],
      gamma = gammas,
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  }

  detail <- lapply(seq_len(nrow(pairs)), function(p) {
    i <- pairs[p, "i"]
    j <- pairs[p, "j"]
    .partgam_one(X[, i], X[, j], score - X[, j], m, G, tG, strata = TRUE)
  })

  out <- data.frame(
    Item1 = items[pairs[, "i"]],
    Item2 = items[pairs[, "j"]],
    gamma = vapply(detail, function(d) d$gamma, numeric(1L)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- cbind(out, do.call(rbind, lapply(detail, .partgam_strata_summary)))
  attr(out, "strata") <- detail
  out
}

#' Partial gamma for one item pair
#'
#' @param x,y Integer response vectors scored from 0.
#' @param z Integer conditioning variable (the rest score).
#' @param m Number of response categories spanning `x` and `y`.
#' @param G,tG The strictly-upper indicator matrix and its transpose.
#' @param strata When `FALSE` (default) the coefficient is returned on its own,
#'   which is the path the bootstrap takes. When `TRUE` the stratum-level
#'   quantities behind it are returned as well.
#' @return With `strata = FALSE`, the partial gamma coefficient, or `NA_real_`
#'   when no pair of observations within a stratum is either concordant or
#'   discordant. With `strata = TRUE`, a list with that value as `gamma` plus
#'   `gamma_k` (per-stratum gamma, `NA_real_` where a stratum yields no
#'   concordant or discordant pair), `weight_k` (that stratum's contribution to
#'   the denominator, \eqn{C_k + D_k}) and `n_k` (stratum size). Strata are in
#'   ascending order of `z`, including any empty intermediate ones.
#' @keywords internal
#' @noRd
.partgam_one <- function(x, y, z, m, G, tG, strata = FALSE) {
  zc <- z - min(z) + 1L
  nz <- max(zc)
  counts <- tabulate(
    (zc - 1L) * m * m + y * m + x + 1L,
    nbins = m * m * nz
  )
  dim(counts) <- c(m, m, nz)

  conc <- 0
  disc <- 0
  if (strata) {
    gamma_k <- rep(NA_real_, nz)
    weight_k <- rep(0, nz)
    n_k <- rep(0L, nz)
  }

  for (s in seq_len(nz)) {
    N <- counts[, , s]
    n_s <- sum(N)
    if (strata) n_k[s] <- n_s
    # A stratum holding fewer than two observations contributes no pairs.
    if (n_s < 2L) next
    A <- G %*% N
    c_s <- sum(N * (A %*% tG))
    d_s <- sum(N * (A %*% G))
    conc <- conc + c_s
    disc <- disc + d_s
    if (strata) {
      weight_k[s] <- c_s + d_s
      # A stratum can hold observations yet no comparable pair, for instance
      # when every respondent in it gave the same answer to one of the items.
      if (c_s + d_s > 0) gamma_k[s] <- (c_s - d_s) / (c_s + d_s)
    }
  }

  total <- conc + disc
  gamma <- if (total == 0) NA_real_ else (conc - disc) / total

  if (!strata) {
    return(gamma)
  }
  list(gamma = gamma, gamma_k = gamma_k, weight_k = weight_k, n_k = n_k)
}

#' Summarise sign homogeneity across strata for one item pair
#'
#' Davis's partial gamma pools concordant and discordant counts over strata, so
#' it is interpretable as a partial correlation only when the stratum-specific
#' associations point the same way. Kreiner (personal communication, 2026)
#' states the condition directly: the stratum gammas need not be equal, but they
#' must be either positive or negative throughout, and if some are negative
#' while others are zero or positive the pooled value is not a meaningful
#' measure of partial correlation.
#'
#' Two cautions on reading the result. Stratum gammas estimated from a handful
#' of respondents change sign readily by chance, so at small sample sizes a
#' heterogeneous verdict is weak evidence of a real violation. And under the
#' null the stratum gammas are zero in the population, so mixed sample signs are
#' expected there and carry no meaning.
#'
#' @param st A list from `.partgam_one(strata = TRUE)`.
#' @return A one-row data.frame with `n_strata` (strata contributing at least
#'   one comparable pair), `n_pos`, `n_neg`, `n_zero`, `homogeneous` (no
#'   stratum positive while another is negative) and `w_opposing`, the share of
#'   the pooled denominator held by strata whose sign opposes the pooled one.
#'   `w_opposing` is `NA_real_` when the pooled gamma is zero or undefined.
#' @keywords internal
#' @noRd
.partgam_strata_summary <- function(st) {
  ok <- !is.na(st$gamma_k)
  g <- st$gamma_k[ok]
  w <- st$weight_k[ok]

  n_pos <- sum(g > 0)
  n_neg <- sum(g < 0)

  pooled_sign <- if (is.na(st$gamma)) NA_real_ else sign(st$gamma)
  w_opposing <- if (is.na(pooled_sign) || pooled_sign == 0 || sum(w) == 0) {
    NA_real_
  } else {
    sum(w[sign(g) == -pooled_sign]) / sum(w)
  }

  data.frame(
    n_strata = length(g),
    n_pos = n_pos,
    n_neg = n_neg,
    n_zero = sum(g == 0),
    homogeneous = !(n_pos > 0 && n_neg > 0),
    w_opposing = w_opposing,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# ---------------------------------------------------------------------------
# Internal: single simulation iteration
# ---------------------------------------------------------------------------

#' Run a single partial gamma LD simulation iteration
#'
#' @param seed Integer seed for reproducibility.
#' @param data_list List produced inside [RMlocdepGammaCutoff()].
#' @return A data.frame with columns `Item1`, `Item2`, and `gamma`, or a
#'   character string on failure.
#' @keywords internal
run_single_partgam_LD_sim <- function(seed, data_list) {
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

  thetas_res <- sample(
    data_list$thetas,
    size = data_list$sample_n,
    replace = TRUE
  )

  tryCatch(
    {
      if (data_list$type == "dichotomous") {
        sim_mat <- psychotools::rrm(
          theta = thetas_res,
          beta = data_list$item_params
        )
        sim_df <- as.data.frame(sim_mat$data)
        colnames(sim_df) <- data_list$item_names

        pos_counts <- colSums(sim_df, na.rm = TRUE)
        if (any(pos_counts < 8L)) {
          return(
            "validation_failed: fewer than 8 positive responses in at least one item"
          )
        }
      } else {
        sim_mat <- sim_partial_score(data_list$deltaslist, thetas_res)
        sim_df <- as.data.frame(sim_mat)
        colnames(sim_df) <- data_list$item_names

        n_cats <- vapply(
          data_list$deltaslist,
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

      # The pair statistic is the larger of the two conditioning directions.
      # Local dependence violates both of the conditional independence
      # hypotheses that the Rasch model implies for a pair (Kreiner &
      # Christensen, 2004), so a pair is tested once, on whichever direction
      # shows the stronger association. Taking the maximum within the same
      # simulated dataset gives the null of that maximum directly, so nothing
      # has to be corrected for having looked at two directions.
      #
      # Computed with the vectorised internal rather than
      # `iarm::partgam_LD()`: the coefficients are identical (see
      # test-ld_partgam_gamma.R), but iarm costs around 300 ms per call against
      # a few milliseconds here, and it is called once per iteration.
      g1 <- .partgam_ld_gamma(sim_df, direction = 1L)
      g2 <- .partgam_ld_gamma(sim_df, direction = 2L)
      # direction 2 enumerates the same pairs with the items swapped, so it is
      # matched on the unordered pair rather than on row order.
      key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "___")
      g1$gamma <- pmax(g1$gamma,
                       g2$gamma[match(key(g1$Item1, g1$Item2),
                                      key(g2$Item1, g2$Item2))])
      g1
    },
    error = function(e) {
      as.character(conditionMessage(e))
    }
  )
}

# ---------------------------------------------------------------------------
# Internal: parallel runner
# ---------------------------------------------------------------------------

#' Run partial gamma LD simulations in parallel using mirai
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param n_cores Number of mirai daemons.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_partgam_LD_sim_parallel <- function(
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
        run_single_partgam_LD_sim(seed, data_list)
      },
      seed = sim_seeds[sim],
      data_list = sim_data_list,
      run_single_partgam_LD_sim = run_single_partgam_LD_sim,
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

#' Run partial gamma LD simulations sequentially
#'
#' @param iterations Number of iterations.
#' @param sim_seeds Integer vector of per-iteration seeds.
#' @param sim_data_list List of data passed to each worker.
#' @param verbose Show progress bar.
#' @return List of raw results (one element per iteration).
#' @keywords internal
run_partgam_LD_sim_sequential <- function(
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
    results[[sim]] <- run_single_partgam_LD_sim(sim_seeds[sim], sim_data_list)
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

#' Plot Distribution of Simulated Partial Gamma LD Values
#'
#' Visualises the distribution of simulation-based partial gamma LD values
#' from \code{\link{RMlocdepGammaCutoff}}, optionally overlaying observed partial
#' gamma values computed from real data via \code{\link[iarm]{partgam_LD}}.
#'
#' Uses `ggdist::stat_dotsinterval()` (when `data` is not supplied) or
#' `ggdist::stat_dots()` (when `data` is supplied) with
#' `point_interval = "median_hdci"` and `.width = c(0.66, 0.95, 0.99)`.
#'
#' @param simfit The return value of \code{\link{RMlocdepGammaCutoff}} (a list with
#'   components `results`, `pair_cutoffs`, `actual_iterations`, `sample_n`, and
#'   `item_names`).
#' @param data Optional. A data.frame or matrix of item responses for computing
#'   and overlaying observed partial gamma values. Items must be scored starting
#'   at 0 (non-negative integers). When provided, the plot includes orange
#'   diamond markers for the observed partial gamma alongside the simulated
#'   distribution, plus segment summaries from the cutoff intervals.
#' @param items Optional character vector of item names to include in the plot.
#'   Only item pairs where **both** items are in this vector will be shown. When
#'   `NULL` (default), all item pairs are plotted.
#' @param n_pairs Optional positive integer. When supplied, only the
#'   `n_pairs` item pairs with the largest absolute partial gamma values
#'   are plotted, sorted by `|gamma|` descending. When `data` is
#'   supplied, the ranking uses the *observed* partial gammas (the
#'   diamonds you actually want to interpret); otherwise it falls back
#'   to the per-pair median of the simulated distributions. Applied
#'   *after* the `items` filter when both are supplied. Values larger
#'   than the number of available pairs are silently capped.
#'
#' @return A `ggplot` object.
#'
#' @details
#' The plot shows one row per item pair (labelled as "Item1 - Item2"). Only
#' direction 1 (rest score = total - Item2) is plotted, matching the
#' convention used in the simulation.
#'
#' When `data` is **not** supplied, the function plots the simulated partial
#' gamma distributions as dot-interval plots using
#' `ggdist::stat_dotsinterval()` with median and Highest Density Continuous
#' Interval (HDCI) summaries.
#'
#' When `data` **is** supplied, the function:
#' \enumerate{
#'   \item Computes observed partial gamma values via
#'     `iarm::partgam_LD()`.
#'   \item Overlays observed gamma values as orange diamond markers on the
#'     simulated distributions.
#'   \item Shows per-pair cutoff intervals (from `simfit$pair_cutoffs`) as
#'     black line segments, with thicker segments for the 66\% interval and
#'     black dots for the median.
#' }
#'
#' The `ggplot2`, `ggdist`, and optionally `iarm` packages must be installed
#' (they are in Suggests, not Imports).
#'
#' @seealso \code{\link{RMlocdepGammaCutoff}}, \code{\link{RMlocdepGamma}}
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("iarm", quietly = TRUE) &&
#'     requireNamespace("ggdist", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:10)
#'
#'   # Run simulation
#'   cutoff_res <- RMlocdepGammaCutoff(sim_data, iterations = 100,
#'                                     parallel = FALSE, seed = 42)
#'
#'   # Simulated distribution only
#'   RMlocdepGammaPlot(cutoff_res)
#'
#'   # With observed partial gamma overlaid
#'   RMlocdepGammaPlot(cutoff_res, data = sim_data)
#'
#'   # Plot only a subset of items
#'   RMlocdepGammaPlot(cutoff_res, data = sim_data,
#'                     items = c("Item1", "Item2", "Item3"))
#' }
#' }
RMlocdepGammaPlot <- function(simfit, data, items = NULL, n_pairs = NULL) {
  # --- Check required packages ------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for RMlocdepGammaPlot() but is not installed.\n",
      "Install it with: install.packages(\"ggplot2\")",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required for RMlocdepGammaPlot() but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")",
      call. = FALSE
    )
  }

  # --- Validate simfit --------------------------------------------------------
  required_names <- c(
    "results",
    "pair_cutoffs",
    "actual_iterations",
    "sample_n",
    "item_names"
  )
  missing_names <- setdiff(required_names, names(simfit))
  if (length(missing_names) > 0L) {
    stop(
      "`simfit` is missing required components: ",
      paste(missing_names, collapse = ", "),
      ".\nExpected the return value of RMlocdepGammaCutoff().",
      call. = FALSE
    )
  }

  results_df <- simfit$results
  pair_cutoffs <- simfit$pair_cutoffs
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

  # --- Validate items parameter -----------------------------------------------
  if (!is.null(items)) {
    unknown_items <- setdiff(items, item_names)
    if (length(unknown_items) > 0L) {
      stop(
        "Unknown item(s) in `items`: ",
        paste(unknown_items, collapse = ", "),
        ".\nAvailable items: ",
        paste(item_names, collapse = ", "),
        call. = FALSE
      )
    }
    if (length(items) < 2L) {
      stop(
        "`items` must contain at least 2 item names to form a pair.",
        call. = FALSE
      )
    }
  }

  # --- Validate n_pairs parameter ---------------------------------------------
  if (!is.null(n_pairs)) {
    if (
      !is.numeric(n_pairs) ||
        length(n_pairs) != 1L ||
        !is.finite(n_pairs) ||
        n_pairs < 1 ||
        n_pairs != as.integer(n_pairs)
    ) {
      stop(
        "`n_pairs` must be a single positive integer or NULL.",
        call. = FALSE
      )
    }
    n_pairs <- as.integer(n_pairs)
  }

  # Create pair labels
  results_df$Pair <- paste(results_df$Item1, "-", results_df$Item2)
  pair_cutoffs$Pair <- paste(pair_cutoffs$Item1, "-", pair_cutoffs$Item2)

  # --- Filter to selected items -----------------------------------------------
  if (!is.null(items)) {
    keep <- results_df$Item1 %in% items & results_df$Item2 %in% items
    results_df <- results_df[keep, , drop = FALSE]
    keep_cut <- pair_cutoffs$Item1 %in% items & pair_cutoffs$Item2 %in% items
    pair_cutoffs <- pair_cutoffs[keep_cut, , drop = FALSE]

    if (nrow(results_df) == 0L) {
      stop("No item pairs remain after filtering by `items`.", call. = FALSE)
    }
  }

  # --- Compute observed partial gammas up-front (when data supplied) ----------
  # We pre-compute here so the `n_pairs` filter below can rank by observed
  # |gamma| when data are available. The result is reused in Case 2.
  observed_df <- NULL
  if (!missing(data)) {
    if (!requireNamespace("iarm", quietly = TRUE)) {
      stop(
        "Package 'iarm' is required to compute observed partial gamma but is not installed.\n",
        "Install it with: install.packages(\"iarm\")",
        call. = FALSE
      )
    }
    validate_response_data(data)

    old_rgl <- getOption("rgl.useNULL")
    options(rgl.useNULL = TRUE)
    on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

    # The simulated distribution is the null of the pair statistic, the larger
    # of the two conditioning directions, so the observed overlay has to be the
    # same quantity rather than one direction's coefficient.
    dc <- data[stats::complete.cases(data), , drop = FALSE]
    o1 <- .partgam_ld_gamma(dc, direction = 1L)
    o2 <- .partgam_ld_gamma(dc, direction = 2L)
    pkey <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "___")

    observed_df <- data.frame(
      Item1 = o1$Item1,
      Item2 = o1$Item2,
      observed_gamma = pmax(
        o1$gamma,
        o2$gamma[match(pkey(o1$Item1, o1$Item2), pkey(o2$Item1, o2$Item2))]
      ),
      stringsAsFactors = FALSE
    )
    observed_df$Pair <- paste(observed_df$Item1, "-", observed_df$Item2)

    if (!is.null(items)) {
      keep_obs <- observed_df$Item1 %in% items & observed_df$Item2 %in% items
      observed_df <- observed_df[keep_obs, , drop = FALSE]
    }
  }

  # --- Apply n_pairs top-N filter ---------------------------------------------
  if (!is.null(n_pairs)) {
    if (!is.null(observed_df)) {
      # Rank by |observed gamma| — the diamonds the user wants to interpret
      ord <- order(abs(observed_df$observed_gamma), decreasing = TRUE)
      keep_n <- min(n_pairs, nrow(observed_df))
      keep_pairs <- observed_df$Pair[ord[seq_len(keep_n)]]
      observed_df <- observed_df[
        observed_df$Pair %in% keep_pairs,
        ,
        drop = FALSE
      ]
    } else {
      # Rank by |median simulated gamma| per pair
      pair_names_all <- unique(results_df$Pair)
      med_g <- vapply(
        pair_names_all,
        function(pp) {
          stats::median(results_df$gamma[results_df$Pair == pp], na.rm = TRUE)
        },
        numeric(1L)
      )
      ord <- order(abs(med_g), decreasing = TRUE)
      keep_n <- min(n_pairs, length(pair_names_all))
      keep_pairs <- pair_names_all[ord[seq_len(keep_n)]]
    }
    results_df <- results_df[results_df$Pair %in% keep_pairs, , drop = FALSE]
    pair_cutoffs <- pair_cutoffs[
      pair_cutoffs$Pair %in% keep_pairs,
      ,
      drop = FALSE
    ]
    # Use the ranked order for y-axis: largest |gamma| at the top
    pair_levels <- rev(keep_pairs)
  } else {
    pair_levels <- rev(unique(results_df$Pair))
  }

  # --- Compute per-pair summary intervals for segment overlays ----------------
  pair_names <- unique(results_df$Pair)
  lo_hi <- do.call(
    rbind,
    lapply(pair_names, function(pair) {
      sub <- results_df[results_df$Pair == pair, ]
      data.frame(
        Pair = pair,
        min_gamma = stats::quantile(sub$gamma, 0.005, na.rm = TRUE),
        max_gamma = stats::quantile(sub$gamma, 0.995, na.rm = TRUE),
        p66lo_gamma = stats::quantile(sub$gamma, 0.167, na.rm = TRUE),
        p66hi_gamma = stats::quantile(sub$gamma, 0.833, na.rm = TRUE),
        median_gamma = stats::median(sub$gamma, na.rm = TRUE),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    })
  )
  rownames(lo_hi) <- NULL

  # --- Case 1: no observed data, show simulation distribution only ------------
  if (missing(data)) {
    results_plot <- data.frame(
      Pair = results_df$Pair,
      Value = results_df$gamma,
      stringsAsFactors = FALSE
    )
    results_plot$Pair <- factor(results_plot$Pair, levels = pair_levels)

    p <- ggplot2::ggplot(
      results_plot,
      ggplot2::aes(
        x = .data$Value,
        y = .data$Pair
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
        y = "Item pair",
        caption = er2_caption(paste0(
          "Results from ",
          actual_iterations,
          " simulated datasets (no true local dependence). ",
          sample_clause,
          " per dataset."
        ))
      ) +
      ggplot2::scale_color_manual(
        values = scales::brewer_pal()(3),
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
  # (`observed_df` was pre-computed above so the n_pairs filter could rank
  # by |observed gamma|.)

  # --- Build plot data --------------------------------------------------------
  gamma_sim <- data.frame(
    Pair = results_df$Pair,
    Value = results_df$gamma,
    stringsAsFactors = FALSE
  )
  gamma_sim <- merge(
    gamma_sim,
    observed_df[, c("Pair", "observed_gamma")],
    by = "Pair",
    sort = FALSE
  )
  gamma_sim$Pair <- factor(gamma_sim$Pair, levels = pair_levels)

  lo_hi$Pair_f <- factor(lo_hi$Pair, levels = pair_levels)

  caption_text <- er2_caption(paste0(
    "Results from ",
    actual_iterations,
    " simulated datasets. ",
    sample_clause,
    " per dataset.\n",
    "Orange diamonds indicate observed partial gamma LD. ",
    "Black dots indicate median gamma from simulations."
  ))

  p <- ggplot2::ggplot(
    gamma_sim,
    ggplot2::aes(
      x = .data$Value,
      y = .data$Pair
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
        x = .data$min_gamma,
        xend = .data$max_gamma,
        y = .data$Pair_f,
        yend = .data$Pair_f
      ),
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x = .data$p66lo_gamma,
        xend = .data$p66hi_gamma,
        y = .data$Pair_f,
        yend = .data$Pair_f
      ),
      color = "black",
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = lo_hi,
      ggplot2::aes(
        x = .data$median_gamma,
        y = .data$Pair_f
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
    ggplot2::labs(
      x = "Partial gamma",
      y = "Item pair",
      caption = caption_text
    ) +
    ggplot2::scale_color_manual(
      values = scales::brewer_pal()(3),
      aesthetics = "slab_fill",
      guide = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
    er2_axis_margins() +
    er2_plot_caption()

  p
}
