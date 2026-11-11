#' DIF analysis via Andersen's likelihood-ratio test
#'
#' Fits a Rasch model by conditional maximum likelihood in the full sample
#' and in each group of an external grouping variable, tests the difference
#' with Andersen's likelihood-ratio test, and reports per-group item
#' locations (or per-group threshold locations) together with their standard
#' errors. \code{level} (item or threshold) and \code{output} (data.frame,
#' kable or ggplot) give six views of the same fits.
#'
#' The Partial Credit Model (PCM) is fitted by default for polytomous data
#' and the dichotomous Rasch Model (RM) is fitted when all responses are
#' 0/1; this can be overridden via \code{model}.
#'
#' @param data A data.frame or matrix of item responses (non-negative
#'   integers, 0-based). One column per item, one row per person. Person
#'   IDs and grouping variables must not be included -- pass the grouping
#'   variable separately via \code{dif_var}.
#' @param dif_var Vector of length \code{nrow(data)} (factor, character,
#'   or numeric) defining the DIF grouping variable. Coerced to factor;
#'   unused levels are dropped. Rows where \code{dif_var} is \code{NA} are
#'   dropped with a message. Must result in at least 2 groups after
#'   cleaning.
#' @param model One of \code{"auto"} (default), \code{"PCM"}, or
#'   \code{"RM"}. \code{"auto"} fits the Rasch model when the data are
#'   dichotomous (max response = 1) and the partial credit model otherwise
#'   (the two coincide for dichotomous data). \code{"RM"} errors on
#'   polytomous data.
#' @param level One of \code{"item"} (default) or \code{"threshold"}.
#'   \code{"item"} reports each item's mean threshold location per group;
#'   \code{"threshold"} reports each individual threshold per group. For
#'   dichotomous (RM) data the two views are equivalent (one threshold per
#'   item).
#' @param output One of \code{"kable"} (default), \code{"dataframe"}, or
#'   \code{"ggplot"}. The data.frame has a logical \code{Flagged} column
#'   (unless \code{cutoff = NULL}), the kable bolds flagged rows, and the
#'   ggplot shows confidence intervals.
#' @param cutoff Numeric or \code{NULL}. Threshold (in logits) for the
#'   \code{Flagged} column: a row is flagged when
#'   \code{MaxDiff > cutoff}, where \code{MaxDiff} is the difference
#'   between the largest and smallest per-group location for that item
#'   (or threshold). Set to \code{NULL} to suppress flagging. Default
#'   \code{0.5}. This default is a convention in applied Rasch work (Jin et
#'   al., 2013; Kim et al., 2017; Pickard et al., 2006), not a value with a
#'   statistical or substantive justification. See Details.
#' @param conf Numeric in (0, 1). Confidence level used for the ggplot
#'   error bars. Default \code{0.95}.
#' @param sort Logical. Sort the rows of the data.frame and kable output by
#'   \code{MaxDiff} (descending). Default \code{FALSE}.
#'
#' @return A data.frame, a \code{knitr_kable} object, or a \code{ggplot}
#'   object, depending on \code{output}.
#'
#'   The data.frame has one row per item (\code{level = "item"}) or per
#'   item x threshold (\code{level = "threshold"}), with columns
#'   \code{Item} (and \code{Threshold} at threshold level), one numeric
#'   column per group level, an \code{All} column for the unsplit fit,
#'   \code{MaxDiff}, \code{Flagged} (when \code{cutoff} is non-\code{NULL}),
#'   and matching \code{SE_*} columns.
#'
#'   The Andersen LR test result is attached as
#'   \code{attr(result, "lr_test")} on the data.frame, in the kable
#'   footnote, and in the ggplot caption (LR \eqn{\chi^2}, df, p-value).
#'
#' @details
#' Each fit is a CML partial credit model from \code{psychotools::pcmodel()}
#' (a dichotomous item is a 2-category PCM). Within each fit the thresholds
#' sum to zero, as in \code{\link{RMitemParameters}}, so each group's
#' locations are centred on that group's own mean threshold. The \code{All}
#' column is the fit to the full dataset. An item's location is the mean of
#' its thresholds, and its SE is the SE of that mean, computed from the
#' thresholds' covariance matrix.
#'
#' The Andersen LR statistic is twice the difference between the sum of the
#' group conditional log-likelihoods and the full-sample one, with df equal
#' to the extra free parameters of the group fits. It equals
#' \code{eRm::LRtest()} on the same data. Respondents with missing responses
#' are kept (CML uses the items they answered).
#'
#' \strong{Response categories.} Every category an item uses must be used in
#' every group, or that item's thresholds cannot be compared across groups.
#' The same holds when a category is used in a group only by respondents
#' with a zero or perfect score there, since CML leaves them out and the
#' threshold cannot be estimated. The function then stops and names the
#' items, categories and groups
#' (\code{eRm::LRtest()} instead left items with an unused category out of
#' the test without saying so, and returned diverging estimates in the
#' second case). With small groups this is common for sparse categories: merge
#' them with an adjacent one, combine groups, or leave the item out.
#'
#' \strong{DIF size and the cutoff.} The \code{Flagged} column applies one
#' fixed cutoff to \code{MaxDiff}, and no single value suits every scale. An
#' appropriate cutoff depends on the number of items, the number of response
#' categories and the type of DIF (Van de Winckel et al., 2022). Treat a flag
#' as a prompt for a closer look, not as a verdict. Judging DIF means
#' weighing its size against its practical consequences (Henninger et al.,
#' 2025). These are best evaluated by simulation: estimate person locations
#' with and without accounting for DIF, for instance by splitting the
#' affected item by group, and check whether the difference changes
#' substantive conclusions (Johansson et al., 2026). For polytomous items,
#' DIF can differ between an item's thresholds, which \code{level =
#' "threshold"} shows (Henninger et al., 2025). \code{MaxDiff} is a point
#' estimate, and with small groups its uncertainty is large (see the
#' \code{SE_*} columns). The LR test addresses a different question: whether
#' any DIF is present, not how large it is.
#'
#' \strong{Real and artificial DIF.} Real DIF in one item induces artificial
#' DIF in the other items, in the opposite direction (Hagquist & Andrich,
#' 2017). Here the mechanism is direct: each group's locations are centred
#' on that group's own mean threshold, so a real shift in one item moves the
#' other items' locations the other way. An item can therefore be flagged
#' only because another item has DIF, and one analysis cannot tell real
#' from artificial DIF. Hagquist and Andrich (2017) recommend resolving DIF
#' items one at a time, starting with the item with the largest DIF (for
#' instance by splitting it into group-specific items), and repeating the
#' analysis before deciding on the next. Do not overinterpret the flags from
#' a single run.
#'
#' @references
#' Andersen, E. B. (1973). A goodness of fit test for the Rasch model.
#' *Psychometrika, 38*(1), 123-140. \doi{10.1007/BF02291180}
#'
#' Hagquist, C., & Andrich, D. (2017). Recent advances in analysis of
#' differential item functioning in health research using the Rasch model.
#' *Health and Quality of Life Outcomes, 15*, 181.
#' \doi{10.1186/s12955-017-0755-0}
#'
#' Henninger, M., Radek, J., Sengewald, M.-A., & Strobl, C. (2025). Partial
#' credit trees meet the partial gamma coefficient for quantifying DIF and
#' DSF in polytomous items. *Behaviormetrika, 52*(2), 221-257.
#' \doi{10.1007/s41237-024-00252-3}
#'
#' Jin, Y., Myers, N. D., Ahn, S., & Penfield, R. D. (2013). A comparison of
#' uniform DIF effect size estimators under the MIMIC and Rasch models.
#' *Educational and Psychological Measurement, 73*(2), 339-358.
#' \doi{10.1177/0013164412462705}
#'
#' Johansson, M., Preuter, M., & Melin, J. (2026). Valid and reliable? Basic
#' criteria and modern methods for assessing psychometric quality. OSF
#' Preprints. \doi{10.31219/osf.io/3htzc_v2}
#'
#' Kim, S., Ouellet, L. J., Mazza, J., & Spaulding, A. C. (2017). Rasch
#' analysis and differential item functioning of a social support measure
#' in jail inmates with HIV infection. *Evaluation & the Health Professions,
#' 40*(1), 33-60. \doi{10.1177/0163278716644954}
#'
#' Pickard, A. S., Dalal, M. R., & Bushnell, D. M. (2006). A comparison of
#' depressive symptoms in stroke and primary care: Applying Rasch models to
#' evaluate the Center for Epidemiologic Studies-Depression Scale. *Value in
#' Health, 9*(1), 59-64. \doi{10.1111/j.1524-4733.2006.00082.x}
#'
#' Strobl, C., Kopf, J., & Zeileis, A. (2015). Rasch trees: A new method for
#' detecting differential item functioning in the Rasch model.
#' *Psychometrika, 80*(2), 289-316. \doi{10.1007/s11336-013-9388-3}
#'
#' Van de Winckel, A., Kozlowski, A. J., Johnston, M. V., Weaver, J.,
#' Grampurohit, N., Terhorst, L., Juengst, S., Ehrlich-Jones, L., Heinemann,
#' A. W., Melvin, J., Sood, P., & Mallinson, T. (2022). Reporting guideline
#' for RULER: Rasch reporting guideline for rehabilitation research:
#' Explanation and elaboration. *Archives of Physical Medicine and
#' Rehabilitation, 103*(7), 1487-1498. \doi{10.1016/j.apmr.2022.03.019}
#'
#' @seealso \code{\link{RMdifGamma}}, with simulation-based cutoffs from
#'   \code{\link{RMdifGammaCutoff}}, and \code{\link{RMdifTree}}, for
#'   continuous DIF variables and interactions between them (Strobl et al.,
#'   2015; Henninger et al., 2025).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("eRm", quietly = TRUE)) { # for the example data
#'   set.seed(1)
#'   data("pcmdat2", package = "eRm")
#'   grp <- factor(sample(c("A", "B"), nrow(pcmdat2), replace = TRUE))
#'
#'   # Default: kable of per-group item locations
#'   RMdifLR(pcmdat2, dif_var = grp)
#'
#'   # ggplot panel of item locations with 95% CIs
#'   RMdifLR(pcmdat2, dif_var = grp, output = "ggplot")
#'
#'   # Threshold-level kable, sorted by MaxDiff
#'   RMdifLR(pcmdat2, dif_var = grp, level = "threshold", sort = TRUE)
#'
#'   # Tidy data.frame for downstream use
#'   df <- RMdifLR(pcmdat2, dif_var = grp, output = "dataframe")
#'   attr(df, "lr_test")
#'   df[df$Flagged, ]
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
RMdifLR <- function(
  data,
  dif_var,
  model = c("auto", "PCM", "RM"),
  level = c("item", "threshold"),
  output = c("kable", "dataframe", "ggplot"),
  cutoff = 0.5,
  conf = 0.95,
  sort = FALSE
) {
  model <- match.arg(model)
  level <- match.arg(level)
  output <- match.arg(output)

  if (output == "ggplot" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for output = \"ggplot\".",
      call. = FALSE
    )
  }
  if (output == "kable" && !requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required for output = \"kable\".", call. = FALSE)
  }
  if (!is.null(cutoff)) {
    if (
      !is.numeric(cutoff) ||
        length(cutoff) != 1L ||
        !is.finite(cutoff) ||
        cutoff < 0
    ) {
      stop(
        "`cutoff` must be a single non-negative numeric value, or NULL.",
        call. = FALSE
      )
    }
  }
  if (
    !is.numeric(conf) ||
      length(conf) != 1L ||
      !is.finite(conf) ||
      conf <= 0 ||
      conf >= 1
  ) {
    stop("`conf` must be a single numeric in (0, 1).", call. = FALSE)
  }

  validate_response_data(data)
  data <- as.data.frame(data)
  # (validate_response_data() also checks for 2+ items and unique names)
  n_total_lr <- nrow(data)

  # ---------------------------------------------------------------------
  # Resolve dif_var and drop NAs (jointly with `data`)
  # ---------------------------------------------------------------------
  if (length(dif_var) != nrow(data)) {
    stop(
      "`dif_var` must have the same length as nrow(data) (",
      nrow(data),
      "). Got ",
      length(dif_var),
      ".",
      call. = FALSE
    )
  }
  dif_factor <- droplevels(as.factor(dif_var))

  na_mask <- is.na(dif_factor)
  if (any(na_mask)) {
    message(sum(na_mask), " row(s) with NA in `dif_var` dropped.")
    data <- data[!na_mask, , drop = FALSE]
    dif_factor <- droplevels(dif_factor[!na_mask])
  }

  # Respondents with no responses at all contribute nothing to CML; drop them jointly with `dif_factor` (the standard
  # .drop_empty_respondents() treatment, inlined to keep the rows aligned
  # with the grouping variable).
  empty <- rowSums(!is.na(as.matrix(data))) == 0L
  if (any(empty)) {
    message(sum(empty), " respondent(s) with no responses dropped.")
    data <- data[!empty, , drop = FALSE]
    dif_factor <- droplevels(dif_factor[!empty])
  }

  # After the drops above: NA in `dif_var`, or in all-NA rows, does not
  # mean that incomplete responses are retained.
  has_na_lr <- anyNA(data)

  if (nlevels(dif_factor) < 2L) {
    stop(
      "`dif_var` must have at least 2 distinct non-missing levels.",
      call. = FALSE
    )
  }
  groups <- levels(dif_factor)
  nr_groups <- length(groups)
  # Group labels become column names next to Item, All, MaxDiff and SE_*.
  clash <- groups[
    groups %in% c("Item", "Threshold", "All", "MaxDiff", "Flagged") |
      startsWith(groups, "SE_") | !nzchar(groups)
  ]
  if (length(clash) > 0L) {
    stop(
      "`dif_var` level(s) ", paste0("'", clash, "'", collapse = ", "),
      " would clash with the output columns (Item, Threshold, All, MaxDiff, ",
      "Flagged, SE_*) or are empty. Rename them, e.g. with levels(dif_var).",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------
  # Pick the model
  # ---------------------------------------------------------------------
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "integer"
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  resolved_model <- if (model == "auto") {
    if (is_polytomous) "PCM" else "RM"
  } else {
    if (model == "RM" && is_polytomous) {
      stop(
        "`model = \"RM\"` requires dichotomous (0/1) data; use ",
        "\"PCM\" or \"auto\" for polytomous responses.",
        call. = FALSE
      )
    }
    model
  }

  # ---------------------------------------------------------------------
  # Fit the full sample and each group by CML, then Andersen's LR test
  # ---------------------------------------------------------------------
  # Every category an item uses in the full data must occur in every group,
  # or the group fits would have different thresholds (psychotools drops an
  # unused category's threshold) and the locations could not be compared.
  .check_response_categories(data) # full data first, the plain message
  # A category only extreme scorers use in the full data is unusable in
  # every group: stop with the full-data message, not the per-group one.
  ec_full <- .extreme_only_categories(data)
  if (nrow(ec_full) > 0L) {
    stop(.category_error(
      paste0(
        .extreme_only_message(ec_full),
        " RMdifLR() stops here, since these thresholds cannot be ",
        "compared across groups."
      ),
      unique(ec_full$Item)
    ))
  }
  .check_group_categories(data, dif_factor)
  item_names <- names(data)
  is_rm <- resolved_model == "RM"
  # "auto" fits raschmodel() for dichotomous data (same estimates as a
  # 2-category pcmodel(), about twice as fast), pcmodel() otherwise.
  fit_lr <- function(d, where) {
    tryCatch(
      # pcmodel() for every fit (a dichotomous item is a 2-category PCM; the
      # raschmodel() Hessian fails on some incomplete data). check = FALSE:
      # the full data and every group were checked above, so any error here
      # comes from the fit itself.
      .cml_fit(d, check = FALSE),
      error = function(e) {
        stop(
          "CML estimation failed in ", where, " (", conditionMessage(e),
          "). Small or sparse groups are the usual cause; inspect the ",
          "counts with RMplotTile(data, group = dif_var).",
          call. = FALSE
        )
      }
    )
  }
  full_max <- apply(as.matrix(data), 2L, max, na.rm = TRUE)
  fit_full <- fit_lr(data, "the full sample")
  fits <- lapply(groups, function(g) {
    fit_lr(data[dif_factor == g, , drop = FALSE], paste0("group '", g, "'"))
  })

  # LR = 2 * (sum of group conditional log-likelihoods - full one), with df
  # the extra free parameters of the group fits (Andersen, 1973).
  ll <- function(f) as.numeric(stats::logLik(f))
  npar <- function(f) length(stats::coef(f))
  lr_stat <- 2 * (sum(vapply(fits, ll, numeric(1L))) - ll(fit_full))
  lr_df <- sum(vapply(fits, npar, integer(1L))) - npar(fit_full)
  lr_p <- stats::pchisq(lr_stat, lr_df, lower.tail = FALSE)

  parts <- c(
    Map(
      function(f, g) {
        lr_fit_summary(f, data[dif_factor == g, , drop = FALSE], g, is_rm,
                       full_max)
      },
      fits, groups
    ),
    list(lr_fit_summary(fit_full, data, "All", is_rm, full_max))
  )
  per_level <- if (level == "item") "item" else "thr"
  long_df <- do.call(rbind, lapply(parts, `[[`, per_level))
  rownames(long_df) <- NULL
  long_df$Item <- factor(long_df$Item, levels = item_names)

  # ---------------------------------------------------------------------
  # Wide table at the requested level
  # ---------------------------------------------------------------------
  if (level == "item") {
    table_df <- build_wide_table(long_df, groups, level = "item")
    plot_df <- long_df
    plot_df$Threshold <- NA_character_
  } else {
    table_df <- build_wide_table(long_df, groups, level = "threshold")
    plot_df <- long_df
  }

  # ---------------------------------------------------------------------
  # Flagged column
  # ---------------------------------------------------------------------
  if (!is.null(cutoff)) {
    table_df$Flagged <- !is.na(table_df$MaxDiff) &
      table_df$MaxDiff > cutoff
  }

  # ---------------------------------------------------------------------
  # Reorder columns: Item [Threshold] groups All MaxDiff [Flagged] SE_*
  # ---------------------------------------------------------------------
  id_cols <- if (level == "item") "Item" else c("Item", "Threshold")
  loc_cols <- c(
    groups,
    "All",
    "MaxDiff",
    if (!is.null(cutoff)) "Flagged" else NULL
  )
  se_cols <- paste0("SE_", c(groups, "All"))
  table_df <- table_df[, c(id_cols, loc_cols, se_cols), drop = FALSE]
  rownames(table_df) <- NULL

  # ---------------------------------------------------------------------
  # LR test summary (Andersen)
  # ---------------------------------------------------------------------
  lr_summary <- list(
    LR = lr_stat,
    df = lr_df,
    p_value = lr_p,
    n_groups = nr_groups,
    groups = groups,
    model = resolved_model,
    n_persons = nrow(data),
    n_items = ncol(data)
  )
  attr(table_df, "lr_test") <- lr_summary

  # ---------------------------------------------------------------------
  # Dispatch on output
  # ---------------------------------------------------------------------
  if (output == "dataframe") {
    if (sort) {
      table_df <- table_df[order(-table_df$MaxDiff), , drop = FALSE]
      rownames(table_df) <- NULL
    }
    return(table_df)
  }

  if (output == "kable") {
    # Display rounding (3 dp; the dataframe output above stays unrounded)
    # (through .round_display(), so test-rounding.R sees these columns)
    num_cols <- c(groups, "All", "MaxDiff", se_cols)
    table_df <- .round_display(
      table_df,
      stats::setNames(as.list(rep(3L, length(num_cols))), num_cols)
    )
    n_clause <- .n_caption(
      lr_summary$n_persons,
      n_total_lr,
      if (has_na_lr) "incomplete responses retained" else character()
    )
    return(render_lr_kable(
      table_df,
      lr_summary,
      level = level,
      cutoff = cutoff,
      sort = sort,
      is_polytomous = is_polytomous,
      n_clause = n_clause
    ))
  }

  # output == "ggplot"
  render_lr_ggplot(
    plot_df,
    lr_summary,
    level = level,
    groups = groups,
    conf = conf,
    is_polytomous = is_polytomous
  )
}

# =====================================================================
# Internal helpers
# =====================================================================

#' Stop unless every category an item uses occurs in every DIF group
#'
#' psychotools fits an item with an unused category in one group with fewer
#' thresholds there (or, for an empty lowest or interior category, shifts or
#' drops them), so its locations could not be compared across groups. The
#' full data have already passed `.check_response_categories()`, so an item's
#' categories there run 0 to its highest observed one. A category used in a
#' group only by respondents with a zero or perfect score there also stops:
#' CML leaves those respondents out, so the category is unused in that
#' group's fit and its threshold diverges.
#'
#' @param data Response data.frame (items from 0).
#' @param dif_factor Factor of the same length, without `NA`.
#' @return Invisibly `TRUE`; stops with an `easyRasch2_category_error`
#'   naming the items, the categories and the groups.
#' @keywords internal
#' @noRd
.check_group_categories <- function(data, dif_factor) {
  data <- as.matrix(data)
  items <- colnames(data)
  problems <- character(0)
  bad <- character(0)
  unused <- matrix(FALSE, ncol(data), nlevels(dif_factor))
  full_max <- apply(data, 2L, max, na.rm = TRUE)
  for (j in seq_len(ncol(data))) {
    x <- data[, j]
    cats <- sort(unique(x[!is.na(x)]))
    for (gi in seq_along(levels(dif_factor))) {
      g <- levels(dif_factor)[gi]
      xg <- x[dif_factor == g]
      miss <- setdiff(cats, xg[!is.na(xg)])
      if (length(miss) > 0L) {
        unused[j, gi] <- TRUE
        bad <- union(bad, items[j])
        problems <- c(problems, paste0(
          items[j], " (",
          if (length(miss) == length(cats)) {
            "no responses"
          } else {
            paste0(
              if (length(miss) == 1L) "category " else "categories ",
              paste(miss, collapse = ", "), " unused"
            )
          },
          " in group '", g, "')"
        ))
      }
    }
  }
  # Extreme-score-only categories per group, for items whose categories are
  # all used there (else the group's item maxima would differ).
  for (gi in seq_along(levels(dif_factor))) {
    g <- levels(dif_factor)[gi]
    # Perfect is judged against the full-data maxima, so an item whose top
    # category is unused in the group does not make others look extreme-only.
    ec <- .extreme_only_categories(
      data[dif_factor == g, , drop = FALSE],
      item_max = full_max
    )
    ec <- ec[!unused[match(ec$Item, items), gi], , drop = FALSE]
    for (it in unique(ec$Item)) {
      k <- ec$Category[ec$Item == it]
      bad <- union(bad, it)
      problems <- c(problems, paste0(
        it, " (", if (length(k) == 1L) "category " else "categories ",
        paste(k, collapse = ", "), " used only by respondents with a zero ",
        "or perfect score in group '", g, "')"
      ))
    }
  }
  if (length(problems) > 0L) {
    msg <- paste0(
      "Item locations cannot be compared across DIF groups for: ",
      paste(problems, collapse = ", "),
      ". Every response category of an item must be used in every group, ",
      "by respondents without a zero or perfect score (conditional maximum ",
      "likelihood leaves those out). ",
      "Merge a sparse category with an adjacent one, combine groups, or ",
      "leave the item out. RMplotTile(data, group = dif_var) shows the ",
      "counts."
    )
    stop(.category_error(msg, bad))
  }
  invisible(TRUE)
}

#' Threshold and item locations with SEs from one CML fit
#'
#' Thresholds come from `psychotools::threshpar()`, which makes them sum to
#' zero within the fit (the identification `RMitemParameters()` uses), with
#' their covariance matrix `V`. An item's location is the mean of its
#' thresholds, and its SE is `sqrt(w' V w)` over that item's block, with
#' `w = 1 / m`, so the covariance between an item's thresholds is included.
#'
#' @param fit A `pcmodel` or `raschmodel` fit (hessian on).
#' @param data The data the fit was made from (the full data or one group).
#' @param full_max The full-data top category per item. The alignment guard
#'   checks every fit against it, independently of the group check.
#' @param dif_group Label for the `DIFgroup` column (`"All"` for the full fit).
#' @param is_rm `TRUE` labels the single threshold `"1"` (dichotomous).
#' @return A list with `thr` (Item, Threshold, DIFgroup, Location, SE) and
#'   `item` (Item, DIFgroup, Location, SE) data.frames.
#' @keywords internal
#' @noRd
lr_fit_summary <- function(fit, data, dif_group, is_rm, full_max) {
  where <- if (dif_group == "All") "the full sample" else paste0("group '", dif_group, "'")
  # Shared extraction: alignment guard, and NA (with a warning) for
  # threshold SEs that cannot be computed.
  ts <- .cml_threshold_summary(fit, data, se = TRUE, label = where,
                               max_cat = full_max)
  item_names <- names(ts$thr_list)
  thr <- ts$thr_list
  V <- ts$vcov_thr
  m <- lengths(thr)
  blocks <- split(seq_len(sum(m)), rep.int(seq_along(m), m))
  thr_df <- data.frame(
    Item = rep(item_names, m),
    Threshold = if (is_rm) "1" else unlist(lapply(m, function(k) paste0("c", seq_len(k)))),
    DIFgroup = dif_group,
    Location = unlist(thr, use.names = FALSE),
    SE = unlist(ts$se_list, use.names = FALSE),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  item_se <- vapply(blocks, function(i) {
    w <- rep(1 / length(i), length(i))
    suppressWarnings(sqrt(sum(w * (V[i, i, drop = FALSE] %*% w))))
  }, numeric(1L))
  # A threshold SE problem has already warned; an item SE can still fail on
  # its own (a negative quadratic form).
  if (any(!is.finite(item_se)) && all(is.finite(thr_df$SE))) {
    warning(
      "Standard errors could not be computed for some item locations in ",
      where, "; they are NA.",
      call. = FALSE
    )
  }
  item_se[!is.finite(item_se)] <- NA_real_
  item_df <- data.frame(
    Item = item_names,
    DIFgroup = dif_group,
    Location = vapply(thr, mean, numeric(1L)),
    SE = unname(item_se),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  list(thr = thr_df, item = item_df)
}

#' Reshape long -> wide on DIFgroup, add MaxDiff, attach SE columns
#'
#' @keywords internal
#' @noRd
build_wide_table <- function(long_df, groups, level) {
  id_cols <- if (level == "item") "Item" else c("Item", "Threshold")

  loc_wide <- stats::reshape(
    long_df[, c(id_cols, "DIFgroup", "Location"), drop = FALSE],
    idvar = id_cols,
    timevar = "DIFgroup",
    direction = "wide"
  )
  names(loc_wide) <- sub("^Location\\.", "", names(loc_wide))

  se_wide <- stats::reshape(
    long_df[, c(id_cols, "DIFgroup", "SE"), drop = FALSE],
    idvar = id_cols,
    timevar = "DIFgroup",
    direction = "wide"
  )
  names(se_wide) <- sub("^SE\\.", "SE_", names(se_wide))

  out <- merge(loc_wide, se_wide, by = id_cols, sort = FALSE)

  # MaxDiff across the per-group columns (excludes "All")
  group_mat <- as.matrix(out[, groups, drop = FALSE])
  out$MaxDiff <- apply(group_mat, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) NA_real_ else max(x) - min(x)
  })

  ord <- if (level == "threshold") {
    # Numeric order, so c10 comes after c9.
    order(out$Item, as.integer(sub("^c", "", out$Threshold)))
  } else {
    order(out$Item)
  }
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Render the wide table as a kable
#'
#' @keywords internal
#' @noRd
render_lr_kable <- function(
  table_df,
  lr_summary,
  level,
  cutoff,
  sort,
  is_polytomous,
  n_clause = NULL
) {
  if (sort) {
    table_df <- table_df[order(-table_df$MaxDiff), , drop = FALSE]
    rownames(table_df) <- NULL
  }

  display <- table_df

  # Bold flagged rows by wrapping each visible cell in **...**
  if ("Flagged" %in% names(display)) {
    flag <- isTRUE_vec(display$Flagged)
    display$Flagged <- ifelse(flag, "yes", "no")
    if (any(flag)) {
      data_cols <- setdiff(
        names(display),
        if (level == "item") "Item" else c("Item", "Threshold")
      )
      for (cc in data_cols) {
        col <- display[[cc]]
        if (is.numeric(col)) {
          col <- format(col, nsmall = 3, trim = TRUE)
        }
        display[[cc]] <- ifelse(
          flag,
          paste0("**", col, "**"),
          as.character(col)
        )
      }
    }
  }

  caption <- paste0(
    if (is_polytomous) "Partial Credit Model" else "Rasch Model",
    " split by DIF variable (",
    lr_summary$n_groups,
    " groups). ",
    "Andersen LR \u03c7\u00b2 = ",
    round(lr_summary$LR, 3),
    ", df = ",
    lr_summary$df,
    ", ",
    .format_p(lr_summary$p_value),
    ". ",
    n_clause,
    ", ",
    lr_summary$n_items,
    " items.",
    if (!is.null(cutoff)) {
      paste0(" **Bold** = MaxDiff > ", cutoff, " logits.")
    } else {
      ""
    }
  )

  knitr::kable(display, format = "pipe", row.names = FALSE, caption = caption)
}

#' Coerce maybe-character flag column back to logical for rendering
#'
#' @keywords internal
#' @noRd
isTRUE_vec <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  vapply(x, isTRUE, logical(1L))
}

#' Render a ggplot panel of per-group locations with confidence bars
#'
#' @keywords internal
#' @noRd
render_lr_ggplot <- function(
  plot_df,
  lr_summary,
  level,
  groups,
  conf,
  is_polytomous
) {
  z <- stats::qnorm(1 - (1 - conf) / 2)

  plot_df$DIFgroup <- factor(plot_df$DIFgroup, levels = c(groups, "All"))
  if (level == "threshold") {
    # Numeric order in the legend, so c10 comes after c9.
    lab <- unique(plot_df$Threshold)
    plot_df$Threshold <- factor(
      plot_df$Threshold,
      levels = lab[order(as.integer(sub("^c", "", lab)))]
    )
  }
  group_only <- plot_df[plot_df$DIFgroup %in% groups, , drop = FALSE]
  all_only <- plot_df[plot_df$DIFgroup == "All", , drop = FALSE]

  group_only$DIFgroup <- droplevels(group_only$DIFgroup)
  diamond_x <- (length(groups) + 1) / 2 + 0.15

  caption <- er2_caption(paste0(
    "Andersen LR \u03c7\u00b2 = ",
    round(lr_summary$LR, 3),
    ", df = ",
    lr_summary$df,
    ", ",
    .format_p(lr_summary$p_value),
    ". Error bars: ",
    round(conf * 100),
    "% CI. ",
    "Black diamond = whole-sample location."
  ))

  if (level == "item") {
    p <- ggplot2::ggplot(
      group_only,
      ggplot2::aes(
        x = .data$DIFgroup,
        y = .data$Location,
        group = .data$Item,
        colour = .data$Item
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data$Location - z * .data$SE,
          ymax = .data$Location + z * .data$SE
        ),
        width = 0.1
      ) +
      ggplot2::geom_point(
        data = all_only,
        ggplot2::aes(x = diamond_x, y = .data$Location),
        shape = 18,
        colour = "black",
        size = 3,
        alpha = 0.6,
        inherit.aes = FALSE
      ) +
      ggplot2::facet_wrap(~Item) +
      ggplot2::labs(
        title = "DIF: item locations by group",
        subtitle = "Item locations are means of threshold locations",
        x = "DIF group",
        y = "Item location (logits)",
        caption = caption
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        panel.spacing = ggplot2::unit(0.7, "cm")
      ) +
      er2_axis_margins() +
      er2_plot_caption()
    return(p)
  }

  # level == "threshold"
  p <- ggplot2::ggplot(
    group_only,
    ggplot2::aes(
      x = .data$DIFgroup,
      y = .data$Location,
      group = .data$Threshold,
      colour = .data$Threshold
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$Location - z * .data$SE,
        ymax = .data$Location + z * .data$SE
      ),
      width = 0.1
    ) +
    ggplot2::geom_errorbar(
      data = all_only,
      ggplot2::aes(
        x = diamond_x,
        ymin = .data$Location - z * .data$SE,
        ymax = .data$Location + z * .data$SE
      ),
      width = 0.1,
      colour = "darkgrey",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = all_only,
      ggplot2::aes(x = diamond_x, y = .data$Location),
      shape = 18,
      colour = "black",
      size = 3,
      alpha = 0.6,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~Item) +
    ggplot2::labs(
      title = "DIF: threshold locations by group",
      x = "DIF group",
      y = "Threshold location (logits)",
      caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.spacing = ggplot2::unit(0.7, "cm")
    ) +
    er2_axis_margins() +
    er2_plot_caption()

  p
}
