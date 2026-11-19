#' Anchored DIF tests for items and thresholds (two groups)
#'
#' Compares two groups' item and threshold locations from conditional
#' maximum likelihood (CML) partial credit fits. The two scales are linked
#' through anchor items chosen from the data, each comparison gets a Wald
#' test, and the p-values are corrected across items or thresholds. The
#' spacing test of thresholds needs no anchor.
#'
#' @param data A data.frame or matrix of item responses (non-negative
#'   integers, 0-based). One column per item, one row per person. Pass the
#'   grouping variable separately via \code{dif_var}.
#' @param dif_var Vector of length \code{nrow(data)} (factor, character or
#'   numeric) with exactly two groups. Rows where it is \code{NA} are dropped
#'   with a message.
#' @param anchor \code{"forward"} (default) or \code{"gini"}. See Details.
#' @param level What is tested: \code{"item"} (default, item locations),
#'   \code{"spacing"} (the distances between an item's thresholds, no anchor
#'   needed) or \code{"threshold"} (each threshold). Each level is its own
#'   family for the multiplicity correction.
#' @param output \code{"kable"} (default), \code{"dataframe"},
#'   \code{"ggplot"} or \code{"list"} (the table together with the
#'   settings and sample sizes, see Value).
#' @param plot_type With \code{output = "ggplot"}: \code{"differences"}
#'   (default) plots the differences at the chosen \code{level} with
#'   confidence intervals. \code{"criterion"} plots the Gini criterion of
#'   Strobl et al. (2021), which is available with either anchor.
#' @param correction Multiplicity correction: \code{"fwer"} (default,
#'   single-step, controls the family-wise error rate), \code{"holm"},
#'   \code{"fdr_bh"} (Benjamini-Hochberg, controls the false discovery rate)
#'   or \code{"none"}.
#' @param alpha A unit is flagged when its corrected p-value is below
#'   \code{alpha}. Default \code{0.05}.
#' @param conf Confidence level of the intervals. Default \code{0.95}. The
#'   intervals are not corrected for multiplicity.
#' @param reference The group subtracted in the differences, one of the two
#'   levels of \code{dif_var}. Default \code{NULL}, the first level.
#' @param sort Logical. Sort the table rows by the size of the difference
#'   (spacing: by the test statistic). Default \code{FALSE}.
#' @param n_seeds Number of random tie-break orderings for the forward
#'   anchor. Default \code{100}. Ignored for \code{anchor = "gini"}.
#' @param seed Integer seed for the tie-break orderings and the Monte-Carlo
#'   single-step adjustment, so that the same data give the same result.
#'   The generator is pinned to R's defaults while it runs, and the caller's
#'   random number stream and \code{RNGkind()} are left unchanged.
#'   \code{NULL} uses the current stream instead, and results then vary
#'   slightly between calls. Default \code{20260113}. See
#'   \link{easyRasch2-reproducibility}.
#'
#' @return A data.frame, a \code{knitr_kable}, a \code{ggplot} or a list,
#'   depending on \code{output}. Nothing is attached as attributes.
#'
#'   The data.frame for \code{level = "item"} has one row per item with
#'   \code{Item}, \code{DIF} (location difference, the other group minus
#'   \code{reference}), \code{SE}, \code{lower} and \code{upper} (confidence
#'   interval), \code{z}, \code{p}, \code{padj}, \code{Flagged}
#'   (\code{padj < alpha}), \code{Anchor} (whether the item is in the
#'   anchor) and, for the forward anchor, \code{Inclusion} (the share of
#'   tie-break orderings whose anchor contains the item). At threshold level
#'   there is one row per threshold, with a \code{Threshold} column
#'   (\code{"c1"}, \code{"c2"}, ...). At spacing level there is one row per
#'   item with at least two thresholds, with \code{Item}, \code{Chisq},
#'   \code{df}, \code{p}, \code{padj}, \code{Flagged} and
#'   \code{MaxSpacingDiff} (the largest absolute difference in the distance
#'   between adjacent thresholds).
#'
#'   With \code{output = "list"}: a list with
#'   \describe{
#'     \item{\code{table}}{the data.frame above.}
#'     \item{\code{settings}}{a one-row data.frame with \code{Anchor}
#'       (method), \code{Anchor_items} (comma-separated, \code{NA} at
#'       spacing level), \code{Untested}
#'       (the item without a location test, \code{NA} if none),
#'       \code{N_seeds} and \code{Seed}, \code{Level},
#'       \code{Correction}, \code{Alpha}, \code{Conf}, \code{Comparison}
#'       (e.g. \code{"B minus A"}), \code{Reference}, \code{Items},
#'       \code{N} (respondents used), \code{N_excluded} (rows dropped for a
#'       missing \code{dif_var} or no responses), \code{N_incomplete}
#'       (respondents used who skipped some items) and one \code{n_<group>}
#'       column per group.}
#'     \item{\code{distances}}{at spacing level, one row per pair of
#'       adjacent thresholds with \code{Item}, \code{Distance},
#'       \code{Difference}, \code{SE}, \code{lower} and \code{upper};
#'       \code{NULL} at the other levels.}
#'     \item{\code{criterion}}{one row per item with \code{Item},
#'       \code{Location_diff} (the item's location difference after
#'       linking) and \code{Gini} (the Gini criterion when the groups are
#'       aligned at that item; Strobl et al., 2021). \code{NULL} at
#'       spacing level, where no anchor is selected.}
#'   }
#'
#' @details
#' \strong{Why an anchor.} The scale of a Rasch model has no natural origin,
#' so two groups' item parameters can only be compared after the scales are
#' linked. Centring each group on its own mean item location makes every
#' item's difference depend on all other items: DIF in some items appears as opposite DIF in
#' the others (Andrich & Hagquist, 2012; Wang, 2004). Linking through anchor
#' items avoids this when the anchor items are free of DIF. Both anchor
#' methods assume that most items are free of DIF, which cannot be tested.
#' When a majority of items share the same DIF, treating them as the
#' reference is an equally valid solution (Strobl et al., 2021).
#'
#' \strong{Forward anchor.} The iterative forward class with the mean test
#' statistic threshold (MTT) selection (Kopf et al., 2015a, 2015b), from
#' \code{psychotools::anchor()}. Each item is tested against several
#' provisional anchors, the items are ranked by how often their statistic
#' exceeds a threshold, and the anchor grows along that ranking while it is
#' shorter than 0.8 times the number of items without a significant test
#' against the current anchor. The ranking criterion is a count, so items
#' often tie, and psychotools orders tied items at random. The same data can
#' then give different anchors and different flags. RMdifAnchor() therefore
#' repeats the selection for \code{n_seeds} random orderings and uses the
#' items that are in the anchor in at least half of them. Each item's
#' inclusion frequency is reported. An item near .5 is in the anchor under
#' about half of the orderings, and more orderings do not change that. A
#' message names such items, since results can then depend on how ties are
#' broken.
#'
#' \strong{Gini anchor.} The anchor point of Strobl et al. (2021): the
#' groups are aligned where the Gini index of the absolute item
#' differences is largest, that is, where a few items differ a lot and the
#' rest overlap. Its optimum always lies on a single item, so the anchor is
#' one item, which is then not tested. It is deterministic in practice, but
#' it flags items without DIF less often than the nominal level (in the
#' simulations below, mostly well under 5% at \code{alpha = 0.05}), so some
#' DIF may be missed, and it has less power than the forward anchor.
#'
#' \strong{Tests.} At item level, each item's location (the mean of its
#' thresholds) is compared after linking. With an anchor of two or more
#' items every item is tested, including the anchor items. Their tests lean
#' towards non-rejection, since they were selected for looking free of DIF.
#' At threshold level each threshold is compared after linking. At spacing
#' level, the differences in the distances between adjacent thresholds are
#' tested jointly per item (Wald \eqn{\chi^2} with one df fewer than the
#' item's thresholds). Shifting a group's scale moves all its thresholds
#' equally, so the distances need no anchor and cannot show artificial DIF.
#' A spacing difference is differential step functioning (DSF; Henninger
#' et al., 2025). Uniform DIF does not change the spacing.
#'
#' \strong{Correction.} \code{"fwer"} is the single-step adjustment
#' (Hothorn et al., 2008), which takes the correlation between the tests
#' into account: max-|z| for items and thresholds, min-p for spacing. It is
#' computed from 10,000 draws from the joint normal of the estimates, so the
#' smallest adjusted p-value is about .0001.
#'
#' \strong{Evidence.} In simulations with two groups (10 items with 4
#' categories, 20 items with 2 or 3 categories, 150 to 600 respondents per
#' group, equal and unequal sizes, uniform DIF and DSF of 0.6 logits), the
#' forward anchor with the single-step correction kept the family-wise error
#' rate close to .05 in most conditions and had the highest power among the
#' anchors that held it. Its family-wise error rate rose above .05 (up to
#' about .10) when several items had small location shifts that the
#' selection missed, and with three of
#' 10 items with DIF at 150 per group. The consensus over orderings kept
#' the error rate and power. The share of datasets where two runs with
#' different seeds flag different thresholds fell from about .037 with a
#' single ordering to about .012 with 100 orderings.
#' More than two groups, other DIF sizes and anchors of fixed length were
#' not studied. The studies are reported in full at
#' \url{https://pgmj.github.io/rasch_book/chapters/dif-anchor-tests.html}.
#'
#' Every response category must be used in both groups, as in
#' \code{\link{RMdifLR}}. Respondents with missing responses are kept.
#'
#' @references
#' Andrich, D., & Hagquist, C. (2012). Real and artificial differential item
#' functioning. *Journal of Educational and Behavioral Statistics, 37*(3),
#' 387-416. \doi{10.3102/1076998611411913}
#'
#' Henninger, M., Radek, J., Sengewald, M.-A., & Strobl, C. (2025). Partial
#' credit trees meet the partial gamma coefficient for quantifying DIF and
#' DSF in polytomous items. *Behaviormetrika, 52*(2), 221-257.
#' \doi{10.1007/s41237-024-00252-3}
#'
#' Hothorn, T., Bretz, F., & Westfall, P. (2008). Simultaneous inference in
#' general parametric models. *Biometrical Journal, 50*(3), 346-363.
#' \doi{10.1002/bimj.200810425}
#'
#' Kopf, J., Zeileis, A., & Strobl, C. (2015a). A framework for anchor
#' methods and an iterative forward approach for DIF detection. *Applied
#' Psychological Measurement, 39*(2), 83-103.
#' \doi{10.1177/0146621614544195}
#'
#' Kopf, J., Zeileis, A., & Strobl, C. (2015b). Anchor selection strategies
#' for DIF analysis: Review, assessment, and new approaches. *Educational
#' and Psychological Measurement, 75*(1), 22-56.
#' \doi{10.1177/0013164414529792}
#'
#' Strobl, C., Kopf, J., Kohler, L., von Oertzen, T., & Zeileis, A. (2021).
#' Anchor point selection: Scale alignment based on an inequality
#' criterion. *Applied Psychological Measurement, 45*(3), 214-230.
#' \doi{10.1177/0146621621990743}
#'
#' Wang, W.-C. (2004). Effects of anchor item methods on the detection of
#' differential item functioning within the family of Rasch models.
#' *Journal of Experimental Education, 72*(3), 221-261.
#' \doi{10.3200/JEXE.72.3.221-261}
#'
#' @seealso \code{\link{RMdifLR}} for Andersen's likelihood-ratio test of any
#'   DIF (any number of groups), \code{\link{RMdifGamma}} and
#'   \code{\link{RMdifTree}}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("eRm", quietly = TRUE)) { # for the example data
#'   set.seed(1)
#'   data("pcmdat2", package = "eRm")
#'   grp <- factor(sample(c("A", "B"), nrow(pcmdat2), replace = TRUE))
#'
#'   # Item locations, forward anchor, single-step correction
#'   RMdifAnchor(pcmdat2, dif_var = grp)
#'
#'   # Spacing of thresholds (no anchor needed)
#'   RMdifAnchor(pcmdat2, dif_var = grp, level = "spacing")
#'
#'   # Thresholds, as a plot of differences
#'   RMdifAnchor(pcmdat2, dif_var = grp, level = "threshold",
#'               output = "ggplot")
#'
#'   # Gini criterion plot (Strobl et al., 2021)
#'   RMdifAnchor(pcmdat2, dif_var = grp, output = "ggplot",
#'               plot_type = "criterion")
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
RMdifAnchor <- function(
  data,
  dif_var,
  anchor = c("forward", "gini"),
  level = c("item", "spacing", "threshold"),
  output = c("kable", "dataframe", "ggplot", "list"),
  plot_type = c("differences", "criterion"),
  correction = c("fwer", "holm", "fdr_bh", "none"),
  alpha = 0.05,
  conf = 0.95,
  reference = NULL,
  sort = FALSE,
  n_seeds = 100,
  seed = 20260113L
) {
  anchor <- match.arg(anchor)
  level <- match.arg(level)
  output <- match.arg(output)
  plot_type <- match.arg(plot_type)
  correction <- match.arg(correction)
  .check_prob_arg(alpha, "alpha")
  .check_prob_arg(conf, "conf")
  if (!is.logical(sort) || length(sort) != 1L || is.na(sort)) {
    stop("`sort` must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.numeric(n_seeds) || length(n_seeds) != 1L || !is.finite(n_seeds) ||
      n_seeds < 1 || n_seeds != round(n_seeds)
  ) {
    stop("`n_seeds` must be a single positive whole number.", call. = FALSE)
  }
  n_seeds <- as.integer(n_seeds)
  if (
    !is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
         seed != round(seed) || abs(seed) > .Machine$integer.max)
  ) {
    stop(
      "`seed` must be NULL or a single whole number between -",
      .Machine$integer.max, " and ", .Machine$integer.max, ".",
      call. = FALSE
    )
  }
  if (output == "ggplot" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for output = \"ggplot\".",
         call. = FALSE)
  }
  if (output == "kable" && !requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required for output = \"kable\".", call. = FALSE)
  }

  validate_response_data(data)
  prep <- .prepare_dif_data(data, dif_var)
  data <- prep$data
  dif_factor <- prep$dif_factor
  n_total <- prep$n_total
  has_na <- prep$has_na
  groups <- levels(dif_factor)
  if (length(groups) != 2L) {
    stop(
      "`dif_var` must have exactly 2 non-missing groups; it has ",
      length(groups), ". The anchored tests have been checked for two ",
      "groups only. For a test of any DIF across more groups use ",
      "RMdifLR(); compare two groups at a time here, or use RMdifTree().",
      call. = FALSE
    )
  }
  if (is.null(reference)) {
    reference <- groups[1L]
  } else if (
    length(reference) != 1L || is.na(reference) ||
      !as.character(reference) %in% groups
  ) {
    stop(
      "`reference` must be one of the groups in `dif_var`: ",
      paste0("'", groups, "'", collapse = ", "), ".",
      call. = FALSE
    )
  }
  reference <- as.character(reference)
  other <- setdiff(groups, reference)
  # Differences are computed as groups[2] - groups[1], then turned around
  # when the reference is the second group.
  sgn <- if (reference == groups[1L]) 1 else -1

  # --- Data checks and fits ---------------------------------------------
  .check_dif_categories(data, dif_factor, "RMdifAnchor")
  items <- names(data)
  full_max <- apply(as.matrix(data), 2L, max, na.rm = TRUE)
  is_polytomous <- any(full_max > 1L)
  sub <- lapply(groups, function(g) data[dif_factor == g, , drop = FALSE])
  fits <- Map(function(d, g) {
    .cml_fit_or_stop(d, paste0("group '", g, "'"))
  }, sub, groups)
  ts <- Map(function(f, d, g) {
    .cml_threshold_summary(f, d, se = TRUE, label = paste0("group '", g, "'"),
                           max_cat = full_max)
  }, fits, sub, groups)
  m <- lengths(ts[[1L]]$thr_list)
  thr_item <- rep(seq_along(items), m)
  thr_k <- unlist(lapply(m, seq_len), use.names = FALSE)
  d_thr <- unlist(ts[[2L]]$thr_list, use.names = FALSE) -
    unlist(ts[[1L]]$thr_list, use.names = FALSE)
  V <- ts[[1L]]$vcov_thr + ts[[2L]]$vcov_thr
  if (anyNA(V)) {
    stop(
      "The threshold covariance matrix could not be computed in one of the ",
      "groups, so the Wald tests are not available. Small or sparse groups ",
      "are the usual cause.",
      call. = FALSE
    )
  }

  run <- function(expr) {
    if (is.null(seed)) expr else .with_local_seed(as.integer(seed), expr)
  }

  # --- Anchor ------------------------------------------------------------
  # The spacing test needs no anchor. Only the criterion plot uses one there.
  use_anchor <- level != "spacing" ||
    (output == "ggplot" && plot_type == "criterion")
  sel <- if (use_anchor) {
    run(.dif_select_anchor(fits[[1L]], fits[[2L]], anchor, n_seeds,
                           length(items)))
  } else {
    list(anchor = integer(0), freq = numeric(length(items)))
  }
  anchor_idx <- sel$anchor
  near <- use_anchor & anchor == "forward" & sel$freq > 0.4 & sel$freq < 0.6
  if (any(near)) {
    one <- sum(near) == 1L
    message(
      "Inclusion in the forward anchor is between .4 and .6 for ",
      paste(items[near], collapse = ", "), ". ",
      if (one) "This item is" else "These items are",
      " in the anchor under about half of the tie-break orderings, so the ",
      "anchor and some flags can depend on how ties are broken. Compare with ",
      "anchor = \"gini\", or rerun with another `seed`."
    )
  }

  # --- Tests ---------------------------------------------------------------
  res <- run(.dif_anchor_tests(d_thr, V, thr_item, thr_k, anchor_idx,
                               length(items), level, correction))
  z_crit <- stats::qnorm(1 - (1 - conf) / 2)
  if (level == "spacing") {
    if (is.null(res)) {
      stop(
        "No item has two or more thresholds, so there is no spacing to ",
        "test. Use level = \"item\".",
        call. = FALSE
      )
    }
    table_df <- data.frame(
      Item = items[res$item_no],
      Chisq = res$stat, df = res$df, p = res$p, padj = res$padj,
      Flagged = !is.na(res$padj) & res$padj < alpha,
      MaxSpacingDiff = res$max_diff,
      stringsAsFactors = FALSE
    )
    spacing_detail <- res$detail
    spacing_detail$Item <- items[spacing_detail$item_no]
    spacing_detail$Difference <- sgn * spacing_detail$est
    spacing_detail$lower <- spacing_detail$Difference - z_crit * spacing_detail$se
    spacing_detail$upper <- spacing_detail$Difference + z_crit * spacing_detail$se
  } else {
    table_df <- data.frame(
      Item = items[res$item_no],
      stringsAsFactors = FALSE
    )
    if (level == "threshold") table_df$Threshold <- paste0("c", res$k)
    # An untested unit (the single Gini anchor item's location) has no
    # estimate of its own: it fixes the scale.
    table_df$DIF <- ifelse(is.na(res$se), NA_real_, sgn * res$est)
    table_df$SE <- res$se
    table_df$lower <- table_df$DIF - z_crit * table_df$SE
    table_df$upper <- table_df$DIF + z_crit * table_df$SE
    table_df$z <- sgn * res$z
    table_df$p <- res$p
    table_df$padj <- res$padj
    table_df$Flagged <- !is.na(table_df$padj) & table_df$padj < alpha
    table_df$Anchor <- res$item_no %in% anchor_idx
    if (anchor == "forward") table_df$Inclusion <- sel$freq[res$item_no]
  }
  rownames(table_df) <- NULL

  # Location differences on the linked scale, for the list output and the
  # criterion plot only.
  need_loc <- (output == "list" && level != "spacing") ||
    (output == "ggplot" && plot_type == "criterion")
  loc <- if (need_loc) {
    .dif_anchor_tests(d_thr, V, thr_item, thr_k, anchor_idx,
                      length(items), "item", "none", locations_only = TRUE)
  }
  untested <- if (length(anchor_idx) == 1L) items[anchor_idx] else character()
  info <- list(
    method = anchor,
    anchor_items = items[anchor_idx],
    n_seeds = if (anchor == "forward") n_seeds else NA_integer_,
    groups = groups,
    reference = reference,
    level = level,
    correction = correction,
    alpha = alpha,
    n_persons = nrow(data),
    n_per_group = stats::setNames(as.integer(table(dif_factor)), groups),
    n_items = length(items)
  )

  if (output %in% c("dataframe", "list") && sort) {
    table_df <- .dif_anchor_sort(table_df, level)
  }
  if (output == "dataframe") {
    return(table_df)
  }
  if (output == "list") {
    settings <- data.frame(
      Anchor = if (!use_anchor) {
        "none (spacing needs no anchor)"
      } else if (anchor == "forward") {
        "forward (MTT), consensus"
      } else {
        "Gini anchor point"
      },
      Anchor_items = if (use_anchor) {
        paste(info$anchor_items, collapse = ", ")
      } else {
        NA_character_
      },
      Untested = if (length(untested)) untested else NA_character_,
      N_seeds = if (use_anchor) info$n_seeds else NA_integer_,
      Seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
      Level = level,
      Correction = correction,
      Alpha = alpha,
      Conf = conf,
      Comparison = paste0(other, " minus ", reference),
      Reference = reference,
      Items = length(items),
      N = nrow(data),
      N_excluded = n_total - nrow(data),
      N_incomplete = sum(!stats::complete.cases(data)),
      stringsAsFactors = FALSE
    )
    for (g in groups) settings[[paste0("n_", g)]] <- info$n_per_group[[g]]
    distances <- if (level == "spacing") {
      d <- spacing_detail[c("Item", "Distance", "Difference", "se", "lower", "upper")]
      names(d)[names(d) == "se"] <- "SE"
      rownames(d) <- NULL
      d
    } else {
      NULL
    }
    criterion <- if (level != "spacing") {
      lsg <- sgn * loc
      data.frame(
        Item = items,
        Location_diff = lsg,
        Gini = vapply(lsg, function(cc) .gini_index(abs(lsg - cc)), numeric(1L)),
        stringsAsFactors = FALSE
      )
    }
    return(list(table = table_df, settings = settings, distances = distances,
                criterion = criterion))
  }
  n_clause <- .n_caption(
    nrow(data), n_total,
    c(if (has_na) "incomplete responses retained",
      paste0(groups[1L], ": ", info$n_per_group[[1L]], ", ", groups[2L], ": ",
             info$n_per_group[[2L]]))
  )
  if (output == "kable") {
    if (sort) table_df <- .dif_anchor_sort(table_df, level)
    return(render_dif_anchor_kable(table_df, info, n_clause, is_polytomous))
  }
  if (plot_type == "criterion") {
    return(render_dif_criterion_ggplot(sgn * loc, items, anchor_idx, info,
                                       n_clause))
  }
  if (level == "spacing") {
    return(render_dif_spacing_ggplot(spacing_detail, table_df, info,
                                     n_clause, conf, other))
  }
  render_dif_anchor_ggplot(table_df, info, n_clause, conf, other)
}

# =====================================================================
# Internal helpers
# =====================================================================

#' Check a single probability-type argument in (0, 1)
#' @keywords internal
#' @noRd
.check_prob_arg <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    stop("`", name, "` must be a single number in (0, 1).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Select the anchor items
#'
#' Forward: `psychotools::anchor()` (iterative forward class, MTT selection)
#' is run `n_seeds` times in the current random stream. Each run orders tied
#' items at random. The consensus anchor is the items in at least half of
#' the runs (if none, the most frequent ones), in item order. Gini: the
#' constant class of length 1 with the Gini criterion (Strobl et al., 2021).
#'
#' @return A list with `anchor` (item indices) and `freq` (inclusion
#'   frequency per item; 1 for the Gini anchor item, 0 otherwise).
#' @keywords internal
#' @noRd
.dif_select_anchor <- function(fit1, fit2, anchor, n_seeds, k) {
  if (anchor == "gini") {
    a <- psychotools::anchor(fit1, fit2, class = "constant",
                             select = "Gini")$anchor_items
    freq <- numeric(k)
    freq[a] <- 1
    return(list(anchor = as.integer(a), freq = freq))
  }
  runs <- lapply(seq_len(n_seeds), function(i) {
    as.integer(psychotools::anchor(fit1, fit2, class = "forward",
                                   select = "MTT")$anchor_items)
  })
  freq <- tabulate(unlist(runs), k) / n_seeds
  if (max(freq) == 0) {
    stop(
      "The forward anchor search returned no anchor items in any of the ",
      n_seeds, " runs. Use anchor = \"gini\".",
      call. = FALSE
    )
  }
  a <- which(freq >= 0.5)
  if (length(a) == 0L) a <- which(freq == max(freq))
  list(anchor = a, freq = freq)
}

#' Draws from a zero-mean normal with covariance V (singular allowed)
#'
#' Uses the symmetric square root Q Lambda^(1/2) Q', which does not change
#' when eigen() returns an eigenvector with the other sign. A root of the
#' form Q Lambda^(1/2) does, and the sign can differ between BLAS/LAPACK
#' builds, so the same seed would give different draws on different
#' machines (the problem with MASS::mvrnorm()). Cholesky is not used because
#' the linked V is singular.
#' @keywords internal
#' @noRd
.mvn_draws <- function(V, B) {
  V <- (V + t(V)) / 2
  e <- eigen(V, symmetric = TRUE)
  root <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(V)) %*%
    t(e$vectors)
  matrix(stats::rnorm(B * nrow(V)), B, nrow(V)) %*% t(root)
}

#' Anchored Wald tests and spacing tests
#'
#' `d` is the threshold differences (group 2 minus group 1) with covariance
#' `V`. Linking shifts group 2 so that the mean location of the anchor items
#' is the same in both groups: L = I - 1 w', with w giving each anchor
#' item's thresholds weight 1 / (|A| m_j). Item locations are the means of
#' the linked threshold differences. With a single anchor item that item's
#' location difference is 0 and is not tested.
#'
#' @return For "item" and "threshold": a data.frame with item_no, k (for
#'   thresholds), est, se, z, p, padj. For "spacing": a list with item_no,
#'   stat, df, p, padj, max_diff and `detail` (one row per adjacent
#'   distance). With `locations_only = TRUE`, the linked location
#'   differences.
#' @keywords internal
#' @noRd
.dif_anchor_tests <- function(d, V, thr_item, thr_k, anchor_idx, k, level,
                              correction, B = 10000L,
                              locations_only = FALSE) {
  n <- length(d)
  adjust <- function(p, stat_draws = NULL, stat_obs = NULL, larger = TRUE) {
    switch(
      correction,
      none = p,
      holm = stats::p.adjust(p, "holm"),
      fdr_bh = stats::p.adjust(p, "BH"),
      # Floored at the raw p, which the Monte-Carlo estimate can fall below
      # by noise (always possible with a single tested unit).
      fwer = pmax(p, vapply(stat_obs, function(s) {
        (1 + sum(if (larger) stat_draws >= s else stat_draws <= s)) / (B + 1)
      }, numeric(1L)))
    )
  }

  if (level == "spacing") {
    poly <- which(tabulate(thr_item, k) >= 2L)
    if (length(poly) == 0L) return(NULL)
    rows <- list()
    blk <- integer(0)
    C <- NULL
    for (i in poly) {
      idx <- which(thr_item == i)
      Ci <- matrix(0, length(idx) - 1L, n)
      for (r in seq_len(length(idx) - 1L)) {
        Ci[r, idx[r]] <- -1
        Ci[r, idx[r + 1L]] <- 1
      }
      C <- rbind(C, Ci)
      blk <- c(blk, rep(i, nrow(Ci)))
    }
    est <- as.numeric(C %*% d)
    S <- C %*% V %*% t(C)
    Sinv <- lapply(poly, function(i) {
      tryCatch(
        solve(S[blk == i, blk == i, drop = FALSE]),
        error = function(e) {
          stop(
            "The covariance of the threshold distances could not be inverted ",
            "for item number ", i, " in `data`, so the spacing tests are ",
            "not available. Small ",
            "or sparse groups are the usual cause; inspect the counts with ",
            "RMplotTile(data, group = dif_var).",
            call. = FALSE
          )
        }
      )
    })
    stat <- vapply(seq_along(poly), function(j) {
      e <- est[blk == poly[j]]
      as.numeric(t(e) %*% Sinv[[j]] %*% e)
    }, numeric(1L))
    df <- vapply(poly, function(i) sum(blk == i), integer(1L))
    p <- stats::pchisq(stat, df, lower.tail = FALSE)
    padj <- if (correction == "fwer") {
      Z <- .mvn_draws(S, B)
      psim <- vapply(seq_along(poly), function(j) {
        Zi <- Z[, blk == poly[j], drop = FALSE]
        stats::pchisq(rowSums((Zi %*% Sinv[[j]]) * Zi), df[j],
                      lower.tail = FALSE)
      }, numeric(B))
      adjust(p, apply(psim, 1L, min), p, larger = FALSE)
    } else {
      adjust(p)
    }
    se <- sqrt(pmax(diag(S), 0))
    detail <- data.frame(
      item_no = blk,
      Distance = paste0("c", sequence(df) + 1L, " - c", sequence(df)),
      est = est, se = se,
      stringsAsFactors = FALSE
    )
    return(list(
      item_no = poly, stat = stat, df = df, p = p, padj = padj,
      max_diff = vapply(poly, function(i) max(abs(est[blk == i])), numeric(1L)),
      detail = detail
    ))
  }

  w <- numeric(n)
  for (j in anchor_idx) {
    idx <- which(thr_item == j)
    w[idx] <- 1 / (length(anchor_idx) * length(idx))
  }
  L <- diag(n) - matrix(1, n, 1L) %*% matrix(w, 1L)
  D <- as.numeric(L %*% d)
  Cv <- L %*% V %*% t(L)
  if (level == "item") {
    M <- t(vapply(seq_len(k), function(i) {
      as.numeric(thr_item == i) / sum(thr_item == i)
    }, numeric(n)))
    D <- as.numeric(M %*% D)
    Cv <- M %*% Cv %*% t(M)
    if (locations_only) return(D)
    unit_item <- seq_len(k)
    unit_k <- rep(NA_integer_, k)
  } else {
    unit_item <- thr_item
    unit_k <- thr_k
  }
  se <- sqrt(pmax(diag(Cv), 0))
  tested <- se > 1e-8
  z <- ifelse(tested, D / se, NA_real_)
  p <- 2 * stats::pnorm(-abs(z))
  padj <- rep(NA_real_, length(D))
  if (any(tested)) {
    padj[tested] <- if (correction == "fwer") {
      Vt <- Cv[tested, tested, drop = FALSE]
      Z <- .mvn_draws(Vt, B)
      maxz <- apply(abs(sweep(Z, 2L, se[tested], "/")), 1L, max)
      adjust(p[tested], maxz, abs(z[tested]))
    } else {
      adjust(p[tested])
    }
  }
  data.frame(
    item_no = unit_item, k = unit_k, est = D,
    se = ifelse(tested, se, NA_real_), z = z, p = p, padj = padj
  )
}

#' Sort the output table by the size of the difference
#' @keywords internal
#' @noRd
.dif_anchor_sort <- function(df, level) {
  key <- if (level == "spacing") -df$Chisq else -abs(df$DIF)
  out <- df[order(key, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' One sentence describing the anchor
#' @keywords internal
#' @noRd
.dif_anchor_sentence <- function(info) {
  if (info$method == "forward") {
    paste0(
      "Anchor: iterative forward method (MTT), items in at least half of ",
      info$n_seeds, " tie-break orderings: ",
      paste(info$anchor_items, collapse = ", "), "."
    )
  } else {
    paste0(
      "Anchor: Gini anchor point (Strobl et al., 2021): ",
      info$anchor_items, ", which fixes the scale",
      if (info$level != "spacing") " and has no location test" else "", "."
    )
  }
}

#' Phrase for the correction
#' @keywords internal
#' @noRd
.dif_correction_phrase <- function(correction, level) {
  switch(
    correction,
    fwer = paste0("single-step correction (",
                  if (level == "spacing") "min-p" else "max-|z|",
                  ", family-wise error rate)"),
    holm = "Holm correction (family-wise error rate)",
    fdr_bh = "Benjamini-Hochberg correction (false discovery rate)",
    none = "no correction for multiple testing"
  )
}

#' Render the RMdifAnchor table as a kable
#' @keywords internal
#' @noRd
render_dif_anchor_kable <- function(table_df, info, n_clause, is_polytomous) {
  level <- info$level
  display <- .round_display(table_df, list(
    DIF = 3L, SE = 3L, lower = 3L, upper = 3L, z = 2L, Chisq = 2L,
    MaxSpacingDiff = 3L, p = 4L, padj = 4L, Inclusion = 2L, df = 0L
  ))
  if ("Anchor" %in% names(display)) {
    display$Anchor <- ifelse(display$Anchor, "yes", "")
  }
  for (cc in intersect(c("p", "padj"), names(display))) {
    display[[cc]] <- .format_p_table(display[[cc]])
  }
  flag <- display$Flagged
  display$Flagged <- ifelse(flag, "yes", "no")
  # Missing cells (the untested Gini anchor item) are blank, flagged or not.
  id_cols <- c("Item", "Threshold")
  for (cc in setdiff(names(display), id_cols)) {
    col <- display[[cc]]
    col <- if (is.numeric(col)) {
      ifelse(is.na(col), "", format(col, trim = TRUE))
    } else {
      as.character(col)
    }
    display[[cc]] <- ifelse(flag & nzchar(col), paste0("**", col, "**"), col)
  }
  what <- switch(
    level,
    item = "item locations",
    threshold = "thresholds",
    spacing = "threshold spacing (distances between adjacent thresholds)"
  )
  other <- setdiff(info$groups, info$reference)
  caption <- paste0(
    "DIF in ", what, ", ",
    if (is_polytomous) "partial credit" else "Rasch", " model, ",
    "'", other, "' minus '", info$reference, "'. ",
    n_clause, ", ", info$n_items, " items. ",
    if (level == "spacing") {
      "No anchor is needed for spacing. Wald chi-square per item, "
    } else {
      paste0(.dif_anchor_sentence(info), " Wald tests, ")
    },
    .dif_correction_phrase(info$correction, level), ". ",
    if (level != "spacing") {
      "Intervals are not corrected for multiple testing. "
    },
    "**Bold** = padj < ", info$alpha, "."
  )
  align <- ifelse(names(display) %in% c("Item", "Threshold", "Flagged", "Anchor"),
                  "l", "r")
  knitr::kable(display, format = "pipe", row.names = FALSE, caption = caption,
               align = align)
}

#' Shared marks for the difference plots (anchor, not flagged, flagged)
#' @keywords internal
#' @noRd
.dif_marks <- function(alpha, with_anchor = TRUE) {
  labs <- c("anchor item", "not flagged", paste0("padj < ", alpha))
  vals <- list(
    colour = c("#0072B2", "black", "#D55E00"),
    shape = c(15, 16, 17),
    linetype = c("dotted", "solid", "longdash")
  )
  keep <- if (with_anchor) 1:3 else 2:3
  list(
    labels = labs[keep],
    colour = stats::setNames(vals$colour[keep], labs[keep]),
    shape = stats::setNames(vals$shape[keep], labs[keep]),
    linetype = stats::setNames(vals$linetype[keep], labs[keep])
  )
}

#' Difference plot for items or thresholds
#' @keywords internal
#' @noRd
render_dif_anchor_ggplot <- function(table_df, info, n_clause, conf, other) {
  d <- table_df[!is.na(table_df$SE), , drop = FALSE]
  is_thr <- info$level == "threshold"
  mk <- .dif_marks(info$alpha)
  flagged <- d$Flagged
  d$Mark <- factor(
    ifelse(flagged, mk$labels[3L], ifelse(d$Anchor, mk$labels[1L], mk$labels[2L])),
    levels = mk$labels
  )
  item_levels <- unique(table_df$Item)
  d$Item <- factor(d$Item, levels = item_levels)
  d$Row <- if (is_thr) {
    # By threshold number, so c10 comes after c9 and not after c1.
    k <- as.integer(sub("^c", "", table_df$Threshold))
    factor(d$Threshold, levels = rev(unique(table_df$Threshold[order(k)])))
  } else {
    factor(as.character(d$Item), levels = rev(item_levels))
  }
  caption <- er2_caption(paste0(
    "Each point is ", if (is_thr) "a threshold's" else "an item's",
    " location in '", other, "' minus '", info$reference, "'. ",
    .dif_anchor_sentence(info), " ",
    "Bars: ", round(conf * 100), "% CI, not corrected. Flags: Wald tests, ",
    .dif_correction_phrase(info$correction, info$level), ". ",
    n_clause, "."
  ))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$DIF, y = .data$Row)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper,
                   colour = .data$Mark, linetype = .data$Mark),
      width = 0.25, orientation = "y"
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Mark, shape = .data$Mark),
                        size = 2.4) +
    ggplot2::scale_colour_manual(values = mk$colour) +
    ggplot2::scale_shape_manual(values = mk$shape) +
    ggplot2::scale_linetype_manual(values = mk$linetype) +
    ggplot2::labs(
      title = paste0("DIF: anchored ", if (is_thr) "threshold" else "item location",
                     " differences"),
      x = "Difference in location (logits)", y = NULL,
      colour = NULL, shape = NULL, linetype = NULL, caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    er2_axis_margins() +
    er2_plot_caption()
  if (is_thr) {
    p <- p +
      ggplot2::facet_grid(Item ~ ., scales = "free_y", space = "free_y",
                          switch = "y") +
      ggplot2::theme(strip.placement = "outside",
                     strip.text.y.left = ggplot2::element_text(angle = 0))
  }
  p
}

#' Plot of differences in adjacent threshold distances
#' @keywords internal
#' @noRd
render_dif_spacing_ggplot <- function(detail, table_df, info, n_clause, conf,
                                      other) {
  mk <- .dif_marks(info$alpha, with_anchor = FALSE)
  flagged_items <- table_df$Item[table_df$Flagged]
  detail$Mark <- factor(
    ifelse(detail$Item %in% flagged_items, mk$labels[2L], mk$labels[1L]),
    levels = mk$labels
  )
  item_levels <- unique(table_df$Item)
  detail$Item <- factor(detail$Item, levels = item_levels)
  detail$Row <- factor(detail$Distance, levels = rev(unique(detail$Distance)))
  caption <- er2_caption(paste0(
    "Each point is the distance between two adjacent thresholds in '", other,
    "' minus the same distance in '", info$reference, "'. No anchor is ",
    "needed. Bars: ", round(conf * 100), "% CI, not corrected. Flags: joint ",
    "Wald test per item, ", .dif_correction_phrase(info$correction, "spacing"),
    ". ", n_clause, "."
  ))
  ggplot2::ggplot(detail, ggplot2::aes(x = .data$Difference, y = .data$Row)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper,
                   colour = .data$Mark, linetype = .data$Mark),
      width = 0.25, orientation = "y"
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Mark, shape = .data$Mark),
                        size = 2.4) +
    ggplot2::scale_colour_manual(values = mk$colour) +
    ggplot2::scale_shape_manual(values = mk$shape) +
    ggplot2::scale_linetype_manual(values = mk$linetype) +
    ggplot2::facet_grid(Item ~ ., scales = "free_y", space = "free_y",
                        switch = "y") +
    ggplot2::labs(
      title = "DSF: differences in threshold spacing",
      x = "Difference in distance (logits)", y = NULL,
      colour = NULL, shape = NULL, linetype = NULL, caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom", strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 0)) +
    er2_axis_margins() +
    er2_plot_caption()
}

#' Gini criterion over candidate anchor points
#'
#' For a shift c of the second group's scale, the criterion is the Gini
#' index of |d_i - c|, where d_i are the linked item location differences
#' (Strobl et al., 2021). Its local maxima lie where c equals some d_i.
#' @keywords internal
#' @noRd
.gini_index <- function(x) {
  k <- length(x)
  s <- sum(x)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  2 * sum(seq_len(k) * sort(x)) / (k * s) - (k + 1) / k
}

#' Plot the Gini criterion
#' @keywords internal
#' @noRd
render_dif_criterion_ggplot <- function(loc, items, anchor_idx, info,
                                        n_clause) {
  rng <- range(loc)
  pad <- 0.1 * max(diff(rng), 0.5)
  grid <- sort(unique(c(seq(rng[1L] - pad, rng[2L] + pad, length.out = 400L),
                        loc)))
  curve <- data.frame(
    Shift = grid,
    Gini = vapply(grid, function(cc) .gini_index(abs(loc - cc)), numeric(1L))
  )
  pts <- data.frame(
    Shift = loc,
    Gini = vapply(loc, function(cc) .gini_index(abs(loc - cc)), numeric(1L)),
    Item = items,
    Anchor = seq_along(items) %in% anchor_idx
  )
  lab_a <- "anchor item"
  lab_o <- "other item"
  pts$Mark <- factor(ifelse(pts$Anchor, lab_a, lab_o), levels = c(lab_a, lab_o))
  best <- items[which.max(pts$Gini)]
  other <- setdiff(info$groups, info$reference)
  caption <- er2_caption(paste0(
    "Gini index of the absolute item location differences ('", other,
    "' minus '", info$reference, "') when the scales are aligned at each ",
    "point on the x-axis (Strobl et al., 2021). Each local maximum is an ",
    "item at which the groups could be aligned. The highest is at ", best,
    ", the Gini anchor point. Several peaks of similar height mean that ",
    "different sets of items agree with each other. The x-axis is relative ",
    "to the current anchor. ", .dif_anchor_sentence(info), " ", n_clause, "."
  ))
  mk_col <- stats::setNames(c("#0072B2", "black"), c(lab_a, lab_o))
  mk_shp <- stats::setNames(c(15, 16), c(lab_a, lab_o))
  ggplot2::ggplot(curve, ggplot2::aes(x = .data$Shift, y = .data$Gini)) +
    ggplot2::geom_line(colour = "grey30") +
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(colour = .data$Mark, shape = .data$Mark),
      size = 2.4
    ) +
    ggplot2::geom_text(
      data = pts, ggplot2::aes(label = .data$Item),
      vjust = -0.9, size = 3, check_overlap = TRUE
    ) +
    ggplot2::scale_colour_manual(values = mk_col) +
    ggplot2::scale_shape_manual(values = mk_shp) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::labs(
      title = "Anchor point criterion (Gini index)",
      x = "Alignment point: item location difference (logits)",
      y = "Gini index", colour = NULL, shape = NULL, caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    er2_axis_margins() +
    er2_plot_caption()
}
