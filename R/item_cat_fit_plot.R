#' Category functioning given the total score
#'
#' Four diagnostic plots of how each item's response categories are used,
#' all built on the probability of each category given the respondent's
#' total score, \eqn{P(X_i = j \mid r)}. Under the partial credit model
#' this probability depends on the item parameters alone, so no person
#' locations are estimated and the comparisons are not affected by their
#' measurement error. Items are shown in rows. [RMitemCatFit()] gives the
#' test that the calibration plot shows.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Respondents with missing responses are left out (complete cases),
#'   since the total score must cover the same items for everyone.
#' @param type Character. Which plot to draw. Default `"residuals"`, or
#'   `"calibration"` when `fit` is given, since that plot shows the test in
#'   [RMitemCatFit()].
#'   * `"residuals"`: standardised residuals of the number of
#'     respondents in each category, within groups of total scores,
#'     against simulated ranges when `fit` is given.
#'   * `"adjacent"`: for each pair of adjacent categories, the share of
#'     respondents in the upper one, against the total score.
#'   * `"restscore"`: the mean and SD of the rest score among respondents
#'     in each category, observed minus expected, against simulated
#'     intervals.
#'   * `"calibration"`: PAV-adjusted calibration of the cumulative events
#'     \eqn{X_i \ge k} against the model's predicted probability.
#' @param fit Optional result of [RMitemCatFit()] with `output = "list"`,
#'   made from the same `data` (checked). Its simulations are used for the
#'   rest-score and calibration plots, so `iterations` and `seed` are
#'   ignored, and give the ranges in the residual plot. Each item's corrected p from
#'   [RMitemCatFit()] is shown under its label, in bold when below .05
#'   (marked with an asterisk when \pkg{ggtext} is not installed), and the
#'   fit's item labels are used unless `item_labels` is given. Not
#'   available with `dif_var`.
#' @param dif_var Optional vector with one value per row of `data` (factor,
#'   character or numeric, treated as groups), at most four groups. Each
#'   group's curves are then drawn against predictions from one fit to all
#'   respondents, with a band showing the range expected without
#'   differential item functioning (DIF). Available for `type =
#'   "calibration"`, `"adjacent"` and `"residuals"`. Rows with a missing
#'   group are left out with a message. See "Groups" in Details.
#' @param grouping How respondents are grouped by total score for
#'   `type = "residuals"`, as in [RMitemICCPlot()]: `"quantile"` (default),
#'   `class_intervals` groups with about equal numbers of respondents,
#'   `"width"`, `class_intervals` intervals of equal width over the
#'   observed total scores, `"score"`, every observed total score its own
#'   group, or `"manual"`, groups set by `score_breaks`. Groups without
#'   respondents are left out. The caption names the grouping used.
#' @param class_intervals Integer, at least 2. Number of groups for
#'   `grouping = "quantile"` and `"width"`. Default `5L`.
#' @param score_breaks Integer vector for `grouping = "manual"`: the total
#'   scores at which a new group starts, in increasing order. Default
#'   `NULL`.
#' @param iterations Integer, at least 20. Number of simulated datasets
#'   for `type = "restscore"` and `type = "calibration"` when `fit` is not
#'   given. Default `200L`. The ends of each pointwise range are extreme
#'   quantiles of the simulated values, so below about 1000 datasets they,
#'   and which points fall outside them, change with the seed (the caption
#'   says so). With `dif_var`, the number of permutations of the group
#'   labels behind the bands, default `1000L`.
#' @param conf_level Numeric in (0, 1). Coverage of the pointwise bands and
#'   intervals, and with `dif_var` of the band for each curve. Default
#'   `0.95`.
#' @param item_labels Optional character vector of unique item labels for
#'   the row strips, one per column of `data`. Default: column names.
#' @param label_wrap Integer. Characters per line in the row strips.
#'   Default `15L`.
#' @param parallel,n_cores Logical / integer. Run the simulations in
#'   parallel via \pkg{mirai} when available. Default sequential. When
#'   `n_cores` is `NULL`, `getOption("mc.cores")` is used. If that is unset
#'   too, a warning is given and the simulations run sequentially. Not used
#'   with `dif_var`, where the permutations run sequentially.
#' @param seed Optional integer for reproducible simulations. See
#'   [easyRasch2-reproducibility].
#' @param output Character. `"ggplot"` (default) or `"dataframe"` for the
#'   data behind the plot.
#'
#' @details
#' Respondents with the lowest or highest possible total score are left
#' out, since their responses are fixed by the score. The item parameters
#' are estimated by conditional maximum likelihood (`psychotools`). As in
#' [RMitemCatFit()], a category that is empty, or used only by respondents
#' with an extreme total score, stops the function with an error that
#' names the item.
#'
#' A plain comparison of observed and expected category counts would show
#' nothing. The category counts are sufficient statistics for the item
#' parameters, so at the CML estimates the expected counts equal the
#' observed ones. Differences show only once respondents are split by
#' total score, which every plot here does.
#'
#' **Residuals.** For each item, category and score group, the observed
#' number of respondents minus the expected number, divided by its
#' standard deviation under the model. This comparison of observed and
#' expected counts within score groups is the one behind Glas's (1988)
#' R1c test. The standard deviation treats the item parameters as known,
#' but they are estimated from the same counts, so the residuals vary less
#' than a standard normal: with 5 score groups about 2.5% of cells fall
#' beyond \eqn{\pm 1.96} when the model holds, fewer with fewer groups.
#' In cells with few expected respondents the residual can only be large
#' in the positive direction. Fixed lines at \eqn{\pm 1.96} would
#' therefore mislead, and none are drawn. With `fit`, each cell gets the
#' pointwise range of the residual in the datasets [RMitemCatFit()]
#' simulated (regenerated for the chosen grouping, without refitting), which accounts for both. The plot has many cells, so some
#' fall outside their range even when the model holds.
#'
#' **Adjacent categories.** Among respondents in category \eqn{j - 1} or
#' \eqn{j} at total score \eqn{r}, the share in \eqn{j}. Under the model
#' this share rises with \eqn{r}. Observed shares that do not rise
#' suggest that respondents do not distinguish the two categories. The
#' band is the pointwise binomial range around the expected share, given
#' the number of respondents in the pair.
#'
#' **Rest score.** The mean and SD of the rest score (total score
#' minus the item's score) among respondents in each category, observed
#' minus expected. A category used by respondents with too low or too
#' high rest scores, or with too wide a spread, stands out. Intervals come
#' from datasets simulated from the model given each respondent's total
#' score, with the item parameters re-estimated in each, so the
#' uncertainty in the estimates is included.
#' An item that discriminates less than the model expects (underfit) gives
#' a rest-score SD above its expectation in every category, and compressed
#' means: low categories above and high categories below their expectation.
#' An item that discriminates more (overfit, often from local dependence)
#' gives the reverse. The SD panel shows this as clearly as the mean panel,
#' and in the middle category only the SD changes. This is a property of the
#' item, not of its categories (see "What a flag means" in
#' [RMitemCatFit()]). A raised SD in a single category points to that
#' category instead, for example a middle category chosen by respondents at
#' all levels.
#'
#' **Calibration.** For each cumulative event \eqn{X_i \ge k}, the model
#' predicts \eqn{p = P(X_i \ge k \mid r)}. The observed events are
#' regressed on \eqn{p} by pool-adjacent-violators (PAV), which gives the
#' conditional event probability, the best non-decreasing fit to the
#' observed events (Dimitriadis et al., 2021). This is the PAV-adjusted
#' reliability diagram recommended for ordinal outcomes by Säilynoja et
#' al. (2025), with the prediction taken from the total score instead of
#' a person estimate. Agreement puts the curve on the diagonal. The band
#' is the pointwise range of the same curve over datasets simulated given
#' the total scores at the estimated item parameters. A curve that is too
#' flat means the event separates respondents less well than the model
#' predicts, and too steep means better. Because CML reproduces each
#' category's count, the curve cannot lie above or below the diagonal
#' overall, only cross it.
#'
#' Without `dif_var`, all bands and intervals are pointwise. They show where
#' to look, not whether an item fits, and with many points some fall outside
#' even when the model holds.
#'
#' **Groups (`dif_var`).** The item parameters are estimated once, from all
#' respondents used, and every group's curves are drawn against the same
#' predictions. Under the Rasch model without DIF, an item's response is
#' independent of the group given the total score (Kreiner & Christensen,
#' 2011), whatever the groups' ability distributions. Curves that separate
#' between groups therefore point to DIF, while curves that depart from the
#' model together point to misfit of the item for everyone.
#'
#' The band is a simultaneous consistency band (Dimitriadis et al., 2021)
#' for each curve with coverage `conf_level` (default 95%): the range the
#' whole curve keeps in that share of permutations of the group labels within total scores, a permutation that
#' keeps the number of respondents of each group at every score (an exact
#' conditional test, Kreiner, 1987). For the residual plot, a "curve" is all
#' cells of one item and group. The band is centred on the pooled data, not
#' on the model, so an item that misfits for everyone moves all groups'
#' bands together. Each band is cut to the values the permutations reached,
#' which matters for cells the model expects to be almost empty. In the
#' residual plot, the y-axis covers the bars and at least -4 to 4, and a
#' band that reaches further is cut at the panel edge. Points (bars for
#' residuals) outside their band are marked. The bands are not corrected across curves, so with many curves
#' some cross their band by chance. In simulations without DIF, at most
#' about 5% of curves crossed their band (dev study in the package
#' repository, `dev/catfit_dif/`).
#'
#' Two conditions limit the reading. First, the independence given the
#' total score holds under the Rasch model, not when items differ in
#' discrimination (Kreiner & Christensen, 2011). Then, if the groups differ
#' in ability, curves can separate without DIF, so check item fit first
#' ([RMitemCatFit()], [RMitemInfit()]). Second, the total score includes the
#' other items: DIF in several items in the same direction shifts the
#' matching, and DIF-free items can then look different between groups.
#' [RMdifAnchor()] (item and threshold locations) and [RMdifGamma()] test
#' DIF.
#'
#' @return
#' * `output = "ggplot"`: a [ggplot2::ggplot] object, one row of panels
#'   per item.
#' * `output = "dataframe"`: the data behind the plot, with columns
#'   `item` and, by `type`:
#'   * `"residuals"`: `score_group`, `category`, `observed`, `expected`,
#'     `z`, `lower`, `upper`, `iterations_used` (the last three `NA`
#'     without `fit`).
#'   * `"adjacent"`: `pair`, `score`, `n_pair`, `observed`, `expected`,
#'     `lower`, `upper`.
#'   * `"restscore"`: `statistic`, `category`, `observed`, `lower`, `upper`,
#'     `iterations_used`.
#'   * `"calibration"`: `event`, `predicted`, `cep`, `n`, `lower`,
#'     `upper`.
#'
#'   With `dif_var`, each type also has `group` (the group with its n),
#'   `perm_mean` (the mean over the permutations), `lower` and `upper`
#'   (the band),
#'   `outside` (the point is outside its band) and `iterations_used` (the
#'   number of permutations). The non-DIF reference columns are replaced:
#'   `adjacent` keeps `expected` (the model's share), `calibration` and
#'   `adjacent` have no binomial or simulated range, and `residuals` has no
#'   cells with zero expected variance.
#'
#' @references
#' Dimitriadis, T., Gneiting, T., & Jordan, A. I. (2021). Stable reliability
#' diagrams for probabilistic classifiers. *Proceedings of the National
#' Academy of Sciences, 118*(8), e2016191118.
#' \doi{10.1073/pnas.2016191118}
#'
#' Glas, C. A. W. (1988). The derivation of some tests for the Rasch model
#' from the multinomial distribution. *Psychometrika, 53*(4), 525--546.
#' \doi{10.1007/BF02294405}
#'
#' Kreiner, S. (1987). Analysis of multidimensional contingency tables by
#' exact conditional tests: Techniques and strategies. *Scandinavian Journal
#' of Statistics, 14*(2), 97--112.
#'
#' Kreiner, S., & Christensen, K. B. (2011). Item screening in graphical
#' loglinear Rasch models. *Psychometrika, 76*(2), 228--256.
#' \doi{10.1007/s11336-011-9203-y}
#'
#' Säilynoja, T., Johnson, A. R., Martin, O. A., & Vehtari, A. (2025).
#' Recommendations for visual predictive checks in Bayesian workflow.
#' *arXiv*. \doi{10.48550/arXiv.2503.01509}
#'
#' @seealso [RMitemCatFit()] for the test, [RMitemCatProb()] for the
#'   model-implied category curves, [RMitemICCPlot()] for conditional item
#'   characteristic curves.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("eRm", quietly = TRUE)) {
#'   data(pcmdat2, package = "eRm")
#'   RMitemCatFitPlot(pcmdat2)
#'   RMitemCatFitPlot(pcmdat2, type = "adjacent")
#'   RMitemCatFitPlot(pcmdat2, type = "calibration", iterations = 100,
#'                      seed = 1)
#'   # By group (here a random split, so no DIF is expected)
#'   set.seed(1)
#'   group <- sample(c("A", "B"), nrow(pcmdat2), replace = TRUE)
#'   RMitemCatFitPlot(pcmdat2, type = "adjacent", dif_var = group, seed = 1)
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
RMitemCatFitPlot <- function(
  data,
  type = c("residuals", "adjacent", "restscore", "calibration"),
  fit = NULL,
  dif_var = NULL,
  grouping = c("quantile", "width", "score", "manual"),
  class_intervals = 5L,
  score_breaks = NULL,
  iterations = 200L,
  conf_level = 0.95,
  item_labels = NULL,
  label_wrap = 15L,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  output = c("ggplot", "dataframe")
) {
  # With `fit`, the default is the plot of the test that `fit` holds.
  type <- if (missing(type) && !is.null(fit)) "calibration" else match.arg(type)
  output <- match.arg(output)
  dif_mode <- !is.null(dif_var)
  if (dif_mode) {
    if (!is.null(fit)) {
      stop("`fit` cannot be combined with `dif_var`: its simulations are for ",
           "item fit in the whole sample.", call. = FALSE)
    }
    if (type == "restscore") {
      stop("type = \"restscore\" is not available with `dif_var`: the rest ",
           "score given a category depends on the group's ability ",
           "distribution, so a group difference there is not DIF.",
           call. = FALSE)
    }
    # Permutations are fast (about 10 s at 1000 for 20 items and 4 groups),
    # and the simultaneous band uses an extreme quantile of each curve's
    # maximum.
    if (missing(iterations)) iterations <- 1000L
  }
  if (output == "ggplot" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for output = \"ggplot\".",
         call. = FALSE)
  }
  g <- .resolve_grouping(grouping, class_intervals, score_breaks)
  grouping <- g$grouping
  class_intervals <- g$class_intervals
  score_breaks <- g$score_breaks
  .check_iterations(iterations, min = 20L)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
        conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.null(fit) && !inherits(fit, "RMitemCatFit")) {
    stop("`fit` must be the result of RMitemCatFit(..., output = \"list\").",
         call. = FALSE)
  }
  if (dif_mode) {
    validate_response_data(data)
    n_dif_missing <- sum(is.na(dif_var))
    dd <- .prepare_dif_data(data, dif_var)
    prep <- .catfit_prepare(dd$data, who = "RMitemCatFitPlot()")
    # .prepare_dif_data() drops respondents with no responses. Without DIF
    # they fall under "complete cases", so the caption says the same here.
    if (nrow(dd$data) < NROW(data) - n_dif_missing &&
          !"complete cases" %in% prep$quals) {
      prep$quals <- c("complete cases", prep$quals)
    }
    item_labels <- .catfit_item_labels(item_labels, prep$item_names)
    if (!is.null(seed)) set.seed(seed)
    return(.catdif_run(
      prep, droplevels(dd$dif_factor[prep$rows]), type, grouping,
      class_intervals, score_breaks, as.integer(iterations), conf_level,
      item_labels, label_wrap, output, n_total = nrow(data),
      n_dif_missing = n_dif_missing
    ))
  }
  prep <- .catfit_prepare(data, who = "RMitemCatFitPlot()")
  if (!is.null(fit) && !identical(prep$X, fit$prep$X)) {
    stop("`fit` was made from different data than `data`.", call. = FALSE)
  }
  if (is.null(item_labels) && !is.null(fit)) item_labels <- fit$item_labels
  item_labels <- .catfit_item_labels(item_labels, prep$item_names)
  X <- prep$X
  r <- prep$r
  P <- prep$P
  if (type == "residuals") {
    K <- nrow(P[[1L]]) - 1L
    if (grouping == "manual" && any(score_breaks > K)) {
      stop("`score_breaks` contains values above the maximum possible total ",
           "score (", K, ").", call. = FALSE)
    }
    bin <- .cicc_class_bins(r, grouping, class_intervals, score_breaks)
    sg <- .score_group_labels(r, bin)
  }
  n_clause <- paste0(.n_caption(prep$n_used, prep$n_total, prep$quals), ".")
  alpha <- 1 - conf_level

  # Simulations: taken from `fit` when given, otherwise run here for the two
  # types that need them.
  S <- NULL
  n_sims <- NA_integer_
  if (type %in% c("restscore", "calibration")) {
    what <- if (type == "restscore") "rest" else "cep"
    if (!is.null(fit)) {
      S <- fit$sims[[what]]
      iterations <- fit$iterations
    } else {
      iterations <- as.integer(iterations)
      if (!is.null(seed)) set.seed(seed)
      S <- .catfit_simulate(prep, iterations, what, parallel, n_cores)[[what]]
    }
    if (is.null(S)) {
      stop("None of the simulated datasets could be refitted.", call. = FALSE)
    }
    n_sims <- ncol(S)
  }

  tab <- switch(
    type,
    residuals = .catplot_residuals(X, r, P, sg),
    adjacent = .catplot_adjacent(X, r, P, alpha),
    restscore = .catplot_restscore(X, r, P, S, alpha),
    calibration = .catplot_calibration(X, r, P, S, alpha)
  )
  # Residual ranges only from the simulations in `fit`: without them the
  # plot draws no reference lines, since +-1.96 does not hold (the
  # residuals are under-dispersed at the estimated item parameters).
  if (type == "residuals") {
    tab$lower <- NA_real_
    tab$upper <- NA_real_
    tab$iterations_used <- NA_integer_
    if (!is.null(fit)) {
      band <- .catplot_residual_band(fit, sg, alpha, parallel, n_cores)
      tab$lower <- band$lower
      tab$upper <- band$upper
      tab$iterations_used <- band$n
      n_sims <- band$n
      iterations <- fit$iterations
    }
  }
  # Row strips: the item label, and the corrected item p when `fit` is given.
  # Labels are wrapped here, so the p stays on a line of its own.
  strip_lab <- vapply(item_labels, function(l) {
    paste(strwrap(l, width = label_wrap), collapse = "\n")
  }, character(1L), USE.NAMES = FALSE)
  # Corrected p < .05 is set in bold when ggtext is installed, and marked
  # with an asterisk otherwise.
  md_strip <- !is.null(fit) && requireNamespace("ggtext", quietly = TRUE)
  if (!is.null(fit)) {
    p_lab <- vapply(fit$table$padj, .format_p, character(1L))
    sig <- !is.na(fit$table$padj) & fit$table$padj < 0.05
    if (md_strip) {
      # Item labels are text, not markdown: escape what gridtext would parse.
      esc <- function(x) {
        x <- gsub("&", "&amp;", x, fixed = TRUE)
        x <- gsub("<", "&lt;", x, fixed = TRUE)
        x <- gsub("\\", "&#92;", x, fixed = TRUE)
        x <- gsub("[", "&#91;", x, fixed = TRUE)
        x <- gsub("*", "&#42;", x, fixed = TRUE)
        x <- gsub("_", "&#95;", x, fixed = TRUE)
        x <- gsub("`", "&#96;", x, fixed = TRUE)
        # A label starting with "1. ", "1) ", "- ", "+ ", "# " or "> " would
        # become a list, heading or quote, which gridtext cannot draw.
        x <- sub("^(\\s*[0-9]+)\\.", "\\1&#46;", x)
        x <- sub("^(\\s*[0-9]+)\\)", "\\1&#41;", x)
        x <- sub("^(\\s*)-", "\\1&#45;", x)
        x <- sub("^(\\s*)\\+", "\\1&#43;", x)
        x <- sub("^(\\s*)#", "\\1&#35;", x)
        sub("^(\\s*)>", "\\1&#62;", x)
      }
      p_lab <- ifelse(sig, paste0("**", p_lab, "**"), p_lab)
      strip_lab <- paste0(gsub("\n", "<br>", esc(strip_lab), fixed = TRUE),
                          "<br>", p_lab)
    } else {
      strip_lab <- paste0(strip_lab, "\n", p_lab, ifelse(sig, " *", ""))
    }
  }
  tab$item <- factor(strip_lab[tab$item], levels = strip_lab)
  tab <- tab[, c("item", setdiff(names(tab), "item"))]
  rownames(tab) <- NULL
  if (output == "dataframe") {
    tab$item <- factor(item_labels[as.integer(tab$item)], levels = item_labels)
    return(tab)
  }

  pct <- paste0(format(100 * conf_level), "%")
  fit_clause <- if (is.null(fit)) "" else paste0(
    " Under each item: its p from RMitemCatFit(), corrected across items ",
    "with the ", .correction_label(fit$correction), ", ",
    if (md_strip) "in bold" else "marked *", " when below .05.")
  # Range ends are extreme order statistics: with 2.5% in each tail, 1000
  # datasets leave 25 values per tail (the convergence point for 95%
  # intervals in Johansson, 2026, for infit).
  range_note <- if (!is.na(n_sims) && n_sims < 1000L) {
    paste0(" With fewer than 1000 simulated datasets, the ends of the ",
           "ranges, and so which points fall outside them, depend on the ",
           "seed.")
  } else {
    ""
  }
  caption <- switch(
    type,
    residuals = paste0(
      "Observed minus expected number of respondents in each category, ",
      "divided by its standard deviation, within groups of total scores. ",
      "Expected counts are from the probability of each category given the ",
      "total score. ",
      if (is.null(fit)) {
        "Give the RMitemCatFit() result as fit for simulated ranges. "
      } else {
        paste0("Grey boxes: pointwise ", pct, " range in ", n_sims,
               " datasets simulated given each respondent's total score, ",
               "with item parameters re-estimated in each. ")
      },
      .cicc_interval_caption(grouping, class_intervals, score_breaks, bin),
      " ", n_clause,
      if (!is.null(fit)) .attrition_clause(n_sims, iterations),
      range_note),
    adjacent = paste(
      "Points: among respondents in either of two adjacent categories at a",
      "total score, the share in the upper one, sized by the number of",
      "respondents in the pair. Line: expected under the model. Band:",
      paste0("pointwise ", pct, " binomial range."), n_clause),
    restscore = paste0(
      "Mean and SD of the rest score among respondents in each category, ",
      "observed minus expected. Bars: pointwise ", pct, " range in ",
      n_sims, " datasets simulated given each respondent's total ",
      "score, with item parameters re-estimated in each. ", n_clause,
      .attrition_clause(n_sims, iterations), range_note),
    calibration = paste0(
      "PAV-adjusted conditional event probability against the predicted ",
      "probability of the event given the total score, sized by the number ",
      "of respondents. Band: pointwise ", pct, " range in ", n_sims,
      " datasets simulated given the total scores. Dotted line: perfect ",
      "calibration. ", n_clause, range_note)
  )
  caption <- paste0(caption, fit_clause)

  strip <- ggplot2::label_value
  p <- switch(
    type,
    residuals = .catplot_residuals_gg(tab, strip),
    adjacent = .catplot_adjacent_gg(tab, strip),
    restscore = .catplot_restscore_gg(tab, strip),
    calibration = .catplot_calibration_gg(tab, strip)
  )
  p +
    ggplot2::labs(caption = er2_caption(caption)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom",
                   strip.text.y = if (md_strip) {
                     ggtext::element_markdown(angle = 0)
                   } else {
                     ggplot2::element_text(angle = 0)
                   }) +
    er2_plot_caption() +
    er2_axis_margins()
}

# ---------------------------------------------------------------------
# Internal: conditional probabilities and sampling
# ---------------------------------------------------------------------

#' P(X_i = j | r) for every item
#'
#' Ratios of elementary symmetric functions: the item's category weight
#' times the ESF of the other items at score r - j, over the ESF of all
#' items at r. Uses `.esf_convolve()` (direct polynomial multiplication,
#' all terms positive).
#' @param thr_list Threshold vectors, one per item.
#' @return List of (K + 1) x (m_i + 1) matrices, rows r = 0..K.
#' @keywords internal
#' @noRd
.cond_cat_probs <- function(thr_list) {
  w <- lapply(thr_list, function(t) exp(-c(0, cumsum(t))))
  g_all <- Reduce(.esf_convolve, w, 1)
  K <- length(g_all) - 1L
  lapply(seq_along(w), function(i) {
    g_wo <- Reduce(.esf_convolve, w[-i], 1)
    m <- length(w[[i]]) - 1L
    P <- matrix(0, K + 1L, m + 1L)
    for (j in 0:m) {
      s <- 0:K - j
      ok <- s >= 0 & s < length(g_wo)
      P[ok, j + 1L] <- w[[i]][j + 1L] * g_wo[s[ok] + 1L] / g_all[ok]
    }
    P
  })
}

#' Sample response patterns given each respondent's total score
#'
#' Vectorised over respondents: items are drawn in turn from their
#' distribution given the remaining score, using suffix ESFs (the same
#' scheme as `.sim_conditional()`, which draws one pattern at a time).
#' @param thr_list Threshold vectors, one per item.
#' @param totals Integer total scores, all strictly between 0 and the
#'   maximum.
#' @return Integer matrix, one row per total.
#' @keywords internal
#' @noRd
.sim_given_totals <- function(thr_list, totals) {
  w <- lapply(thr_list, function(t) exp(-c(0, cumsum(t))))
  k <- length(w)
  suffix <- vector("list", k + 1L)
  suffix[[k + 1L]] <- 1
  for (i in k:1) suffix[[i]] <- .esf_convolve(suffix[[i + 1L]], w[[i]])
  n <- length(totals)
  out <- matrix(0L, n, k)
  s <- totals
  for (i in seq_len(k)) {
    m <- length(w[[i]]) - 1L
    g <- suffix[[i + 1L]]
    cp <- matrix(0, n, m + 1L)
    for (j in 0:m) {
      idx <- s - j
      ok <- idx >= 0 & idx < length(g)
      cp[ok, j + 1L] <- w[[i]][j + 1L] * g[idx[ok] + 1L]
    }
    for (j in seq_len(m)) cp[, j + 1L] <- cp[, j + 1L] + cp[, j]
    u <- stats::runif(n) * cp[, m + 1L]
    x <- as.integer(rowSums(u > cp[, seq_len(m), drop = FALSE]))
    out[, i] <- x
    s <- s - x
  }
  out
}

# ---------------------------------------------------------------------
# Internal: tables behind the four plots
# ---------------------------------------------------------------------

#' Observed and expected counts per item, score and category
#' @keywords internal
#' @noRd
.cat_score_table <- function(X, r, P) {
  K <- nrow(P[[1L]]) - 1L
  n_r <- tabulate(r + 1L, nbins = K + 1L)
  do.call(rbind, lapply(seq_along(P), function(i) {
    m <- ncol(P[[i]]) - 1L
    O <- matrix(table(factor(r, levels = 0:K), factor(X[, i], levels = 0:m)),
                K + 1L)
    d <- expand.grid(score = 0:K, category = 0:m)
    d$item <- i
    d$n_r <- n_r[d$score + 1L]
    d$observed <- O[cbind(d$score + 1L, d$category + 1L)]
    d$p <- P[[i]][cbind(d$score + 1L, d$category + 1L)]
    d$expected <- d$n_r * d$p
    d[d$n_r > 0, ]
  }))
}

#' Score-group factor labelled by the total scores it holds
#'
#' Groups from `.cicc_class_bins()` (shared with RMitemICCPlot()), relabelled
#' "Score s" or "Scores a-b" by the scores observed in each, with empty
#' groups dropped.
#' @return Factor of group labels, one per element of `r`.
#' @keywords internal
#' @noRd
.score_group_labels <- function(r, bin) {
  bin <- droplevels(bin)
  # Labels built without names: tapply() names would be re-encoded, which
  # warns for the en dash in a non-UTF-8 locale.
  lab <- vapply(seq_len(nlevels(bin)), function(j) {
    s <- r[as.integer(bin) == j]
    if (min(s) == max(s)) paste("Score", s[1L])
    else paste0("Scores ", min(s), "\u2013", max(s))
  }, character(1L))
  factor(lab[as.integer(bin)], levels = lab)
}

#' @keywords internal
#' @noRd
.catplot_residuals <- function(X, r, P, sg) {
  s <- .residual_sums(X, r, P, as.integer(sg))
  G <- nlevels(sg)
  do.call(rbind, lapply(seq_along(s), function(i) {
    m1 <- ncol(s[[i]]$O)
    data.frame(
      item = i,
      score_group = factor(rep(levels(sg), each = m1), levels = levels(sg)),
      category = rep(seq_len(m1) - 1L, times = G),
      observed = as.vector(t(s[[i]]$O)),
      expected = as.vector(t(s[[i]]$E)),
      z = .residual_z(s[[i]])
    )
  }))
}

#' Observed, expected and variance sums per item, score group and category
#'
#' Counts within each total score, summed over the scores in each group.
#' The variance is binomial within each score, treating the item
#' parameters as known.
#' @param g Integer score group (1..G) of each respondent.
#' @return List over items of list(O, E, V), each a G x (m_i + 1) matrix.
#' @keywords internal
#' @noRd
.residual_sums <- function(X, r, P, g) {
  K <- nrow(P[[1L]]) - 1L
  n_r <- tabulate(r + 1L, nbins = K + 1L)
  ok <- n_r > 0L
  g_score <- integer(K + 1L)
  g_score[r + 1L] <- g
  lapply(seq_along(P), function(i) {
    m1 <- ncol(P[[i]])
    O <- matrix(tabulate(r * m1 + X[, i] + 1L, nbins = (K + 1L) * m1),
                K + 1L, m1, byrow = TRUE)
    E <- n_r * P[[i]]
    V <- E * (1 - P[[i]])
    gs <- g_score[ok]
    list(O = rowsum(O[ok, , drop = FALSE], gs, reorder = TRUE),
         E = rowsum(E[ok, , drop = FALSE], gs, reorder = TRUE),
         V = rowsum(V[ok, , drop = FALSE], gs, reorder = TRUE))
  })
}

#' Standardised residuals from one item's sums, row-major (group, category)
#' @keywords internal
#' @noRd
.residual_z <- function(s) {
  z <- ifelse(s$V > 0, (s$O - s$E) / sqrt(s$V), NA_real_)
  as.vector(t(z))
}

#' Pointwise range of the residuals in the datasets RMitemCatFit() simulated
#'
#' Each refitted dataset is regenerated from its stored seed (sampled given
#' the observed total scores, so the score groups are the same as in the
#' data) and its residuals computed at its re-estimated thresholds. The
#' caller's random number state is restored.
#' @return list(lower, upper, n), in the row order of .catplot_residuals().
#' @keywords internal
#' @noRd
.catplot_residual_band <- function(fit, sg, alpha, parallel, n_cores) {
  sims <- fit$sims
  prep <- fit$prep
  if (is.null(sims$seeds) || is.null(sims$thr)) {
    stop("`fit` holds no simulations for type = \"residuals\". Rerun ",
         "RMitemCatFit() with this version of easyRasch2.", call. = FALSE)
  }
  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv())) {
        rm(".Random.seed", envir = globalenv())
      }
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  g <- as.integer(sg)
  thr <- prep$thr
  r <- prep$r
  one <- function(b) {
    # The same generator and seed as in .catfit_simulate().
    set.seed(sims$seeds[b], kind = "Mersenne-Twister",
             normal.kind = "Inversion", sample.kind = "Rejection")
    Xb <- .sim_given_totals(thr, r)
    s <- .residual_sums(Xb, r, .cond_cat_probs(sims$thr[[b]]), g)
    unlist(lapply(s, .residual_z))
  }
  Z <- do.call(cbind, .maybe_parallel_lapply(seq_along(sims$seeds), one,
                                             parallel, n_cores))
  q <- function(p) apply(Z, 1L, function(z) {
    if (all(is.na(z))) NA_real_ else stats::quantile(z, p, na.rm = TRUE,
                                                      names = FALSE)
  })
  list(lower = q(alpha / 2), upper = q(1 - alpha / 2), n = ncol(Z))
}

#' @keywords internal
#' @noRd
.catplot_adjacent <- function(X, r, P, alpha) {
  tab <- .cat_score_table(X, r, P)
  do.call(rbind, lapply(split(tab, tab$item), function(d) {
    m <- max(d$category)
    do.call(rbind, lapply(seq_len(m), function(j) {
      lo <- d[d$category == j - 1L, c("score", "observed", "expected")]
      hi <- d[d$category == j, c("score", "observed", "expected")]
      a <- merge(lo, hi, by = "score", suffixes = c("_lo", "_hi"))
      a <- a[a$expected_lo + a$expected_hi > 0, ]
      a$n_pair <- a$observed_lo + a$observed_hi
      a <- a[a$n_pair > 0, ]
      p_exp <- a$expected_hi / (a$expected_lo + a$expected_hi)
      data.frame(
        item = d$item[1L],
        pair = paste0(j, " vs ", j - 1L),
        score = a$score,
        n_pair = a$n_pair,
        observed = a$observed_hi / a$n_pair,
        expected = p_exp,
        lower = stats::qbinom(alpha / 2, a$n_pair, p_exp) / a$n_pair,
        upper = stats::qbinom(1 - alpha / 2, a$n_pair, p_exp) / a$n_pair
      )
    }))
  }))
}

#' Rest-score mean and SD per item and category, observed minus expected
#'
#' Expected values weight each rest score by the expected count at that
#' total score. A category no respondent used gives NA.
#' @return data.frame with item, statistic, category, diff.
#' @keywords internal
#' @noRd
.restscore_stats <- function(X, r, P) {
  tab <- .cat_score_table(X, r, P)
  tab$rest <- tab$score - tab$category
  tab <- tab[tab$rest >= 0, ]
  wm <- function(w, x) sum(w * x) / sum(w)
  ws <- function(w, x) sqrt(sum(w * (x - wm(w, x))^2) / sum(w))
  out <- do.call(rbind, lapply(split(tab, list(tab$item, tab$category),
                                     drop = TRUE), function(d) {
    o_ok <- sum(d$observed) > 0
    data.frame(
      item = d$item[1L], category = d$category[1L],
      statistic = c("mean", "sd"),
      diff = c(if (o_ok) wm(d$observed, d$rest) - wm(d$expected, d$rest)
               else NA_real_,
               if (o_ok) ws(d$observed, d$rest) - ws(d$expected, d$rest)
               else NA_real_)
    )
  }))
  out[order(out$item, out$statistic, out$category), ]
}

#' @keywords internal
#' @noRd
.catplot_restscore <- function(X, r, P, S, alpha) {
  obs <- .restscore_stats(X, r, P)
  data.frame(
    item = obs$item,
    statistic = obs$statistic,
    category = obs$category,
    observed = obs$diff,
    lower = apply(S, 1L, stats::quantile, alpha / 2, na.rm = TRUE),
    upper = apply(S, 1L, stats::quantile, 1 - alpha / 2, na.rm = TRUE),
    iterations_used = ncol(S)
  )
}

#' Predicted P(X_i >= k | r) for each respondent, item and event
#' @return List over items of n x m_i matrices.
#' @keywords internal
#' @noRd
.event_forecasts <- function(P, r) {
  lapply(P, function(Pi) {
    m <- ncol(Pi) - 1L
    E <- vapply(seq_len(m), function(k) {
      rowSums(Pi[, (k + 1L):(m + 1L), drop = FALSE])
    }, numeric(nrow(Pi)))
    matrix(E, ncol = m)[r + 1L, , drop = FALSE]
  })
}

#' Weighted PAV of y on p, ties in p pooled first
#'
#' Returns the conditional event probability (CEP) per distinct value of
#' p. Agrees with `stats::isoreg()` on the tie-pooled data.
#' @return data.frame with p, cep, n.
#' @keywords internal
#' @noRd
.pav <- function(p, y) {
  u <- sort(unique(p))
  idx <- match(p, u)
  n <- tabulate(idx, length(u))
  v <- as.numeric(rowsum(y, idx, reorder = TRUE)[, 1L] / n)
  # Stack of blocks: value, weight, number of distinct p values.
  bv <- bw <- numeric(0)
  bl <- integer(0)
  for (j in seq_along(u)) {
    bv <- c(bv, v[j]); bw <- c(bw, n[j]); bl <- c(bl, 1L)
    while (length(bv) > 1L && bv[length(bv) - 1L] > bv[length(bv)]) {
      t <- length(bv)
      bv[t - 1L] <- (bv[t - 1L] * bw[t - 1L] + bv[t] * bw[t]) /
        (bw[t - 1L] + bw[t])
      bw[t - 1L] <- bw[t - 1L] + bw[t]
      bl[t - 1L] <- bl[t - 1L] + bl[t]
      bv <- bv[-t]; bw <- bw[-t]; bl <- bl[-t]
    }
  }
  data.frame(p = u, cep = rep(bv, bl), n = n)
}

#' CEPs of one simulated dataset at the observed forecasts
#' @return Numeric vector, CEPs stacked over items and events in the order
#'   of `.catplot_calibration()`.
#' @keywords internal
#' @noRd
.pav_replicate <- function(Xb, ev) {
  unlist(lapply(seq_along(ev), function(i) {
    lapply(seq_len(ncol(ev[[i]])), function(k) {
      .pav(ev[[i]][, k], as.numeric(Xb[, i] >= k))$cep
    })
  }))
}

#' @keywords internal
#' @noRd
.catplot_calibration <- function(X, r, P, S, alpha) {
  ev <- .event_forecasts(P, r)
  obs <- do.call(rbind, lapply(seq_along(ev), function(i) {
    do.call(rbind, lapply(seq_len(ncol(ev[[i]])), function(k) {
      d <- .pav(ev[[i]][, k], as.numeric(X[, i] >= k))
      data.frame(item = i, event = paste0("X \u2265 ", k),
                 predicted = d$p, cep = d$cep, n = d$n)
    }))
  }))
  obs$lower <- apply(S, 1L, stats::quantile, alpha / 2)
  obs$upper <- apply(S, 1L, stats::quantile, 1 - alpha / 2)
  obs
}

# ---------------------------------------------------------------------
# Internal: plots. Facet labeller and theme are added by the caller.
# ---------------------------------------------------------------------

#' @keywords internal
#' @noRd
.catplot_residuals_gg <- function(tab, strip) {
  tab$sign <- factor(ifelse(tab$z >= 0, "More than expected",
                            "Fewer than expected"),
                     levels = c("More than expected", "Fewer than expected"))
  has_band <- any(!is.na(tab$lower))
  p <-ggplot2::ggplot(tab, ggplot2::aes(factor(.data$category), .data$z,
                                         fill = .data$sign))
  if (has_band) {
    p <- p + ggplot2::geom_tile(
      ggplot2::aes(x = factor(.data$category),
                   y = (.data$lower + .data$upper) / 2,
                   height = .data$upper - .data$lower),
      inherit.aes = FALSE, fill = "grey85", width = 0.85, na.rm = TRUE
    )
  }
  p +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_col(width = if (has_band) 0.5 else 0.7, colour = "grey30",
                      linewidth = 0.2, na.rm = TRUE) +
    ggplot2::scale_fill_manual(NULL, values = c(
      "More than expected" = "#e34948", "Fewer than expected" = "#0d366b"
    ), drop = FALSE) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$score_group)) +
    ggplot2::labs(x = "Category", y = "Standardised residual")
}

#' @keywords internal
#' @noRd
.catplot_adjacent_gg <- function(tab, strip) {
  tab$pair <- paste("Categories", tab$pair)
  ggplot2::ggplot(tab, ggplot2::aes(.data$score)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "grey85") +
    ggplot2::geom_line(ggplot2::aes(y = .data$expected,
                                    linetype = "Expected"),
                       colour = "#1c5cab", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(y = .data$observed, size = .data$n_pair,
                                     shape = "Observed")) +
    ggplot2::scale_size_area("Respondents in the pair", max_size = 3) +
    ggplot2::scale_shape_manual(NULL, values = 1) +
    ggplot2::scale_linetype_manual(NULL, values = "solid") +
    ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$pair)) +
    ggplot2::labs(x = "Total score", y = "Share in the upper category")
}

#' @keywords internal
#' @noRd
.catplot_restscore_gg <- function(tab, strip) {
  tab$statistic <- factor(c(mean = "Rest-score mean",
                            sd = "Rest-score SD")[tab$statistic],
                          levels = c("Rest-score mean", "Rest-score SD"))
  tab$outside <- factor(
    ifelse(tab$observed < tab$lower | tab$observed > tab$upper,
           "Outside the range", "Inside the range"),
    levels = c("Inside the range", "Outside the range")
  )
  ggplot2::ggplot(tab, ggplot2::aes(factor(.data$category))) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40",
                        linetype = "dotted") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower,
                                         ymax = .data$upper),
                            colour = "grey85", linewidth = 2.5,
                            na.rm = TRUE) +
    # Outside points differ in both fill (shape) and colour, so they stand
    # out without colour vision. Light bars keep the vermillion visible.
    ggplot2::geom_point(ggplot2::aes(y = .data$observed,
                                     shape = .data$outside,
                                     colour = .data$outside),
                        size = 2.2, na.rm = TRUE) +
    ggplot2::scale_shape_manual(NULL, values = c(1, 16), drop = FALSE) +
    ggplot2::scale_colour_manual(NULL, values = c("black", "#D55E00"),
                                 drop = FALSE) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$statistic),
                        scales = "free_y") +
    ggplot2::labs(x = "Category", y = "Observed minus expected")
}

#' @keywords internal
#' @noRd
.catplot_calibration_gg <- function(tab, strip) {
  ggplot2::ggplot(tab, ggplot2::aes(.data$predicted)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40",
                         linetype = "dotted") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "grey85") +
    ggplot2::geom_step(ggplot2::aes(y = .data$cep, linetype = "PAV-adjusted"),
                       colour = "#1c5cab", linewidth = 0.6,
                       direction = "mid") +
    ggplot2::geom_point(ggplot2::aes(y = .data$cep, size = .data$n,
                                     shape = "PAV-adjusted"),
                        colour = "#1c5cab") +
    ggplot2::scale_size_area("Respondents", max_size = 2.5) +
    ggplot2::scale_shape_manual(NULL, values = 16) +
    ggplot2::scale_linetype_manual(NULL, values = "solid") +
    ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1),
                                labels = c("0", ".5", "1")) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1),
                                labels = c("0", ".5", "1")) +
    # Not coord_equal(): with one row per item, square panels leave most of
    # the width empty. The dotted diagonal still marks perfect calibration.
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$event)) +
    ggplot2::labs(x = "Predicted probability given the total score",
                  y = "Conditional event probability")
}
