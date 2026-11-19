#' Practical impact of DIF on group comparisons and person scores
#'
#' Estimates how much uniform DIF in the given items biases a comparison of
#' two groups and the scores of persons, so that DIF flagged by
#' \code{\link{RMdifAnchor}} or \code{\link{RMdifTree}} can be judged by its
#' size. The items are given a focal-group shift in a split model, and the
#' estimates when the items are kept (scored with pooled parameters) are
#' compared with the split model's. Everything is reported on several
#' scales, without a threshold or a verdict.
#'
#' @param data A data.frame or matrix of item responses (non-negative
#'   integers, 0-based). One column per item, one row per person. Pass the
#'   grouping variable separately via \code{dif_var}.
#' @param dif_var Vector of length \code{nrow(data)} (factor, character or
#'   numeric) with exactly two groups. Rows where it is \code{NA} are dropped
#'   with a message.
#' @param items The items with DIF, by column name or position. At least
#'   one item must stay common to both groups.
#' @param reference The reference group, one of the two levels of
#'   \code{dif_var}. Default \code{NULL}, the first level. The other group
#'   is the focal group, and differences are focal minus reference.
#' @param cut Optional cut-off for the classification rates and the local
#'   bias at the cut. On the logit scale (\code{cut_scale = "theta"}) or a
#'   raw sum score (\code{cut_scale = "sum"}, a whole number from 1 to the
#'   maximum score). A person is classified at or above the cut when their
#'   WLE (or raw score) is at or above it. Default \code{NULL}, no rates.
#' @param cut_scale \code{"theta"} (default) or \code{"sum"}.
#' @param target Optional target population for a projection: a list with
#'   \code{n} (group sizes, reference then focal) and optionally
#'   \code{mean} and \code{sd} (the latent mean and SD of each group, on the
#'   scale described in Details; default the estimates from the data).
#'   Elements may also be named by the group labels. Default \code{NULL},
#'   no projection.
#' @param conf Confidence level of the intervals. Default \code{0.95}.
#' @param table With \code{output = "kable"} or \code{"dataframe"}, the
#'   table to return: \code{"impact"} (default, the group comparison),
#'   \code{"dif"} (the estimated shifts), \code{"items"} (the impact of
#'   each item, with two or more items), \code{"local"} (the local bias),
#'   \code{"rates"} (needs \code{cut}) or \code{"projection"} (needs
#'   \code{target}).
#' @param output \code{"kable"} (default), \code{"dataframe"},
#'   \code{"ggplot"} (the local bias curve) or \code{"list"} (all tables,
#'   see Value).
#'
#' @return A data.frame, a \code{knitr_kable}, a \code{ggplot} or a list,
#'   depending on \code{output}. Nothing is attached as attributes.
#'
#'   With \code{output = "list"}: a list with
#'   \describe{
#'     \item{\code{impact}}{one row per estimator (\code{Estimator}:
#'       latent mean, WLE mean, sum score; \code{Units}) with
#'       \code{Delta_keep} and \code{Delta_split} (the group difference,
#'       focal minus reference, when the items are kept and when they are
#'       split), \code{SE_Delta} (the sampling SE of \code{Delta_split} in
#'       these data), \code{Impact} (\code{Delta_keep - Delta_split}, the
#'       bias in the group difference from keeping the items), \code{SE},
#'       \code{lower} and \code{upper}, \code{Pct_of_Delta}
#'       (\code{100 * Impact / abs(Delta_split)}), \code{SMD}
#'       (\code{Impact} divided by the estimator's SD under split, the
#'       square root of the average of the two group variances) and
#'       \code{Impact_in_SE}
#'       (\code{Impact / SE_Delta}).}
#'     \item{\code{dif}}{one row per item in \code{items} with \code{Item},
#'       \code{DIF} (the focal shift in logits, positive when the item is
#'       harder for the focal group), \code{SE}, \code{lower} and
#'       \code{upper}.}
#'     \item{\code{items}}{with two or more items, the impact of each:
#'       one row per item and estimator with \code{Item},
#'       \code{Estimator}, \code{Units}, \code{Impact} (the change in the
#'       group difference when this item is kept and the others stay
#'       split), \code{SE}, \code{lower} and \code{upper}. \code{NULL}
#'       with one item.}
#'     \item{\code{local}}{the local bias: one row per summary with
#'       \code{Quantity}, \code{Estimate}, \code{SE}, \code{lower},
#'       \code{upper} and \code{Theta} (where it occurs). The rows are the
#'       local bias where its absolute value is largest, in logits and in
#'       SEM units, and, with a \code{cut}, the local bias at the cut in
#'       logits and in SEM units. All four are signed.}
#'     \item{\code{rates}}{with a \code{cut}: one row per population
#'       (\code{"sample"}, and \code{"target"} with a \code{target}),
#'       group, strategy (\code{"keep"}, \code{"split"},
#'       \code{"keep - split"}) and rate, with \code{Estimate}, \code{SE},
#'       \code{lower} and \code{upper}. The rates are shares of the group:
#'       \code{"at or above"} (classified at or above the cut),
#'       \code{"false negative"} (latent value at or above the cut,
#'       classified below) and \code{"false positive"}. \code{NULL}
#'       without a cut.}
#'     \item{\code{rule}}{with a \code{cut}: the classification rule, one
#'       row per group and strategy with \code{Lowest_raw_score} (the
#'       lowest raw score classified at or above the cut; the maximum
#'       score plus 1 when none is) and, for a cut on the logit scale,
#'       \code{WLE} (that score's WLE).}
#'     \item{\code{projection}}{with a \code{target}: one row per
#'       estimator with \code{Impact}, \code{SE}, \code{lower},
#'       \code{upper}, \code{SE_Delta} (the SE of the group difference in
#'       a study of the target's size) and \code{Impact_in_SE}.}
#'     \item{\code{curve}}{the local bias over the central 95% of the
#'       latent distribution: \code{Theta}, \code{b}, \code{SE},
#'       \code{lower}, \code{upper} and \code{SEM} (the reference group's
#'       standard error of measurement at \code{Theta}).}
#'     \item{\code{settings}}{a one-row data.frame with \code{Reference},
#'       \code{Focal}, \code{Comparison}, \code{Items_split},
#'       \code{Common_items}, \code{Cut}, \code{Cut_scale}, \code{Conf},
#'       the latent means and SDs (\code{Mean_ref}, \code{SD_ref},
#'       \code{Mean_focal}, \code{SD_focal}), \code{Range_lo} and
#'       \code{Range_hi} (the central 95% of the latent distribution),
#'       \code{Theta_to_pooled} (see Details), the target (\code{Target_n_ref},
#'       \code{Target_n_focal}, \code{Target_mean_ref},
#'       \code{Target_sd_ref}, \code{Target_mean_focal},
#'       \code{Target_sd_focal}, \code{NA} without one), \code{Items},
#'       \code{N} (respondents used), \code{N_excluded} (rows dropped for a
#'       missing \code{dif_var}, no responses or incomplete responses),
#'       \code{N_incomplete} (of these, rows dropped for incomplete
#'       responses) and one \code{n_<group>} column per group.}
#'   }
#'
#' @details
#' \strong{Models.} Under \emph{keep}, everyone is scored with one set of
#' item parameters, a CML partial credit (or Rasch) fit to both groups.
#' Under \emph{split}, each item in \code{items} has its thresholds shifted
#' by the same amount in the focal group (uniform DIF), the other items are
#' common to both groups and link them, and the model is fitted by CML.
#' Three estimators of the group difference are compared: the latent mean
#' (normal latent distribution per group, fitted by marginal ML with the
#' item parameters fixed), the mean WLE (Warm, 1989), and the raw sum score
#' (under split, the focal group's sum score is equated to the reference
#' group's items by true-score equating). Removing the items gives about
#' the same group difference as splitting them, so it is not reported.
#'
#' \strong{Size of the impact.} A DIF test that is significant says that
#' DIF exists, not that it matters (Meade, 2010). The impact is the bias in
#' the group difference from keeping the items, reported in the
#' estimator's own units, as a percentage of the group difference, as a
#' standardised mean difference (the SD is the square root of the average
#' of the two group variances, so it does not depend on the ratio of the
#' group sizes), and in units of the SE of the group
#' difference in these data. The last one grows with the sample size for
#' the same DIF: it shows whether a confidence interval for the group
#' difference can be trusted, not whether the DIF matters. When the items
#' were flagged in the same data, an impact interval that excludes zero
#' adds nothing to the DIF test. Its size is what to judge.
#'
#' \strong{Local bias.} For a person at latent value \eqn{\theta}, the
#' local bias \eqn{b(\theta)} is how far a focal person lands from a
#' reference person with the same \eqn{\theta} when both are scored with
#' the reference group's parameters (the horizontal distance between the
#' two test characteristic curves). It is negative when the DIF makes the
#' items harder for the focal group. It is reported where its absolute
#' value is largest over the central 95% of the latent distribution of
#' both groups together, with its sign, in logits and in units of the
#' reference group's standard error of measurement (SEM). With DIF items in opposite
#' directions the group-mean bias can cancel while persons in parts of the
#' scale are still shifted.
#'
#' \strong{Rates at a cut.} With a \code{cut}, persons are classified by
#' their WLE (or raw score). Under split, the focal group's WLE comes from
#' its own item parameters, and its raw score is equated. The rule (which
#' raw scores are classified at or above the cut, per group and strategy)
#' is taken from the fit and held fixed (an equated focal score counts as
#' at or above a sum-score cut when it is at or above \code{cut - 0.5}), and its rates are computed from
#' the split model and each group's estimated latent distribution. The
#' rates therefore describe the rule these data give, which can differ
#' from the rule of the population by one raw score, and the rule is
#' reported with them. The difference keep - split is zero whenever the
#' two strategies classify the same raw scores. On the sum-score scale,
#' the true status is a latent value at or above the point where the
#' reference group's expected score is \code{cut - 0.5}.
#'
#' \strong{Projection.} A study that uses the scale in another population,
#' or with other group sizes, gets the same item parameters and DIF but a
#' different bias in the group difference, a different SE, and different
#' rates. With \code{target}, these are computed for the target
#' population from the estimated split model (population CML for the keep
#' parameters). Their uncertainty comes from the present data and does
#' not shrink with the target size. For the rates it includes the
#' uncertainty of the estimated latent means and SDs when \code{target}
#' leaves them at their defaults. The rates in the target population are
#' those of the rule fitted to these data, as when a scoring table from one
#' study is used in another population, while the projected impact refits
#' the pooled model in the target population. With a sum-score cut the keep
#' rule is the cut itself, so only a \code{theta} cut is affected. The
#' projection assumes that the item parameters and the DIF are the same in
#' the target population.
#'
#' \strong{Scale.} Latent values (\code{cut}, the latent means and SDs, the
#' range and location of the local bias, \code{target} means) are on the
#' logit scale of the split model with the reference group's thresholds
#' centred to mean zero. \code{\link{RMpersonParameters}} and
#' \code{\link{RMitemParameters}} centre all thresholds of the pooled
#' model instead. Linked through the common items, a value on this
#' function's scale plus \code{Theta_to_pooled} (in the settings) is the
#' same latent value on that scale.
#'
#' \strong{Intervals.} The impact, the local bias and the projected impact
#' have delta-method intervals on the estimated shifts. The impact of each
#' item has the SE \eqn{|\mathrm{impact}_j / \mathrm{DIF}_j| \times
#' SE(\mathrm{DIF}_j)}. The per-item impacts need not add up to the joint
#' impact (in the simulations their sum was somewhat larger). The rates have delta-method intervals on all
#' parameters of the split model and the latent means and SDs. In
#' simulations with 5 and 10 four-category items, one or two DIF items of
#' 0.5 or 0.8 logits and 250 to 800 persons per group, these intervals
#' covered close to 95% when the true DIF items were split (most cells
#' .92 to .97). The interval for keep - split rates was the least
#' accurate, down to .84 in a few cells with a theta cut when the two
#' rules differed. The results hold for the items given: an item with DIF
#' that is not in \code{items} biases every estimate, and with few items a
#' falsely flagged item can cancel part of a real impact. The studies are
#' reported in full at
#' \url{https://pgmj.github.io/rasch_book/chapters/dif-impact.html}.
#'
#' \strong{Data.} Respondents with missing responses are dropped (with a
#' message), since the latent means are fitted to the raw-score
#' distribution. Every category an item uses must be used in the data,
#' but not in both groups: the shift model shares an item's thresholds
#' across groups.
#'
#' @references
#' Meade, A. W. (2010). A taxonomy of effect size measures for the
#' differential functioning of items and scales. *Journal of Applied
#' Psychology, 95*(4), 728-743. \doi{10.1037/a0018966}
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. *Psychometrika, 54*(3), 427-450.
#' \doi{10.1007/BF02294627}
#'
#' @seealso \code{\link{RMdifAnchor}} to find the items with DIF,
#'   \code{\link{RMdifLR}} and \code{\link{RMdifTree}}.
#'
#' @examples
#' \donttest{
#' data(phq9)
#' items <- phq9[, 1:9]
#'
#' # q9 is the item RMdifAnchor() flags for gender. Impact of keeping it on
#' # the comparison of men and women (rows with missing gender are dropped)
#' RMdifImpact(items, dif_var = phq9$gender, items = "q9")
#'
#' # Local bias, and rates at a sum score of 10
#' RMdifImpact(items, dif_var = phq9$gender, items = "q9", table = "local")
#' RMdifImpact(items, dif_var = phq9$gender, items = "q9", cut = 10,
#'             cut_scale = "sum", table = "rates")
#'
#' # Two items, with the impact of each
#' RMdifImpact(items, dif_var = phq9$gender, items = c("q2", "q9"),
#'             table = "items")
#'
#' # Projection to a survey of 2 000 women and 2 000 men
#' RMdifImpact(items, dif_var = phq9$gender, items = "q9",
#'             target = list(n = c(2000, 2000)), table = "projection")
#'
#' # Local bias curve
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   RMdifImpact(items, dif_var = phq9$gender, items = "q9",
#'               output = "ggplot")
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
RMdifImpact <- function(
  data,
  dif_var,
  items,
  reference = NULL,
  cut = NULL,
  cut_scale = c("theta", "sum"),
  target = NULL,
  conf = 0.95,
  table = c("impact", "dif", "items", "local", "rates", "projection"),
  output = c("kable", "dataframe", "ggplot", "list")
) {
  cut_scale <- match.arg(cut_scale)
  table <- match.arg(table)
  output <- match.arg(output)
  .check_prob_arg(conf, "conf")
  if (missing(items) || is.null(items) || length(items) == 0L) {
    stop("`items` must name at least one item with DIF.", call. = FALSE)
  }
  if (output == "ggplot" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for output = \"ggplot\".",
         call. = FALSE)
  }
  if (output == "kable" && !requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required for output = \"kable\".", call. = FALSE)
  }

  # --- Data ----------------------------------------------------------------
  validate_response_data(data)
  prep <- .prepare_dif_data(data, dif_var)
  data <- prep$data
  dif_factor <- prep$dif_factor
  n_total <- prep$n_total
  groups <- levels(dif_factor)
  if (length(groups) != 2L) {
    stop(
      "`dif_var` must have exactly 2 non-missing groups; it has ",
      length(groups), ".",
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
  focal_lab <- setdiff(groups, reference)

  # Categories in the full data, then in the complete cases used here.
  .check_complete_case_categories(data, who = "RMdifImpact()")
  complete <- stats::complete.cases(data)
  n_incomplete <- sum(!complete)
  if (!all(complete)) {
    message(sum(!complete), " respondent(s) with missing responses dropped ",
            "(RMdifImpact() uses complete cases).")
    data <- data[complete, , drop = FALSE]
    dif_factor <- dif_factor[complete]
  }
  item_names <- names(data)
  items <- .dimp_resolve_items(items, item_names)
  if (length(items) >= length(item_names)) {
    stop(
      "At least one item must stay common to both groups, to link them. ",
      "`items` names all ", length(item_names), " items.",
      call. = FALSE
    )
  }
  X <- as.matrix(data)
  storage.mode(X) <- "integer"
  ec <- .extreme_only_categories(X)
  if (nrow(ec) > 0L) {
    stop(.category_error(
      paste0(.extreme_only_message(ec),
             if (n_incomplete > 0L) {
               paste0(" (In the ", nrow(X), " respondents with complete ",
                      "responses, which RMdifImpact() uses.)")
             },
             " RMdifImpact() stops here."),
      unique(ec$Item)
    ))
  }
  focal <- dif_factor == focal_lab
  n <- c(sum(!focal), sum(focal))
  if (any(n < 2L)) {
    stop("Each group needs at least 2 respondents with complete responses.",
         call. = FALSE)
  }
  .dimp_check_shifts(X, focal, items, c(reference, focal_lab))
  R <- sum(apply(X, 2L, max))

  # --- Cut and target --------------------------------------------------------
  if (!is.null(cut)) {
    if (!is.numeric(cut) || length(cut) != 1L || !is.finite(cut)) {
      stop("`cut` must be a single finite number, or NULL.", call. = FALSE)
    }
    if (cut_scale == "sum" && (cut != round(cut) || cut < 1 || cut > R)) {
      stop("With cut_scale = \"sum\", `cut` must be a whole number from 1 ",
           "to the maximum score (", R, ").", call. = FALSE)
    }
  }
  target <- .dimp_check_target(target, reference, focal_lab)
  shows_table <- output %in% c("kable", "dataframe")
  if (shows_table && table == "rates" && is.null(cut)) {
    stop("table = \"rates\" needs a `cut`.", call. = FALSE)
  }
  if (shows_table && table == "projection" && is.null(target)) {
    stop("table = \"projection\" needs a `target`.", call. = FALSE)
  }
  if (shows_table && table == "items" && length(items) < 2L) {
    stop("table = \"items\" needs two or more items in `items`.",
         call. = FALSE)
  }

  # --- Estimation ------------------------------------------------------------
  est <- .dimp_estimate(X, focal, items, cut, cut_scale, target)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  res <- .dimp_tables(est, z, item_names, reference, focal_lab, cut,
                      cut_scale)

  settings <- data.frame(
    Reference = reference,
    Focal = focal_lab,
    Comparison = paste0(focal_lab, " minus ", reference),
    Items_split = paste(item_names[est$items], collapse = ", "),
    Common_items = length(item_names) - length(est$items),
    Cut = if (is.null(cut)) NA_real_ else cut,
    Cut_scale = if (is.null(cut)) NA_character_ else cut_scale,
    Conf = conf,
    Mean_ref = est$lat_mu[1L], SD_ref = est$lat_sd[1L],
    Mean_focal = est$lat_mu[2L], SD_focal = est$lat_sd[2L],
    Range_lo = est$range[1L], Range_hi = est$range[2L],
    Theta_to_pooled = est$theta_to_pooled,
    Target_n_ref = if (is.null(est$target)) NA_real_ else est$target$n[1L],
    Target_n_focal = if (is.null(est$target)) NA_real_ else est$target$n[2L],
    Target_mean_ref = if (is.null(est$target)) NA_real_ else est$target$mean[1L],
    Target_sd_ref = if (is.null(est$target)) NA_real_ else est$target$sd[1L],
    Target_mean_focal = if (is.null(est$target)) NA_real_ else est$target$mean[2L],
    Target_sd_focal = if (is.null(est$target)) NA_real_ else est$target$sd[2L],
    Items = length(item_names),
    N = nrow(X),
    N_excluded = n_total - nrow(X),
    N_incomplete = n_incomplete,
    stringsAsFactors = FALSE
  )
  settings[[paste0("n_", reference)]] <- n[1L]
  settings[[paste0("n_", focal_lab)]] <- n[2L]

  if (output == "list") {
    return(c(res, list(settings = settings)))
  }
  if (output == "ggplot") {
    return(render_dif_impact_ggplot(res$curve, est, settings, conf,
                                    n_total))
  }
  df <- res[[table]]
  if (output == "dataframe") {
    return(df)
  }
  render_dif_impact_kable(df, table, settings, conf, cut, cut_scale, n_total)
}

# =====================================================================
# Argument helpers
# =====================================================================

#' Item names or positions to sorted column indices
#' @keywords internal
#' @noRd
.dimp_resolve_items <- function(items, item_names) {
  if (is.character(items)) {
    bad <- setdiff(items, item_names)
    if (length(bad)) {
      stop("`items` not found in `data`: ", paste(bad, collapse = ", "), ".",
           call. = FALSE)
    }
    idx <- match(items, item_names)
  } else if (is.numeric(items) && all(is.finite(items)) &&
               all(items == round(items))) {
    idx <- as.integer(items)
    if (any(idx < 1L | idx > length(item_names))) {
      stop("`items` positions must be between 1 and ", length(item_names),
           ".", call. = FALSE)
    }
  } else {
    stop("`items` must be item names or column positions.", call. = FALSE)
  }
  if (anyDuplicated(idx)) {
    stop("`items` lists an item more than once.", call. = FALSE)
  }
  sort(idx)
}

#' Check and order the target population (reference first)
#' @keywords internal
#' @noRd
.dimp_check_target <- function(target, reference, focal_lab) {
  if (is.null(target)) return(NULL)
  if (!is.list(target) || is.null(target$n)) {
    stop("`target` must be a list with at least `n` (group sizes).",
         call. = FALSE)
  }
  bad <- setdiff(names(target), c("n", "mean", "sd"))
  if (length(bad)) {
    stop("Unknown element(s) in `target`: ", paste(bad, collapse = ", "),
         ". Use `n`, `mean` and `sd`.", call. = FALSE)
  }
  ord <- function(v, what) {
    if (is.null(v)) return(NULL)
    if (!is.numeric(v) || length(v) != 2L || any(!is.finite(v))) {
      stop("`target$", what, "` must be two finite numbers (reference, ",
           "focal).", call. = FALSE)
    }
    if (!is.null(names(v)) && all(nzchar(names(v)))) {
      if (!setequal(names(v), c(reference, focal_lab))) {
        stop("The names of `target$", what, "` must be the groups: '",
             reference, "' and '", focal_lab, "'.", call. = FALSE)
      }
      v <- v[c(reference, focal_lab)]
    }
    unname(v)
  }
  out <- list(n = ord(target$n, "n"), mean = ord(target$mean, "mean"),
              sd = ord(target$sd, "sd"))
  if (any(out$n < 2)) {
    stop("`target$n` must be at least 2 in each group.", call. = FALSE)
  }
  if (!is.null(out$sd) && any(out$sd <= 0)) {
    stop("`target$sd` must be positive.", call. = FALSE)
  }
  out
}

# =====================================================================
# Estimation
# =====================================================================

#' All estimates and their gradients-based SEs
#'
#' @param X Complete integer response matrix (items from 0).
#' @param focal Logical, focal group membership.
#' @param items Sorted column indices of the split items.
#' @param target Checked target list or `NULL`.
#' @return A list of raw results for `.dimp_tables()`.
#' @keywords internal
#' @noRd
.dimp_estimate <- function(X, focal, items, cut, cut_scale, target,
                           h = 1e-3) {
  nm <- colnames(X)
  n <- c(sum(!focal), sum(focal))
  m <- apply(X, 2L, max)
  R <- sum(m)
  r <- rowSums(X)
  pr <- list(tabulate(r[!focal] + 1L, R + 1L) / n[1L],
             tabulate(r[focal] + 1L, R + 1L) / n[2L])
  common <- setdiff(seq_along(nm), items)

  # Keep: pooled CML (the data were checked by the caller).
  pooled <- .cml_fit(X, model = "pcm", hessian = FALSE, check = FALSE)
  thr_keep <- .center_thresholds(
    .cml_threshold_summary(pooled, X, se = FALSE)$thr_list
  )
  keep <- .dimp_compare(list(ref = thr_keep, focal = thr_keep), pr, n,
                        equate = FALSE)

  # Split: shift CML, started from the pooled fit with zero shifts.
  stats <- .dimp_data_stats(X, focal, m)
  prep <- .dimp_prep(m, items, stats)
  nd <- length(items)
  start <- c(unname(stats::coef(pooled)), rep(0, nd))
  if (length(start) != prep$P - 1L + nd) start <- NULL
  fit <- .dimp_fit(prep, start = start, hessian = TRUE)
  if (fit$convergence != 0L) {
    warning("The split model did not converge (optim code ",
            fit$convergence, "). Results may be unreliable.", call. = FALSE)
  }
  if (is.null(fit$vcov)) {
    stop("The covariance matrix of the split model could not be computed. ",
         "Small or sparse groups are the usual cause.", call. = FALSE)
  }
  P <- prep$P
  idx <- P - 1L + seq_len(nd)
  V <- fit$vcov
  Vd <- V[idx, idx, drop = FALSE]
  par_hat <- fit$par
  np <- length(par_hat)
  d_hat <- par_hat[idx]
  thr <- .dimp_par_to_thr(par_hat, prep, nm)
  sp <- .dimp_compare(thr, pr, n, equate = TRUE)
  if (anyNA(c(keep$group[, "se"], sp$group[, "se"]))) {
    stop("The latent distribution of a group could not be estimated ",
         "(its information matrix is not positive definite).",
         call. = FALSE)
  }
  not_conv <- character()
  if (any(vapply(c(keep$lat, sp$lat), `[[`, integer(1L), "convergence") != 0L)) {
    not_conv <- c(not_conv, "the latent mean fits")
  }
  lat_par <- lapply(sp$lat, `[[`, "par")
  lat_mu <- unname(vapply(sp$lat, `[[`, numeric(1L), "mu"))
  lat_sd <- unname(vapply(sp$lat, `[[`, numeric(1L), "sd"))
  unit <- function(q) h * (seq_len(nd) == q)
  shifted <- function(q, sgn) .dimp_with_shifts(thr, items, d_hat + sgn * unit(q))

  # Group comparison and the impact (delta method on the shifts).
  impact <- keep$group[, "delta"] - sp$group[, "delta"]
  Gi <- matrix(0, 3L, nd)
  for (q in seq_len(nd)) {
    dp <- .dimp_compare(shifted(q, 1), pr, n, TRUE, lat_par,
                        hessian = FALSE)$group[, "delta"]
    dm <- .dimp_compare(shifted(q, -1), pr, n, TRUE, lat_par,
                        hessian = FALSE)$group[, "delta"]
    Gi[, q] <- -(dp - dm) / (2 * h)
  }
  se_impact <- sqrt(pmax(rowSums((Gi %*% Vd) * Gi), 0))

  # Per-item impacts: keep item q common, the other DIF items split.
  item_rows <- NULL
  if (nd >= 2L) {
    item_rows <- lapply(seq_len(nd), function(q) {
      prep_q <- .dimp_prep(m, items[-q], stats)
      fq <- .dimp_fit(prep_q, start = c(par_hat[seq_len(P - 1L)], d_hat[-q]),
                      hessian = FALSE)
      if (fq$convergence != 0L) {
        not_conv <<- union(not_conv, "the per-item refits")
      }
      cq <- .dimp_compare(.dimp_par_to_thr(fq$par, prep_q, nm), pr, n, TRUE,
                          lat_par, hessian = FALSE)$group[, "delta"]
      imp_q <- cq - sp$group[, "delta"]
      list(impact = imp_q, se = abs(imp_q / d_hat[q]) * sqrt(Vd[q, q]))
    })
  }

  # Classification parts (needed for b at a sum-score cut too).
  pt_hat <- NULL
  theta_cut <- NULL
  if (!is.null(cut)) {
    pt_hat <- .dimp_rate_parts(thr, items,
                               if (cut_scale == "sum") cut else NULL)
    theta_cut <- if (cut_scale == "theta") cut else pt_hat$theta_c_sum
  }

  # Local bias over the central 95% of the latent mixture. The curve and
  # the table use the same grid, so the plotted maximum is the tabled one.
  rng <- .dimp_range(lat_mu, lat_sd, n[2L] / sum(n))
  grid <- seq(rng[1L], rng[2L], length.out = 401L)
  bc <- .dimp_b_at(thr, grid)
  where <- .dimp_local_where(bc, grid)
  lb <- .dimp_local_at(thr, where, theta_cut)
  Gl <- matrix(0, length(lb), nd)
  Gc <- matrix(0, length(grid), nd)
  for (q in seq_len(nd)) {
    Gl[, q] <- (.dimp_local_at(shifted(q, 1), where, theta_cut) -
                  .dimp_local_at(shifted(q, -1), where, theta_cut)) / (2 * h)
    Gc[, q] <- (.dimp_b_at(shifted(q, 1), grid)$b -
                  .dimp_b_at(shifted(q, -1), grid)$b) / (2 * h)
  }
  se_lb <- sqrt(pmax(rowSums((Gl %*% Vd) * Gl), 0))
  names(se_lb) <- names(lb)
  se_curve <- sqrt(pmax(rowSums((Gc %*% Vd) * Gc), 0))

  # Scale offset to the pooled model (through the common items).
  loc <- function(tl) vapply(tl, mean, numeric(1L))
  link <- mean(loc(thr$ref)[common]) - mean(loc(thr_keep)[common])

  # Target population: means and SDs default to the estimates.
  target_default <- c(mean = FALSE, sd = FALSE)
  if (!is.null(target)) {
    target_default <- c(mean = is.null(target$mean), sd = is.null(target$sd))
    if (is.null(target$mean)) target$mean <- lat_mu
    if (is.null(target$sd)) target$sd <- lat_sd
  }

  # Rates of the fitted rule, delta method on all split-model parameters.
  rates <- NULL
  if (!is.null(cut)) {
    wk <- .dimp_wle_by_score(thr_keep) + link
    rule <- .dimp_rule(pt_hat, wk, cut_scale, cut)
    th_c <- function(pt) if (cut_scale == "theta") cut else pt$theta_c_sum
    ms_hat <- c(lat_par[[1L]], lat_par[[2L]])  # mu, log sd per group
    Vms <- matrix(0, 4L, 4L)
    Vms[1:2, 1:2] <- sp$lat[[1L]]$vcov1 / n[1L]
    Vms[3:4, 3:4] <- sp$lat[[2L]]$vcov1 / n[2L]
    rs <- function(pt, ms) {
      .dimp_fixed_rates(pt, rule, ms[c(1L, 3L)], exp(ms[c(2L, 4L)]), th_c(pt))
    }
    # Target means and SDs left at their defaults are the estimates, and
    # carry their uncertainty as in the sample rates.
    rp <- function(pt, ms) {
      mu <- if (target_default[["mean"]]) ms[c(1L, 3L)] else target$mean
      sd <- if (target_default[["sd"]]) exp(ms[c(2L, 4L)]) else target$sd
      .dimp_fixed_rates(pt, rule, mu, sd, th_c(pt))
    }
    est_s <- rs(pt_hat, ms_hat)
    est_p <- if (!is.null(target)) rp(pt_hat, ms_hat)
    cut_raw <- if (cut_scale == "sum") cut else NULL
    # (mu, log sd) follow the item parameters (forward differences).
    Jms <- matrix(0, 4L, np)
    for (p in seq_len(np)) {
      tp <- .dimp_par_to_thr(par_hat + h * (seq_len(np) == p), prep, nm)
      Jms[, p] <- (c(.dimp_mml_fit(tp$ref, pr[[1L]], lat_par[[1L]],
                                   hessian = FALSE)$par,
                     .dimp_mml_fit(tp$focal, pr[[2L]], lat_par[[2L]],
                                   hessian = FALSE)$par) - ms_hat) / h
    }
    Gs <- matrix(0, length(est_s), np)
    Gp <- if (!is.null(target)) matrix(0, length(est_p), np)
    for (p in seq_len(np)) {
      e <- h * (seq_len(np) == p)
      pt_u <- .dimp_rate_parts(.dimp_par_to_thr(par_hat + e, prep, nm), items,
                               cut_raw, rule = FALSE, grid = pt_hat$grid)
      pt_l <- .dimp_rate_parts(.dimp_par_to_thr(par_hat - e, prep, nm), items,
                               cut_raw, rule = FALSE, grid = pt_hat$grid)
      Gs[, p] <- (rs(pt_u, ms_hat + h * Jms[, p]) -
                    rs(pt_l, ms_hat - h * Jms[, p])) / (2 * h)
      if (!is.null(target)) {
        Gp[, p] <- (rp(pt_u, ms_hat + h * Jms[, p]) -
                      rp(pt_l, ms_hat - h * Jms[, p])) / (2 * h)
      }
    }
    Gm <- matrix(0, length(est_s), 4L)
    Gmp <- if (!is.null(target)) matrix(0, length(est_p), 4L)
    for (k in 1:4) {
      e <- h * (1:4 == k)
      Gm[, k] <- (rs(pt_hat, ms_hat + e) - rs(pt_hat, ms_hat - e)) / (2 * h)
      if (!is.null(target)) {
        Gmp[, k] <- (rp(pt_hat, ms_hat + e) - rp(pt_hat, ms_hat - e)) / (2 * h)
      }
    }
    vq <- function(G, Vx) rowSums((G %*% Vx) * G)
    se_s <- sqrt(pmax(vq(Gs, V) + vq(Gm, Vms), 0))
    rates <- list(sample = list(est = est_s, se = se_s))
    if (!is.null(target)) {
      rates$target <- list(est = est_p,
                           se = sqrt(pmax(vq(Gp, V) + vq(Gmp, Vms), 0)))
    }
    rates$rule <- .dimp_rule_cuts(rule)
    rates$wle <- list(keep = wk, ref = pt_hat$wle_ref,
                      focal = pt_hat$wle_focal)
  }

  # Projection of the impact (delta method on the shifts).
  projection <- NULL
  if (!is.null(target)) {
    pj <- .dimp_project(thr, target)
    if (!pj$convergence) not_conv <- c(not_conv, "the projection")
    if (anyNA(pj$se_target)) {
      warning("The SE of the latent mean difference in the target ",
              "population could not be computed (information matrix not ",
              "positive definite), so its SE_Delta and Impact_in_SE are NA.",
              call. = FALSE)
    }
    proj_impact <- function(tq) {
      .dimp_project(tq, target, start = pj$start, se = FALSE)$impact
    }
    Gpj <- matrix(0, 3L, nd)
    for (q in seq_len(nd)) {
      Gpj[, q] <- (proj_impact(shifted(q, 1)) -
                     proj_impact(shifted(q, -1))) / (2 * h)
    }
    projection <- list(impact = pj$impact, se_target = pj$se_target,
                       se = sqrt(pmax(rowSums((Gpj %*% Vd) * Gpj), 0)))
  }

  if (length(not_conv)) {
    warning("Some fits did not converge (", paste(not_conv, collapse = ", "),
            "). Results may be unreliable.", call. = FALSE)
  }

  list(items = items, d_hat = d_hat, se_d = sqrt(diag(Vd)),
       keep = keep$group, split = sp$group, impact = impact,
       se_impact = se_impact, item_rows = item_rows,
       local = lb, se_local = se_lb, where = where, theta_cut = theta_cut,
       curve = data.frame(Theta = grid, b = bc$b, SE = se_curve,
                          SEM = bc$sem),
       rates = rates, projection = projection, target = target,
       lat_mu = lat_mu, lat_sd = lat_sd, range = rng,
       theta_to_pooled = -link, R = R)
}

# =====================================================================
# Tables
# =====================================================================

.dimp_est_labels <- c(latent = "Latent mean", wle = "WLE mean",
                      sum = "Sum score")
.dimp_units <- c(latent = "logits", wle = "logits", sum = "raw points")

#' Build the output tables from `.dimp_estimate()`
#' @keywords internal
#' @noRd
.dimp_tables <- function(est, z, item_names, reference, focal_lab, cut,
                         cut_scale) {
  e3 <- c("latent", "wle", "sum")
  ci <- function(x, se) list(lower = x - z * se, upper = x + z * se)

  ci_i <- ci(est$impact[e3], est$se_impact)
  impact <- data.frame(
    Estimator = unname(.dimp_est_labels[e3]),
    Units = unname(.dimp_units[e3]),
    Delta_keep = unname(est$keep[e3, "delta"]),
    Delta_split = unname(est$split[e3, "delta"]),
    SE_Delta = unname(est$split[e3, "se"]),
    Impact = unname(est$impact[e3]),
    SE = est$se_impact,
    lower = unname(ci_i$lower),
    upper = unname(ci_i$upper),
    Pct_of_Delta = unname(100 * est$impact[e3] / abs(est$split[e3, "delta"])),
    SMD = unname(est$impact[e3] / est$split[e3, "sd_pooled"]),
    Impact_in_SE = unname(est$impact[e3] / est$split[e3, "se"]),
    stringsAsFactors = FALSE
  )

  ci_d <- ci(est$d_hat, est$se_d)
  dif <- data.frame(Item = item_names[est$items], DIF = unname(est$d_hat),
                    SE = unname(est$se_d), lower = unname(ci_d$lower),
                    upper = unname(ci_d$upper), stringsAsFactors = FALSE)

  items_tab <- NULL
  if (!is.null(est$item_rows)) {
    items_tab <- do.call(rbind, lapply(seq_along(est$item_rows), function(q) {
      x <- est$item_rows[[q]]
      cq <- ci(x$impact[e3], x$se[e3])
      data.frame(Item = item_names[est$items[q]],
                 Estimator = unname(.dimp_est_labels[e3]),
                 Units = unname(.dimp_units[e3]),
                 Impact = unname(x$impact[e3]), SE = unname(x$se[e3]),
                 lower = unname(cq$lower), upper = unname(cq$upper),
                 stringsAsFactors = FALSE)
    }))
  }

  lq <- c(b_max = "b where absolute b is largest (logits)",
          b_sem_max = "b / SEM where its absolute value is largest",
          b_cut = "b at the cut (logits)",
          b_sem_cut = "b / SEM at the cut")
  th <- c(b_max = est$where[["theta_b"]],
          b_sem_max = est$where[["theta_b_sem"]],
          b_cut = if (is.null(est$theta_cut)) NA_real_ else est$theta_cut,
          b_sem_cut = if (is.null(est$theta_cut)) NA_real_ else est$theta_cut)
  keep_q <- intersect(names(lq), names(est$local))
  ci_l <- ci(est$local[keep_q], est$se_local[keep_q])
  local <- data.frame(Quantity = unname(lq[keep_q]),
                      Estimate = unname(est$local[keep_q]),
                      SE = unname(est$se_local[keep_q]),
                      lower = unname(ci_l$lower), upper = unname(ci_l$upper),
                      Theta = unname(th[keep_q]), stringsAsFactors = FALSE)

  rates <- NULL
  rule <- NULL
  if (!is.null(est$rates)) {
    grp_lab <- c(ref = reference, focal = focal_lab)
    st_lab <- c(keep = "keep", split = "split", diff = "keep - split")
    rate_lab <- c(prev = "at or above", fn = "false negative",
                  fp = "false positive")
    pops <- intersect(c("sample", "target"), names(est$rates))
    rates <- do.call(rbind, lapply(pops, function(pp) {
      x <- est$rates[[pp]]
      key <- do.call(rbind, strsplit(names(x$est), ".", fixed = TRUE))
      cr <- ci(x$est, x$se)
      data.frame(Population = pp, Group = unname(grp_lab[key[, 2L]]),
                 Strategy = unname(st_lab[key[, 1L]]),
                 Rate = unname(rate_lab[key[, 3L]]),
                 Estimate = unname(x$est), SE = unname(x$se),
                 lower = unname(cr$lower), upper = unname(cr$upper),
                 stringsAsFactors = FALSE)
    }))
    o <- order(match(rates$Population, pops),
               match(rates$Group, grp_lab),
               match(rates$Strategy, st_lab),
               match(rates$Rate, rate_lab))
    rates <- rates[o, , drop = FALSE]
    rownames(rates) <- NULL
    key <- do.call(rbind, strsplit(names(est$rates$rule), ".", fixed = TRUE))
    rule <- data.frame(Group = unname(grp_lab[key[, 2L]]),
                       Strategy = key[, 1L],
                       Lowest_raw_score = unname(est$rates$rule),
                       stringsAsFactors = FALSE)
    if (cut_scale == "theta") {
      tabs <- list(keep.ref = est$rates$wle$keep,
                   keep.focal = est$rates$wle$keep,
                   split.ref = est$rates$wle$ref,
                   split.focal = est$rates$wle$focal)
      rule$WLE <- vapply(seq_len(nrow(rule)), function(i) {
        s <- est$rates$rule[[i]]
        tb <- tabs[[names(est$rates$rule)[i]]]
        if (s + 1L <= length(tb)) tb[s + 1L] else NA_real_
      }, numeric(1L))
    }
  }

  projection <- NULL
  if (!is.null(est$projection)) {
    pj <- est$projection
    cp <- ci(pj$impact[e3], pj$se)
    projection <- data.frame(
      Estimator = unname(.dimp_est_labels[e3]),
      Units = unname(.dimp_units[e3]),
      Impact = unname(pj$impact[e3]), SE = pj$se,
      lower = unname(cp$lower), upper = unname(cp$upper),
      SE_Delta = unname(pj$se_target[e3]),
      Impact_in_SE = unname(pj$impact[e3] / pj$se_target[e3]),
      stringsAsFactors = FALSE
    )
  }

  curve <- est$curve
  curve$lower <- curve$b - z * curve$SE
  curve$upper <- curve$b + z * curve$SE
  curve <- curve[c("Theta", "b", "SE", "lower", "upper", "SEM")]

  list(impact = impact, dif = dif, items = items_tab, local = local,
       rates = rates, rule = rule, projection = projection, curve = curve)
}

# =====================================================================
# Rendering
# =====================================================================

#' Caption clause for the sample
#' @keywords internal
#' @noRd
.dimp_n_clause <- function(settings, n_total) {
  .n_caption(
    settings$N, n_total,
    c(if (settings$N_incomplete > 0L) "complete cases",
      paste0(settings$Reference, ": ", settings[[paste0("n_", settings$Reference)]],
             ", ", settings$Focal, ": ",
             settings[[paste0("n_", settings$Focal)]]))
  )
}

#' Render one RMdifImpact table as a kable
#' @keywords internal
#' @noRd
render_dif_impact_kable <- function(df, table, settings, conf, cut, cut_scale,
                                    n_total) {
  display <- .round_display(df, list(
    Delta_keep = 3L, Delta_split = 3L, SE_Delta = 3L, Impact = 3L, SE = 3L,
    lower = 3L, upper = 3L, Pct_of_Delta = 1L, SMD = 3L, Impact_in_SE = 2L,
    DIF = 3L, Estimate = 3L, Theta = 2L
  ))
  pct <- paste0(round(conf * 100), "%")
  f2 <- function(x) formatC(x, format = "f", digits = 2)
  comp <- paste0("'", settings$Focal, "' minus '", settings$Reference, "'")
  n_clause <- .dimp_n_clause(settings, n_total)
  items_clause <- paste0("Split: ", settings$Items_split, " (",
                         settings$Common_items, " common items). ")
  cut_text <- if (!is.null(cut)) {
    if (cut_scale == "theta") paste0("theta = ", cut) else
      paste0("a raw score of ", cut)
  }
  caption <- switch(
    table,
    impact = paste0(
      "Impact of keeping the DIF items on the group difference (", comp,
      "). Impact = Delta_keep - Delta_split, ", pct, " CI by the delta ",
      "method. Pct_of_Delta: impact as a percentage of |Delta_split|. SMD: ",
      "impact / SD (square root of the average of the two group ",
      "variances). Impact_in_SE: impact / SE of ",
      "Delta_split in these data. No threshold is implied. ", items_clause,
      n_clause, "."
    ),
    dif = paste0(
      "Estimated DIF (uniform focal shift in logits, ", comp, "), with ", pct,
      " CI. Positive: the item is harder for '", settings$Focal, "'. ",
      items_clause, n_clause, "."
    ),
    items = paste0(
      "Impact of each item: change in the group difference (", comp,
      ") when this item is kept and the other DIF items stay split, ", pct,
      " CI. ", items_clause, n_clause, "."
    ),
    local = paste0(
      "Local bias b(theta): the shift of a focal person relative to a ",
      "reference person with the same latent value, when both are scored ",
      "with the reference parameters. Signed values where the absolute ",
      "value is largest over the central 95% ",
      "of the latent distribution (theta ", f2(settings$Range_lo),
      " to ", f2(settings$Range_hi), "), ", pct, " CI. SEM: the ",
      "reference group's standard error of measurement.",
      if (!is.null(cut)) paste0(" Cut at ", cut_text, "."), " ",
      items_clause, n_clause, "."
    ),
    rates = paste0(
      "Classification at ", cut_text, ", as shares of each group, for the ",
      "rule these data give (lowest raw score classified at or above the ",
      "cut, per group and strategy, in the list output). ", pct, " CI by ",
      "the delta method.",
      if (any(df$Population == "target")) {
        paste0(" Target: n = ", settings$Target_n_ref, " + ",
               settings$Target_n_focal, ".")
      }, " ", items_clause, n_clause, "."
    ),
    projection = paste0(
      "Projected impact on the group difference (", comp, ") in a target ",
      "population of n = ", settings$Target_n_ref, " ('", settings$Reference,
      "', latent mean ", f2(settings$Target_mean_ref), ", SD ",
      f2(settings$Target_sd_ref), ") and ", settings$Target_n_focal,
      " ('", settings$Focal, "', latent mean ",
      f2(settings$Target_mean_focal), ", SD ", f2(settings$Target_sd_focal),
      "). SE_Delta: SE of the ",
      "group difference at that size. ", pct, " CI from the present data. ",
      items_clause, n_clause, "."
    )
  )
  align <- ifelse(vapply(display, is.numeric, logical(1L)), "r", "l")
  knitr::kable(display, format = "pipe", row.names = FALSE, caption = caption,
               align = align)
}

#' Local bias curve
#' @keywords internal
#' @noRd
render_dif_impact_ggplot <- function(curve, est, settings, conf, n_total) {
  sem <- rbind(
    data.frame(Theta = curve$Theta, y = curve$SEM),
    data.frame(Theta = curve$Theta, y = -curve$SEM)
  )
  sem$side <- rep(c("upper", "lower"), each = nrow(curve))
  imax <- which.min(abs(curve$Theta - est$where[["theta_b"]]))
  caption <- er2_caption(paste0(
    "b(theta): the shift of a '", settings$Focal, "' person relative to a '",
    settings$Reference, "' person with the same latent value, when both are ",
    "scored with the reference parameters. Band: ", round(conf * 100),
    "% pointwise CI (delta method). Dashed lines: plus and minus the ",
    "reference group's SEM. Range: central 95% of the latent distribution. ",
    "Split: ", settings$Items_split, ". ",
    .dimp_n_clause(settings, n_total), "."
  ))
  p <- ggplot2::ggplot(curve, ggplot2::aes(x = .data$Theta)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_line(data = sem,
                       ggplot2::aes(y = .data$y, group = .data$side),
                       linetype = "dashed", colour = "grey45") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "#0072B2", alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = .data$b), colour = "#0072B2",
                       linewidth = 0.8) +
    ggplot2::geom_point(data = curve[imax, , drop = FALSE],
                        ggplot2::aes(y = .data$b), colour = "#0072B2",
                        size = 2.4)
  if (!is.null(est$theta_cut) && is.finite(est$theta_cut)) {
    p <- p + ggplot2::geom_vline(xintercept = est$theta_cut,
                                 linetype = "dotted", colour = "grey20")
  }
  p +
    ggplot2::labs(
      title = "Local bias from keeping the DIF items",
      x = "Latent value (logits)", y = "b (logits)", caption = caption
    ) +
    ggplot2::theme_bw() +
    er2_axis_margins() +
    er2_plot_caption()
}
