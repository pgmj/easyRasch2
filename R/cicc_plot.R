#' Conditional item characteristic curves
#'
#' Plots one panel per item with the model-expected item score given the
#' total score, \eqn{E[X_i \mid R = r]}, and the observed mean item score in
#' groups of total scores, with confidence intervals. With a `dif_var`, the
#' observed means are shown per group, which makes each panel a visual DIF
#' check, and the partial gamma DIF coefficient is shown per item.
#'
#' The approach follows Buchardt, Christensen and Jensen (2023) and their
#' `RASCHplot` package, implemented on easyRasch2's CML engine: CML item
#' thresholds from `psychotools`, and the exact distribution of each item
#' score given the total score from elementary symmetric functions.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. One column per item, one row per person.
#' @param dif_var Optional vector of length `nrow(data)` (a factor, or
#'   coercible to one) defining groups for DIF. The observed means are then
#'   shown per group, and each panel shows the item's partial gamma from
#'   `iarm::partgam_DIF()`, as in [RMdifGamma()]. Default `NULL`.
#' @param grouping How the observed means are grouped along the total score
#'   (see Class intervals in Details): `"quantile"` (default),
#'   `class_intervals` groups with about equal numbers of respondents,
#'   `"width"`, `class_intervals` intervals of equal width over the
#'   observed total scores, `"score"`, every observed total score its own
#'   group, or `"manual"`, groups set by `score_breaks`. `"cut"` is an old
#'   name for `"quantile"`, still accepted.
#' @param class_intervals Integer, at least 2. Number of groups for
#'   `grouping = "quantile"` and `"width"`. Default `4`. Ignored by
#'   `"score"` and `"manual"`.
#' @param score_breaks Integer vector for `grouping = "manual"`: the total
#'   scores at which a new group starts, in increasing order. For example,
#'   `score_breaks = c(3, 6, 8)` on a 0 to 9 scale gives the groups 0-2,
#'   3-5, 6-7 and 8-9. It corresponds to `lower.groups` in RASCHplot's
#'   `CICCplot()`, without the leading zero. Default `NULL`.
#' @param ci Logical. Draw confidence intervals on the observed group means.
#'   Default `TRUE`.
#' @param error_band Logical. Draw a band around the expected curve showing
#'   where the observed mean at each total score should fall under the
#'   model (see Details). It complements `ci` and does not replace it.
#'   Default `FALSE`.
#' @param conf_level Numeric in (0, 1). Confidence level of the error bars
#'   and the band. Default `0.95`.
#' @param min_n Integer. A cell (group by score interval) needs at least
#'   this many respondents to be shown with its interval. Default `8`.
#' @param items Optional character or integer vector of the items to plot.
#'   The model is always fitted on all items. Default `NULL`, all items.
#' @param output `"patchwork"` (default) for one combined figure, or
#'   `"list"` for a named list of one `ggplot` per item.
#' @param method Old name for `grouping`, still accepted. Default `NULL`.
#'
#' @return A `patchwork` figure, or with `output = "list"` a named list of
#'   `ggplot` objects. With `dif_var`, the partial gamma table is also
#'   attached as `attr(., "dif_gamma")`, for convenience. It is the result
#'   of `RMdifGamma(data, dif_var, output = "dataframe")` on the same
#'   respondents, so call [RMdifGamma()] to get it as a table, with p-values
#'   corrected across items and, with a cutoff, a calibrated test.
#'
#' @details
#' \strong{Data.} Only complete cases are used, and with `dif_var` also only
#' rows where it is not missing. A message says how many rows were dropped.
#' An unused response category stops the
#' function and names the item. A category used only by respondents with
#' a zero or perfect total score gives a warning, since its threshold
#' cannot be estimated well.
#'
#' \strong{Conditioning.} The expected curve is the conditional expectation
#' of the item score given the total score (Andersen, 1995), which accounts
#' for the item being part of the total. It is not a marginal ICC.
#'
#' \strong{Class intervals.} Following Buchardt et al. (2023), the observed
#' item means can be shown for each total score or for groups of total
#' scores. All grouping is on the total score, never on an estimated person
#' location. `grouping = "quantile"` splits the observed total scores into
#' `class_intervals` groups of about equal numbers of respondents, with
#' common boundaries so that DIF groups share the x-axis. This is the
#' approach of `iarm::ICCplot()` (`Hmisc::cut2(g = ...)`). `grouping =
#' "width"` splits the observed score range into intervals of equal width,
#' which reads naturally on the score axis but can leave sparse intervals,
#' which `min_n` removes. `grouping = "score"` shows every observed total
#' score, the finest resolution, sensible for short scales. `grouping =
#' "manual"` places the boundaries where `score_breaks` says. Each group
#' gives one observed point at its mean total score. If the scores are too
#' sparse to form `class_intervals` distinct groups, `"quantile"` and
#' `"width"` fall back to one point per score.
#'
#' The caption names the grouping that was used, which is not always the
#' one requested. Quantile grouping forms fewer groups than asked for when
#' total scores tie at the boundaries, equal-width and manual grouping can
#' define intervals that no respondent falls into, and both binned methods
#' fall back to score level when the scores are too sparse. The caption
#' reports the number of groups formed, how many contain respondents when
#' that is fewer, and any fall back.
#'
#' \strong{Confidence intervals.} The error bars use the normal
#' approximation \eqn{\bar{x}_l \pm z \sqrt{\mathrm{var}(x_l) / n_l}}
#' within each cell, limited to the item's score range. Cells with fewer
#' than `min_n` respondents are dropped. With `dif_var`, this shows sparse
#' group differences as uncertain.
#'
#' \strong{The model band (`error_band`).} The band is
#' \eqn{E_{ri} \pm z \sqrt{V_{ri}/n_r}}, where \eqn{E_{ri}} and
#' \eqn{V_{ri}} are the model's conditional mean and variance of the item
#' score at total score \eqn{r}, and \eqn{n_r} is the observed number of
#' respondents with that total score. If the model holds, this is where
#' the observed mean item score at that total score should fall, given how
#' many respondents have it. An observed point outside the band shows
#' local misfit. The band uses the standard error of the mean, not the
#' much wider SD of single responses, which would flag nothing. Its width
#' therefore follows the data. It is narrow where many respondents are
#' (a strict check) and wide where few are (a lenient one), with gaps where
#' \eqn{n_r = 0}, and it closes on the curve at total score 0 and at the
#' maximum, where the item score is fixed. The band shows the model's
#' uncertainty of the observed mean per total score, and the error bars
#' the sampling uncertainty of the observed group means.
#'
#' \strong{DIF magnitude.} The partial gamma (Bjorner et al., 1998) is the
#' association between the item score and the group given the rest score,
#' a non-parametric effect size, the same statistic as in [RMdifGamma()].
#'
#' @references
#' Andersen, E. B. (1995). Polytomous Rasch models and their estimation. In
#' G. H. Fischer & I. W. Molenaar (Eds.), *Rasch Models: Foundations, Recent
#' Developments, and Applications* (pp. 271-291). Springer. (Conditional
#' expectation of item scores given the total score, formula 15.22.)
#'
#' Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
#' (1998). Differential item functioning in the Danish translation of the
#' SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189-1202.
#'
#' Buchardt, A.-S., Christensen, K. B., & Jensen, S. N. (2023). Visualizing
#' Rasch item fit using conditional item characteristic curves in R.
#' *Psychological Test and Assessment Modeling, 65*(2), 206-219.
#'
#' @seealso [RMitemInfit()] and [RMitemRestscore()] for item fit tests,
#'   [RMdifGamma()] and [RMdifAnchor()] for DIF tests,
#'   [RMitemCatFitPlot()] for category-level plots given the total score.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE) &&
#'     requireNamespace("iarm", quietly = TRUE)) {
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:2, 200 * 4, replace = TRUE), nrow = 200, ncol = 4)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:4)
#'
#'   # Default: quantile groups (~equal numbers of respondents per group)
#'   RMitemICCPlot(sim_data)
#'
#'   # Equal-width intervals on the total-score scale, with the model band
#'   RMitemICCPlot(sim_data, grouping = "width", error_band = TRUE)
#'
#'   # Every total score as its own group (short scales)
#'   RMitemICCPlot(sim_data, grouping = "score")
#'
#'   # Manual grouping: new groups start at total scores 3, 5, and 7
#'   RMitemICCPlot(sim_data, grouping = "manual", score_breaks = c(3, 5, 7))
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
RMitemICCPlot <- function(data,
                          dif_var         = NULL,
                          grouping        = c("quantile", "width", "score",
                                              "manual"),
                          class_intervals = 4,
                          score_breaks    = NULL,
                          ci              = TRUE,
                          error_band      = FALSE,
                          conf_level      = 0.95,
                          min_n           = 8L,
                          items           = NULL,
                          output          = c("patchwork", "list"),
                          method          = NULL) {

  # `method` is the pre-1.5.0 name of `grouping`, kept as a silent alias so
  # existing code (including the jamovi module) keeps working.
  if (!is.null(method)) {
    if (!missing(grouping)) {
      stop("Give `grouping` only. `method` is its old name.", call. = FALSE)
    }
    grouping <- method
  }
  g <- .resolve_grouping(grouping, class_intervals, score_breaks)
  method <- g$grouping
  class_intervals <- g$class_intervals
  score_breaks <- g$score_breaks
  output <- match.arg(output)

  validate_response_data(data) # also checks for 2+ items
  data <- as.data.frame(data)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(min_n) || length(min_n) != 1L || min_n < 1) {
    stop("`min_n` must be a single positive integer.", call. = FALSE)
  }
  min_n <- as.integer(min_n)

  # --- DIF variable + joint NA handling --------------------------------------
  dif_mode   <- !is.null(dif_var)
  dif_factor <- NULL
  if (dif_mode) {
    if (length(dif_var) != nrow(data)) {
      stop("`dif_var` must have the same length as `nrow(data)` (",
           nrow(data), "). Got ", length(dif_var), ".", call. = FALSE)
    }
    dif_factor <- droplevels(as.factor(dif_var))
  }
  keep <- stats::complete.cases(data)
  if (dif_mode) keep <- keep & !is.na(dif_factor)
  if (!all(keep)) {
    message(
      sum(!keep), " respondent(s) with missing responses",
      if (dif_mode) " or missing `dif_var`",
      " dropped (RMitemICCPlot() uses complete cases)."
    )
  }
  data <- data[keep, , drop = FALSE]
  if (dif_mode) dif_factor <- droplevels(dif_factor[keep])
  if (nrow(data) == 0L) stop("No complete cases in `data`.", call. = FALSE)
  if (dif_mode && nlevels(dif_factor) < 2L) {
    stop("`dif_var` must have at least 2 distinct non-missing levels ",
         "after NA filtering.", call. = FALSE)
  }

  # --- Which items to render (model is always fitted on all) -----------------
  item_names <- names(data)
  show_items <- if (is.null(items)) item_names else {
    if (is.numeric(items)) {
      if (any(items < 1L | items > ncol(data))) {
        stop("`items` indices out of range.", call. = FALSE)
      }
      item_names[as.integer(items)]
    } else {
      unknown <- setdiff(items, item_names)
      if (length(unknown) > 0L) {
        stop("Unknown item(s) in `items`: ", paste(unknown, collapse = ", "),
             call. = FALSE)
      }
      items
    }
  }

  # --- Plotting packages (checked after input validation) --------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for RMitemICCPlot(). ",
         "Install with: install.packages(\"ggplot2\")", call. = FALSE)
  }
  if (output == "patchwork" && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for output = \"patchwork\". ",
         "Install with: install.packages(\"patchwork\")", call. = FALSE)
  }

  # --- Model curve + conditional variance (shared CML engine) ----------------
  data_mat <- as.matrix(data)
  thr_list <- .fit_cml_thresholds(data_mat)
  .sparsity_warning(data) # after the fit: category-check stops come first
  mom      <- .cond_moments(thr_list)              # $E, $V by total score; $max_score
  max_sc   <- mom$max_score
  steps    <- vapply(thr_list, length, integer(1L))   # per-item max score
  rs         <- rowSums(data_mat)
  n_by_score <- tabulate(rs + 1L, nbins = max_sc + 1L)  # respondents per total score
  z          <- stats::qnorm(1 - (1 - conf_level) / 2)

  # --- Common class-interval assignment --------------------------------------
  if (method == "manual" && any(score_breaks > max_sc)) {
    stop(
      "`score_breaks` contains values above the maximum possible total ",
      "score (", max_sc, ").",
      call. = FALSE
    )
  }
  bin <- .cicc_class_bins(rs, method, class_intervals, score_breaks)
  grp <- if (dif_mode) dif_factor else factor(rep("all", nrow(data_mat)))

  # Dodge width for per-group points/bars: scales with the number of DIF groups,
  # but capped at the gap between adjacent class-interval centres so groups stay
  # within their interval (more groups -> wider spread; 0 when not in DIF mode).
  n_grp       <- nlevels(grp)
  bin_centres <- sort(unique(as.numeric(tapply(rs, bin, mean))))
  bin_spacing <- if (length(bin_centres) > 1L) stats::median(diff(bin_centres)) else 1
  dodge_width <- if (dif_mode && n_grp > 1L) min(n_grp * 1, 0.9 * bin_spacing) else 0

  # --- Partial-gamma DIF magnitude (per item; engine-invariant via iarm) -----
  gamma_tab <- NULL
  if (dif_mode) {
    gamma_tab <- tryCatch(
      RMdifGamma(data, dif_var = dif_factor, output = "dataframe"),
      error = function(e) NULL
    )
  }

  # --- Build per-item panels -------------------------------------------------
  panels <- lapply(show_items, function(nm) {
    i     <- match(nm, item_names)
    curve <- .cicc_curve_df(mom, i, steps[i], z, n_by_score)
    obs   <- .cicc_obs_df(data_mat[, i], rs, grp, bin, min_n, z, steps[i])
    gl    <- .cicc_gamma_label(gamma_tab, nm)
    .cicc_panel(curve, obs, nm, steps[i], error_band, ci, dif_mode, gl, dodge_width)
  })
  names(panels) <- show_items

  if (output == "list") {
    if (!is.null(gamma_tab)) attr(panels, "dif_gamma") <- gamma_tab
    return(panels)
  }

  out <- patchwork::wrap_plots(panels, axes = "collect", guides = "collect") +
    patchwork::plot_annotation(
      title   = "Conditional Item Characteristic Curves",
      caption = er2_caption(.cicc_caption(
        dif_mode,
        conf_level,
        .cicc_interval_caption(method, class_intervals, score_breaks, bin)
      )),
      theme   = er2_plot_caption()
    )
  if (!is.null(gamma_tab)) attr(out, "dif_gamma") <- gamma_tab
  out
}

# ===========================================================================
# Internal helpers
# ===========================================================================

#' Validate the class-interval arguments shared by RMitemICCPlot() and
#' RMitemCatFitPlot()
#' @return list(grouping, class_intervals, score_breaks), normalised.
#' @keywords internal
#' @noRd
.resolve_grouping <- function(grouping, class_intervals, score_breaks) {
  grouping <- match.arg(grouping[1L],
                        c("quantile", "width", "score", "manual", "cut"))
  # "cut" is the pre-1.1.0 name of the quantile grouping, kept as a silent
  # legacy alias.
  if (grouping == "cut") grouping <- "quantile"
  if (!is.null(score_breaks)) {
    if (grouping != "manual") {
      stop("`score_breaks` is only used with grouping = 'manual'.",
           call. = FALSE)
    }
    if (!is.numeric(score_breaks) || length(score_breaks) < 1L ||
        anyNA(score_breaks) || any(score_breaks != as.integer(score_breaks)) ||
        is.unsorted(score_breaks, strictly = TRUE) || any(score_breaks < 1)) {
      stop(
        "`score_breaks` must be a strictly increasing vector of positive ",
        "integers (each value starts a new total-score group).",
        call. = FALSE
      )
    }
    score_breaks <- as.integer(score_breaks)
  } else if (grouping == "manual") {
    stop(
      "grouping = 'manual' requires `score_breaks` (the total scores at ",
      "which a new group starts, e.g. c(3, 6, 8)).",
      call. = FALSE
    )
  }
  if (grouping %in% c("quantile", "width")) {
    if (!is.numeric(class_intervals) || length(class_intervals) != 1L ||
        !is.finite(class_intervals) ||
        class_intervals != round(class_intervals) || class_intervals < 2L) {
      stop("`class_intervals` must be a single integer >= 2.", call. = FALSE)
    }
    class_intervals <- as.integer(class_intervals)
  }
  list(grouping = grouping, class_intervals = class_intervals,
       score_breaks = score_breaks)
}

#' Assign respondents to total-score class intervals
#' @keywords internal
#' @noRd
.cicc_class_bins <- function(rs, method, class_intervals,
                             score_breaks = NULL) {
  if (method == "score") return(factor(rs))

  if (method == "manual") {
    # Each value in score_breaks starts a new group: breaks c(3, 6, 8) on a
    # 0-9 scale give the groups 0-2, 3-5, 6-7, 8-9. Groups left empty by the
    # data simply contribute no observed point (min_n prunes them anyway).
    lows   <- c(0L, score_breaks)
    highs  <- c(score_breaks - 1L, max(rs, score_breaks))
    labels <- ifelse(lows == highs, as.character(lows),
                     paste0(lows, "-", highs))
    return(cut(rs,
               breaks = c(-0.5, score_breaks - 0.5, max(rs, score_breaks) + 0.5),
               labels = labels, include.lowest = TRUE))
  }

  # The two binned methods fall back to score level when the data cannot
  # support the requested grouping. The fallback is flagged on the result so
  # the figure caption can report what was actually used rather than what was
  # asked for.
  fell_back <- function(x) {
    attr(x, "fallback") <- TRUE
    x
  }

  if (method == "width") {
    # Equal-width intervals over the observed total-score range.
    breaks <- unique(seq(min(rs) - 0.5, max(rs) + 0.5,
                         length.out = class_intervals + 1L))
    if (length(breaks) < 3L) {
      # Too few distinct total scores to form bins: fall back to score level.
      return(fell_back(factor(rs)))
    }
    return(cut(rs, breaks = breaks, include.lowest = TRUE))
  }

  # method == "quantile": equal-count bins over the observed total scores.
  qs <- unique(stats::quantile(
    rs, probs = seq(0, 1, length.out = class_intervals + 1L), na.rm = TRUE))
  if (length(qs) < 3L) {
    # Too few distinct total scores to form bins: fall back to score level.
    return(fell_back(factor(rs)))
  }
  cut(rs, breaks = qs, include.lowest = TRUE)
}

#' Caption clause naming the class-interval grouping actually used
#'
#' Reports the realised grouping rather than the requested one. `"quantile"`
#' can form fewer groups than asked for when total scores tie at the quantile
#' boundaries, and both binned methods fall back to one group per total score
#' when there are too few distinct scores to bin at all. Either way the plot
#' no longer shows what the arguments imply, so the caption says so.
#'
#' @param method Normalised method (`"cut"` already resolved to `"quantile"`).
#' @param class_intervals The requested number of groups.
#' @param score_breaks Breaks for `method = "manual"`, else `NULL`.
#' @param bin The factor returned by `.cicc_class_bins()`.
#' @return A single sentence.
#' @keywords internal
#' @noRd
.cicc_interval_caption <- function(method, class_intervals, score_breaks, bin) {
  n_groups <- nlevels(bin)
  # Equal-width and manual grouping can define intervals that no respondent
  # falls into, and those contribute no marker, so the caption reports both
  # counts rather than implying more markers than the figure shows.
  n_used <- sum(table(bin) > 0L)
  empty_clause <- if (n_used < n_groups) {
    sprintf(" %d of them contain respondents.", n_used)
  } else {
    ""
  }

  if (isTRUE(attr(bin, "fallback"))) {
    return(sprintf(paste0(
      "Class intervals: %s grouping was requested, but there are too few ",
      "distinct total scores to form it, so every observed total score is ",
      "its own group (%d groups)."), method, n_groups))
  }
  body <- switch(
    method,
    score = sprintf(
      "Class intervals: every observed total score is its own group (%d groups).",
      n_groups),
    manual = sprintf(
      "Class intervals: manual, a new group starting at total score%s %s (%d groups).",
      if (length(score_breaks) > 1L) "s" else "",
      paste(score_breaks, collapse = ", "), n_groups),
    width = sprintf(
      "Class intervals: %d equal-width groups over the observed total-score range.",
      n_groups),
    quantile = if (n_groups < class_intervals) {
      sprintf(paste0(
        "Class intervals: %d groups of approximately equal size (%d were ",
        "requested, fewer were formed because total scores tie at the ",
        "quantile boundaries)."), n_groups, class_intervals)
    } else {
      sprintf(
        "Class intervals: %d groups of approximately equal size (quantile).",
        n_groups)
    }
  )
  paste0(body, empty_clause)
}

#' Expected conditional item score (and SE-of-mean band) per total score
#'
#' The optional band is the model 95\% interval for the *observed mean* at each
#' total score, \eqn{E \pm z \sqrt{\mathrm{Var}/n_r}}, with `n_r` the number of
#' respondents at that total score (as in RASCHplot's `CICCplot`) -- not the
#' much wider individual-response SD. It is `NA` where no respondent has that
#' score and collapses to the point estimate at the degenerate extremes
#' (`Var = 0`).
#' @keywords internal
#' @noRd
.cicc_curve_df <- function(mom, i, max_y, z, n_by_score) {
  max_sc   <- mom$max_score
  expected <- numeric(max_sc + 1L)
  v        <- numeric(max_sc + 1L)
  expected[1L] <- 0; expected[max_sc + 1L] <- max_y      # degenerate extremes
  v[1L] <- 0; v[max_sc + 1L] <- 0
  if (max_sc >= 2L) {
    interior <- 2L:max_sc                                # total scores 1..max-1
    expected[interior] <- mom$E[interior, i]
    v[interior]        <- mom$V[interior, i]
  }
  # Standard error of the conditional *mean* at each total score (variance
  # divided by the count at that score), not the individual-response SD.
  se <- sqrt(pmax(v, 0) / n_by_score)
  se[!is.finite(se)] <- NA_real_                         # no respondents at score
  se[v == 0]         <- 0                                # deterministic extremes
  data.frame(
    score    = 0:max_sc,
    expected = expected,
    band_lo  = pmax(0,     expected - z * se),
    band_hi  = pmin(max_y, expected + z * se),
    stringsAsFactors = FALSE
  )
}

#' Observed class-interval means + CIs for one item (per group)
#' @keywords internal
#' @noRd
.cicc_obs_df <- function(y, rs, grp, bin, min_n, z, max_y) {
  # Common x per bin (shared across groups) so per-group points/error bars at
  # the same interval can be separated cleanly with position_dodge().
  bin_x <- tapply(rs, bin, mean)
  key   <- interaction(grp, bin, drop = TRUE)
  parts <- split(seq_along(y), key)
  rows  <- lapply(parts, function(idx) {
    if (length(idx) < min_n) return(NULL)
    yy <- y[idx]
    m  <- mean(yy)
    se <- sqrt(stats::var(yy) / length(idx))
    data.frame(
      grp  = grp[idx][1L],
      x    = as.numeric(bin_x[as.character(bin[idx][1L])]),
      mean = m,
      lo   = max(0,     m - z * se),
      hi   = min(max_y, m + z * se),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

#' Per-item partial-gamma DIF annotation string (or NULL)
#' @keywords internal
#' @noRd
.cicc_gamma_label <- function(gamma_tab, nm) {
  if (is.null(gamma_tab)) return(NULL)
  row <- gamma_tab[gamma_tab$Item == nm, , drop = FALSE]
  if (nrow(row) == 0L) return(NULL)
  row <- row[1L, ]
  # Short label for the in-plot annotation; `p` is the BH-adjusted DIF p-value.
  # \u03b3 is the Greek small letter gamma (escaped to keep R source ASCII).
  sprintf("\u03b3 = %.2f [%.2f, %.2f]\np = %.3f",
          row$gamma, row$lower, row$upper, row$padj_bh)
}

#' Composite caption explaining the markers and the DIF annotation
#' @keywords internal
#' @noRd
.cicc_caption <- function(dif_mode, conf_level, interval_clause = NULL) {
  pct <- round(conf_level * 100)
  g   <- intToUtf8(947L)   # Greek small letter gamma (keeps R source ASCII)
  cap <- sprintf(paste0(
    "Black line: model-expected conditional item score. Diamonds: observed mean ",
    "item score per total-score group with %d%% confidence-interval error ",
    "bars."), pct)
  if (!is.null(interval_clause)) {
    cap <- paste(cap, interval_clause)
  }
  if (dif_mode) {
    cap <- paste0(cap, sprintf(paste0(
      " %s is Goodman-Kruskal's gamma (partial gamma DIF coefficient). The ",
      "bracketed values are its %d%% confidence interval and p is the ",
      "Benjamini-Hochberg adjusted p-value."), g, pct))
  }
  cap
}

#' Build one cICC ggplot panel
#' @keywords internal
#' @noRd
.cicc_panel <- function(curve, obs, item_name, max_y, error_band, ci, dif_mode,
                        gamma_label, dodge_width) {
  curve$item <- item_name                     # single facet -> grey strip label
  dodge      <- ggplot2::position_dodge(width = dodge_width)

  p <- ggplot2::ggplot()
  if (error_band) {
    p <- p + ggplot2::geom_ribbon(
      data = curve,
      ggplot2::aes(x = .data$score, ymin = .data$band_lo, ymax = .data$band_hi),
      fill = "grey85", alpha = 0.6
    )
  }
  p <- p + ggplot2::geom_line(
    data = curve, ggplot2::aes(x = .data$score, y = .data$expected),
    colour = "black", linewidth = 0.8
  )
  if (!is.null(obs)) {
    obs$item <- item_name
    if (dif_mode) {
      if (ci) {
        p <- p + ggplot2::geom_errorbar(
          data = obs,
          ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi,
                       colour = .data$grp),
          width = 0, linewidth = 0.75, position = dodge
        )
      }
      p <- p + ggplot2::geom_point(
        data = obs,
        ggplot2::aes(x = .data$x, y = .data$mean, colour = .data$grp),
        shape = 18, size = 2.9, position = dodge
      )
    } else {
      if (ci) {
        p <- p + ggplot2::geom_errorbar(
          data = obs, ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
          width = 0, linewidth = 0.75, colour = "sienna2"
        )
      }
      p <- p + ggplot2::geom_point(
        data = obs, ggplot2::aes(x = .data$x, y = .data$mean),
        shape = 18, size = 2.9, colour = "sienna2"
      )
    }
  }
  # DIF magnitude annotated top-left inside the panel (like iarm's cICC).
  if (!is.null(gamma_label)) {
    p <- p + ggplot2::annotate(
      "text", x = 0, y = max_y, label = gamma_label,
      hjust = 0, vjust = 1, size = 3, colour = "grey20"
    )
  }
  xmax  <- max(curve$score)
  xstep <- if (xmax <= 15) 1L else if (xmax <= 40) 3L else 5L
  p +
    ggplot2::facet_wrap(~ item) +
    ggplot2::scale_x_continuous(breaks = seq(0L, xmax, by = xstep)) +
    ggplot2::labs(x = "Total score", y = "Item score", colour = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold")
    ) +
    er2_axis_margins()
}
