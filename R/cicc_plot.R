#' Conditional Item Characteristic Curves
#'
#' Faceted panel of Conditional Item Characteristic Curves (cICCs), one per
#' item. Each panel shows the model-expected conditional item score
#' \eqn{E[X_i \mid R = r]} against the total score \eqn{r}, together with the
#' average observed item score within class intervals of the total score and
#' their confidence intervals. When a `dif_var` is supplied, observed averages
#' are computed separately per group, turning each panel into a visual DIF
#' check, and the partial-gamma DIF magnitude is reported per item.
#'
#' The model curve and its conditional variance come from the shared CML
#' engine (CML item thresholds via `psychotools` and the exact conditional
#' distribution of each item score given the total score, via elementary
#' symmetric functions). The conditional-ICC approach follows Buchardt,
#' Christensen & Jensen (2023) and their `RASCHplot` package; the implementation
#' here is native to easyRasch2's CML/WLE engine.
#'
#' @param data A data.frame or matrix of item responses (non-negative
#'   integers, 0-based). One column per item, one row per person.
#' @param dif_var Optional vector of length `nrow(data)` (factor or coercible
#'   to factor) defining a DIF grouping variable. When supplied, observed
#'   averages are shown separately per group and the partial-gamma DIF
#'   coefficient (`iarm::partgam_DIF()`, as in [RMdifGamma()]) is annotated
#'   per panel. Default `NULL` (no DIF).
#' @param method One of `"cut"` (default) or `"score"`. `"cut"` groups
#'   respondents into `class_intervals` total-score bins (equal-count
#'   quantile bins, common across groups); `"score"` uses every observed
#'   total score as its own group.
#' @param class_intervals Integer >= 2. Number of class intervals when
#'   `method = "cut"`. Default `4`. Ignored when `method = "score"`.
#' @param ci Logical. Draw confidence intervals (error bars) on the observed
#'   class-interval means. Default `TRUE`.
#' @param error_band Logical. Add a shaded band for the model's interval of the
#'   *observed mean* at each total score, \eqn{E \pm z \sqrt{\mathrm{Var}/n_r}}
#'   (`n_r` = number of respondents at that total score), as in RASCHplot's
#'   `CICCplot`. Default `FALSE`.
#' @param conf_level Numeric in (0, 1). Confidence level for the observed
#'   error bars and the model band. Default `0.95`.
#' @param min_n Integer. A group-by-interval cell needs at least this many
#'   respondents to contribute an observed point + CI; sparser cells are
#'   skipped. Default `8` (the package's per-cell stability floor).
#' @param items Optional character or integer vector selecting which items to
#'   plot. The model is always fitted on all items; only the rendering is
#'   filtered. Default `NULL` (all items).
#' @param output One of `"patchwork"` (default) -- composite patchwork figure
#'   -- or `"list"` -- a named list of per-item `ggplot` objects.
#'
#' @return Either a `patchwork`/`ggplot` composite (default) or a named list
#'   of per-item `ggplot` objects (`output = "list"`). In DIF mode the per-item
#'   partial-gamma table is attached as `attr(., "dif_gamma")`.
#'
#' @details
#' \strong{Conditioning.} The expected curve is the exact conditional
#' expectation of the item score given the total score (it accounts for the
#' item being part of the total), not a marginal ICC.
#'
#' \strong{Class intervals.} With `method = "cut"`, total scores are split into
#' `class_intervals` equal-count bins using common boundaries (so groups share
#' the x-axis); each bin contributes one observed point at its mean total
#' score. With `method = "score"`, every observed total score is a point. If
#' the score distribution is too sparse to form `class_intervals` distinct
#' bins, the function falls back to score-level points.
#'
#' \strong{Confidence intervals.} Observed error bars use the normal
#' approximation \eqn{\bar{x}_l \pm z \sqrt{\mathrm{var}(x_l) / n_l}} within
#' each (group, interval) cell, clamped to the item's score range; cells with
#' fewer than `min_n` respondents are dropped. In DIF mode this makes sparse
#' group differences visibly uncertain rather than over-interpreted.
#'
#' \strong{DIF magnitude.} The annotated partial gamma (Bjorner et al., 1998)
#' is the association between the item score and the group conditional on the
#' rest score -- a non-parametric effect size, complementary to the
#' total-score-conditional visual. It is the same statistic as [RMdifGamma()].
#'
#' @references
#' Buchardt, A.-S., Christensen, K. B., & Jensen, S. N. (2023). Visualizing
#' Rasch item fit using conditional item characteristic curves in R.
#' *Psychological Test and Assessment Modeling, 65*(2), 206-219.
#'
#' Andersen, E. B. (1995). Polytomous Rasch models and their estimation. In
#' G. H. Fischer & I. W. Molenaar (Eds.), *Rasch Models: Foundations, Recent
#' Developments, and Applications* (pp. 271-291). Springer. (Conditional
#' expectation of item scores given the total score; formula 15.22.)
#'
#' Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
#' (1998). Differential item functioning in the Danish translation of the
#' SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189-1202.
#'
#' @seealso \code{\link{RMdifGamma}}, \code{\link{RMitemRestscore}}
#'
#' @importFrom rlang .data
#' @export
RMitemICCPlot <- function(data,
                          dif_var         = NULL,
                          method          = c("cut", "score"),
                          class_intervals = 4,
                          ci              = TRUE,
                          error_band      = FALSE,
                          conf_level      = 0.95,
                          min_n           = 8L,
                          items           = NULL,
                          output          = c("patchwork", "list")) {

  method <- match.arg(method)
  output <- match.arg(output)

  validate_response_data(data)
  data <- as.data.frame(data)
  if (ncol(data) < 2L) {
    stop("RMitemICCPlot() requires at least 2 items.", call. = FALSE)
  }
  if (!is.numeric(class_intervals) || length(class_intervals) != 1L ||
      !is.finite(class_intervals) || class_intervals < 2L) {
    stop("`class_intervals` must be a single integer >= 2.", call. = FALSE)
  }
  class_intervals <- as.integer(class_intervals)
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
  .sparsity_warning(data, is_poly = max(data_mat, na.rm = TRUE) > 1L)
  thr_list <- .fit_cml_thresholds(data_mat)
  mom      <- .cond_moments(thr_list)              # $E, $V by total score; $max_score
  max_sc   <- mom$max_score
  steps    <- vapply(thr_list, length, integer(1L))   # per-item max score
  rs         <- rowSums(data_mat)
  n_by_score <- tabulate(rs + 1L, nbins = max_sc + 1L)  # respondents per total score
  z          <- stats::qnorm(1 - (1 - conf_level) / 2)

  # --- Common class-interval assignment --------------------------------------
  bin <- .cicc_class_bins(rs, method, class_intervals)
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
      caption = er2_caption(.cicc_caption(dif_mode, conf_level)),
      theme   = er2_plot_caption()
    )
  if (!is.null(gamma_tab)) attr(out, "dif_gamma") <- gamma_tab
  out
}

# ===========================================================================
# Internal helpers
# ===========================================================================

#' Assign respondents to total-score class intervals
#' @keywords internal
#' @noRd
.cicc_class_bins <- function(rs, method, class_intervals) {
  if (method == "score") return(factor(rs))
  qs <- unique(stats::quantile(
    rs, probs = seq(0, 1, length.out = class_intervals + 1L), na.rm = TRUE))
  if (length(qs) < 3L) {
    # Too few distinct total scores to form bins: fall back to score level.
    return(factor(rs))
  }
  cut(rs, breaks = qs, include.lowest = TRUE)
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
.cicc_caption <- function(dif_mode, conf_level) {
  pct <- round(conf_level * 100)
  g   <- intToUtf8(947L)   # Greek small letter gamma (keeps R source ASCII)
  cap <- sprintf(paste0(
    "Black line: model-expected conditional item score. Diamonds: observed mean ",
    "item score per total-score group with %d%% confidence-interval error ",
    "bars."), pct)
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
          width = 0, linewidth = 0.6, position = dodge
        )
      }
      p <- p + ggplot2::geom_point(
        data = obs,
        ggplot2::aes(x = .data$x, y = .data$mean, colour = .data$grp),
        shape = 18, size = 2.6, position = dodge
      )
    } else {
      if (ci) {
        p <- p + ggplot2::geom_errorbar(
          data = obs, ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
          width = 0, linewidth = 0.6, colour = "sienna2"
        )
      }
      p <- p + ggplot2::geom_point(
        data = obs, ggplot2::aes(x = .data$x, y = .data$mean),
        shape = 18, size = 2.6, colour = "sienna2"
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
  xstep <- if (xmax <= 15) 1L else if (xmax <= 40) 2L else 5L
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
