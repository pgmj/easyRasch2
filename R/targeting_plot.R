#' Person-Item Targeting Plot (Wright Map)
#'
#' Produces a three-panel targeting plot with a shared logit scale x-axis:
#' \enumerate{
#'   \item \strong{Top}: Histogram of person location estimates, with a
#'     reference line for the mean (or median) and shading for ±1 SD
#'     (or ±1 MAD).
#'   \item \strong{Middle}: Inverted histogram of item threshold locations,
#'     with the same summary annotations.
#'   \item \strong{Bottom}: one bar per item, either partitioned into
#'     response-category bands (`panel = "categories"`, the default) or drawn
#'     as a dot-and-whisker plot of the individual thresholds
#'     (`panel = "thresholds"`).
#' }
#'
#' Together, the top and middle panels form a back-to-back histogram that
#' makes it easy to assess whether the test is well-targeted to the sample.
#' The bottom panel places the items on the same scale, so the category bands
#' show which response is the most likely one at the locations where the
#' persons actually sit.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed.
#' @param panel Character string selecting the bottom panel. `"categories"`
#'   (the default) draws each item as a bar partitioned into response-category
#'   bands, with the threshold estimates and their confidence intervals below
#'   it. `"thresholds"` draws the dot-and-whisker plot of item thresholds that
#'   was the only option before version 1.3.0.
#' @param robust Logical. If `FALSE` (the default), histogram annotations use
#'   mean ± SD. If `TRUE`, median ± MAD is used instead.
#' @param sort_items Character string controlling item ordering on the y-axis
#'   of the bottom panel. `"data"` (the default) preserves the column order in
#'   `data` (first item at top). `"location"` sorts items by their average
#'   threshold location (easiest at top, hardest at bottom).
#' @param bins Integer. Number of bins for both histograms. Default is number of
#'   unique scores divided by 2 (rounded up), but no less than 11.
#' @param xlim Numeric vector of length 2. Initial lower and upper limits for
#'   the shared x-axis. Automatically expanded if any person or item threshold
#'   values fall outside these limits.
#' @param ci_level Numeric. Confidence level for the item threshold error bars.
#'   Default is `0.95` (95% CI). Set to `NULL` to hide error bars.
#' @param category_labels Optional character vector of labels for the response
#'   categories, in ascending order and one per category. Used for the legend
#'   of the `"categories"` panel. Default `NULL` uses the category scores.
#' @param person_fill Fill colour for the person histogram. Default
#'   `"#0072B2"` (blue).
#' @param threshold_fill Fill colour for the item threshold histogram, and for
#'   the dot-and-whisker panel. Default `"#D55E00"` (vermillion).
#' @param viridis_option Character. Viridis palette option for the category
#'   bands. Default `"G"` (mako).
#' @param viridis_begin,viridis_end Numeric in \eqn{[0, 1]}. Start and end
#'   points of the viridis palette for the category bands. Defaults `0.9` and
#'   `0.2`, which runs the palette from light to dark so that higher
#'   categories are darker.
#' @param row_gap Numeric. Vertical spacing between item rows in the
#'   `"categories"` panel. Default `NULL` uses `1`, widened to `1.18` when at
#'   least one category collapses, so that its label has room above the bar.
#' @param height_ratios Numeric vector of length 3 specifying the relative
#'   heights of the top (person), middle (threshold), and bottom (dot-whisker)
#'   panels. Default `c(3, 2, 5)`.
#' @param output Character string. `"patchwork"` (the default) returns the
#'   combined patchwork plot. `"list"` returns a named list of the three
#'   ggplot objects (`p1`, `p2`, `p3`) for further customisation.
#'
#' @return
#' * If `output = "patchwork"`: a `patchwork` object (combined `ggplot`).
#' * If `output = "list"`: a named list with elements `p1` (person histogram),
#'   `p2` (threshold histogram), and `p3` (the bottom panel selected by
#'   `panel`).
#'
#' @details
#' \strong{Estimation method selection.}
#' The function checks whether any item response category has fewer than 3
#' observations. If all categories have at least 3 responses, item threshold
#' locations and their standard errors are estimated via Conditional Maximum
#' Likelihood (CML) using `psychotools::pcmodel()` (a dichotomous item is a
#' 2-category PCM). If any category has fewer than 3 responses, the function
#' falls back to Marginal Maximum Likelihood (MML) estimation via
#' `mirt::mirt()` with `itemtype = "Rasch"` and `SE = TRUE`, which is more
#' numerically stable under sparse-category conditions. A message is emitted
#' when the MML fallback is used.
#'
#' In both cases, item threshold locations are centered (shifted so the grand
#' mean of all thresholds equals zero).
#'
#' \strong{Person estimates} are obtained by Warm's weighted likelihood (WLE)
#' from the fitted item thresholds, consistent with the rest of the package.
#' WLE is finite at extreme scores, so all-zero and perfect responders are
#' located rather than dropped.
#'
#' \strong{Confidence intervals} for item thresholds are based on Wald-type
#' intervals: threshold estimate ± z × SE, where z is the standard normal
#' quantile corresponding to `ci_level`.
#'
#' \strong{Category bands.} With `panel = "categories"`, each band spans the
#' locations at which its response category is the most likely response. When
#' an item's thresholds are ordered these boundaries are the Andrich
#' thresholds themselves. When they are not, the disordered run is pooled by
#' averaging and the categories it skips over, which are never the most likely
#' response at any location, collapse to a red tick labelled with the category
#' number. Red arrows below the bar give the size of each threshold reversal in
#' logits. Ordered thresholds therefore leave no red marks at all.
#'
#' The two outer bands are open-ended and fade towards the panel edge, since
#' the lowest and highest categories have no outer boundary.
#'
#' The `ggplot2` and `patchwork` packages must be installed (they are in
#' Suggests, not Imports).
#'
#' @references
#' Wright, B. D. & Stone, M. H. (1979). *Best Test Design*. MESA Press.
#'
#' @seealso [psychotools::pcmodel()], [mirt::mirt()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   # Polytomous example
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:3, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:8)
#'
#'   # Default: category bands, mean/SD, data order, 95% CI
#'   RMtargeting(sim_data)
#'
#'   # Category bands with labels
#'   RMtargeting(sim_data, category_labels = c("Never", "Sometimes",
#'                                             "Often", "Always"))
#'
#'   # The dot-and-whisker panel
#'   RMtargeting(sim_data, panel = "thresholds")
#'
#'   # Robust (median/MAD), sorted by location, 84% CI
#'   RMtargeting(sim_data, robust = TRUE, sort_items = "location",
#'               ci_level = 0.84)
#'
#'   # Get list of sub-plots for customisation
#'   plots <- RMtargeting(sim_data, output = "list")
#'   plots$p1 + ggplot2::ggtitle("My custom title")
#'
#'   # Dichotomous example
#'   sim_bin <- as.data.frame(
#'     matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#'   )
#'   colnames(sim_bin) <- paste0("Item", 1:10)
#'   RMtargeting(sim_bin)
#' }
#' }
RMtargeting <- function(
  data,
  panel = c("categories", "thresholds"),
  robust = FALSE,
  sort_items = c("data", "location"),
  bins,
  xlim = c(-4, 4),
  ci_level = 0.95,
  category_labels = NULL,
  person_fill = "#0072B2",
  threshold_fill = "#D55E00",
  viridis_option = "G",
  viridis_begin = 0.9,
  viridis_end = 0.2,
  row_gap = NULL,
  height_ratios = c(3, 2, 5),
  output = "patchwork"
) {
  # --- Check required packages ------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for RMtargeting() but is not installed.\n",
      "Install it with: install.packages(\"ggplot2\")",
      call. = FALSE
    )
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required for RMtargeting() but is not installed.\n",
      "Install it with: install.packages(\"patchwork\")",
      call. = FALSE
    )
  }

  # if no manual value was selected for number of bins, use the number of
  # unique scores divided by 2, but no less than 11. ggplot2 requires a whole
  # number, so round up (an odd maximum score would otherwise give e.g. 13.5).
  if (missing(bins)) {
    bins <- ceiling(max(rowSums(data, na.rm = TRUE)) / 2)
    if (bins < 11) {
      bins <- 11
    }
  }

  panel <- match.arg(panel)
  sort_items <- match.arg(sort_items)
  output <- match.arg(output, c("patchwork", "list"))

  validate_response_data(data)

  data <- as.data.frame(data)
  n_total <- nrow(data)
  has_na <- anyNA(data)
  data <- .drop_empty_respondents(data)
  n_used <- nrow(data)

  data_mat <- as.matrix(data)
  max_score <- max(data_mat, na.rm = TRUE)
  is_dicho <- max_score == 1L
  item_names <- names(data)
  n_items <- ncol(data)

  # --- Check for sparse categories -------------------------------------------
  sparse <- .has_sparse_categories(data, min_n = 3L)

  # --- Estimate item thresholds -----------------------------------------------
  if (sparse) {
    message(
      "Some response categories have fewer than 3 observations. ",
      "Using MML estimation (mirt) for item thresholds."
    )
    thresh_info <- .estimate_thresholds_mml(data, is_dicho)
  } else {
    thresh_info <- .estimate_thresholds_cml(data, is_dicho)
  }

  # data.frame: Item, Threshold, Location, SE
  item_thresholds <- thresh_info$thresholds

  # --- Person estimates: WLE on the fitted item thresholds --------------------
  # Reconstruct the per-item threshold list from the (centred) locations and
  # estimate person locations by Warm's WLE, consistent with the rest of the
  # package and finite at extreme scores. Persons and items share a scale by
  # construction, so the targeting overlap is unaffected by the centring.
  thr_list <- lapply(
    split(
      item_thresholds$Location,
      factor(item_thresholds$Item, levels = names(data))
    ),
    as.numeric
  )
  person_theta <- .estimate_thetas(
    as.matrix(data),
    thr_list,
    method = "WLE"
  )$theta
  person_theta <- person_theta[is.finite(person_theta)]

  # --- Compute CI bounds for thresholds ---------------------------------------
  show_ci <- !is.null(ci_level) && all(!is.na(item_thresholds$SE))
  if (show_ci) {
    z_val <- stats::qnorm(1 - (1 - ci_level) / 2)
    item_thresholds$CI_low <- item_thresholds$Location -
      z_val * item_thresholds$SE
    item_thresholds$CI_high <- item_thresholds$Location +
      z_val * item_thresholds$SE
  }

  # --- Auto-expand xlim -------------------------------------------------------
  all_values <- c(person_theta, item_thresholds$Location)
  if (show_ci) {
    all_values <- c(all_values, item_thresholds$CI_low, item_thresholds$CI_high)
  }
  if (max(all_values, na.rm = TRUE) > xlim[2]) {
    xlim[2] <- ceiling(max(all_values, na.rm = TRUE))
  }
  if (min(all_values, na.rm = TRUE) < xlim[1]) {
    xlim[1] <- floor(min(all_values, na.rm = TRUE))
  }

  # --- Summary statistics -----------------------------------------------------
  if (robust) {
    p_center <- stats::median(person_theta, na.rm = TRUE)
    p_spread <- stats::mad(person_theta, na.rm = TRUE)
    t_center <- stats::median(item_thresholds$Location, na.rm = TRUE)
    t_spread <- stats::mad(item_thresholds$Location, na.rm = TRUE)
    center_label <- "Median"
    spread_label <- "MAD"
  } else {
    p_center <- mean(person_theta, na.rm = TRUE)
    p_spread <- stats::sd(person_theta, na.rm = TRUE)
    t_center <- mean(item_thresholds$Location, na.rm = TRUE)
    t_spread <- stats::sd(item_thresholds$Location, na.rm = TRUE)
    center_label <- "Mean"
    spread_label <- "SD"
  }

  # --- Item ordering for bottom panel -----------------------------------------
  if (sort_items == "location") {
    item_means <- stats::aggregate(
      Location ~ Item,
      data = item_thresholds,
      FUN = mean,
      na.rm = TRUE
    )
    # Easiest at top = lowest location at top.
    # ggplot places first factor level at bottom, so reverse.
    item_order <- rev(item_means$Item[order(item_means$Location)])
  } else {
    # Data order: first column at top -> reverse for ggplot
    item_order <- rev(item_names)
  }

  # ============================================================================
  # TOP PANEL: Person histogram
  # ============================================================================
  person_df <- data.frame(theta = person_theta)

  p1 <- ggplot2::ggplot(person_df, ggplot2::aes(x = .data$theta)) +
    ggplot2::geom_histogram(
      bins = bins,
      fill = person_fill,
      colour = "white",
      alpha = 0.85
    ) +
    ggplot2::annotate(
      "rect",
      xmin = p_center - p_spread,
      xmax = p_center + p_spread,
      ymin = -Inf,
      ymax = Inf,
      fill = person_fill,
      alpha = 0.12
    ) +
    ggplot2::geom_vline(
      xintercept = p_center,
      linewidth = 0.8,
      linetype = "dashed",
      colour = "grey20"
    ) +
    ggplot2::annotate(
      "text",
      x = p_center,
      y = Inf,
      vjust = -0.5,
      label = paste0(
        center_label,
        " = ",
        round(p_center, 2),
        ", ",
        spread_label,
        " = ",
        round(p_spread, 2)
      ),
      size = 3.2,
      colour = "grey20"
    ) +
    ggplot2::scale_y_continuous(
      breaks = function(lim) {
        seq(0, floor(lim[2]), by = max(1, round(lim[2] / 6)))
      }
    ) +
    ggplot2::scale_x_continuous(breaks = seq(xlim[1], xlim[2], by = 1)) +
    ggplot2::coord_cartesian(xlim = xlim, clip = "off") +
    ggplot2::labs(x = NULL, y = "Persons") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 5, 0, 5)
    ) +
    er2_axis_margins()

  # ============================================================================
  # MIDDLE PANEL: Inverted threshold histogram
  # ============================================================================
  thresh_hist_df <- data.frame(location = item_thresholds$Location)

  p2 <- ggplot2::ggplot(thresh_hist_df, ggplot2::aes(x = .data$location)) +
    ggplot2::geom_histogram(
      bins = bins,
      fill = threshold_fill,
      colour = "white",
      alpha = 0.85
    ) +
    ggplot2::annotate(
      "rect",
      xmin = t_center - t_spread,
      xmax = t_center + t_spread,
      ymin = -Inf,
      ymax = Inf,
      fill = threshold_fill,
      alpha = 0.12
    ) +
    ggplot2::geom_vline(
      xintercept = t_center,
      linewidth = 0.8,
      linetype = "dashed",
      colour = "grey20"
    ) +
    ggplot2::annotate(
      "text",
      x = t_center,
      y = -Inf,
      vjust = 1.5,
      label = paste0(
        center_label,
        " = ",
        round(t_center, 2),
        ", ",
        spread_label,
        " = ",
        round(t_spread, 2)
      ),
      size = 3.2,
      colour = "grey20"
    ) +
    ggplot2::scale_y_reverse(
      minor_breaks = NULL,
      breaks = function(lim) {
        max_val <- abs(floor(lim[1]))
        seq(0, max_val, by = max(1, round(max_val / 4)))
      }
    ) +
    ggplot2::scale_x_continuous(breaks = seq(xlim[1], xlim[2], by = 1)) +
    ggplot2::coord_cartesian(xlim = xlim, clip = "off") +
    ggplot2::labs(x = NULL, y = "Thresholds") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 5, 0, 5)
    ) +
    er2_axis_margins()

  # ============================================================================
  # BOTTOM PANEL: Item threshold dot-whisker plot
  # ============================================================================
  item_thresholds$Item <- factor(item_thresholds$Item, levels = item_order)

  # Build CI caption fragment
  ci_caption <- ""
  if (show_ci) {
    ci_caption <- paste0(
      "Error bars show ",
      round(ci_level * 100),
      "% confidence intervals."
    )
  }

  caption_body <- paste0(
    "Person location ",
    tolower(center_label),
    ": ",
    round(p_center, 2),
    " (",
    spread_label,
    " ",
    round(p_spread, 2),
    "). Item threshold location ",
    tolower(center_label),
    ": ",
    round(t_center, 2),
    " (",
    spread_label,
    " ",
    round(t_spread, 2),
    "). ",
    .n_caption(
      n_used,
      n_total,
      if (has_na) "incomplete responses retained" else character()
    ),
    "."
  )

  # The threshold panel keeps its own wording; the category panel describes
  # its extra marks instead, and says the same thing about the intervals.
  caption_text <- er2_caption(paste0(caption_body, "\n", ci_caption))

  if (panel == "categories") {
    # Integer threshold index, needed by the band helpers
    item_thresholds$k <- as.integer(sub("^T", "", item_thresholds$Threshold))
    n_cat <- max(item_thresholds$k) + 1L

    if (is.null(category_labels)) {
      category_labels <- as.character(seq_len(n_cat) - 1L)
    } else if (length(category_labels) != n_cat) {
      stop(
        "`category_labels` must have one label per response category.\n",
        "The data have ", n_cat, " categories and ",
        length(category_labels), " labels were supplied.",
        call. = FALSE
      )
    }

    fills <- scales::viridis_pal(
      option = viridis_option,
      begin = viridis_begin,
      end = viridis_end
    )(n_cat)

    bands_preview <- .category_bands(
      item_thresholds, as.character(item_order), xlim
    )
    # The collapsed-category label needs headroom above the bar, so open the
    # row spacing only when there is something to label.
    if (is.null(row_gap)) {
      row_gap <- if (is.null(bands_preview$drops)) 1 else 1.18
    }

    band_caption <- paste(
      "Bands span the locations at which each response category is the most",
      "likely response. Outer bands fade because they are open-ended."
    )
    if (show_ci) {
      band_caption <- paste(
        band_caption,
        sprintf(paste(
          "Points and intervals below each bar are the Andrich thresholds with",
          "%d%% confidence intervals, coloured by the category entered."
        ), round(ci_level * 100))
      )
    }
    if (!is.null(bands_preview$drops)) {
      band_caption <- paste(
        band_caption,
        "A red tick marks a category that is never the most likely response,",
        "labelled with its number."
      )
    }
    if (!is.null(.threshold_reversals(item_thresholds, as.character(item_order)))) {
      band_caption <- paste(
        band_caption,
        "Red arrows give the size in logits of each threshold reversal."
      )
    }

    p3 <- .targeting_category_panel(
      thr = item_thresholds,
      item_order = as.character(item_order),
      xlim = xlim,
      category_labels = category_labels,
      fills = fills,
      show_ci = show_ci,
      row_gap = row_gap,
      caption_text = er2_caption(paste(caption_body, band_caption))
    )
  } else if (is_dicho) {
    p3 <- ggplot2::ggplot(
      item_thresholds,
      ggplot2::aes(x = .data$Location, y = .data$Item)
    ) +
      ggplot2::geom_point(size = 3, colour = threshold_fill)

    if (show_ci) {
      p3 <- p3 +
        ggplot2::geom_errorbar(
          ggplot2::aes(xmin = .data$CI_low, xmax = .data$CI_high),
          width = 0.25,
          linewidth = 0.5,
          colour = threshold_fill,
          orientation = "y"
        )
    }

    p3 <- p3 +
      ggplot2::coord_cartesian(xlim = xlim) +
      ggplot2::labs(
        x = "Location (logit scale)",
        y = NULL,
        caption = caption_text
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 5, 5, 5)
      ) +
      er2_axis_margins() +
      er2_plot_caption()
  } else {
    p3 <- ggplot2::ggplot(
      item_thresholds,
      ggplot2::aes(
        x = .data$Location,
        y = .data$Item,
        colour = .data$Threshold
      )
    )

    if (show_ci) {
      p3 <- p3 +
        ggplot2::geom_errorbar(
          ggplot2::aes(xmin = .data$CI_low, xmax = .data$CI_high),
          width = 0.25,
          linewidth = 0.5,
          position = ggplot2::position_dodge(width = 0.4),
          orientation = "y"
        )
    }

    p3 <- p3 +
      ggplot2::geom_point(
        size = 2.5,
        position = ggplot2::position_dodge(width = 0.4)
      ) +
      ggplot2::scale_colour_viridis_d(end = 0.9) +
      ggplot2::scale_x_continuous(breaks = seq(xlim[1], xlim[2], by = 1)) +
      ggplot2::coord_cartesian(xlim = xlim) +
      ggplot2::labs(
        x = "Location (logit scale)",
        y = NULL,
        colour = "Threshold",
        caption = caption_text
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.margin = ggplot2::margin(0, 5, 5, 5)
      ) +
      er2_axis_margins() +
      er2_plot_caption()
  }

  # --- Return -----------------------------------------------------------------
  if (output == "list") {
    return(list(p1 = p1, p2 = p2, p3 = p3))
  }

  p1 / p2 / p3 + patchwork::plot_layout(heights = height_ratios)
}


# ── Internal helpers ─────────────────────────────────────────────────────────

#' Check for sparse response categories
#'
#' Returns `TRUE` if any item has a response category (within the item's
#' observed range 0..max) with fewer than `min_n` responses.
#'
#' @param data A data.frame of item responses.
#' @param min_n Minimum count threshold.
#' @return Logical scalar.
#' @noRd
.has_sparse_categories <- function(data, min_n = 3L) {
  for (j in seq_len(ncol(data))) {
    col_vals <- data[[j]]
    col_vals <- col_vals[!is.na(col_vals)]
    if (length(col_vals) == 0L) {
      next
    }
    max_cat <- max(col_vals)
    counts <- tabulate(col_vals + 1L, nbins = max_cat + 1L)
    if (any(counts < min_n)) return(TRUE)
  }
  FALSE
}


#' Estimate item thresholds via CML (eRm)
#'
#' Returns a list with `$thresholds` (data.frame: Item, Threshold, Location,
#' SE) and `$erm_out` (the fitted eRm model object, needed for person
#' parameter estimation).
#'
#' @param data A data.frame of item responses.
#' @param is_dicho Logical. `TRUE` for dichotomous data.
#' @return Named list.
#' @noRd
.estimate_thresholds_cml <- function(data, is_dicho) {
  # CML Andrich thresholds + SEs via psychotools (a dichotomous item is a
  # 2-category PCM), on the grand-mean-zero scale. Locations match the previous
  # eRm values; SEs differ slightly (psychotools vs eRm vcov). Person locations
  # are estimated downstream by WLE from these thresholds.
  fit <- psychotools::pcmodel(data)
  tp <- psychotools::threshpar(fit, vcov = TRUE)
  se_all <- sqrt(diag(attr(tp, "vcov")))
  loc <- unlist(lapply(tp, as.numeric), use.names = FALSE)
  loc <- loc - mean(loc, na.rm = TRUE)

  thresh_df <- data.frame(
    Item = rep(names(tp), lengths(tp)),
    Threshold = paste0("T", unlist(lapply(tp, seq_along), use.names = FALSE)),
    Location = loc,
    SE = as.numeric(se_all),
    stringsAsFactors = FALSE
  )
  rownames(thresh_df) <- NULL
  list(thresholds = thresh_df, erm_out = NULL)
}


#' Estimate item thresholds via MML (mirt) for sparse-category data
#'
#' Uses `mirt::mirt()` with `itemtype = "Rasch"` and `SE = TRUE` for threshold
#' estimation (more stable with sparse categories). Person locations are
#' estimated downstream by WLE from these thresholds.
#'
#' @param data A data.frame of item responses.
#' @param is_dicho Logical. `TRUE` for dichotomous data.
#' @return Named list with `$thresholds` and `$erm_out`.
#' @noRd
.estimate_thresholds_mml <- function(data, is_dicho) {
  mirt_out <- mirt::mirt(
    data,
    model = 1,
    itemtype = "Rasch",
    SE = TRUE,
    verbose = FALSE
  )

  # Extract per-item coefficients with SEs in IRT parameterisation
  item_coefs <- mirt::coef(
    mirt_out,
    IRTpars = TRUE,
    printSE = TRUE,
    simplify = FALSE
  )
  # item_coefs is a list: one element per item + "GroupPars"
  # Each item element is a matrix with rows: "par", "SE" and cols: "a", "b" (or "b1","b2",...)
  n_items <- ncol(data)

  thresh_list <- vector("list", n_items)

  for (i in seq_len(n_items)) {
    item_mat <- item_coefs[[i]]
    b_cols <- grep("^b\\d*$", colnames(item_mat), value = TRUE)
    if (length(b_cols) == 0L) {
      b_cols <- "b"
    }
    b_vals <- item_mat["par", b_cols]
    b_ses <- item_mat["SE", b_cols]

    thresh_list[[i]] <- data.frame(
      Item = names(data)[i],
      Threshold = paste0("T", seq_along(b_cols)),
      Location = as.numeric(b_vals),
      SE = as.numeric(b_ses),
      stringsAsFactors = FALSE
    )
  }

  thresh_df <- do.call(rbind, thresh_list)
  rownames(thresh_df) <- NULL

  # Center thresholds (SEs are invariant to location shift)
  grand_mean <- mean(thresh_df$Location, na.rm = TRUE)
  thresh_df$Location <- thresh_df$Location - grand_mean

  # Person locations are estimated downstream by WLE from these thresholds.
  list(thresholds = thresh_df, erm_out = NULL)
}


#' Modal-category boundaries for one item
#'
#' For a PCM item with Andrich thresholds \eqn{\tau_1 \dots \tau_m},
#' \eqn{\log P(X = k \mid \theta)} is \eqn{k\theta - S_k} up to a constant,
#' with \eqn{S_k = \sum_{j \le k} \tau_j} and \eqn{S_0 = 0}. Category `k` is
#' the most likely response where that line is the upper envelope of the
#' \eqn{m + 1} lines, which happens when the point \eqn{(k, S_k)} lies on the
#' lower convex hull of \eqn{\{(k, S_k)\}}. The boundary between two
#' consecutive hull vertices \eqn{k < l} falls at \eqn{(S_l - S_k)/(l - k)},
#' the mean of the thresholds they span.
#'
#' So when the thresholds are ordered the boundaries are the thresholds
#' themselves, and when they are not, a disordered run is pooled by averaging
#' and the categories it skips over are those that are never the most likely
#' response at any location. Closed form, no grid search.
#'
#' @param tau Numeric vector of Andrich thresholds for one item.
#' @return A list with `cats` (categories that are modal somewhere),
#'   `dropped` (categories that never are) and `bounds` (the boundary
#'   locations between consecutive elements of `cats`).
#' @noRd
.modal_boundaries <- function(tau) {
  s_cum <- c(0, cumsum(tau))
  k_idx <- seq_along(s_cum) - 1L
  hull <- 1L
  for (i in 2:length(s_cum)) {
    while (length(hull) >= 2L) {
      a <- hull[length(hull) - 1L]
      b <- hull[length(hull)]
      above <- (s_cum[b] - s_cum[a]) * (k_idx[i] - k_idx[a]) >=
        (s_cum[i] - s_cum[a]) * (k_idx[b] - k_idx[a])
      if (above) hull <- hull[-length(hull)] else break
    }
    hull <- c(hull, i)
  }
  vert <- k_idx[hull]
  list(
    cats    = vert,
    dropped = setdiff(k_idx, vert),
    bounds  = diff(s_cum[hull]) / diff(k_idx[hull])
  )
}


#' Category bands and collapsed categories for every item
#'
#' @param thr A data.frame with `Item`, `k` (integer threshold index) and
#'   `Location`.
#' @param item_levels Character vector of item names, in plotting order
#'   (first level at the bottom of the panel).
#' @param xlim Numeric length 2. Used to close the two open-ended outer bands.
#' @return A list with `bands` (Item, cat, xmin, xmax) and `drops`
#'   (Item, cat, x), the latter `NULL` when no category collapses.
#' @noRd
.category_bands <- function(thr, item_levels, xlim) {
  taus <- split(thr$Location[order(thr$Item, thr$k)],
                factor(thr$Item[order(thr$Item, thr$k)], levels = item_levels))
  bands <- list()
  drops <- list()
  for (it in item_levels) {
    tau <- as.numeric(taus[[it]])
    mb <- .modal_boundaries(tau)
    bands[[it]] <- data.frame(
      Item = it,
      cat  = mb$cats,
      xmin = c(xlim[1], mb$bounds),
      xmax = c(mb$bounds, xlim[2]),
      stringsAsFactors = FALSE
    )
    if (length(mb$dropped)) {
      s_cum <- c(0, cumsum(tau))
      drops[[it]] <- data.frame(
        Item = it,
        cat  = mb$dropped,
        x    = vapply(mb$dropped, function(cc) {
          lo <- max(mb$cats[mb$cats < cc])
          hi <- min(mb$cats[mb$cats > cc])
          (s_cum[hi + 1L] - s_cum[lo + 1L]) / (hi - lo)
        }, numeric(1)),
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    bands = do.call(rbind, bands),
    drops = if (length(drops)) do.call(rbind, drops) else NULL
  )
}


#' Reversed adjacent threshold pairs
#'
#' @inheritParams .category_bands
#' @return A data.frame with `Item`, `k`, `lo`, `hi` and `gap`, or `NULL`
#'   when every item has ordered thresholds.
#' @noRd
.threshold_reversals <- function(thr, item_levels) {
  out <- lapply(item_levels, function(it) {
    tau <- thr$Location[thr$Item == it][order(thr$k[thr$Item == it])]
    d <- which(diff(tau) < 0)
    if (!length(d)) return(NULL)
    data.frame(
      Item = it, k = d, lo = tau[d + 1L], hi = tau[d], gap = tau[d] - tau[d + 1L],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out) || !nrow(out)) NULL else out
}


#' Darken pale fills until they read as a thin line on white
#'
#' The lightest end of a sequential palette disappears when drawn as a 1 px
#' error bar. Colours already below the luminance floor are returned unchanged.
#'
#' @param cols Character vector of colours.
#' @param floor Numeric. Maximum perceived luminance on the 0-255 scale.
#' @return A character vector of colours, same length as `cols`.
#' @noRd
.ci_contrast <- function(cols, floor = 130) {
  vapply(cols, function(cl) {
    v <- as.numeric(grDevices::col2rgb(cl)[, 1])
    while (sum(v * c(0.299, 0.587, 0.114)) > floor) v <- v * 0.78
    grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
  }, character(1), USE.NAMES = FALSE)
}


#' Pick readable label colour for text drawn on a filled band
#'
#' @param fills Character vector of band fill colours.
#' @return A character vector of `"grey15"` / `"white"`.
#' @noRd
.band_label_colour <- function(fills) {
  lum <- apply(grDevices::col2rgb(fills), 2, function(v) {
    sum(v * c(0.299, 0.587, 0.114))
  })
  ifelse(lum > 145, "grey15", "white")
}


#' Build the response-category band panel
#'
#' @param thr A data.frame with `Item`, `k`, `Location`, `SE` and, when
#'   `show_ci` is `TRUE`, `CI_low` and `CI_high`.
#' @param item_order Character vector of item names, first level at the bottom.
#' @param xlim Numeric length 2.
#' @param category_labels Character vector of category labels, length `ncat`.
#' @param fills Character vector of band fill colours, length `ncat`.
#' @param show_ci Logical. Draw the threshold error-bar row.
#' @param row_gap Numeric. Vertical spacing between item rows.
#' @param caption_text Character. Pre-built caption.
#' @return A `ggplot` object.
#' @noRd
.targeting_category_panel <- function(
  thr,
  item_order,
  xlim,
  category_labels,
  fills,
  show_ci,
  row_gap,
  caption_text
) {
  rev_col <- "#B2182B"
  bar_h <- 0.52
  hh <- bar_h / 2
  yof <- function(x) match(as.character(x), item_order) * row_gap

  bd <- .category_bands(thr, item_order, xlim)
  b <- bd$bands
  b$y <- yof(b$Item)
  b$category <- factor(category_labels[b$cat + 1L], levels = category_labels)
  b$fill <- fills[b$cat + 1L]

  # Open-ended outer bands fade towards the panel edge, so the bar does not
  # read as though the extreme categories stopped at the axis limits.
  is_out <- b$cat == 0L | b$cat == max(b$cat)
  n_step <- 40L
  fade <- do.call(rbind, lapply(which(is_out), function(i) {
    r <- b[i, ]
    br <- seq(r$xmin, r$xmax, length.out = n_step + 1L)
    data.frame(
      y = r$y,
      category = r$category,
      alpha = if (r$cat == 0L) {
        seq(0.22, 1, length.out = n_step)
      } else {
        seq(1, 0.22, length.out = n_step)
      },
      xmin = br[-length(br)],
      xmax = br[-1],
      stringsAsFactors = FALSE
    )
  }))

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = fade,
      ggplot2::aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$y - hh, ymax = .data$y + hh,
        fill = .data$category, alpha = .data$alpha
      )
    ) +
    ggplot2::geom_rect(
      data = b[!is_out, ],
      ggplot2::aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$y - hh, ymax = .data$y + hh,
        fill = .data$category
      )
    )

  # Hairline at every band join
  joins <- b[b$xmin > xlim[1], ]
  p <- p +
    ggplot2::geom_segment(
      data = joins,
      ggplot2::aes(
        x = .data$xmin, xend = .data$xmin,
        y = .data$y - hh, yend = .data$y + hh
      ),
      colour = "white",
      linewidth = 0.6,
      show.legend = FALSE
    )

  # Category number inside the band, where the band is wide enough to hold it
  lab <- b
  lab$xc <- (pmax(lab$xmin, xlim[1]) + pmin(lab$xmax, xlim[2])) / 2
  keep <- (pmin(lab$xmax, xlim[2]) - pmax(lab$xmin, xlim[1])) >
    0.30 * diff(xlim) / 8
  lab <- lab[keep, ]
  if (nrow(lab)) {
    p <- p +
      ggplot2::geom_text(
        data = lab,
        ggplot2::aes(x = .data$xc, y = .data$y, label = .data$cat),
        colour = .band_label_colour(lab$fill),
        size = 3.1,
        show.legend = FALSE
      )
  }

  # Categories that are never the most likely response collapse to a tick.
  # Red, to tie the tick to the reversal span that explains it, with a white
  # halo so it stays visible against the dark end of the palette.
  drops <- bd$drops
  if (!is.null(drops)) {
    drops <- drops[order(drops$Item, drops$cat), ]
    key <- paste(drops$Item, round(drops$x, 8))
    drops$label <- vapply(
      key,
      function(kk) paste(drops$cat[key == kk], collapse = ","),
      character(1),
      USE.NAMES = FALSE
    )
    drops <- drops[!duplicated(key), ]
    drops$y <- yof(drops$Item)
    p <- p +
      ggplot2::geom_segment(
        data = drops,
        ggplot2::aes(
          x = .data$x, xend = .data$x,
          y = .data$y - hh - 0.06, yend = .data$y + hh + 0.06
        ),
        colour = "white", linewidth = 2.8, lineend = "butt", show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        data = drops,
        ggplot2::aes(
          x = .data$x, xend = .data$x,
          y = .data$y - hh - 0.04, yend = .data$y + hh + 0.04
        ),
        colour = rev_col, linewidth = 1.7, lineend = "butt", show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = drops,
        ggplot2::aes(x = .data$x, y = .data$y + hh + 0.14, label = .data$label),
        colour = rev_col, size = 2.7, fontface = "bold",
        hjust = 0.5, vjust = 0, show.legend = FALSE
      )
  }

  # Threshold estimates below the bar, coloured by the category the threshold
  # gives entry to, so intervals that overlap can still be told apart. A white
  # halo separates intervals from each other.
  if (show_ci) {
    ci <- thr
    ci$y <- yof(ci$Item) - hh - 0.17
    ci$col <- .ci_contrast(fills)[ci$k + 1L]
    p <- p +
      ggplot2::geom_errorbar(
        data = ci,
        ggplot2::aes(y = .data$y, xmin = .data$CI_low, xmax = .data$CI_high),
        width = 0.17, linewidth = 2.1, colour = "white",
        orientation = "y", show.legend = FALSE
      ) +
      ggplot2::geom_errorbar(
        data = ci,
        ggplot2::aes(y = .data$y, xmin = .data$CI_low, xmax = .data$CI_high),
        width = 0.17, linewidth = 0.8, colour = ci$col,
        orientation = "y", show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = ci,
        ggplot2::aes(x = .data$Location, y = .data$y),
        colour = ci$col, size = 1.15, show.legend = FALSE
      )
  }

  # Reversal spans, stacked under the threshold row so several on one item
  # cannot overprint each other
  revs <- .threshold_reversals(thr, item_order)
  n_rev <- 0L
  if (!is.null(revs)) {
    revs <- revs[order(revs$Item, revs$k), ]
    revs$row <- stats::ave(revs$k, revs$Item, FUN = seq_along)
    n_rev <- max(revs$row)
    revs$y <- yof(revs$Item) - hh - 0.32 - (revs$row - 1L) * 0.15
    p <- p +
      ggplot2::geom_segment(
        data = revs,
        ggplot2::aes(x = .data$lo, xend = .data$hi, y = .data$y, yend = .data$y),
        colour = rev_col, linewidth = 0.8, show.legend = FALSE,
        arrow = ggplot2::arrow(
          ends = "both", length = ggplot2::unit(3, "pt"), type = "closed"
        )
      ) +
      ggplot2::geom_text(
        data = revs,
        ggplot2::aes(
          x = .data$hi, y = .data$y, label = sprintf("%.2f", .data$gap)
        ),
        colour = rev_col, size = 2.5, hjust = -0.3, vjust = 0.42,
        show.legend = FALSE
      )
  }

  y_bottom <- max(0.55, 0.32 + 0.15 * n_rev + 0.15)

  p +
    ggplot2::scale_fill_manual(
      values = stats::setNames(fills, category_labels),
      name = NULL,
      drop = FALSE,
      guide = ggplot2::guide_legend(
        nrow = 1,
        override.aes = list(alpha = 1)
      )
    ) +
    ggplot2::scale_alpha_identity(guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = seq_along(item_order) * row_gap,
      labels = item_order,
      expand = ggplot2::expansion(add = c(y_bottom, max(0.55, hh + 0.34)))
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(xlim[1], xlim[2], by = 1),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = "Location (logit scale)", y = NULL, caption = caption_text) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 5, 5, 5)
    ) +
    er2_axis_margins() +
    er2_plot_caption()
}
