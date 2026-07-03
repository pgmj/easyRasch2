#' Person-Item Targeting Plot (Wright Map)
#'
#' Produces a three-panel targeting plot with a shared logit scale x-axis:
#' \enumerate{
#'   \item \strong{Top}: Histogram of person location estimates, with a
#'     reference line for the mean (or median) and shading for ±1 SD
#'     (or ±1 MAD).
#'   \item \strong{Middle}: Inverted histogram of item threshold locations,
#'     with the same summary annotations.
#'   \item \strong{Bottom}: Dot-and-whisker plot of individual item thresholds
#'     with confidence intervals based on threshold standard errors.
#' }
#'
#' Together, the top and middle panels form a back-to-back histogram that
#' makes it easy to assess whether the test is well-targeted to the sample.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed.
#' @param robust Logical. If `FALSE` (the default), histogram annotations use
#'   mean ± SD. If `TRUE`, median ± MAD is used instead.
#' @param sort_items Character string controlling item ordering on the y-axis
#'   of the bottom panel. `"data"` (the default) preserves the column order in
#'   `data` (first item at top). `"location"` sorts items by their average
#'   threshold location (easiest at top, hardest at bottom).
#' @param bins Integer. Number of bins for both histograms. Default is number of
#'   unique scores plus one, but no less than 15.
#' @param xlim Numeric vector of length 2. Initial lower and upper limits for
#'   the shared x-axis. Automatically expanded if any person or item threshold
#'   values fall outside these limits.
#' @param ci_level Numeric. Confidence level for the item threshold error bars.
#'   Default is `0.95` (95% CI). Set to `NULL` to hide error bars.
#' @param person_fill Fill colour for the person histogram. Default
#'   `"#0072B2"` (blue).
#' @param threshold_fill Fill colour for the item threshold histogram. Default
#'   `"#D55E00"` (vermillion).
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
#'   `p2` (threshold histogram), and `p3` (item threshold dot-whisker plot).
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
#'   # Default: mean/SD, data order, 95% CI
#'   RMtargeting(sim_data)
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
  robust = FALSE,
  sort_items = c("data", "location"),
  bins,
  xlim = c(-4, 4),
  ci_level = 0.95,
  person_fill = "#0072B2",
  threshold_fill = "#D55E00",
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
  # unique scores plus one, but no less than 15
  if (missing(bins)) {
    bins <- max(rowSums(data, na.rm = TRUE)) + 1
    if (bins < 15) {
      bins <- 15
    }
  }

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
      breaks = function(lim) {
        max_val <- abs(floor(lim[1]))
        seq(0, max_val, by = max(1, round(max_val / 4)))
      }
    ) +
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

  caption_text <- er2_caption(paste0(
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
    ".\n",
    ci_caption
  ))

  if (is_dicho) {
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
