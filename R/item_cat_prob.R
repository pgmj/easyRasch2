#' Item Category Probability Curves
#'
#' Plots model-implied response-category probability curves for each item
#' as a function of the latent trait \eqn{\theta}. Polytomous items are
#' fitted with the Partial Credit Model via `eRm::PCM()`; dichotomous
#' items are fitted with the Rasch model via `eRm::RM()`. Each item gets
#' its own facet panel, with one curve per response category coloured
#' from low to high using the viridis palette. Comparable in scope to
#' `eRm::plotICC()` and `mirt`'s trace plots, with a `ggplot2` /
#' viridis output and optional descriptive labels for items and
#' categories.
#'
#' @param data A data.frame or matrix of item responses. Items must be
#'   scored starting at 0 (non-negative integers). Missing values
#'   (`NA`) are allowed --- the model fit handles them.
#' @param item_labels Optional character vector of descriptive item
#'   labels (facet strip titles). Must be the same length as
#'   `ncol(data)`. If `NULL` (the default), column names are used.
#' @param category_labels Optional character vector of labels for the
#'   response categories (legend). Must be the same length as the
#'   number of categories spanning from 0 to the maximum observed
#'   value. If `NULL` (the default), numeric category values are used.
#' @param theta_range Numeric length 2. Range of the latent trait
#'   \eqn{\theta} (logits) plotted along the x-axis. When `NULL`
#'   (default) the range is `c(-4, 4)` for `label_curves = "legend"`
#'   and a wider `c(-5, 5)` for `label_curves = "path"` --- the
#'   wider range gives `geomtextpath` more horizontal room to place
#'   per-curve labels at their peaks. Pass an explicit
#'   `c(lo, hi)` to override.
#' @param n_points Integer. Number of evenly-spaced \eqn{\theta} values
#'   at which the curves are evaluated. Default `200L`.
#' @param viridis_option Character. Viridis palette identifier. One of
#'   `"A"` through `"H"`. Default `"D"` (viridis green).
#' @param viridis_end Numeric in (0, 1]. Upper end of the viridis
#'   palette range; lower values keep the palette inside its mid-tones
#'   (avoids the very bright yellow at `1.0`). Default `0.95`.
#' @param facet_ncol Optional integer. Number of columns in the facet
#'   layout. Default `NULL` (`ggplot2::facet_wrap()` auto-layout).
#' @param label_wrap Integer. Characters per line for facet-strip
#'   label wrapping. Default `25L`.
#' @param line_width Numeric. Line width for the probability curves.
#'   Default `0.9`.
#' @param font Character. Font family for all text. Default `"sans"`.
#' @param output Character. Either `"ggplot"` (default) for a
#'   `ggplot` object, or `"dataframe"` for the underlying long-format
#'   probability table.
#' @param label_curves Character. How response categories are
#'   identified. `"legend"` (default) draws lines and a separate
#'   colour legend with one swatch per category --- the standard
#'   multi-item presentation. `"path"` uses
#'   `geomtextpath::geom_textpath()` to write each category's label
#'   *along* its own curve (a la classic IRT trace plots) and
#'   suppresses the legend; valid for **a single item at a time**
#'   because narrow facets do not give path labels enough horizontal
#'   room to render legibly. The model is still fit on the full
#'   multi-item `data` (PCM thresholds require multi-item input for
#'   CML estimation); the `item` argument then selects which item's
#'   curves to plot. Requires the optional `geomtextpath` package.
#' @param item Character or integer. Used only when
#'   `label_curves = "path"`. Either an item name (matching a column
#'   in `data`) or a column index, identifying which item's category
#'   curves to plot. Required for path mode.
#' @param text_size Numeric. Used only when `label_curves = "path"`.
#'   Size of the path labels in mm. Default `4`.
#'
#' @return
#' * If `output = "ggplot"`: a [ggplot2::ggplot] object, one facet per
#'   item.
#' * If `output = "dataframe"`: a long-format data.frame with columns
#'   `Item` (factor in column order), `Category` (integer), and
#'   `Theta`, `Probability` (numeric), one row per item × category ×
#'   theta gridpoint.
#'
#' @details
#' For each polytomous item *i* with response categories
#' \eqn{0, 1, \ldots, K_i} and threshold parameters
#' \eqn{\delta_{i,1}, \ldots, \delta_{i,K_i}} from
#' `eRm::thresholds()`, the PCM category probability is
#'
#' \deqn{P(X_i = k \mid \theta) = \frac{\exp(\sum_{j=1}^{k}
#'   (\theta - \delta_{i,j}))}{\sum_{k'=0}^{K_i} \exp(\sum_{j=1}^{k'}
#'   (\theta - \delta_{i,j}))},}
#'
#' with the empty sum (when \eqn{k = 0}) taken as zero. For
#' dichotomous items the function fits a Rasch model and treats
#' the item difficulty \eqn{\delta_i = -\beta_i} as the single
#' threshold, recovering the standard two-category logistic ICC.
#'
#' The colour mapping uses `scale_color_viridis_c()` against the
#' integer category value, so the natural ordering of response
#' categories is preserved visually --- low categories at one end of
#' the palette, high categories at the other. When `category_labels`
#' is provided, the legend uses those labels (e.g., "Never" /
#' "Sometimes" / "Often") while the colour mapping stays on the
#' integer category value.
#'
#' Items with fewer response categories than the maximum (e.g., an
#' otherwise four-category scale with one three-category item)
#' contribute only the categories they actually have to their own
#' facet --- the y-axis still spans \eqn{[0, 1]}.
#'
#' @seealso [RMitemICCPlot()] for conditional ICCs binned by total
#'   score, [RMitemHierarchy()] for item threshold locations on the
#'   logit scale.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("eRm", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(pcmdat2, package = "eRm")
#'
#'   # Default plot
#'   RMitemCatProb(pcmdat2)
#'
#'   # Custom item and category labels
#'   RMitemCatProb(
#'     pcmdat2,
#'     item_labels     = c("Mood", "Sleep", "Appetite", "Energy"),
#'     category_labels = c("Never", "Sometimes", "Often")
#'   )
#'
#'   # Underlying probability data
#'   df <- RMitemCatProb(pcmdat2, output = "dataframe")
#'   head(df)
#'
#'   # Single-item plot with labels written along each curve
#'   # (classic IRT trace-plot style). Model is still fit on all four
#'   # items; `item` picks which one to plot.
#'   if (requireNamespace("geomtextpath", quietly = TRUE)) {
#'     RMitemCatProb(
#'       pcmdat2,
#'       category_labels = c("Never", "Sometimes", "Often"),
#'       label_curves    = "path",
#'       item            = "I1"
#'     )
#'   }
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
RMitemCatProb <- function(
    data,
    item_labels     = NULL,
    category_labels = NULL,
    theta_range     = NULL,
    n_points        = 200L,
    viridis_option  = "D",
    viridis_end     = 0.95,
    facet_ncol      = NULL,
    label_wrap      = 25L,
    line_width      = 0.9,
    font            = "sans",
    output          = c("ggplot", "dataframe"),
    label_curves    = c("legend", "path"),
    item            = NULL,
    text_size       = 4
) {

  output       <- match.arg(output)
  label_curves <- match.arg(label_curves)

  # Resolve default theta_range based on label_curves: wider for path mode
  # so geom_textpath has enough horizontal room to place labels at peaks.
  if (is.null(theta_range)) {
    theta_range <- if (label_curves == "path") c(-5, 5) else c(-4, 4)
  }

  if (output == "ggplot" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for output = \"ggplot\".",
         call. = FALSE)
  }
  if (!requireNamespace("eRm", quietly = TRUE)) {
    stop("Package 'eRm' is required for RMitemCatProb().",
         call. = FALSE)
  }
  if (label_curves == "path" &&
      !requireNamespace("geomtextpath", quietly = TRUE)) {
    stop(
      "Package 'geomtextpath' is required when `label_curves = \"path\"`.\n",
      "Install it with: install.packages(\"geomtextpath\")",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------
  validate_response_data(data)

  # PCM / RM require multi-item data for CML estimation, so `data` must
  # always have >= 2 items. Path mode uses the `item` argument to filter
  # which item's curves to plot AFTER the model is fit on all items.
  if (ncol(data) < 2L) {
    stop("`data` must contain at least 2 columns (items). Found ",
         ncol(data), " column(s).", call. = FALSE)
  }

  if (!is.numeric(text_size) || length(text_size) != 1L ||
      !is.finite(text_size) || text_size <= 0) {
    stop("`text_size` must be a single positive number.", call. = FALSE)
  }

  if (!is.numeric(theta_range) || length(theta_range) != 2L ||
      theta_range[1L] >= theta_range[2L]) {
    stop("`theta_range` must be a numeric vector of length 2 with ",
         "theta_range[1] < theta_range[2].", call. = FALSE)
  }

  if (!is.numeric(n_points) || length(n_points) != 1L ||
      n_points < 2L || n_points != as.integer(n_points)) {
    stop("`n_points` must be a single integer >= 2.", call. = FALSE)
  }
  n_points <- as.integer(n_points)

  if (!is.numeric(viridis_end) || length(viridis_end) != 1L ||
      viridis_end <= 0 || viridis_end > 1) {
    stop("`viridis_end` must be a single number in (0, 1].",
         call. = FALSE)
  }

  item_names <- colnames(data)
  if (is.null(item_names)) item_names <- paste0("V", seq_len(ncol(data)))
  n_items <- ncol(data)

  if (!is.null(item_labels) && length(item_labels) != n_items) {
    stop("`item_labels` must have the same length as ncol(data) (",
         n_items, "). Got ", length(item_labels), ".", call. = FALSE)
  }

  # Resolve `item` for path mode. Accept either a single item name or a
  # column index; produce a single character item name.
  selected_item <- NULL
  if (label_curves == "path") {
    if (is.null(item)) {
      stop(
        "`label_curves = \"path\"` requires an `item` argument (a single ",
        "item name or column index) identifying which item's category ",
        "curves to plot. Available items: ",
        paste(item_names, collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (length(item) != 1L) {
      stop("`item` must identify exactly one item (length 1).",
           call. = FALSE)
    }
    if (is.numeric(item)) {
      if (!is.finite(item) || item != as.integer(item) ||
          item < 1L || item > n_items) {
        stop("`item` index must be an integer in 1:", n_items, ".",
             call. = FALSE)
      }
      selected_item <- item_names[as.integer(item)]
    } else if (is.character(item)) {
      if (!(item %in% item_names)) {
        stop("`item` \"", item, "\" not found in data. ",
             "Available items: ", paste(item_names, collapse = ", "), ".",
             call. = FALSE)
      }
      selected_item <- item
    } else {
      stop("`item` must be a character item name or a numeric index.",
           call. = FALSE)
    }
  }

  data_mat <- as.matrix(data)
  max_score <- max(data_mat, na.rm = TRUE)
  if (!is.finite(max_score) || max_score < 1L) {
    stop("`data` must contain at least one item with maximum score >= 1.",
         call. = FALSE)
  }
  n_categories <- max_score + 1L

  if (!is.null(category_labels) &&
      length(category_labels) != n_categories) {
    stop("`category_labels` must have the same length as the number of ",
         "response categories (0 to ", max_score, " = ", n_categories,
         " categories). Got ", length(category_labels), ".",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------
  # Fit model + extract per-item threshold vectors
  # ---------------------------------------------------------------------
  is_polytomous <- max_score > 1L

  # Per-item thresholds: a list, one numeric vector per item.
  if (is_polytomous) {
    pcm_fit <- eRm::PCM(data)
    thresh_obj <- eRm::thresholds(pcm_fit)
    thresh_mat <- thresh_obj$threshtable[[1L]]
    # Drop the "Location" column (always first); keep only the actual
    # threshold columns (Threshold 1, Threshold 2, ...).
    if ("Location" %in% colnames(thresh_mat)) {
      thresh_mat <- thresh_mat[, -1L, drop = FALSE]
    }
    item_thresholds <- lapply(seq_len(nrow(thresh_mat)), function(i) {
      v <- as.numeric(thresh_mat[i, ])
      v[!is.na(v)]
    })
    names(item_thresholds) <- rownames(thresh_mat)
  } else {
    rm_fit <- eRm::RM(data)
    # Item difficulty delta_i = -beta_i (eRm parametrises beta as
    # easiness; we want difficulty for the threshold convention).
    deltas <- as.numeric(-rm_fit$betapar)
    item_thresholds <- lapply(deltas, function(d) d)
    names(item_thresholds) <- item_names
  }

  # Make sure the order matches the input data column order
  item_thresholds <- item_thresholds[item_names]

  # ---------------------------------------------------------------------
  # Compute category probability curves on the theta grid
  # ---------------------------------------------------------------------
  theta <- seq(theta_range[1L], theta_range[2L], length.out = n_points)

  # For one item: returns a length(theta) x (K+1) probability matrix.
  pcm_probs <- function(thresholds, theta) {
    K <- length(thresholds)
    # Cumulative sums of (theta - delta_j) for j = 0, 1, ..., K
    numerators <- matrix(0, nrow = length(theta), ncol = K + 1L)
    for (k in seq_len(K)) {
      numerators[, k + 1L] <- numerators[, k] + (theta - thresholds[k])
    }
    exp_num <- exp(numerators)
    exp_num / rowSums(exp_num)
  }

  # Build long-format probability data.frame
  per_item_dfs <- lapply(seq_along(item_thresholds), function(i) {
    item   <- item_names[i]
    thr    <- item_thresholds[[i]]
    probs  <- pcm_probs(thr, theta)
    K_i    <- length(thr)
    out <- data.frame(
      Item        = rep(item, times = length(theta) * (K_i + 1L)),
      Category    = rep(0L:K_i, each = length(theta)),
      Theta       = rep(theta, times = K_i + 1L),
      Probability = as.numeric(probs),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    out
  })
  plot_df <- do.call(rbind, per_item_dfs)
  rownames(plot_df) <- NULL

  # Apply item labels to facet titles (preserve input column order)
  facet_labels <- if (!is.null(item_labels)) {
    paste0(item_names, " - ", item_labels)
  } else {
    item_names
  }
  facet_lookup <- stats::setNames(facet_labels, item_names)
  plot_df$Item <- factor(facet_lookup[plot_df$Item],
                         levels = facet_labels)

  # In path mode, subset to the chosen item AFTER the model has been
  # fit on all items (so thresholds are CML-estimated against the full
  # multi-item dataset).
  if (label_curves == "path") {
    facet_label_selected <- facet_lookup[[selected_item]]
    plot_df <- plot_df[as.character(plot_df$Item) == facet_label_selected,
                       , drop = FALSE]
    rownames(plot_df) <- NULL
  }

  # ---------------------------------------------------------------------
  # Return early for dataframe output
  # ---------------------------------------------------------------------
  if (output == "dataframe") {
    return(plot_df)
  }

  # ---------------------------------------------------------------------
  # Build ggplot
  # ---------------------------------------------------------------------
  cat_breaks <- 0L:max_score
  cat_labels <- if (!is.null(category_labels)) {
    category_labels
  } else {
    as.character(cat_breaks)
  }

  # Attach a human-readable label per row for use in path mode.
  plot_df$CategoryLabel <- cat_labels[plot_df$Category + 1L]

  # Common scale + axis components for both modes
  base_scales <- list(
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ),
    ggplot2::scale_x_continuous(
      limits = theta_range,
      breaks = seq(ceiling(theta_range[1L]),
                   floor(theta_range[2L]), by = 1)
    ),
    ggplot2::labs(
      x = expression(paste("Latent trait ", theta, " (logits)")),
      y = "Category probability"
    ),
    ggplot2::theme_bw(base_family = font),
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95",
                                               colour = NA),
      strip.text       = ggplot2::element_text(face = "bold")
    ),
    er2_axis_margins()
  )

  if (label_curves == "legend") {

    # -----------------------------------------------------------------
    # Mode 1: multi-item faceted plot with a separate colour legend
    # -----------------------------------------------------------------
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x     = .data$Theta,
        y     = .data$Probability,
        group = .data$Category,
        color = .data$Category
      )
    ) +
      ggplot2::geom_line(linewidth = line_width) +
      ggplot2::scale_color_viridis_c(
        name    = "Response\ncategory",
        breaks  = cat_breaks,
        labels  = cat_labels,
        limits  = c(0, max_score),
        option  = viridis_option,
        end     = viridis_end,
        # guide_legend gives one swatch per category (not a continuous bar),
        # so all `n_categories` labels are visible. ggplot's default places
        # the smallest break at the top of a vertical legend, which is the
        # natural reading order for ordinal "least -> most" categories
        # (e.g. PHQ-9: "Not at all" -> "Nearly every day").
        guide   = ggplot2::guide_legend(
          override.aes = list(linewidth = 2)
        )
      ) +
      ggplot2::facet_wrap(
        ~ Item,
        ncol     = facet_ncol,
        labeller = ggplot2::labeller(
          Item = ggplot2::label_wrap_gen(label_wrap)
        )
      ) +
      base_scales

  } else {

    # -----------------------------------------------------------------
    # Mode 2: single-item plot with labels along each curve.
    #
    # For each category we compute the *modal theta* (the theta value
    # at which that category has its highest probability) and convert
    # that to a numeric `hjust` fraction along the path. This puts the
    # label exactly where the curve peaks:
    #   - Monotone-decreasing extreme (Category 0): peak at the left
    #     edge of the plot.
    #   - Bell-shaped middle categories: peak at the modal theta
    #     between the surrounding thresholds.
    #   - Monotone-increasing extreme (Category K): peak at the right
    #     edge of the plot.
    #
    # Why not `hjust = "auto"`? geomtextpath's "auto" mode picks the
    # *smoothest* path segment, which for sharp bell curves means the
    # long, flat near-zero tails --- so labels get pushed to the plot
    # edges, not the peaks. Computing modal theta directly is robust.
    #
    # A small margin (3% of the plotted range) clamps `hjust` away
    # from 0 and 1 so the extreme-category labels don't get clipped
    # at the plot edges.
    # -----------------------------------------------------------------
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x     = .data$Theta,
        y     = .data$Probability,
        group = .data$Category,
        color = .data$Category,
        label = .data$CategoryLabel
      )
    ) +
      ggplot2::scale_color_viridis_c(
        breaks = cat_breaks,
        limits = c(0, max_score),
        option = viridis_option,
        end    = viridis_end,
        guide  = "none"
      ) +
      ggplot2::ggtitle(as.character(plot_df$Item[1L])) +
      base_scales

    edge_margin   <- 0.03 * diff(theta_range)
    hjust_lo_clip <- edge_margin / diff(theta_range)
    hjust_hi_clip <- 1 - hjust_lo_clip

    for (cat in sort(unique(plot_df$Category))) {
      sub_cat     <- plot_df[plot_df$Category == cat, , drop = FALSE]
      modal_theta <- sub_cat$Theta[which.max(sub_cat$Probability)]
      hjust_val   <- (modal_theta - theta_range[1L]) / diff(theta_range)
      hjust_val   <- max(hjust_lo_clip, min(hjust_hi_clip, hjust_val))
      p <- p + geomtextpath::geom_textpath(
        data      = sub_cat,
        hjust     = hjust_val,
        vjust     = -0.08,
        linewidth = line_width,
        size      = text_size
      )
    }
  }

  p
}
