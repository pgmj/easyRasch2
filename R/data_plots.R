#' Item Response Distribution Bar Chart
#'
#' Creates a faceted bar chart showing the response distribution for each
#' item, with counts and percentages displayed on each bar. Each item gets
#' its own panel, with response categories on the x-axis and percentage of
#' responses on the y-axis. This is a descriptive data visualization tool
#' intended for use before model fitting.
#'
#' @param data A data.frame in wide format containing only the item
#'   response columns. Each column is one item, each row is one person.
#'   All columns must be numeric (integer-valued). Response categories
#'   may be coded starting from 0 or 1. Do not include person IDs,
#'   grouping variables, or other non-item columns.
#' @param item_labels An optional character vector of descriptive labels
#'   for the items (facet strips). Must be the same length as
#'   `ncol(data)`. If `NULL` (the default), column names are used.
#'   Labels are displayed as `"column_name - label"`.
#' @param category_labels An optional character vector of labels for the
#'   response categories (x-axis). Must be the same length as the number
#'   of response categories spanning from the minimum to the maximum
#'   observed value. If `NULL` (the default), numeric category values
#'   are used.
#' @param ncol Integer. Number of columns in the faceted layout.
#'   Default is `1L`.
#' @param label_wrap Integer. Number of characters per line in facet
#'   strip labels before wrapping. Default is `25L`.
#' @param text_y Numeric. Vertical position (in percent units) for the
#'   count labels on each bar. Adjust upward if bars are tall. Default
#'   is `6`.
#' @param viridis_option Character. Viridis palette option for the
#'   count-text colour. One of `"A"` through `"H"`. Default is `"A"`.
#' @param viridis_end Numeric in \eqn{[0, 1]}. End point of the viridis
#'   colour scale for count text. Adjust if text is hard to read against
#'   the bar colours. Default is `0.9`.
#' @param font Character. Font family for all text. Default is
#'   `"sans"`.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @details
#' Each item is displayed as a separate facet panel with the item label
#' in the strip on the left side. Bars are coloured by response category
#' using the viridis palette. Each bar shows the count (`n = X`) as text.
#'
#' **Input requirements:**
#' \itemize{
#'   \item All columns must be numeric (integer-valued).
#'   \item The data frame must contain at least 2 columns (items)
#'     and at least 1 row (person).
#' }
#'
#' @seealso [RMplotStackedbar()], [RMplotTile()]
#'
#' @examples
#' if (requireNamespace("eRm", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(pcmdat2, package = "eRm")
#'
#'   # Basic response distribution plot
#'   RMplotBar(pcmdat2)
#'
#'   # With custom item labels
#'   RMplotBar(
#'     pcmdat2,
#'     item_labels = c("Mood", "Sleep", "Appetite", "Energy")
#'   )
#'
#'   # Two-column layout with wrapped labels
#'   RMplotBar(
#'     pcmdat2,
#'     item_labels = c(
#'       "General mood and emotional wellbeing",
#'       "Quality of sleep at night",
#'       "Appetite and eating habits",
#'       "Overall energy level during the day"
#'     ),
#'     ncol = 2, label_wrap = 20
#'   )
#'
#'   # With custom category labels
#'   RMplotBar(
#'     pcmdat2,
#'     category_labels = c("Never", "Sometimes", "Often")
#'   )
#' }
#'
#' @importFrom rlang .data
#' @export
RMplotBar <- function(
    data,
    item_labels     = NULL,
    category_labels = NULL,
    ncol            = 1L,
    label_wrap      = 25L,
    text_y          = 6,
    viridis_option  = "A",
    viridis_end     = 0.9,
    font            = "sans"
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for RMplotBar().",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame in wide format ",
         "(one column per item, one row per person).", call. = FALSE)
  }
  if (ncol(data) < 2L) {
    stop("`data` must contain at least 2 columns (items). Found ",
         ncol(data), " column(s).", call. = FALSE)
  }
  if (nrow(data) < 1L) {
    stop("`data` must contain at least 1 row (person). Found 0 rows.",
         call. = FALSE)
  }

  non_numeric <- names(data)[!vapply(data, is.numeric, logical(1L))]
  if (length(non_numeric) > 0L) {
    stop("All columns must be numeric. Non-numeric column(s): ",
         paste(non_numeric, collapse = ", "),
         ". Remove non-item columns (e.g., person IDs, grouping ",
         "variables) before passing to RMplotBar().",
         call. = FALSE)
  }

  all_vals <- unlist(data, use.names = FALSE)
  all_vals <- all_vals[!is.na(all_vals)]
  if (length(all_vals) == 0L) {
    stop("`data` contains only NA values.", call. = FALSE)
  }
  if (!all(all_vals == round(all_vals))) {
    stop("All response values must be integers. Found non-integer values.",
         call. = FALSE)
  }

  item_names <- names(data)

  if (!is.null(item_labels) && length(item_labels) != length(item_names)) {
    stop("`item_labels` must have the same length as the number of ",
         "columns in `data`. Expected ", length(item_names),
         ", got ", length(item_labels), ".", call. = FALSE)
  }

  min_val        <- min(all_vals)
  max_val        <- max(all_vals)
  all_categories <- seq(min_val, max_val)
  n_categories   <- length(all_categories)

  if (!is.null(category_labels) &&
      length(category_labels) != n_categories) {
    stop("`category_labels` must have the same length as the number of ",
         "response categories (", min_val, " to ", max_val, " = ",
         n_categories, " categories). Got ", length(category_labels), ".",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------
  # Prepare data
  # ---------------------------------------------------------------------
  facet_labels <- if (!is.null(item_labels)) {
    paste0(item_names, " - ", item_labels)
  } else {
    item_names
  }

  count_list <- vector("list", length(item_names))
  for (i in seq_along(item_names)) {
    item <- item_names[i]
    item_vals <- data[[item]]
    item_vals <- item_vals[!is.na(item_vals)]
    tab   <- table(factor(item_vals, levels = all_categories))
    total <- sum(tab)

    count_list[[i]] <- data.frame(
      item_name   = item,
      facet_label = facet_labels[i],
      category    = as.integer(names(tab)),
      n           = as.integer(tab),
      percent     = if (total > 0L) {
        round(as.integer(tab) / total * 100, 1)
      } else {
        rep(0, length(tab))
      },
      stringsAsFactors = FALSE
    )
  }
  plot_df <- do.call(rbind, count_list)
  rownames(plot_df) <- NULL

  plot_df$facet_label  <- factor(plot_df$facet_label, levels = facet_labels)
  plot_df$category_fct <- factor(plot_df$category,    levels = all_categories)

  if (!is.null(category_labels)) {
    levels(plot_df$category_fct) <- category_labels
  }

  # ---------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$category_fct, y = .data$percent)
  ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$category_fct)
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0("n = ", .data$n),
        color = .data$category_fct
      ),
      y = text_y
    ) +
    ggplot2::facet_wrap(
      ~ facet_label,
      ncol           = ncol,
      strip.position = "left",
      labeller       = ggplot2::labeller(
        facet_label = ggplot2::label_wrap_gen(label_wrap)
      )
    ) +
    ggplot2::scale_fill_viridis_d(guide = "none") +
    ggplot2::scale_color_viridis_d(
      guide     = "none",
      direction = -1,
      option    = viridis_option,
      end       = viridis_end
    ) +
    ggplot2::scale_y_continuous(position = "right") +
    ggplot2::labs(
      x = "Response category",
      y = "% of responses"
    ) +
    ggplot2::theme_bw(base_family = font) +
    ggplot2::theme(
      legend.position   = "none",
      strip.text.y.left = ggplot2::element_text(angle = 0)
    )
}


#' Stacked Bar Chart of Item Response Distributions
#'
#' Creates a horizontal stacked bar chart showing the response
#' distribution for all items. Each bar represents one item, with segments
#' coloured by response category. Counts are displayed as text labels
#' within each segment. This is a descriptive data visualization tool
#' intended for use before model fitting.
#'
#' @param data A data.frame in wide format containing only the item
#'   response columns. Each column is one item, each row is one person.
#'   All columns must be numeric (integer-valued). Response categories
#'   may be coded starting from 0 or 1. Do not include person IDs,
#'   grouping variables, or other non-item columns.
#' @param item_labels An optional character vector of descriptive labels
#'   for the items (y-axis). Must be the same length as `ncol(data)`.
#'   If `NULL` (the default), column names are used.
#' @param category_labels An optional character vector of labels for the
#'   response categories (legend). Must be the same length as the number
#'   of response categories spanning from the minimum to the maximum
#'   observed value, ordered from lowest to highest category. If `NULL`
#'   (the default), numeric category values are used.
#' @param show_n Logical. If `TRUE` (the default), the count of responses
#'   is displayed as a text label inside each bar segment.
#' @param show_percent Logical. If `TRUE`, the percentage of responses is
#'   displayed instead of (or in addition to) counts. Default is `FALSE`.
#' @param text_color Character. Colour for the count/percentage labels.
#'   Default is `"sienna1"`.
#' @param text_size Numeric. Size of the count/percentage labels.
#'   Default is `3`.
#' @param min_label_n Integer. Minimum count required for a label to be
#'   displayed within a bar segment. Segments with fewer responses are
#'   left unlabelled to avoid clutter. Default is `0L` (all segments
#'   labelled).
#' @param viridis_option Character. Viridis palette option. One of
#'   `"A"` through `"H"`. Default is `"D"`.
#' @param viridis_end Numeric in \eqn{[0, 1]}. End point of the viridis
#'   colour scale. Default is `0.99`.
#' @param title Character. Plot title. Default is `"Item responses"`.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @details
#' Items are displayed on the y-axis in the same order as the columns
#' in `data` (first column at the top). Each bar is divided into
#' segments representing response categories, with the lowest category
#' on the left and the highest on the right. The total bar length
#' equals the number of non-missing responses for that item.
#'
#' Categories with zero responses still appear in the legend but produce
#' no visible bar segment, which helps identify gaps in the response
#' distribution.
#'
#' **Input requirements:**
#' \itemize{
#'   \item All columns must be numeric (integer-valued).
#'   \item The data frame must contain at least 2 columns (items)
#'     and at least 1 row (person).
#' }
#'
#' @seealso [RMplotBar()], [RMplotTile()]
#'
#' @examples
#' if (requireNamespace("eRm", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(pcmdat2, package = "eRm")
#'
#'   # Basic stacked bar chart
#'   RMplotStackedbar(pcmdat2)
#'
#'   # With custom item and category labels
#'   RMplotStackedbar(
#'     pcmdat2,
#'     item_labels     = c("Mood", "Sleep", "Appetite", "Energy"),
#'     category_labels = c("Never", "Sometimes", "Often")
#'   )
#'
#'   # Show percentages, suppress small segments
#'   RMplotStackedbar(
#'     pcmdat2,
#'     show_percent = TRUE,
#'     show_n       = FALSE,
#'     min_label_n  = 5
#'   )
#' }
#'
#' @importFrom rlang .data
#' @export
RMplotStackedbar <- function(
    data,
    item_labels     = NULL,
    category_labels = NULL,
    show_n          = TRUE,
    show_percent    = FALSE,
    text_color      = "sienna1",
    text_size       = 3,
    min_label_n     = 0L,
    viridis_option  = "D",
    viridis_end     = 0.99,
    title           = "Item responses"
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for RMplotStackedbar().",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame in wide format ",
         "(one column per item, one row per person).", call. = FALSE)
  }
  if (ncol(data) < 2L) {
    stop("`data` must contain at least 2 columns (items). Found ",
         ncol(data), " column(s).", call. = FALSE)
  }
  if (nrow(data) < 1L) {
    stop("`data` must contain at least 1 row (person). Found 0 rows.",
         call. = FALSE)
  }

  non_numeric <- names(data)[!vapply(data, is.numeric, logical(1L))]
  if (length(non_numeric) > 0L) {
    stop("All columns must be numeric. Non-numeric column(s): ",
         paste(non_numeric, collapse = ", "),
         ". Remove non-item columns (e.g., person IDs, grouping ",
         "variables) before passing to RMplotStackedbar().",
         call. = FALSE)
  }

  all_vals <- unlist(data, use.names = FALSE)
  all_vals <- all_vals[!is.na(all_vals)]
  if (length(all_vals) == 0L) {
    stop("`data` contains only NA values.", call. = FALSE)
  }
  if (!all(all_vals == round(all_vals))) {
    stop("All response values must be integers. Found non-integer values.",
         call. = FALSE)
  }

  item_names <- names(data)

  if (!is.null(item_labels) && length(item_labels) != length(item_names)) {
    stop("`item_labels` must have the same length as the number of ",
         "columns in `data`. Expected ", length(item_names),
         ", got ", length(item_labels), ".", call. = FALSE)
  }

  min_val        <- min(all_vals)
  max_val        <- max(all_vals)
  all_categories <- seq(min_val, max_val)
  n_categories   <- length(all_categories)

  if (!is.null(category_labels) &&
      length(category_labels) != n_categories) {
    stop("`category_labels` must have the same length as the number of ",
         "response categories (", min_val, " to ", max_val, " = ",
         n_categories, " categories). Got ", length(category_labels), ".",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------
  # Prepare data
  # ---------------------------------------------------------------------
  y_labels <- if (!is.null(item_labels)) {
    stats::setNames(item_labels, item_names)
  } else {
    stats::setNames(item_names, item_names)
  }

  count_list <- vector("list", length(item_names))
  for (i in seq_along(item_names)) {
    item <- item_names[i]
    item_vals <- data[[item]]
    item_vals <- item_vals[!is.na(item_vals)]
    tab   <- table(factor(item_vals, levels = all_categories))
    total <- sum(tab)

    count_list[[i]] <- data.frame(
      item_name = item,
      category  = as.integer(names(tab)),
      n         = as.integer(tab),
      percent   = if (total > 0L) {
        round(as.integer(tab) / total * 100, 1)
      } else {
        rep(0, length(tab))
      },
      stringsAsFactors = FALSE
    )
  }
  plot_df <- do.call(rbind, count_list)
  rownames(plot_df) <- NULL

  # Items: reversed so the first column appears at the top of the y-axis
  plot_df$item_fct <- factor(
    plot_df$item_name,
    levels = rev(item_names),
    labels = rev(y_labels[item_names])
  )

  # Categories: reversed for stacking order (lowest category on the left)
  plot_df$category_fct <- factor(
    plot_df$category,
    levels = rev(all_categories)
  )

  # Build text labels
  plot_df$label <- ""
  show_mask <- plot_df$n >= min_label_n

  if (show_n && show_percent) {
    plot_df$label[show_mask] <- paste0(
      plot_df$n[show_mask], " (", plot_df$percent[show_mask], "%)"
    )
  } else if (show_n) {
    plot_df$label[show_mask] <- as.character(plot_df$n[show_mask])
  } else if (show_percent) {
    plot_df$label[show_mask] <- paste0(plot_df$percent[show_mask], "%")
  }

  # Suppress labels for zero-count segments
  plot_df$label[plot_df$n == 0L] <- ""

  # Legend labels: lowest -> highest (reverse of stacking factor)
  legend_labels <- if (!is.null(category_labels)) {
    stats::setNames(rev(category_labels),
                    rev(as.character(all_categories)))
  } else {
    stats::setNames(rev(as.character(all_categories)),
                    rev(as.character(all_categories)))
  }

  # ---------------------------------------------------------------------
  # Build plot
  # ---------------------------------------------------------------------
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x    = .data$n,
      y    = .data$item_fct,
      fill = .data$category_fct
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_viridis_d(
      "Response\ncategory",
      direction = -1,
      option    = viridis_option,
      end       = viridis_end,
      labels    = legend_labels
    ) +
    ggplot2::labs(
      title = title,
      x     = "Number of responses",
      y     = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")

  if (show_n || show_percent) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$label),
        color    = text_color,
        size     = text_size,
        position = ggplot2::position_stack(vjust = 0.5)
      )
  }

  p
}
