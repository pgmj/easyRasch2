#' Tile Plot of Item Response Distributions
#'
#' Creates a tile (heat map) plot showing the distribution of responses
#' across all items and response categories. Each cell displays the count
#' (or percentage) of responses, with optional conditional highlighting
#' for cells with low counts. Optional faceting by a grouping variable
#' is provided for inspecting subgroup response distributions before DIF
#' analyses -- particularly useful for spotting empty categories or
#' under-represented subgroups before fitting Rasch models per group.
#'
#' Adapted from `easyRaschBayes::plot_tile()` and extended with the
#' `group` parameter for faceted display.
#'
#' @param data A data.frame in wide format containing only the item
#'   response columns. Each column is one item, each row is one person.
#'   All columns must be numeric (integer-valued). Response categories
#'   may be coded starting from 0 or 1. Do not include person IDs,
#'   grouping variables, or other non-item columns -- supply the grouping
#'   variable separately via `group`.
#' @param group Optional vector of length `nrow(data)` (factor, character,
#'   or numeric) defining a grouping variable. When provided, the plot is
#'   faceted by group and counts / percentages are computed within each
#'   group. Default `NULL` (no faceting). Persons with `NA` group are
#'   excluded.
#' @param cutoff Integer. Cells with counts below this value are
#'   highlighted (when `highlight = TRUE`). Default `10`.
#' @param highlight Logical. If `TRUE` (default), cell labels with counts
#'   below `cutoff` are displayed in red. This includes empty cells
#'   (`n = 0`), useful for identifying gaps in the response distribution.
#' @param percent Logical. If `TRUE`, cell labels show percentages
#'   instead of raw counts. Percentages are computed within item (and
#'   within group, when `group` is supplied). Default `FALSE`.
#' @param text_color Character. Colour for non-highlighted cell labels.
#'   Default `"orange"`.
#' @param item_labels Optional character vector of descriptive labels
#'   for the items (y-axis), same length as `ncol(data)`. Default `NULL`
#'   uses the column names.
#' @param category_labels Optional character vector of labels for the
#'   response categories (x-axis), same length as the number of
#'   categories spanning `min` to `max` observed value. Default `NULL`
#'   uses the numeric category values.
#' @param group_labels Optional character vector of length
#'   `nlevels(as.factor(group))` to override the displayed facet labels.
#'   Order corresponds to `levels(as.factor(group))`.
#' @param facet_ncol Integer or `NULL`. Number of columns in the facet
#'   grid when `group` is supplied. Default `NULL` (ggplot2 chooses).
#' @param output Character. `"ggplot"` (default) returns the plot;
#'   `"dataframe"` returns the underlying per-cell counts (one row per
#'   item x category, plus group when supplied).
#'
#' @return Either a `ggplot` object or a data.frame, depending on
#'   `output`.
#'
#' @details
#' Items are placed on the y-axis (in the same order as the columns of
#' `data`, top to bottom) and response categories on the x-axis. Cell
#' shading represents the count of responses (darker = more responses).
#' Categories with zero responses are explicitly shown (`n = 0`), which
#' helps identify gaps in the response distribution -- one of the primary
#' purposes of the plot, especially before DIF analyses where
#' under-represented categories within a subgroup can break model
#' fitting on that subgroup.
#'
#' When `group` is supplied, percentages and the highlight cutoff are
#' applied within each group, so a cell labelled "5" in the group-A
#' facet contains the count for group A only.
#'
#' @examples
#' \donttest{
#' data("pcmdat2", package = "eRm")
#'
#' # Basic tile plot
#' RMplotTile(pcmdat2)
#'
#' # With percentages
#' RMplotTile(pcmdat2, percent = TRUE)
#'
#' # Faceted by an external grouping variable
#' set.seed(1)
#' grp <- sample(c("A", "B"), nrow(pcmdat2), replace = TRUE)
#' RMplotTile(pcmdat2, group = grp)
#'
#' # With custom labels and tighter cutoff
#' RMplotTile(pcmdat2,
#'            group = grp,
#'            group_labels = c("Female", "Male"),
#'            cutoff = 5,
#'            facet_ncol = 2)
#'
#' # Underlying counts as a data.frame
#' RMplotTile(pcmdat2, group = grp, output = "dataframe")
#' }
#'
#' @importFrom rlang .data
#' @export
RMplotTile <- function(
    data,
    group           = NULL,
    cutoff          = 10,
    highlight       = TRUE,
    percent         = FALSE,
    text_color      = "orange",
    item_labels     = NULL,
    category_labels = NULL,
    group_labels    = NULL,
    facet_ncol      = NULL,
    output          = c("ggplot", "dataframe")
) {

  output <- match.arg(output)

  if (output == "ggplot" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for output = \"ggplot\".",
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
         "variables) before passing to RMplotTile() -- pass groups via ",
         "the `group` argument instead.", call. = FALSE)
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

  if (!is.null(item_labels) && length(item_labels) != ncol(data)) {
    stop("`item_labels` must have the same length as the number of ",
         "columns in `data`. Expected ", ncol(data), ", got ",
         length(item_labels), ".", call. = FALSE)
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

  # Group validation
  group_factor <- NULL
  if (!is.null(group)) {
    if (length(group) != nrow(data)) {
      stop("`group` must have the same length as nrow(data) (",
           nrow(data), "). Got ", length(group), ".", call. = FALSE)
    }
    group_factor <- as.factor(group)
    if (nlevels(droplevels(group_factor)) < 2L) {
      stop("`group` must have at least 2 distinct non-missing values.",
           call. = FALSE)
    }
    if (!is.null(group_labels)) {
      if (length(group_labels) != nlevels(group_factor)) {
        stop("`group_labels` must have the same length as the number of ",
             "group levels (", nlevels(group_factor), "). Got ",
             length(group_labels), ".", call. = FALSE)
      }
      levels(group_factor) <- as.character(group_labels)
    }
  }

  # ---------------------------------------------------------------------
  # Prepare data
  # ---------------------------------------------------------------------
  item_names <- names(data)
  label_map <- if (!is.null(item_labels)) {
    stats::setNames(item_labels, item_names)
  } else {
    stats::setNames(item_names, item_names)
  }

  # Per-item x category counts; optionally per group too
  build_counts <- function(sub_data, group_label = NULL) {
    rows <- lapply(item_names, function(item) {
      vals <- sub_data[[item]]
      vals <- vals[!is.na(vals)]
      tab  <- table(factor(vals, levels = all_categories))
      data.frame(
        item     = item,
        category = as.integer(names(tab)),
        n        = as.integer(tab),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
    if (!is.null(group_label)) out$group <- group_label
    rownames(out) <- NULL
    out
  }

  if (is.null(group_factor)) {
    count_df <- build_counts(data)
    agg_keys <- "item"
  } else {
    parts <- lapply(levels(group_factor), function(g) {
      mask <- !is.na(group_factor) & group_factor == g
      build_counts(data[mask, , drop = FALSE], group_label = g)
    })
    count_df <- do.call(rbind, parts)
    rownames(count_df) <- NULL
    # Preserve factor ordering for downstream facet order
    count_df$group <- factor(count_df$group, levels = levels(group_factor))
    agg_keys <- c("item", "group")
  }

  # Per-(item [x group]) totals -> percentages
  totals <- stats::aggregate(
    count_df$n,
    by  = count_df[, agg_keys, drop = FALSE],
    FUN = sum
  )
  colnames(totals)[ncol(totals)] <- "total"
  count_df <- merge(count_df, totals, by = agg_keys, sort = FALSE)
  count_df$percentage <- round(count_df$n / count_df$total * 100, 1)

  # Item labels with column-order ordering (top item = first column)
  count_df$item_label <- label_map[count_df$item]
  count_df$item_label <- factor(
    count_df$item_label,
    levels = rev(label_map[item_names])
  )

  # ---------------------------------------------------------------------
  # Output: dataframe
  # ---------------------------------------------------------------------
  if (output == "dataframe") {
    keep <- c(if (!is.null(group_factor)) "group" else NULL,
              "item", "item_label", "category", "n", "percentage")
    out_df <- count_df[, keep, drop = FALSE]
    rownames(out_df) <- NULL
    return(out_df)
  }

  # ---------------------------------------------------------------------
  # Output: ggplot
  # ---------------------------------------------------------------------
  p <- ggplot2::ggplot(
    count_df,
    ggplot2::aes(x = .data$category,
                 y = .data$item_label,
                 fill = .data$n)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(
      expression(italic(n)),
      limits = c(0, NA)
    ) +
    ggplot2::labs(y = "Items") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8),
      panel.grid  = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12)),
      panel.spacing = ggplot2::unit(0.7, "cm")
    )

  # X-axis: numeric breaks or category labels
  if (!is.null(category_labels)) {
    p <- p + ggplot2::scale_x_continuous(
      "Response category",
      expand = c(0, 0),
      breaks = all_categories,
      labels = category_labels
    )
  } else {
    p <- p + ggplot2::scale_x_continuous(
      "Response category",
      expand = c(0, 0),
      breaks = all_categories
    )
  }

  # Cell labels (count or percent), with optional <-cutoff highlighting
  if (percent) {
    if (highlight) {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(
            label = paste0(.data$percentage, "%"),
            color = ifelse(.data$n < cutoff, "red", text_color)
          )
        ) +
        ggplot2::guides(color = "none") +
        ggplot2::scale_color_identity()
    } else {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(label = paste0(.data$percentage, "%")),
          color = text_color
        )
    }
  } else {
    if (highlight) {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(
            label = .data$n,
            color = ifelse(.data$n < cutoff, "red", text_color)
          )
        ) +
        ggplot2::guides(color = "none") +
        ggplot2::scale_color_identity()
    } else {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(label = .data$n),
          color = text_color
        )
    }
  }

  # Faceting
  if (!is.null(group_factor)) {
    p <- p + ggplot2::facet_wrap(~ group, ncol = facet_ncol)
  }

  p
}
