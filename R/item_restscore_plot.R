#' Plot the simulated item-restscore differences
#'
#' Plots each item's simulated observed minus expected item-restscore gamma
#' from [RMitemRestscoreCutoff()], and optionally the observed differences
#' in your data.
#'
#' @param simfit The result of [RMitemRestscoreCutoff()].
#' @param data Optional. The item responses the simulations were made from,
#'   items scored from 0. The observed differences are then computed with
#'   [RMitemRestscore()] and drawn as orange diamonds.
#'
#' @return A `ggplot` object.
#'
#' @details
#' The x-axis is the observed minus expected gamma, the statistic
#' [RMitemRestscore()] tests, with a dashed line at zero. Positive values
#' indicate overfit, negative values underfit. In small samples the
#' simulated differences are generally not centred on zero, which is one
#' reason the asymptotic test is miscalibrated (see
#' [RMitemRestscoreCutoff()]).
#'
#' Each item's simulated values are drawn as dots, with the median, an
#' inner 66% interval and an outer interval. The outer interval is the one
#' stored in `simfit` and tabulated by [RMitemRestscore()]: the HDCI of
#' width `hdci_width`, or the 2.5th to 97.5th percentiles with
#' `cutoff_method = "quantile"`. The inner interval uses the same method.
#' Without `data` the intervals are drawn with
#' `ggdist::stat_dotsinterval()`. With `data` they are black segments with
#' a black dot for the median, and the observed value is the orange
#' diamond.
#'
#' The intervals are per item and not corrected for testing several items.
#' With the full cutoff object, [RMitemRestscore()] flags items by adjusted
#' p-values instead, so an item just outside its interval need not be
#' flagged there. The plot shows no p-values.
#'
#' The \pkg{ggplot2} and \pkg{ggdist} packages must be installed (they are
#' in Suggests).
#'
#' @seealso [RMitemRestscoreCutoff()], [RMitemRestscore()]
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("iarm", quietly = TRUE) &&
#'     requireNamespace("ggdist", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:10)
#'
#'   cutoff_res <- RMitemRestscoreCutoff(sim_data, iterations = 100,
#'                                       parallel = FALSE, seed = 42)
#'
#'   # Simulated differences only
#'   RMitemRestscorePlot(cutoff_res)
#'
#'   # With the observed differences overlaid
#'   RMitemRestscorePlot(cutoff_res, data = sim_data)
#' }
#' }
RMitemRestscorePlot <- function(simfit, data) {
  # --- Check required packages ------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for RMitemRestscorePlot() but is not installed.\n",
      "Install it with: install.packages(\"ggplot2\")",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required for RMitemRestscorePlot() but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")",
      call. = FALSE
    )
  }

  # --- Validate simfit --------------------------------------------------------
  required_names <- c(
    "results",
    "item_cutoffs",
    "actual_iterations",
    "sample_n",
    "item_names"
  )
  missing_names <- setdiff(required_names, names(simfit))
  if (length(missing_names) > 0L || !"Difference" %in% names(simfit$results)) {
    stop(
      "`simfit` is not an RMitemRestscoreCutoff() result",
      if (length(missing_names) > 0L) {
        paste0(" (missing: ", paste(missing_names, collapse = ", "), ")")
      },
      ".",
      call. = FALSE
    )
  }

  results_df <- simfit$results
  actual_iterations <- simfit$actual_iterations
  sample_n <- simfit$sample_n
  item_names <- simfit$item_names

  sample_clause <- .n_caption(
    sample_n,
    if (is.null(simfit$sample_n_total)) sample_n else simfit$sample_n_total,
    if (isTRUE(simfit$sample_has_na)) "complete cases" else character()
  )

  # Item factor levels (reversed for plotting top-to-bottom)
  item_levels <- rev(item_names)

  # Intervals. The outer interval is the one stored in `simfit` (the one
  # RMitemRestscore() tabulates) and the inner 66% interval uses the same
  # method, so the plot and the table describe the same interval.
  spec <- .sim_interval_spec(simfit)
  outer_width <- spec$width
  interval_note <- .sim_interval_caption(spec, "item", "RMitemRestscore()")

  sim_df <- data.frame(
    Item = factor(results_df$Item, levels = item_levels),
    Value = results_df$Difference,
    stringsAsFactors = FALSE
  )
  sim_df <- sim_df[is.finite(sim_df$Value), ]

  x_lab <- "Observed minus expected item-restscore gamma"
  zero_line <- ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  )
  fill_scale <- ggplot2::scale_color_manual(
    values = scales::brewer_pal()(3)[-1],
    aesthetics = "slab_fill",
    guide = "none"
  )

  # --- Case 1: no observed data, show simulation distribution only ------------
  if (missing(data)) {
    p <- ggplot2::ggplot(
      sim_df,
      ggplot2::aes(x = .data$Value, y = .data$Item)
    ) +
      zero_line +
      ggdist::stat_dotsinterval(
        ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
        quantiles = actual_iterations,
        point_interval = .sim_point_interval(spec),
        layout = "weave",
        slab_color = NA,
        .width = c(0.66, outer_width)
      ) +
      ggplot2::labs(
        x = x_lab,
        y = "Item",
        caption = er2_caption(paste0(
          "Results from ",
          actual_iterations,
          " simulated datasets. ",
          sample_clause,
          " per dataset.\n",
          interval_note
        ))
      ) +
      fill_scale +
      ggplot2::theme_minimal() +
      er2_axis_margins() +
      er2_plot_caption()

    return(p)
  }

  # --- Case 2: observed data supplied -----------------------------------------
  validate_response_data(data)
  .check_cutoff_sample(
    sample_n,
    nrow(stats::na.omit(as.data.frame(data))),
    "RMitemRestscoreCutoff()"
  )
  observed_df <- RMitemRestscore(data, output = "dataframe")
  if (!setequal(observed_df$Item, item_names)) {
    stop(
      "Item names in `data` do not match the items in `simfit`.",
      call. = FALSE
    )
  }
  observed_df <- data.frame(
    Item = factor(observed_df$Item, levels = item_levels),
    observed = observed_df$Difference,
    stringsAsFactors = FALSE
  )

  cut_idx <- match(item_names, simfit$item_cutoffs$Item)
  lo_hi <- .sim_plot_intervals(
    sim_df$Value, as.character(sim_df$Item), item_names,
    simfit$item_cutoffs$diff_low[cut_idx],
    simfit$item_cutoffs$diff_high[cut_idx],
    spec$method
  )
  names(lo_hi)[names(lo_hi) == "key"] <- "Item"
  names(lo_hi)[names(lo_hi) == "lo66"] <- "p66lo"
  names(lo_hi)[names(lo_hi) == "hi66"] <- "p66hi"
  lo_hi$Item_f <- factor(lo_hi$Item, levels = item_levels)

  ggplot2::ggplot(
    sim_df,
    ggplot2::aes(x = .data$Value, y = .data$Item)
  ) +
    zero_line +
    ggdist::stat_dots(
      ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
      quantiles = actual_iterations,
      point_interval = .sim_point_interval(spec),
      layout = "weave",
      slab_color = NA,
      .width = c(0.66, outer_width)
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x = .data$lo,
        xend = .data$hi,
        y = .data$Item_f,
        yend = .data$Item_f
      ),
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x = .data$p66lo,
        xend = .data$p66hi,
        y = .data$Item_f,
        yend = .data$Item_f
      ),
      color = "black",
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = lo_hi,
      ggplot2::aes(x = .data$median, y = .data$Item_f),
      size = 3.6
    ) +
    ggplot2::geom_point(
      data = observed_df,
      ggplot2::aes(x = .data$observed, y = .data$Item),
      color = "sienna2",
      shape = 18,
      position = ggplot2::position_nudge(y = -0.1),
      size = 4
    ) +
    ggplot2::labs(
      x = x_lab,
      y = "Item",
      caption = er2_caption(paste0(
        "Results from ",
        actual_iterations,
        " simulated datasets. ",
        sample_clause,
        " per dataset.\n",
        "Orange diamonds indicate the observed difference. ",
        "Black dots indicate the median from simulations.\n",
        interval_note
      ))
    ) +
    fill_scale +
    ggplot2::theme_minimal() +
    er2_axis_margins() +
    er2_plot_caption()
}
