#' Plot the Simulated Item-Restscore Null Distribution
#'
#' Visualises the per-item null distribution of the observed minus expected
#' item-restscore gamma from \code{\link{RMitemRestscoreCutoff}}, optionally
#' overlaying the observed differences from the original data.
#'
#' Uses `ggdist::stat_dotsinterval()` (when `data` is not supplied) or
#' `ggdist::stat_dots()` (when `data` is supplied) with
#' `point_interval = "median_hdci"`. The outer `.width` follows
#' `simfit$hdci_width`, so the shaded interval matches the one
#' [RMitemRestscore()] tabulates.
#'
#' @param simfit The return value of \code{\link{RMitemRestscoreCutoff}} (a
#'   list with components `results`, `item_cutoffs`, `actual_iterations`,
#'   `sample_n`, and `item_names`).
#' @param data Optional. A data.frame or matrix of item responses for computing
#'   and overlaying the observed item-restscore differences. Items must be
#'   scored starting at 0 (non-negative integers). When provided, the plot
#'   includes orange diamond markers for the observed difference alongside the
#'   simulated distribution, plus segment summaries of the intervals.
#'
#' @return A `ggplot` object.
#'
#' @details
#' The x-axis is the difference between observed and model-expected gamma, the
#' statistic \code{\link{RMitemRestscore}} tests. A dashed line marks zero.
#' The simulated distributions are generally not centred on zero in small
#' samples, which is one of the reasons the asymptotic test is miscalibrated
#' (see \code{\link{RMitemRestscoreCutoff}}). Positive values indicate
#' over-discrimination (overfit), negative values under-discrimination
#' (underfit).
#'
#' When `data` **is** supplied, the observed differences are computed with
#' \code{\link{RMitemRestscore}} and overlaid as orange diamonds, with per-item
#' intervals drawn as black line segments (thicker for the 66% range) and
#' black dots for the simulated median.
#'
#' The `ggplot2` and `ggdist` packages must be installed (they are in
#' Suggests, not Imports).
#'
#' @seealso \code{\link{RMitemRestscoreCutoff}}, \code{\link{RMitemRestscore}}
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
#'   # Simulated distribution only
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

  # Outer interval width follows the cutoff object, so the plot and the table
  # describe the same interval. `cutoff_method = "quantile"` fixes it at the
  # 2.5th/97.5th percentiles.
  outer_width <- if (
    identical(simfit$cutoff_method, "hdci") && !is.null(simfit$hdci_width)
  ) {
    simfit$hdci_width
  } else {
    0.95
  }
  outer_lo <- (1 - outer_width) / 2
  outer_hi <- 1 - outer_lo

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
        point_interval = "median_hdci",
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
          " per dataset."
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

  lo_hi <- do.call(
    rbind,
    lapply(item_names, function(item) {
      v <- sim_df$Value[sim_df$Item == item]
      data.frame(
        Item = item,
        lo = stats::quantile(v, outer_lo, names = FALSE),
        hi = stats::quantile(v, outer_hi, names = FALSE),
        p66lo = stats::quantile(v, 0.167, names = FALSE),
        p66hi = stats::quantile(v, 0.833, names = FALSE),
        median = stats::median(v),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    })
  )
  lo_hi$Item_f <- factor(lo_hi$Item, levels = item_levels)

  ggplot2::ggplot(
    sim_df,
    ggplot2::aes(x = .data$Value, y = .data$Item)
  ) +
    zero_line +
    ggdist::stat_dots(
      ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
      quantiles = actual_iterations,
      layout = "weave",
      slab_color = NA,
      .width = c(0.666, outer_width)
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
        "Orange dots indicate the observed difference. ",
        "Black dots indicate the median from simulations."
      ))
    ) +
    fill_scale +
    ggplot2::theme_minimal() +
    er2_axis_margins() +
    er2_plot_caption()
}
