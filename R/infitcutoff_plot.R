#' Plot the simulated infit and outfit MSQ values
#'
#' Plots the simulated conditional item fit values from
#' [RMitemInfitCutoff()] (or [RMitemInfitCutoffMI()]) per item, and
#' optionally the observed values in your data.
#'
#' @param simfit The result of [RMitemInfitCutoff()] or
#'   [RMitemInfitCutoffMI()].
#' @param data Optional. The item responses the simulations were made from,
#'   items scored from 0. The observed conditional infit and outfit MSQ are
#'   then computed (on complete cases) and drawn as orange diamonds.
#' @param statistic Character. `"infit"`, `"outfit"`, or `"both"` side
#'   by side. With `data` the default is `"infit"`, and `"both"` needs
#'   \pkg{patchwork}. Without `data` the default is `"both"`, drawn as two
#'   facets.
#'
#' @return A `ggplot` object, or a `patchwork` object with
#'   `statistic = "both"` and `data`.
#'
#' @details
#' Each item's simulated values are drawn as dots, with the median, an
#' inner 66% interval and an outer interval. The outer interval is the one
#' stored in `simfit` and tabulated by [RMitemInfit()]: the HDCI of width
#' `hdci_width`, or the 2.5th to 97.5th percentiles with
#' `cutoff_method = "quantile"`. The inner interval uses the same method.
#'
#' **Without `data`** the intervals are drawn with
#' `ggdist::stat_dotsinterval()`.
#'
#' **With `data`**, the model is fitted by CML (`psychotools::pcmodel()`)
#' and the observed conditional MSQ computed with `iarm::out_infit()`. The
#' intervals are black segments with a black dot for the median, and the
#' observed value is the orange diamond.
#'
#' The intervals are per item and not corrected for testing several items.
#' With the full cutoff object, [RMitemInfit()] flags items by adjusted
#' p-values instead, so an item just outside its interval need not be
#' flagged there. The plot shows no p-values.
#'
#' The \pkg{ggplot2} and \pkg{ggdist} packages must be installed, and
#' \pkg{patchwork} for `statistic = "both"` (all are in Suggests).
#'
#' @seealso [RMitemInfitCutoff()], [RMitemInfit()]
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
#'   # Run simulation
#'   cutoff_res <- RMitemInfitCutoff(sim_data, iterations = 100,
#'                                   parallel = FALSE, seed = 42)
#'
#'   # Simulated values only (infit and outfit side by side)
#'   RMitemInfitPlot(cutoff_res)
#'
#'   # With observed fit overlaid (infit only, the default)
#'   RMitemInfitPlot(cutoff_res, data = sim_data)
#'
#'   # Both infit and outfit panels side by side
#'   if (requireNamespace("patchwork", quietly = TRUE)) {
#'     RMitemInfitPlot(cutoff_res, data = sim_data, statistic = "both")
#'   }
#' }
#' }
RMitemInfitPlot <- function(simfit, data, statistic = "infit") {
  # --- Check required packages ------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for RMitemInfitPlot() but is not installed.\n",
      "Install it with: install.packages(\"ggplot2\")",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required for RMitemInfitPlot() but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")",
      call. = FALSE
    )
  }

  # Without `data` the default shows both panels, as before `statistic` was
  # honoured there. An explicit `statistic` is followed in both cases.
  statistic <- if (missing(statistic) && missing(data)) {
    "both"
  } else {
    match.arg(statistic, c("infit", "outfit", "both"))
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
  if (length(missing_names) > 0L) {
    stop(
      "`simfit` is missing required components: ",
      paste(missing_names, collapse = ", "),
      ".\nExpected the return value of RMitemInfitCutoff().",
      call. = FALSE
    )
  }

  results_df <- simfit$results
  item_cutoffs <- simfit$item_cutoffs
  actual_iterations <- simfit$actual_iterations
  sample_n <- simfit$sample_n
  item_names <- simfit$item_names

  # Standard sample-size clause for the (complete-case) simulation sample.
  # `sample_n_total` / `sample_has_na` are absent in cutoff objects made by
  # older versions, so fall back to the plain count.
  sample_clause <- .n_caption(
    sample_n,
    if (is.null(simfit$sample_n_total)) sample_n else simfit$sample_n_total,
    if (isTRUE(simfit$sample_has_na)) "complete cases" else character()
  )

  # Item factor levels (reversed for plotting top-to-bottom)
  item_levels <- rev(item_names)

  # Intervals. The outer segment is the interval stored in `simfit` (the
  # one RMitemInfit() tabulates) and the inner 66% segment uses the same
  # method, so the plot and the table describe the same interval.
  spec <- .sim_interval_spec(simfit)
  outer_width <- spec$width
  cut_idx <- match(item_names, item_cutoffs$Item)
  iv_infit <- .sim_plot_intervals(
    results_df$InfitMSQ, results_df$Item, item_names,
    item_cutoffs$infit_low[cut_idx], item_cutoffs$infit_high[cut_idx],
    spec$method
  )
  iv_outfit <- .sim_plot_intervals(
    results_df$OutfitMSQ, results_df$Item, item_names,
    item_cutoffs$outfit_low[cut_idx], item_cutoffs$outfit_high[cut_idx],
    spec$method
  )
  lo_hi <- data.frame(
    Item = item_names,
    min_infit_msq = iv_infit$lo,
    max_infit_msq = iv_infit$hi,
    p66lo_infit_msq = iv_infit$lo66,
    p66hi_infit_msq = iv_infit$hi66,
    median_infit = iv_infit$median,
    min_outfit_msq = iv_outfit$lo,
    max_outfit_msq = iv_outfit$hi,
    p66lo_outfit_msq = iv_outfit$lo66,
    p66hi_outfit_msq = iv_outfit$hi66,
    median_outfit = iv_outfit$median,
    stringsAsFactors = FALSE
  )
  interval_note <- .sim_interval_caption(spec, "item", "RMitemInfit()")

  # --- Case 1: no observed data, show simulation distribution only ------------
  if (missing(data)) {
    # Pivot results to long format
    stat_cols <- switch(
      statistic,
      infit = "InfitMSQ",
      outfit = "OutfitMSQ",
      both = c("InfitMSQ", "OutfitMSQ")
    )
    results_long <- do.call(rbind, lapply(stat_cols, function(sc) {
      data.frame(
        Item = results_df$Item,
        statistic = sc,
        Value = results_df[[sc]],
        stringsAsFactors = FALSE
      )
    }))
    results_long$Item <- factor(results_long$Item, levels = item_levels)

    p <- ggplot2::ggplot(
      results_long,
      ggplot2::aes(
        x = .data$Value,
        y = .data$Item
      )
    ) +
      ggdist::stat_dotsinterval(
        ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
        quantiles = actual_iterations,
        point_interval = .sim_point_interval(spec),
        layout = "weave",
        slab_color = NA,
        .width = c(0.66, outer_width)
      ) +
      ggplot2::labs(
        x = "Conditional MSQ",
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
      ggplot2::scale_color_manual(
        values = scales::brewer_pal()(3)[-1],
        aesthetics = "slab_fill",
        guide = "none"
      ) +
      ggplot2::facet_wrap(~statistic, ncol = length(stat_cols)) +
      ggplot2::scale_x_continuous(
        breaks = seq(0.5, 1.5, 0.1),
        minor_breaks = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
      er2_axis_margins() +
      er2_plot_caption()

    return(p)
  }

  # --- Case 2: observed data supplied -----------------------------------------
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop(
      "Package 'iarm' is required to compute observed item fit but is not installed.\n",
      "Install it with: install.packages(\"iarm\")",
      call. = FALSE
    )
  }
  validate_response_data(data)

  # Respondents with no responses at all break the CML fit below
  # (psychotools errors on all-NA rows); drop them, as in RMitemInfit().
  data <- .drop_empty_respondents(data)
  # iarm's conditional infit refits on complete cases, so check them too.
  .check_complete_case_categories(data)
  .check_cutoff_sample(
    simfit$sample_n,
    nrow(stats::na.omit(data)),
    "RMitemInfitCutoff()"
  )

  # rgl workaround
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  # CML fit via psychotools (a dichotomous item is a 2-category PCM);
  # iarm::out_infit() accepts pcmodel objects and matches the former
  # eRm::RM()/PCM() route to ~1e-6 (same pairing as RMitemInfitCutoff()).
  model_fit <- .cml_fit(data) # shared category check, then pcmodel()
  cfit <- iarm::out_infit(model_fit)

  observed_df <- data.frame(
    Item = names(data),
    observed_infit = cfit$Infit,
    observed_outfit = cfit$Outfit,
    stringsAsFactors = FALSE
  )

  # --- Build infit data -------------------------------------------------------
  infit_sim <- data.frame(
    Item = results_df$Item,
    Value = results_df$InfitMSQ,
    stringsAsFactors = FALSE
  )
  infit_sim <- merge(
    infit_sim,
    observed_df[, c("Item", "observed_infit")],
    by = "Item",
    sort = FALSE
  )
  infit_sim$Item <- factor(infit_sim$Item, levels = item_levels)

  lo_hi$Item_f <- factor(lo_hi$Item, levels = item_levels)

  caption_text <- er2_caption(paste0(
    "Results from ",
    actual_iterations,
    " simulated datasets. ",
    sample_clause,
    " per dataset.\n",
    "Orange diamonds indicate observed conditional item fit. ",
    "Black dots indicate median fit from simulations.\n",
    interval_note
  ))

  infit_p <- ggplot2::ggplot(
    infit_sim,
    ggplot2::aes(
      x = .data$Value,
      y = .data$Item
    )
  ) +
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
        x = .data$min_infit_msq,
        xend = .data$max_infit_msq,
        y = .data$Item_f,
        yend = .data$Item_f
      ),
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x = .data$p66lo_infit_msq,
        xend = .data$p66hi_infit_msq,
        y = .data$Item_f,
        yend = .data$Item_f
      ),
      color = "black",
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = lo_hi,
      ggplot2::aes(
        x = .data$median_infit,
        y = .data$Item_f
      ),
      size = 3.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$observed_infit),
      color = "sienna2",
      shape = 18,
      position = ggplot2::position_nudge(y = -0.1),
      size = 4
    ) +
    ggplot2::labs(
      x = "Conditional Infit MSQ",
      y = "Item"
    ) +
    ggplot2::scale_color_manual(
      values = scales::brewer_pal()(3)[-1],
      aesthetics = "slab_fill",
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0.5, 1.5, 0.1),
      minor_breaks = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
    er2_axis_margins() +
    er2_plot_caption()

  # --- Build outfit data ------------------------------------------------------
  outfit_sim <- data.frame(
    Item = results_df$Item,
    Value = results_df$OutfitMSQ,
    stringsAsFactors = FALSE
  )
  outfit_sim <- merge(
    outfit_sim,
    observed_df[, c("Item", "observed_outfit")],
    by = "Item",
    sort = FALSE
  )
  outfit_sim$Item <- factor(outfit_sim$Item, levels = item_levels)

  outfit_p <- ggplot2::ggplot(
    outfit_sim,
    ggplot2::aes(
      x = .data$Value,
      y = .data$Item
    )
  ) +
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
        x = .data$min_outfit_msq,
        xend = .data$max_outfit_msq,
        y = .data$Item_f,
        yend = .data$Item_f
      ),
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x = .data$p66lo_outfit_msq,
        xend = .data$p66hi_outfit_msq,
        y = .data$Item_f,
        yend = .data$Item_f
      ),
      color = "black",
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = lo_hi,
      ggplot2::aes(
        x = .data$median_outfit,
        y = .data$Item_f
      ),
      size = 3.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$observed_outfit),
      color = "sienna2",
      shape = 18,
      position = ggplot2::position_nudge(y = -0.1),
      size = 4
    ) +
    ggplot2::labs(
      x = "Conditional Outfit MSQ",
      y = "Item",
      caption = caption_text
    ) +
    ggplot2::scale_color_manual(
      values = scales::brewer_pal()(3)[-1],
      aesthetics = "slab_fill",
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0.5, 1.5, 0.1),
      minor_breaks = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
    er2_axis_margins() +
    er2_plot_caption()

  # --- Return the requested panel(s) ------------------------------------------
  if (statistic == "both") {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop(
        "Package 'patchwork' is required for statistic = \"both\" but is not installed.\n",
        "Install it with: install.packages(\"patchwork\")",
        call. = FALSE
      )
    }
    return(infit_p + outfit_p)
  } else if (statistic == "outfit") {
    return(outfit_p)
  } else {
    # default: infit only, add caption to infit panel
    infit_p <- infit_p +
      ggplot2::labs(caption = caption_text)
    return(infit_p)
  }
}

#' @rdname RMitemInfitPlot
#' @param ... Arguments passed on to \code{RMitemInfitPlot()}.
#' @details
#' `RMitemInfitCutoffPlot()` is a deprecated alias for `RMitemInfitPlot()`,
#' kept for code written for easyRasch2 0.8.0. It warns and calls
#' `RMitemInfitPlot()`.
#' @export
RMitemInfitCutoffPlot <- function(...) {
  .Deprecated("RMitemInfitPlot")
  RMitemInfitPlot(...)
}
