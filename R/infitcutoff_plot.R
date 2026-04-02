#' Plot Distribution of Simulated Infit and Outfit MSQ Values
#'
#' Visualises the distribution of simulation-based conditional item fit values
#' from \code{\link{RMinfitcutoff}}, optionally overlaying observed item fit
#' from the original data.
#'
#' Uses `ggdist::stat_dotsinterval()` (when `data` is not supplied) or
#' `ggdist::stat_dots()` (when `data` is supplied) with
#' `point_interval = "median_hdci"` and `.width = c(0.66, 0.999)`.
#'
#' @param simfit The return value of \code{\link{RMinfitcutoff}} (a list with
#'   components `results`, `item_cutoffs`, `actual_iterations`, `sample_n`, and
#'   `item_names`).
#' @param data Optional. A data.frame or matrix of item responses for computing
#'   and overlaying observed conditional item fit values. Items must be scored
#'   starting at 0 (non-negative integers). When provided, the plot includes
#'   orange diamond markers for the observed infit/outfit MSQ alongside the
#'   simulated distribution, plus segment summaries from the cutoff intervals.
#' @param output Character string. Either `"infit"` (default) to show only the
#'   infit panel, `"outfit"` to show only the outfit panel, or `"both"` to show
#'   infit and outfit side by side (requires the `patchwork` package when
#'   `data` is supplied).
#'
#' @return A `ggplot` object (or a `patchwork` object when `output = "both"`
#'   and `data` is supplied).
#'
#' @details
#' When `data` is **not** supplied, the function plots the simulated MSQ
#' distributions as dot-interval plots using `ggdist::stat_dotsinterval()` with
#' median and Highest Density Continuous Interval (HDCI) summaries, faceted by
#' statistic (InfitMSQ / OutfitMSQ).
#'
#' When `data` **is** supplied, the function:
#' \enumerate{
#'   \item Fits a Rasch model (`eRm::RM()` for dichotomous data or
#'     `eRm::PCM()` for polytomous data) and computes observed conditional
#'     infit and outfit MSQ via `iarm::out_infit()`.
#'   \item Overlays observed fit values as orange diamond markers on the
#'     simulated distributions.
#'   \item Shows per-item cutoff intervals (from `simfit$item_cutoffs`) as
#'     black line segments, with thicker segments for the 66% HDCI range and
#'     black dots for the median.
#' }
#'
#' The `ggplot2`, `ggdist`, and optionally `patchwork` packages must be
#' installed (they are in Suggests, not Imports).
#'
#' @seealso \code{\link{RMinfitcutoff}}, \code{\link{RMiteminfit}}
#' 
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#' )
#' colnames(sim_data) <- paste0("Item", 1:10)
#'
#' # Run simulation
#' cutoff_res <- RMinfitcutoff(sim_data, iterations = 100, parallel = FALSE,
#'                             seed = 42)
#'
#' # Simulated distribution only (infit + outfit faceted)
#' RMinfitcutoffPlot(cutoff_res)
#'
#' # With observed fit overlaid (infit only, the default)
#' RMinfitcutoffPlot(cutoff_res, data = sim_data)
#'
#' # Both infit and outfit panels side by side
#' RMinfitcutoffPlot(cutoff_res, data = sim_data, output = "both")
#' }
RMinfitcutoffPlot <- function(simfit, data, output = "infit") {
  
  # --- Check required packages ------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for RMinfitcutoffPlot() but is not installed.\n",
      "Install it with: install.packages(\"ggplot2\")",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop(
      "Package 'ggdist' is required for RMinfitcutoffPlot() but is not installed.\n",
      "Install it with: install.packages(\"ggdist\")",
      call. = FALSE
    )
  }
  
  output <- match.arg(output, c("infit", "outfit", "both"))
  
  # --- Validate simfit --------------------------------------------------------
  required_names <- c("results", "item_cutoffs", "actual_iterations",
                      "sample_n", "item_names")
  missing_names <- setdiff(required_names, names(simfit))
  if (length(missing_names) > 0L) {
    stop(
      "`simfit` is missing required components: ",
      paste(missing_names, collapse = ", "),
      ".\nExpected the return value of RMinfitcutoff().",
      call. = FALSE
    )
  }
  
  results_df       <- simfit$results
  item_cutoffs     <- simfit$item_cutoffs
  actual_iterations <- simfit$actual_iterations
  sample_n         <- simfit$sample_n
  item_names       <- simfit$item_names
  
  # Item factor levels (reversed for plotting top-to-bottom)
  item_levels <- rev(item_names)
  
  # --- Compute per-item summary intervals for segment overlays ----------------
  lo_hi <- do.call(rbind, lapply(item_names, function(item) {
    sub <- results_df[results_df$Item == item, ]
    data.frame(
      Item              = item,
      min_infit_msq     = stats::quantile(sub$InfitMSQ,  0.001, na.rm = TRUE),
      max_infit_msq     = stats::quantile(sub$InfitMSQ,  0.999, na.rm = TRUE),
      p66lo_infit_msq   = stats::quantile(sub$InfitMSQ,  0.167, na.rm = TRUE),
      p66hi_infit_msq   = stats::quantile(sub$InfitMSQ,  0.833, na.rm = TRUE),
      median_infit      = stats::median(sub$InfitMSQ, na.rm = TRUE),
      min_outfit_msq    = stats::quantile(sub$OutfitMSQ, 0.001, na.rm = TRUE),
      max_outfit_msq    = stats::quantile(sub$OutfitMSQ, 0.999, na.rm = TRUE),
      p66lo_outfit_msq  = stats::quantile(sub$OutfitMSQ, 0.167, na.rm = TRUE),
      p66hi_outfit_msq  = stats::quantile(sub$OutfitMSQ, 0.833, na.rm = TRUE),
      median_outfit     = stats::median(sub$OutfitMSQ, na.rm = TRUE),
      stringsAsFactors  = FALSE,
      row.names         = NULL
    )
  }))
  rownames(lo_hi) <- NULL
  
  # --- Case 1: no observed data, show simulation distribution only ------------
  if (missing(data)) {
    
    # Pivot results to long format
    infit_long <- data.frame(
      Item      = results_df$Item,
      statistic = "InfitMSQ",
      Value     = results_df$InfitMSQ,
      stringsAsFactors = FALSE
    )
    outfit_long <- data.frame(
      Item      = results_df$Item,
      statistic = "OutfitMSQ",
      Value     = results_df$OutfitMSQ,
      stringsAsFactors = FALSE
    )
    results_long <- rbind(infit_long, outfit_long)
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
        point_interval = "median_hdci",
        layout = "weave",
        slab_color = NA,
        .width = c(0.66, 0.999)
      ) +
      ggplot2::labs(
        x = "Conditional MSQ",
        y = "Item",
        caption = paste0(
          "Note: Results from ", actual_iterations,
          " simulated datasets with ", sample_n, " respondents."
        )
      ) +
      ggplot2::scale_color_manual(
        values = scales::brewer_pal()(3)[-1],
        aesthetics = "slab_fill",
        guide = "none"
      ) +
      ggplot2::facet_wrap(~statistic, ncol = 2) +
      ggplot2::scale_x_continuous(
        breaks = seq(0.5, 1.5, 0.1),
        minor_breaks = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"))
    
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
  
  data_mat <- as.matrix(data)
  
  if (max(data_mat, na.rm = TRUE) == 1L) {
    erm_out <- eRm::RM(data)
  } else {
    erm_out <- eRm::PCM(data)
  }
  
  # rgl workaround
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)
  
  cfit <- iarm::out_infit(erm_out)
  
  observed_df <- data.frame(
    Item      = names(data),
    observed_infit  = cfit$Infit,
    observed_outfit = cfit$Outfit,
    stringsAsFactors = FALSE
  )
  
  # --- Build infit data -------------------------------------------------------
  infit_sim <- data.frame(
    Item      = results_df$Item,
    Value     = results_df$InfitMSQ,
    stringsAsFactors = FALSE
  )
  infit_sim <- merge(infit_sim, observed_df[, c("Item", "observed_infit")],
                     by = "Item", sort = FALSE)
  infit_sim$Item <- factor(infit_sim$Item, levels = item_levels)
  
  lo_hi$Item_f <- factor(lo_hi$Item, levels = item_levels)
  
  caption_text <- paste0(
    "Note: Results from ", actual_iterations,
    " simulated datasets with ", sample_n, " respondents.\n",
    "Orange dots indicate observed conditional item fit. ",
    "Black dots indicate median fit from simulations."
  )
  
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
      layout = "weave",
      slab_color = NA,
      .width = c(0.666, 0.999)
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
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"))
  
  # --- Build outfit data ------------------------------------------------------
  outfit_sim <- data.frame(
    Item      = results_df$Item,
    Value     = results_df$OutfitMSQ,
    stringsAsFactors = FALSE
  )
  outfit_sim <- merge(outfit_sim, observed_df[, c("Item", "observed_outfit")],
                      by = "Item", sort = FALSE)
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
      layout = "weave",
      slab_color = NA,
      .width = c(0.666, 0.999)
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
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"))
  
  # --- Return the requested panel(s) ------------------------------------------
  if (output == "both") {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop(
        "Package 'patchwork' is required for output = \"both\" but is not installed.\n",
        "Install it with: install.packages(\"patchwork\")",
        call. = FALSE
      )
    }
    return(infit_p + outfit_p)
  } else if (output == "outfit") {
    return(outfit_p)
  } else {
    # default: infit only, add caption to infit panel
    infit_p <- infit_p +
      ggplot2::labs(caption = caption_text)
    return(infit_p)
  }
}