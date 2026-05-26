#' Plot Distribution of Simulated Q3 Residual Correlations
#'
#' Visualises the distribution of simulation-based Yen's Q3 residual
#' correlations per item pair from \code{\link{RMlocdepQ3Cutoff}},
#' optionally overlaying observed Q3 values computed from real data via
#' \code{mirt::residuals(..., type = "Q3")}.
#'
#' Uses `ggdist::stat_dotsinterval()` (when `data` is not supplied) or
#' `ggdist::stat_dots()` (when `data` is supplied) with
#' `point_interval = "median_hdci"` and `.width = c(0.66, 0.95, 0.99)`.
#'
#' @param simfit The return value of \code{\link{RMlocdepQ3Cutoff}} (a list
#'   with components `pair_results`, `pair_cutoffs`, `actual_iterations`,
#'   `sample_n`, and `item_names`).
#' @param data Optional. A data.frame or matrix of item responses for
#'   computing and overlaying observed Q3 values. Items must be scored
#'   starting at 0 (non-negative integers). When provided, the plot
#'   includes orange diamond markers for the observed Q3 alongside the
#'   simulated distribution, plus segment summaries from the cutoff
#'   intervals.
#' @param items Optional character vector of item names to include in the
#'   plot. Only item pairs where **both** items are in this vector will be
#'   shown. When `NULL` (default), all item pairs are plotted.
#' @param n_pairs Optional positive integer. When supplied, only the
#'   `n_pairs` item pairs with the largest **deviation** from the simulated
#'   null are plotted, sorted by `|observed Q3 - median(simulated Q3 per
#'   pair)|` descending when `data` is supplied, or by
#'   `|median(simulated Q3 per pair)|` otherwise. Applied *after* the
#'   `items` filter when both are supplied. Values larger than the number
#'   of available pairs are silently capped.
#'
#' @return A `ggplot` object.
#'
#' @details
#' The plot shows one row per item pair (labelled as "Item1 - Item2"). Only
#' the upper triangle of the Q3 matrix is plotted (pairs are unordered
#' under symmetric Q3, unlike partial gamma which is direction-dependent).
#'
#' When `data` is **not** supplied, the function plots the simulated Q3
#' distributions as dot-interval plots using `ggdist::stat_dotsinterval()`
#' with median and Highest Density Continuous Interval (HDCI) summaries.
#'
#' When `data` **is** supplied, the function:
#' \enumerate{
#'   \item Fits a Rasch model to `data` via `mirt::mirt()` and extracts
#'     observed Q3 residual correlations.
#'   \item Overlays observed Q3 values as orange diamond markers on the
#'     simulated distributions.
#'   \item Shows per-pair cutoff intervals (from `simfit$pair_cutoffs`)
#'     as black line segments, with thicker segments for the 66\%
#'     interval and black dots for the median.
#' }
#'
#' The `ggplot2`, `ggdist`, `mirt`, and `scales` packages must be
#' installed (most are in Suggests, not Imports).
#'
#' @seealso \code{\link{RMlocdepQ3}}, \code{\link{RMlocdepQ3Cutoff}},
#'   \code{\link{RMlocdepGammaPlot}}
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
#' cutoff_res <- RMlocdepQ3Cutoff(sim_data, iterations = 100,
#'                                parallel = FALSE, seed = 42)
#'
#' # Simulated distribution only
#' RMlocdepQ3Plot(cutoff_res)
#'
#' # With observed Q3 overlaid
#' RMlocdepQ3Plot(cutoff_res, data = sim_data)
#'
#' # Top 10 pairs by departure from null
#' RMlocdepQ3Plot(cutoff_res, data = sim_data, n_pairs = 10)
#' }
RMlocdepQ3Plot <- function(simfit, data, items = NULL, n_pairs = NULL) {

  # --- Check required packages ------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for RMlocdepQ3Plot().", call. = FALSE)
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required for RMlocdepQ3Plot().", call. = FALSE)
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required for RMlocdepQ3Plot().", call. = FALSE)
  }

  # --- Validate simfit --------------------------------------------------------
  required_names <- c("pair_results", "pair_cutoffs", "actual_iterations",
                      "sample_n", "item_names")
  missing_names <- setdiff(required_names, names(simfit))
  if (length(missing_names) > 0L) {
    stop(
      "`simfit` is missing required components: ",
      paste(missing_names, collapse = ", "),
      ".\nExpected the return value of RMlocdepQ3Cutoff() (this requires ",
      "easyRasch2 >= 0.7; rerun the cutoff function if you have an older ",
      "cached result).",
      call. = FALSE
    )
  }

  results_df        <- simfit$pair_results
  pair_cutoffs      <- simfit$pair_cutoffs
  actual_iterations <- simfit$actual_iterations
  sample_n          <- simfit$sample_n
  item_names        <- simfit$item_names

  # --- Validate items parameter -----------------------------------------------
  if (!is.null(items)) {
    unknown_items <- setdiff(items, item_names)
    if (length(unknown_items) > 0L) {
      stop(
        "Unknown item(s) in `items`: ",
        paste(unknown_items, collapse = ", "),
        ".\nAvailable items: ", paste(item_names, collapse = ", "),
        call. = FALSE
      )
    }
    if (length(items) < 2L) {
      stop("`items` must contain at least 2 item names to form a pair.",
           call. = FALSE)
    }
  }

  # --- Validate n_pairs parameter ---------------------------------------------
  if (!is.null(n_pairs)) {
    if (!is.numeric(n_pairs) || length(n_pairs) != 1L ||
        !is.finite(n_pairs) || n_pairs < 1 ||
        n_pairs != as.integer(n_pairs)) {
      stop("`n_pairs` must be a single positive integer or NULL.",
           call. = FALSE)
    }
    n_pairs <- as.integer(n_pairs)
  }

  # --- Create pair labels -----------------------------------------------------
  results_df$Pair   <- paste(results_df$Item1, "-", results_df$Item2)
  pair_cutoffs$Pair <- paste(pair_cutoffs$Item1, "-", pair_cutoffs$Item2)

  # --- Filter to selected items -----------------------------------------------
  if (!is.null(items)) {
    keep <- results_df$Item1 %in% items & results_df$Item2 %in% items
    results_df <- results_df[keep, , drop = FALSE]
    keep_cut <- pair_cutoffs$Item1 %in% items & pair_cutoffs$Item2 %in% items
    pair_cutoffs <- pair_cutoffs[keep_cut, , drop = FALSE]

    if (nrow(results_df) == 0L) {
      stop("No item pairs remain after filtering by `items`.", call. = FALSE)
    }
  }

  # --- Compute observed Q3 (when data supplied) -------------------------------
  observed_df <- NULL
  if (!missing(data)) {
    if (!requireNamespace("mirt", quietly = TRUE)) {
      stop(
        "Package 'mirt' is required to compute observed Q3 but is not ",
        "installed.\nInstall it with: install.packages(\"mirt\")",
        call. = FALSE
      )
    }
    validate_response_data(data)

    mirt_fit <- mirt::mirt(
      data,
      model      = 1L,
      itemtype   = "Rasch",
      verbose    = FALSE,
      accelerate = "squarem"
    )
    q3_mat <- mirt::residuals(mirt_fit, type = "Q3", digits = 4,
                              verbose = FALSE)
    diag(q3_mat) <- NA
    item_names_q3 <- colnames(q3_mat)
    if (is.null(item_names_q3)) item_names_q3 <- item_names
    upper_idx <- which(upper.tri(q3_mat), arr.ind = TRUE)
    observed_df <- data.frame(
      Item1       = item_names_q3[upper_idx[, "row"]],
      Item2       = item_names_q3[upper_idx[, "col"]],
      observed_Q3 = q3_mat[upper_idx],
      stringsAsFactors = FALSE
    )
    observed_df$Pair <- paste(observed_df$Item1, "-", observed_df$Item2)

    if (!is.null(items)) {
      keep_obs <- observed_df$Item1 %in% items & observed_df$Item2 %in% items
      observed_df <- observed_df[keep_obs, , drop = FALSE]
    }
  }

  # --- Apply n_pairs filter (rank by |observed - median(sim)| or |median|) ---
  if (!is.null(n_pairs)) {
    # Per-pair simulated median (used by either ranking metric below)
    pair_names_all <- unique(results_df$Pair)
    med_sim <- vapply(pair_names_all, function(pp) {
      stats::median(results_df$Q3[results_df$Pair == pp], na.rm = TRUE)
    }, numeric(1L))
    names(med_sim) <- pair_names_all

    if (!is.null(observed_df)) {
      # Rank by |observed Q3 − median(simulated Q3 per pair)|
      obs_lookup <- stats::setNames(observed_df$observed_Q3, observed_df$Pair)
      deviation  <- abs(obs_lookup[pair_names_all] - med_sim[pair_names_all])
    } else {
      # Rank by |median simulated Q3 per pair|
      deviation <- abs(med_sim[pair_names_all])
    }
    ord <- order(deviation, decreasing = TRUE)
    keep_n     <- min(n_pairs, length(pair_names_all))
    keep_pairs <- pair_names_all[ord[seq_len(keep_n)]]

    results_df   <- results_df[results_df$Pair %in% keep_pairs, , drop = FALSE]
    pair_cutoffs <- pair_cutoffs[pair_cutoffs$Pair %in% keep_pairs, ,
                                 drop = FALSE]
    if (!is.null(observed_df)) {
      observed_df <- observed_df[observed_df$Pair %in% keep_pairs, ,
                                 drop = FALSE]
    }
    # Y-axis order: largest deviation at the top
    pair_levels <- rev(keep_pairs)
  } else {
    pair_levels <- rev(unique(results_df$Pair))
  }

  # --- Per-pair summary intervals for segment overlays ------------------------
  pair_names <- unique(results_df$Pair)
  lo_hi <- do.call(rbind, lapply(pair_names, function(pair) {
    sub <- results_df[results_df$Pair == pair, ]
    data.frame(
      Pair         = pair,
      min_Q3       = stats::quantile(sub$Q3, 0.005, na.rm = TRUE),
      max_Q3       = stats::quantile(sub$Q3, 0.995, na.rm = TRUE),
      p66lo_Q3     = stats::quantile(sub$Q3, 0.167, na.rm = TRUE),
      p66hi_Q3     = stats::quantile(sub$Q3, 0.833, na.rm = TRUE),
      median_Q3    = stats::median(sub$Q3, na.rm = TRUE),
      stringsAsFactors = FALSE,
      row.names    = NULL
    )
  }))
  rownames(lo_hi) <- NULL

  # --- Case 1: no observed data, show simulation distribution only ------------
  if (missing(data)) {

    results_plot <- data.frame(
      Pair  = results_df$Pair,
      Value = results_df$Q3,
      stringsAsFactors = FALSE
    )
    results_plot$Pair <- factor(results_plot$Pair, levels = pair_levels)

    p <- ggplot2::ggplot(
      results_plot,
      ggplot2::aes(x = .data$Value, y = .data$Pair)
    ) +
      ggdist::stat_dotsinterval(
        ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
        quantiles      = actual_iterations,
        point_interval = "median_hdci",
        layout         = "weave",
        slab_color     = NA,
        .width         = c(0.66, 0.95, 0.99)
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype   = "dashed",
        color      = "grey50",
        linewidth  = 0.4
      ) +
      ggplot2::labs(
        x       = "Q3 residual correlation",
        y       = "Item pair",
        caption = paste0(
          "Note: Results from ", actual_iterations,
          " simulated datasets with ", sample_n,
          " respondents (no true local dependence)."
        )
      ) +
      ggplot2::scale_color_manual(
        values     = scales::brewer_pal()(3),
        aesthetics = "slab_fill",
        guide      = "none"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"))

    return(p)
  }

  # --- Case 2: observed data supplied -----------------------------------------
  Q3_sim <- data.frame(
    Pair  = results_df$Pair,
    Value = results_df$Q3,
    stringsAsFactors = FALSE
  )
  Q3_sim <- merge(Q3_sim, observed_df[, c("Pair", "observed_Q3")],
                  by = "Pair", sort = FALSE)
  Q3_sim$Pair <- factor(Q3_sim$Pair, levels = pair_levels)

  lo_hi$Pair_f <- factor(lo_hi$Pair, levels = pair_levels)

  caption_text <- paste0(
    "Note: Results from ", actual_iterations,
    " simulated datasets with ", sample_n, " respondents.\n",
    "Orange diamonds indicate observed Q3 residual correlation. ",
    "Black dots indicate median Q3 from simulations."
  )

  p <- ggplot2::ggplot(
    Q3_sim,
    ggplot2::aes(x = .data$Value, y = .data$Pair)
  ) +
    ggdist::stat_dots(
      ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
      quantiles  = actual_iterations,
      layout     = "weave",
      slab_color = NA,
      .width     = c(0.66, 0.95, 0.99)
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x    = .data$min_Q3,
        xend = .data$max_Q3,
        y    = .data$Pair_f,
        yend = .data$Pair_f
      ),
      color     = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x    = .data$p66lo_Q3,
        xend = .data$p66hi_Q3,
        y    = .data$Pair_f,
        yend = .data$Pair_f
      ),
      color     = "black",
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = lo_hi,
      ggplot2::aes(
        x = .data$median_Q3,
        y = .data$Pair_f
      ),
      size = 3.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$observed_Q3),
      color    = "sienna2",
      shape    = 18,
      position = ggplot2::position_nudge(y = -0.1),
      size     = 4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      color      = "grey50",
      linewidth  = 0.4
    ) +
    ggplot2::labs(
      x       = "Q3 residual correlation",
      y       = "Item pair",
      caption = caption_text
    ) +
    ggplot2::scale_color_manual(
      values     = scales::brewer_pal()(3),
      aesthetics = "slab_fill",
      guide      = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm"))

  p
}
