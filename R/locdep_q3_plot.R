#' Plot the simulated Q3 residual correlations
#'
#' Plots the simulated \eqn{Q_3} of each item pair from [RMlocdepQ3Cutoff()],
#' and optionally the observed \eqn{Q_3} in your data and a heatmap of the
#' observed matrix.
#'
#' @param simfit The result of [RMlocdepQ3Cutoff()].
#' @param data Optional. The item responses the simulations were made from,
#'   items scored from 0. The observed \eqn{Q_3} is then computed with the
#'   estimator of `simfit` and drawn as orange diamonds, and the heatmap is
#'   returned.
#' @param items Optional character vector of item names. Only pairs with
#'   both items in it are shown. Default `NULL`, all pairs.
#' @param n_pairs Optional positive integer. Only the `n_pairs` pairs with
#'   the largest departure are shown: `|observed - simulated median|` with
#'   `data`, `|simulated median|` without. Applied after `items`. Values
#'   above the number of pairs show all pairs.
#'
#' @return A named list of two `ggplot` objects, with the same structure as
#'   the tables of [RMlocdepQ3()]:
#'   \describe{
#'     \item{`pairs`}{one row per item pair, labelled "Item1 - Item2".}
#'     \item{`matrix`}{a heatmap of the lower triangle of the observed
#'       \eqn{Q_3} matrix, with pairs above the global dynamic cutoff
#'       outlined. It needs `data`, so it is `NULL` (with a message)
#'       without it. With `items`, it shows those items only. `n_pairs`
#'       does not apply to it.}
#'   }
#'
#' @details
#' \eqn{Q_3} is symmetric, so each pair is shown once (unlike partial
#' gamma, which depends on the direction).
#'
#' Each pair's simulated values are drawn as dots, with the median, an
#' inner 66% interval and an outer interval. The outer interval is the one
#' in `simfit$pair_cutoffs`, which [RMlocdepQ3()] tabulates: the HDCI of
#' width `hdci_width`, or the 2.5th to 97.5th percentiles with
#' `cutoff_method = "quantile"`. The inner interval uses the same method.
#' Without `data` the intervals are drawn with
#' `ggdist::stat_dotsinterval()`. With `data` they are black segments with
#' a black dot for the median, and the observed value is an orange diamond.
#'
#' The intervals are per pair and not corrected for testing many pairs.
#' With the full cutoff object, [RMlocdepQ3()] flags pairs by adjusted
#' p-values instead, so a pair just outside its interval need not be
#' flagged there. The plot shows no p-values.
#'
#' The \pkg{ggplot2}, \pkg{ggdist} and \pkg{scales} packages must be
#' installed, and \pkg{mirt} for `estimator = "MML"`.
#'
#' @seealso [RMlocdepQ3()], [RMlocdepQ3Cutoff()], [RMlocdepGammaPlot()]
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("ggdist", quietly = TRUE)) {
#'   set.seed(42)
#'   sim_data <- as.data.frame(
#'     matrix(sample(0:1, 200 * 10, replace = TRUE), nrow = 200, ncol = 10)
#'   )
#'   colnames(sim_data) <- paste0("Item", 1:10)
#'
#'   # Run the simulation (use 1000 or more iterations in real analyses)
#'   cutoff_res <- RMlocdepQ3Cutoff(sim_data, iterations = 50,
#'                                  parallel = FALSE, seed = 42)
#'
#'   # Simulated distribution only
#'   RMlocdepQ3Plot(cutoff_res)
#'
#'   # With observed Q3 overlaid
#'   RMlocdepQ3Plot(cutoff_res, data = sim_data)
#'
#'   # Top 10 pairs by departure from null
#'   RMlocdepQ3Plot(cutoff_res, data = sim_data, n_pairs = 10)
#' }
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
  required_names <- c(
    "pair_results",
    "pair_cutoffs",
    "actual_iterations",
    "sample_n",
    "item_names"
  )
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

  results_df <- simfit$pair_results
  pair_cutoffs <- simfit$pair_cutoffs
  actual_iterations <- simfit$actual_iterations
  sample_n <- simfit$sample_n
  item_names <- simfit$item_names
  # Match the estimator used to build the cut-off (older objects predate it).
  estimator <- if (is.null(simfit$estimator)) "CML" else simfit$estimator

  # Standard sample-size clause. The Q3 cutoff retains incomplete responses
  # (matching RMlocdepQ3); `sample_n_total` / `sample_has_na` are absent in
  # cutoff objects made by older versions, so fall back to the plain count.
  sample_clause <- .n_caption(
    sample_n,
    if (is.null(simfit$sample_n_total)) sample_n else simfit$sample_n_total,
    if (isTRUE(simfit$sample_has_na)) {
      "incomplete responses retained"
    } else {
      character()
    }
  )

  # --- Validate items parameter -----------------------------------------------
  if (!is.null(items)) {
    unknown_items <- setdiff(items, item_names)
    if (length(unknown_items) > 0L) {
      stop(
        "Unknown item(s) in `items`: ",
        paste(unknown_items, collapse = ", "),
        ".\nAvailable items: ",
        paste(item_names, collapse = ", "),
        call. = FALSE
      )
    }
    if (length(items) < 2L) {
      stop(
        "`items` must contain at least 2 item names to form a pair.",
        call. = FALSE
      )
    }
  }

  # --- Validate n_pairs parameter ---------------------------------------------
  if (!is.null(n_pairs)) {
    if (
      !is.numeric(n_pairs) ||
        length(n_pairs) != 1L ||
        !is.finite(n_pairs) ||
        n_pairs < 1 ||
        n_pairs != as.integer(n_pairs)
    ) {
      stop(
        "`n_pairs` must be a single positive integer or NULL.",
        call. = FALSE
      )
    }
    n_pairs <- as.integer(n_pairs)
  }

  # --- Create pair labels -----------------------------------------------------
  results_df$Pair <- paste(results_df$Item1, "-", results_df$Item2)
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
    if (estimator == "MML" && !requireNamespace("mirt", quietly = TRUE)) {
      stop(
        "Package 'mirt' is required to compute observed Q3 with ",
        "estimator = \"MML\" but is not installed.\n",
        "Install it with: install.packages(\"mirt\")",
        call. = FALSE
      )
    }
    validate_response_data(data)

    # Respondents with no responses at all break the CML fit inside
    # .q3_residual_matrix() (psychotools errors on all-NA rows); drop them,
    # as RMlocdepQ3() and RMlocdepQ3Cutoff() do.
    data <- .drop_empty_respondents(data)
    .check_cutoff_sample(
      simfit$sample_n,
      nrow(data),
      "RMlocdepQ3Cutoff()",
      policy = "respondents with at least one response"
    )

    # Same estimator as the simulated cut-off (stored in simfit$estimator).
    q3_mat <- .q3_residual_matrix(data, estimator = estimator)
    item_names_q3 <- colnames(q3_mat)
    if (is.null(item_names_q3)) {
      item_names_q3 <- item_names
    }
    upper_idx <- which(upper.tri(q3_mat), arr.ind = TRUE)
    observed_df <- data.frame(
      Item1 = item_names_q3[upper_idx[, "row"]],
      Item2 = item_names_q3[upper_idx[, "col"]],
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
    med_sim <- vapply(
      pair_names_all,
      function(pp) {
        stats::median(results_df$Q3[results_df$Pair == pp], na.rm = TRUE)
      },
      numeric(1L)
    )
    names(med_sim) <- pair_names_all

    if (!is.null(observed_df)) {
      # Rank by |observed Q3 − median(simulated Q3 per pair)|
      obs_lookup <- stats::setNames(observed_df$observed_Q3, observed_df$Pair)
      deviation <- abs(obs_lookup[pair_names_all] - med_sim[pair_names_all])
    } else {
      # Rank by |median simulated Q3 per pair|
      deviation <- abs(med_sim[pair_names_all])
    }
    ord <- order(deviation, decreasing = TRUE)
    keep_n <- min(n_pairs, length(pair_names_all))
    keep_pairs <- pair_names_all[ord[seq_len(keep_n)]]

    results_df <- results_df[results_df$Pair %in% keep_pairs, , drop = FALSE]
    pair_cutoffs <- pair_cutoffs[
      pair_cutoffs$Pair %in% keep_pairs,
      ,
      drop = FALSE
    ]
    if (!is.null(observed_df)) {
      observed_df <- observed_df[
        observed_df$Pair %in% keep_pairs,
        ,
        drop = FALSE
      ]
    }
    # Y-axis order: largest deviation at the top
    pair_levels <- rev(keep_pairs)
  } else {
    pair_levels <- rev(unique(results_df$Pair))
  }

  # --- Per-pair summary intervals for segment overlays ------------------------
  # The outer interval is the one stored in `simfit$pair_cutoffs` (the one
  # RMlocdepQ3() tabulates) and the inner 66% interval uses the same method.
  spec <- .sim_interval_spec(simfit)
  interval_note <- .sim_interval_caption(spec, "pair", "RMlocdepQ3()")
  pair_names <- unique(results_df$Pair)
  cut_idx <- match(pair_names, pair_cutoffs$Pair)
  iv <- .sim_plot_intervals(
    results_df$Q3, results_df$Pair, pair_names,
    pair_cutoffs$Q3_low[cut_idx], pair_cutoffs$Q3_high[cut_idx],
    spec$method
  )
  lo_hi <- data.frame(
    Pair = pair_names,
    min_Q3 = iv$lo,
    max_Q3 = iv$hi,
    p66lo_Q3 = iv$lo66,
    p66hi_Q3 = iv$hi66,
    median_Q3 = iv$median,
    stringsAsFactors = FALSE
  )

  # --- Case 1: no observed data, show simulation distribution only ------------
  if (missing(data)) {
    results_plot <- data.frame(
      Pair = results_df$Pair,
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
        quantiles = actual_iterations,
        point_interval = .sim_point_interval(spec),
        layout = "weave",
        slab_color = NA,
        .width = c(0.66, spec$width)
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "grey50",
        linewidth = 0.4
      ) +
      ggplot2::labs(
        x = "Q3 residual correlation",
        y = "Item pair",
        caption = er2_caption(paste0(
          "Results from ",
          actual_iterations,
          " simulated datasets (no true local dependence). ",
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
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
      er2_axis_margins() +
      er2_plot_caption()

    # The $matrix tile heatmap needs the observed Q3 matrix, which requires
    # `data`; without it only the per-pair distribution ($pairs) is available.
    message(
      "The Q3 tile heatmap ($matrix) requires `data`; ",
      "returning $matrix = NULL."
    )
    return(list(matrix = NULL, pairs = p))
  }

  # --- Case 2: observed data supplied -----------------------------------------
  Q3_sim <- data.frame(
    Pair = results_df$Pair,
    Value = results_df$Q3,
    stringsAsFactors = FALSE
  )
  Q3_sim <- merge(
    Q3_sim,
    observed_df[, c("Pair", "observed_Q3")],
    by = "Pair",
    sort = FALSE
  )
  Q3_sim$Pair <- factor(Q3_sim$Pair, levels = pair_levels)

  lo_hi$Pair_f <- factor(lo_hi$Pair, levels = pair_levels)

  caption_text <- er2_caption(paste0(
    "Results from ",
    actual_iterations,
    " simulated datasets. ",
    sample_clause,
    " per dataset.\n",
    "Orange diamonds indicate observed Q3 residual correlation. ",
    "Black dots indicate median Q3 from simulations.\n",
    interval_note
  ))

  p <- ggplot2::ggplot(
    Q3_sim,
    ggplot2::aes(x = .data$Value, y = .data$Pair)
  ) +
    ggdist::stat_dots(
      ggplot2::aes(slab_fill = ggplot2::after_stat(.data$level)),
      quantiles = actual_iterations,
      point_interval = .sim_point_interval(spec),
      layout = "weave",
      slab_color = NA,
      .width = c(0.66, spec$width)
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x = .data$min_Q3,
        xend = .data$max_Q3,
        y = .data$Pair_f,
        yend = .data$Pair_f
      ),
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = lo_hi,
      ggplot2::aes(
        x = .data$p66lo_Q3,
        xend = .data$p66hi_Q3,
        y = .data$Pair_f,
        yend = .data$Pair_f
      ),
      color = "black",
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
      color = "sienna2",
      shape = 18,
      position = ggplot2::position_nudge(y = -0.1),
      size = 4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.4
    ) +
    ggplot2::labs(
      x = "Q3 residual correlation",
      y = "Item pair",
      caption = caption_text
    ) +
    ggplot2::scale_color_manual(
      values = scales::brewer_pal()(3),
      aesthetics = "slab_fill",
      guide = "none"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.7, "cm")) +
    er2_axis_margins() +
    er2_plot_caption()

  # --- $matrix: tile heatmap of the observed Q3 matrix ------------------------
  # Same heatmap earlier versions returned from RMlocdepQ3()$plot, now here so
  # the table and plot outputs share the $matrix / $pairs structure. The global
  # dynamic cut-off mirrors RMlocdepQ3(): mean off-diagonal Q3 + suggested
  # cut-off. The `items` filter subsets the matrix; `n_pairs` does not apply.
  q3_tile <- if (is.null(items)) q3_mat else q3_mat[items, items, drop = FALSE]
  dyn_cutoff <- if (!is.null(simfit$suggested_cutoff)) {
    mean(q3_tile, na.rm = TRUE) + as.numeric(simfit$suggested_cutoff)
  } else {
    NULL
  }
  matrix_plot <- .q3_tile_plot(
    q3_tile,
    dyn_cutoff = dyn_cutoff,
    estimator = estimator,
    actual_iterations = actual_iterations
  )

  list(matrix = matrix_plot, pairs = p)
}
