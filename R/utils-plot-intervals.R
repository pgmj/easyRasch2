# Interval segments for the simulation plots.
#
# The plots draw, per item or pair, an inner 66% segment, an outer segment
# and the median of the simulated values. The outer segment is the interval
# stored in the cutoff object (the one the table reports), and the inner
# segment uses the same method, so plot and table describe the same
# interval. Before 2026-11-20 the plots drew fixed quantile intervals (help-
# page review, findings group A).

#' How a cutoff object computed its intervals
#'
#' @param simfit A cutoff object with `cutoff_method` and `hdci_width`.
#' @return A list with `method` (`"hdci"` or `"quantile"`), `width` (the
#'   outer width, 0.95 for the 2.5th/97.5th percentile interval) and
#'   `label` for captions.
#' @noRd
.sim_interval_spec <- function(simfit) {
  if (identical(simfit$cutoff_method, "hdci") && !is.null(simfit$hdci_width)) {
    width <- simfit$hdci_width
    list(
      method = "hdci",
      width = width,
      label = paste0(format(100 * width), "% HDCI")
    )
  } else {
    list(
      method = "quantile",
      width = 0.95,
      label = "central 95% interval"
    )
  }
}

#' Inner interval and median of simulated values, outer from the cutoffs
#'
#' @param values Numeric simulated values.
#' @param groups Item or pair label of each value.
#' @param keys Labels to summarise, in output order.
#' @param low,high Stored outer bounds, matched to `keys`.
#' @param method `"hdci"` or `"quantile"`, for the inner 66% interval.
#' @return data.frame with `key`, `lo`, `hi`, `lo66`, `hi66`, `median`.
#' @noRd
.sim_plot_intervals <- function(values, groups, keys, low, high, method) {
  rows <- lapply(seq_along(keys), function(i) {
    v <- values[groups == keys[i]]
    v <- v[is.finite(v)]
    inner <- if (length(v) == 0L) {
      c(NA_real_, NA_real_)
    } else if (identical(method, "hdci")) {
      ggdist::hdci(v, .width = 0.66)[1L, ]
    } else {
      stats::quantile(v, c(0.17, 0.83), names = FALSE)
    }
    data.frame(
      key = keys[i],
      lo = as.numeric(low[i]),
      hi = as.numeric(high[i]),
      lo66 = as.numeric(inner[1L]),
      hi66 = as.numeric(inner[2L]),
      median = if (length(v)) stats::median(v) else NA_real_,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

#' Caption sentence describing the segments
#'
#' @param spec Result of `.sim_interval_spec()`.
#' @param unit `"item"` or `"pair"`.
#' @param table_fun Name of the table function, for example
#'   `"RMitemInfit()"`, or `NULL` when the table flags against the interval.
#' @noRd
.sim_interval_caption <- function(spec, unit = "item", table_fun = NULL) {
  out <- paste0(
    "Intervals: 66% and ", spec$label, " of the simulated values per ",
    unit, ", not corrected for testing several ", unit, "s."
  )
  if (!is.null(table_fun)) {
    out <- paste0(
      out, " ", table_fun, " flags by adjusted p-values by default."
    )
  }
  out
}

#' Point-interval function for ggdist dot shading
#' @noRd
.sim_point_interval <- function(spec) {
  if (identical(spec$method, "hdci")) "median_hdci" else "median_qi"
}
