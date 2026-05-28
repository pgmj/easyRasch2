#' Item-Threshold Hierarchy Plot
#'
#' Visualises item and threshold locations on the logit scale for a
#' Partial Credit Model. Items are sorted by their (mean-threshold)
#' location; each item shows its location as a black diamond and its
#' individual thresholds as coloured dots with confidence-interval error
#' bars.
#'
#' Threshold locations are centred at the grand mean of all thresholds
#' across items, so the dashed reference line at 0 represents the mean
#' threshold location across the scale.
#'
#' Confidence intervals around thresholds are 84% by default
#' (`sem_multiplier = 1.405`), following Payton, Greenstone, & Schenker
#' (2003) -- non-overlap of 84% intervals approximately corresponds to a
#' two-sample significance test at \eqn{\alpha = 0.05}. Use
#' `sem_multiplier = 1.96` for 95% intervals.
#'
#' @param data A data.frame or matrix of polytomous item responses
#'   (non-negative integers, 0-based, max value > 1). One column per
#'   item, one row per person.
#' @param show_numbers Logical. When `TRUE` (default), prints the
#'   numerical item location and threshold values next to the points on
#'   the plot. When `FALSE`, only the threshold labels (`T1`, `T2`, ...)
#'   are shown.
#' @param sem_multiplier Numeric multiplier for the threshold SE used to
#'   draw the error bars. Default `1.405` (84% CI). Common alternatives:
#'   `1.96` (95% CI), `2.576` (99% CI).
#' @param item_labels Optional character vector of length `ncol(data)`
#'   providing descriptive labels for items. Default `NULL` uses the
#'   column names of `data`. Labels are appended to the column names on
#'   the y-axis (`name - label`) and wrapped at 36 characters via base-R
#'   `strwrap()`.
#' @param output One of `"ggplot"` (default) -- a faceted hierarchy
#'   figure -- or `"dataframe"` -- the long-format underlying data with
#'   one row per (item × threshold).
#'
#' @return Either a `ggplot` (default) or a data.frame with columns
#'   `Item`, `ItemLabel`, `Threshold`, `ThresholdLocation`,
#'   `ThresholdSE`, and `ItemLocation` (the per-item mean of the centred
#'   thresholds).
#'
#' @details
#' \strong{Polytomous only.} Dichotomous items have a single threshold
#' that coincides with the item location; the hierarchy plot is
#' visually degenerate in that case. For dichotomous data use
#' [RMtargeting()] or [RMscoreSE()] instead.
#'
#' \strong{Centring convention.} The PCM thresholds returned by
#' [eRm::thresholds()] are shifted so that their grand mean is zero;
#' each item's location is then the mean of its centred thresholds.
#' The dashed horizontal reference line on the plot marks this zero --
#' i.e., the average threshold across all items.
#'
#' @references
#' Payton, M. E., Greenstone, M. H., & Schenker, N. (2003). Overlapping
#' confidence intervals or standard error intervals: What do they mean
#' in terms of statistical significance? \emph{Journal of Insect Science,
#' 3}(34), 1-6. \doi{10.1093/jis/3.1.34}
#'
#' @seealso [RMtargeting()], [RMscoreSE()], [RMitemICCPlot()]
#'
#' @examples
#' \donttest{
#' data("pcmdat2", package = "eRm")
#' RMitemHierarchy(pcmdat2)
#'
#' # 95% CI instead of 84%
#' RMitemHierarchy(pcmdat2, sem_multiplier = 1.96)
#'
#' # Underlying data.frame
#' RMitemHierarchy(pcmdat2, output = "dataframe")
#' }
#'
#' @importFrom rlang .data
#' @export
RMitemHierarchy <- function(data,
                            show_numbers   = TRUE,
                            sem_multiplier = 1.405,
                            item_labels    = NULL,
                            output         = c("ggplot", "dataframe")) {

  output <- match.arg(output)

  if (output == "ggplot" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for output = \"ggplot\". ",
         "Install with: install.packages(\"ggplot2\")",
         call. = FALSE)
  }

  # --- Data validation ----------------------------------------------------
  validate_response_data(data)
  data <- as.data.frame(data)
  if (ncol(data) < 2L) {
    stop("RMitemHierarchy() requires at least 2 items.", call. = FALSE)
  }

  data_mat <- as.matrix(data)
  if (max(data_mat, na.rm = TRUE) <= 1L) {
    stop("RMitemHierarchy() requires polytomous data (max response > 1). ",
         "For dichotomous items the threshold coincides with the item ",
         "location, so the hierarchy plot is degenerate. Use ",
         "RMtargeting() or RMscoreSE() instead.",
         call. = FALSE)
  }

  if (!is.numeric(sem_multiplier) || length(sem_multiplier) != 1L ||
      !is.finite(sem_multiplier) || sem_multiplier <= 0) {
    stop("`sem_multiplier` must be a positive numeric value. Common ",
         "choices: 1.405 (84% CI), 1.96 (95%), 2.576 (99%).",
         call. = FALSE)
  }

  if (is.null(item_labels)) {
    item_labels <- names(data)
  } else {
    if (length(item_labels) != ncol(data)) {
      stop("`item_labels` must have the same length as ncol(data) (",
           ncol(data), "). Got ", length(item_labels), ".", call. = FALSE)
    }
    item_labels <- as.character(item_labels)
  }

  # --- Fit PCM and extract centred thresholds + SEs -----------------------
  # rgl workaround for any iarm-related dependency chains
  old_rgl <- getOption("rgl.useNULL")
  options(rgl.useNULL = TRUE)
  on.exit(options(rgl.useNULL = old_rgl), add = TRUE)

  erm_fit <- eRm::PCM(data)
  thr_obj <- eRm::thresholds(erm_fit)

  # threshpar / se.thresh: named numeric vectors in the same item-threshold
  # order, e.g. "thresh beta I1.c1", "thresh beta I1.c2", ...
  loc_named <- thr_obj$threshpar
  se_named  <- thr_obj$se.thresh
  if (is.list(loc_named)) loc_named <- unlist(loc_named, use.names = TRUE)
  if (is.list(se_named))  se_named  <- unlist(se_named,  use.names = TRUE)

  parsed <- parse_threshold_names(names(loc_named), names(data))

  # Centre thresholds at the grand mean (matching the original
  # easyRasch::RIitemHierarchy convention)
  thresh_raw     <- as.numeric(loc_named)
  thresh_centred <- thresh_raw - mean(thresh_raw, na.rm = TRUE)

  long_df <- data.frame(
    Item              = parsed$Item,
    ItemLabel         = item_labels[match(parsed$Item, names(data))],
    Threshold         = paste0("T", as.integer(
                          sub("[^0-9]+", "", parsed$Threshold)
                        )),
    ThresholdLocation = thresh_centred,
    ThresholdSE       = as.numeric(se_named),
    stringsAsFactors  = FALSE
  )

  # Per-item Location = mean of that item's centred thresholds
  item_loc <- stats::aggregate(
    ThresholdLocation ~ Item,
    data = long_df,
    FUN  = function(x) mean(x, na.rm = TRUE)
  )
  names(item_loc)[2L] <- "ItemLocation"
  long_df <- merge(long_df, item_loc, by = "Item", sort = FALSE)

  # Order items by item location (lowest first, so they appear at the
  # bottom of the y-axis after coord_flip)
  item_order  <- item_loc$Item[order(item_loc$ItemLocation)]
  label_order <- item_labels[match(item_order, names(data))]
  long_df$Item <- factor(long_df$Item, levels = item_order)

  rownames(long_df) <- NULL

  if (output == "dataframe") return(long_df)

  # --- Plot ---------------------------------------------------------------
  # Axis labels: "I1 - label" wrapped at 36 chars (base R; matches the
  # original's str_wrap usage).
  axis_labels <- vapply(
    paste0(item_order, " - ", label_order),
    function(s) paste(strwrap(s, width = 36L), collapse = "\n"),
    character(1L)
  )

  # `scale_x_discrete(labels = ...)` does a named lookup against the
  # factor levels (= `item_order`), so we must name the vector
  # accordingly. If we don't, the level names slip through unchanged
  # and the `item_labels` argument silently has no effect.
  names(axis_labels) <- item_order

  caption_base <-
    "Note. Item locations are indicated by black diamond shapes"
  caption_thr <-
    "Item threshold locations are indicated by coloured dots and labels."
  caption_ci  <- sprintf(
    "Horizontal error bars indicate +/-%.3f x SE around each threshold.",
    sem_multiplier
  )
  caption_text <- paste(
    if (show_numbers) paste0(caption_base, " and black text.")
    else              paste0(caption_base, "."),
    caption_thr, caption_ci,
    sep = "\n"
  )

  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = .data$Item,
                                    colour = .data$Threshold)) +
    # Item location: black diamond
    ggplot2::geom_point(
      ggplot2::aes(y = .data$ItemLocation),
      size = 4, shape = 18, colour = "black"
    ) +
    # Threshold locations: coloured dots + 84%-style CI error bars
    ggplot2::geom_point(ggplot2::aes(y = .data$ThresholdLocation)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$ThresholdLocation - sem_multiplier * .data$ThresholdSE,
        ymax = .data$ThresholdLocation + sem_multiplier * .data$ThresholdSE
      ),
      width = 0.11
    ) +
    # Threshold labels (T1, T2, ...) below each coloured dot
    ggplot2::geom_text(
      ggplot2::aes(y = .data$ThresholdLocation, label = .data$Threshold),
      size = 3, hjust = 0.5, vjust = 1.5, show.legend = FALSE
    )

  if (show_numbers) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(y = .data$ItemLocation,
                     label = sprintf("%.2f", .data$ItemLocation)),
        hjust = 0.5, vjust = -1.3,
        colour = "black", size = 3, show.legend = FALSE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(y = .data$ThresholdLocation,
                     label = sprintf("%.2f", .data$ThresholdLocation)),
        hjust = 0.5, vjust = -1.1,
        size = 3, show.legend = FALSE
      )
  }

  p +
    ggplot2::geom_hline(
      yintercept = mean(long_df$ItemLocation, na.rm = TRUE),
      linetype = "dashed", colour = "darkgrey"
    ) +
    ggplot2::geom_rug(
      ggplot2::aes(y = .data$ThresholdLocation),
      colour = "darkgrey", sides = "l",
      length = ggplot2::unit(0.02, "npc")
    ) +
    ggplot2::scale_x_discrete(labels = axis_labels) +
    ggplot2::coord_flip() +
    ggplot2::scale_colour_viridis_d(guide = "none", option = "H") +
    ggplot2::labs(y = "Logit location",
                  x = NULL,
                  caption = caption_text) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      plot.caption    = ggplot2::element_text(hjust = 0, face = "italic",
                                              size = 9)
    ) +
    er2_axis_margins()
}
