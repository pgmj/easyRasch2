#' DIF analysis via Andersen's likelihood-ratio test
#'
#' Splits a Rasch model by an external grouping variable using
#' \code{eRm::LRtest()} and reports per-group item locations (or per-group
#' threshold locations) together with their standard errors. A single
#' function replaces the four legacy helpers (`RIdifTableLR`,
#' `RIdifThreshTblLR`, `RIdifFigureLR`, `RIdifThreshFigLR`) by exposing the
#' two underlying axes -- \code{level} (item or threshold) and \code{output}
#' (data.frame, kable, or ggplot) -- as arguments. The same data
#' preparation pipeline feeds all six combinations.
#'
#' The Partial Credit Model (PCM) is fitted by default for polytomous data
#' and the dichotomous Rasch Model (RM) is fitted when all responses are
#' 0/1; this can be overridden via \code{model}.
#'
#' @param data A data.frame or matrix of item responses (non-negative
#'   integers, 0-based). One column per item, one row per person. Person
#'   IDs and grouping variables must not be included -- pass the grouping
#'   variable separately via \code{dif_var}.
#' @param dif_var Vector of length \code{nrow(data)} (factor, character,
#'   or numeric) defining the DIF grouping variable. Coerced to factor;
#'   unused levels are dropped. Rows where \code{dif_var} is \code{NA} are
#'   dropped with a message. Must result in at least 2 groups after
#'   cleaning.
#' @param model One of \code{"auto"} (default), \code{"PCM"}, or
#'   \code{"RM"}. \code{"auto"} fits \code{eRm::RM()} when the data are
#'   dichotomous (max response = 1) and \code{eRm::PCM()} otherwise.
#'   \code{"RM"} errors on polytomous data.
#' @param level One of \code{"item"} (default) or \code{"threshold"}.
#'   \code{"item"} reports each item's mean threshold location per group;
#'   \code{"threshold"} reports each individual threshold per group. For
#'   dichotomous (RM) data the two views are equivalent (one threshold per
#'   item).
#' @param output One of \code{"kable"} (default), \code{"dataframe"}, or
#'   \code{"ggplot"}. The data.frame view always carries a
#'   \code{Flagged} logical column; the kable view bolds flagged rows; the
#'   ggplot view shows confidence intervals.
#' @param cutoff Numeric or \code{NULL}. Threshold (in logits) for the
#'   \code{Flagged} column: a row is flagged when
#'   \code{MaxDiff > cutoff}, where \code{MaxDiff} is the difference
#'   between the largest and smallest per-group location for that item
#'   (or threshold). Set to \code{NULL} to suppress flagging. Default
#'   \code{0.5}.
#' @param conf Numeric in (0, 1). Confidence level used for the ggplot
#'   error bars. Default \code{0.95}.
#' @param sort Logical. \code{kable} output only: sort rows by
#'   \code{MaxDiff} (descending). Default \code{FALSE}.
#'
#' @return A data.frame, a \code{knitr_kable} object, or a \code{ggplot}
#'   object, depending on \code{output}.
#'
#'   The data.frame has one row per item (\code{level = "item"}) or per
#'   item x threshold (\code{level = "threshold"}), with columns
#'   \code{Item} (and \code{Threshold} at threshold level), one numeric
#'   column per group level, an \code{All} column for the unsplit fit,
#'   \code{MaxDiff}, \code{Flagged} (when \code{cutoff} is non-\code{NULL}),
#'   and matching \code{SE_*} columns.
#'
#'   The Andersen LR test result is attached as
#'   \code{attr(result, "lr_test")} on the data.frame, in the kable
#'   footnote, and in the ggplot caption (LR \eqn{\chi^2}, df, p-value).
#'
#' @details
#' For the data.frame and kable outputs, locations are reported on the
#' centred eRm parameterisation returned by \code{eRm::thresholds()}.
#' Per-group fits come from \code{eRm::LRtest(..., splitcr = dif_var)};
#' the unsplit fit (\code{All} column) is the model fitted to the full
#' dataset. The Andersen LR statistic, df, and p-value reported as the
#' \code{lr_test} attribute / caption come directly from
#' \code{LRtest()}'s return value.
#'
#' \code{cell_spec()}-style HTML cell colouring used in the legacy
#' easyRasch package has been dropped in favour of a logical
#' \code{Flagged} column (and bold rendering in the kable output), so the
#' kable renders correctly in HTML, LaTeX, and pipe/markdown.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("eRm", quietly = TRUE)) {
#'   set.seed(1)
#'   data("pcmdat2", package = "eRm")
#'   grp <- factor(sample(c("A", "B"), nrow(pcmdat2), replace = TRUE))
#'
#'   # Default: kable of per-group item locations
#'   RMdifLR(pcmdat2, dif_var = grp)
#'
#'   # ggplot panel of item locations with 95% CIs
#'   RMdifLR(pcmdat2, dif_var = grp, output = "ggplot")
#'
#'   # Threshold-level kable, sorted by MaxDiff
#'   RMdifLR(pcmdat2, dif_var = grp, level = "threshold", sort = TRUE)
#'
#'   # Tidy data.frame for downstream use
#'   df <- RMdifLR(pcmdat2, dif_var = grp, output = "dataframe")
#'   attr(df, "lr_test")
#'   df[df$Flagged, ]
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
RMdifLR <- function(
  data,
  dif_var,
  model = c("auto", "PCM", "RM"),
  level = c("item", "threshold"),
  output = c("kable", "dataframe", "ggplot"),
  cutoff = 0.5,
  conf = 0.95,
  sort = FALSE
) {
  model <- match.arg(model)
  level <- match.arg(level)
  output <- match.arg(output)

  if (output == "ggplot" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for output = \"ggplot\".",
      call. = FALSE
    )
  }
  if (output == "kable" && !requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required for output = \"kable\".", call. = FALSE)
  }
  if (!requireNamespace("eRm", quietly = TRUE)) {
    stop("Package 'eRm' is required for RMdifLR().", call. = FALSE)
  }

  if (!is.null(cutoff)) {
    if (
      !is.numeric(cutoff) ||
        length(cutoff) != 1L ||
        !is.finite(cutoff) ||
        cutoff < 0
    ) {
      stop(
        "`cutoff` must be a single non-negative numeric value, or NULL.",
        call. = FALSE
      )
    }
  }
  if (
    !is.numeric(conf) ||
      length(conf) != 1L ||
      !is.finite(conf) ||
      conf <= 0 ||
      conf >= 1
  ) {
    stop("`conf` must be a single numeric in (0, 1).", call. = FALSE)
  }

  validate_response_data(data)
  data <- as.data.frame(data)
  n_total_lr <- nrow(data)
  has_na_lr <- anyNA(data) || anyNA(dif_var)

  # ---------------------------------------------------------------------
  # Resolve dif_var and drop NAs (jointly with `data`)
  # ---------------------------------------------------------------------
  if (length(dif_var) != nrow(data)) {
    stop(
      "`dif_var` must have the same length as nrow(data) (",
      nrow(data),
      "). Got ",
      length(dif_var),
      ".",
      call. = FALSE
    )
  }
  dif_factor <- droplevels(as.factor(dif_var))

  na_mask <- is.na(dif_factor)
  if (any(na_mask)) {
    message(sum(na_mask), " row(s) with NA in `dif_var` dropped.")
    data <- data[!na_mask, , drop = FALSE]
    dif_factor <- droplevels(dif_factor[!na_mask])
  }

  if (nlevels(dif_factor) < 2L) {
    stop(
      "`dif_var` must have at least 2 distinct non-missing levels.",
      call. = FALSE
    )
  }
  groups <- levels(dif_factor)
  nr_groups <- length(groups)

  # ---------------------------------------------------------------------
  # Pick the model
  # ---------------------------------------------------------------------
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "integer"
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L

  resolved_model <- if (model == "auto") {
    if (is_polytomous) "PCM" else "RM"
  } else {
    if (model == "RM" && is_polytomous) {
      stop(
        "`model = \"RM\"` requires dichotomous (0/1) data; use ",
        "\"PCM\" or \"auto\" for polytomous responses.",
        call. = FALSE
      )
    }
    model
  }

  # ---------------------------------------------------------------------
  # Fit, run LRtest, extract per-group + overall thresholds and SEs
  # ---------------------------------------------------------------------
  fit_full <- if (resolved_model == "PCM") eRm::PCM(data) else eRm::RM(data)

  lrt <- tryCatch(
    eRm::LRtest(fit_full, splitcr = dif_factor),
    error = function(e) {
      stop(
        "eRm::LRtest() failed: ",
        conditionMessage(e),
        "\nThis often indicates an empty response category in one ",
        "subgroup. Inspect with RMplotTile(data, group = dif_var).",
        call. = FALSE
      )
    }
  )

  per_group <- extract_lr_locations(lrt, groups, names(data))
  overall <- extract_overall_locations(fit_full, names(data))

  # Combine into a long-form table: Item, Threshold, DIFgroup, Location, SE
  long_df <- rbind(per_group, overall)
  long_df$Item <- factor(long_df$Item, levels = names(data))

  # ---------------------------------------------------------------------
  # Aggregate to the requested level
  # ---------------------------------------------------------------------
  if (level == "item") {
    agg <- stats::aggregate(
      cbind(Location = long_df$Location, SE = long_df$SE),
      by = list(Item = long_df$Item, DIFgroup = long_df$DIFgroup),
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    agg$Item <- factor(agg$Item, levels = names(data))
    table_df <- build_wide_table(agg, groups, level = "item")
    plot_df <- agg
    plot_df$Threshold <- NA_character_
  } else {
    table_df <- build_wide_table(long_df, groups, level = "threshold")
    plot_df <- long_df
  }

  # ---------------------------------------------------------------------
  # Flagged column
  # ---------------------------------------------------------------------
  if (!is.null(cutoff)) {
    table_df$Flagged <- !is.na(table_df$MaxDiff) &
      table_df$MaxDiff > cutoff
  }

  # ---------------------------------------------------------------------
  # Reorder columns: Item [Threshold] groups All MaxDiff [Flagged] SE_*
  # ---------------------------------------------------------------------
  id_cols <- if (level == "item") "Item" else c("Item", "Threshold")
  loc_cols <- c(
    groups,
    "All",
    "MaxDiff",
    if (!is.null(cutoff)) "Flagged" else NULL
  )
  se_cols <- paste0("SE_", c(groups, "All"))
  table_df <- table_df[, c(id_cols, loc_cols, se_cols), drop = FALSE]
  rownames(table_df) <- NULL

  # ---------------------------------------------------------------------
  # Round numeric columns (3 dp) -- match other easyRasch2 functions
  # ---------------------------------------------------------------------
  num_cols <- c(groups, "All", "MaxDiff", se_cols)
  for (cc in num_cols) {
    if (is.numeric(table_df[[cc]])) {
      table_df[[cc]] <- round(table_df[[cc]], 3)
    }
  }

  # ---------------------------------------------------------------------
  # LR test summary (Andersen)
  # ---------------------------------------------------------------------
  lr_summary <- list(
    LR = unname(lrt$LR),
    df = unname(lrt$df),
    p_value = unname(lrt$pvalue),
    n_groups = nr_groups,
    groups = groups,
    model = resolved_model,
    n_persons = nrow(data),
    n_items = ncol(data)
  )
  attr(table_df, "lr_test") <- lr_summary

  # ---------------------------------------------------------------------
  # Dispatch on output
  # ---------------------------------------------------------------------
  if (output == "dataframe") {
    if (sort) {
      table_df <- table_df[order(-table_df$MaxDiff), , drop = FALSE]
      rownames(table_df) <- NULL
    }
    return(table_df)
  }

  if (output == "kable") {
    n_clause <- .n_caption(
      lr_summary$n_persons,
      n_total_lr,
      if (has_na_lr) "complete cases" else character()
    )
    return(render_lr_kable(
      table_df,
      lr_summary,
      level = level,
      cutoff = cutoff,
      sort = sort,
      is_polytomous = is_polytomous,
      n_clause = n_clause
    ))
  }

  # output == "ggplot"
  render_lr_ggplot(
    plot_df,
    lr_summary,
    level = level,
    groups = groups,
    conf = conf,
    is_polytomous = is_polytomous
  )
}

# =====================================================================
# Internal helpers
# =====================================================================

#' Extract per-group threshold locations + SEs from an LRtest fit
#'
#' Returns a long-format data.frame with columns
#' Item, Threshold, DIFgroup, Location, SE.
#'
#' @keywords internal
#' @noRd
extract_lr_locations <- function(lrt, groups, item_names) {
  parts <- lapply(seq_along(groups), function(g) {
    fit_g <- lrt$fitobj[[g]]
    one_fit_long(fit_g, dif_group = groups[g], item_names = item_names)
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' Extract threshold locations + SEs from the unsplit fit
#'
#' @keywords internal
#' @noRd
extract_overall_locations <- function(fit, item_names) {
  one_fit_long(fit, dif_group = "All", item_names = item_names)
}

#' Pull thresholds + SEs out of one eRm `Rm` fit, as a long data.frame
#'
#' Both \code{eRm::PCM()} and \code{eRm::RM()} fits are supported via
#' \code{eRm::thresholds()}, which returns the centred threshold matrix
#' and matching SEs. For RM fits there is exactly one threshold per item.
#'
#' @keywords internal
#' @noRd
one_fit_long <- function(fit, dif_group, item_names) {
  is_rm <- isTRUE(fit$model == "RM")

  if (is_rm) {
    # eRm::thresholds() is polytomous-only. For RM the "threshold" is
    # the item difficulty, available as -betapar with matching se.beta.
    beta <- fit$betapar
    se <- fit$se.beta
    raw_names <- sub("^beta\\s+", "", names(beta))
    if (!all(raw_names %in% item_names)) {
      stop(
        "Could not match RM item parameter names to data columns: ",
        paste(utils::head(setdiff(raw_names, item_names), 5L), collapse = ", "),
        call. = FALSE
      )
    }
    return(data.frame(
      Item = raw_names,
      Threshold = "1",
      DIFgroup = dif_group,
      Location = unname(-as.numeric(beta)),
      SE = unname(as.numeric(se)),
      stringsAsFactors = FALSE
    ))
  }

  thr <- eRm::thresholds(fit)
  loc <- thr$threshpar
  se <- thr$se.thresh
  if (is.list(loc)) {
    loc <- unlist(loc, use.names = TRUE)
  }
  if (is.list(se)) {
    se <- unlist(se, use.names = TRUE)
  }

  parsed <- parse_threshold_names(names(loc), item_names)

  data.frame(
    Item = parsed$Item,
    Threshold = parsed$Threshold,
    DIFgroup = dif_group,
    Location = unname(as.numeric(loc)),
    SE = unname(as.numeric(se)),
    stringsAsFactors = FALSE
  )
}

#' Parse "beta Item.Tk" style names into Item / Threshold columns
#'
#' Tolerates name shapes seen across eRm versions:
#'   "beta I1.c1", "beta I1.T1", "I1.c1", "I1".
#'
#' @keywords internal
#' @noRd
parse_threshold_names <- function(nm, item_names) {
  bare <- sub("^(thresh\\s+)?beta\\s+", "", nm)

  # Match the longest item name as a prefix -- needed because item names
  # themselves may legitimately contain dots ("Q1.a").
  matched_item <- vapply(
    bare,
    function(b) {
      hits <- item_names[
        startsWith(b, item_names) |
          startsWith(b, paste0(item_names, "."))
      ]
      if (length(hits) == 0L) {
        NA_character_
      } else {
        hits[which.max(nchar(hits))]
      }
    },
    character(1L)
  )

  thr_part <- mapply(
    function(b, it) {
      if (is.na(it)) {
        return(NA_character_)
      }
      rest <- substr(b, nchar(it) + 1L, nchar(b))
      rest <- sub("^\\.", "", rest)
      if (!nzchar(rest)) "1" else rest
    },
    bare,
    matched_item,
    USE.NAMES = FALSE
  )

  if (any(is.na(matched_item))) {
    bad <- nm[is.na(matched_item)]
    stop(
      "Could not match threshold parameter name(s) to items: ",
      paste(utils::head(bad, 5L), collapse = ", "),
      if (length(bad) > 5L) ", ..." else "",
      call. = FALSE
    )
  }

  list(Item = unname(matched_item), Threshold = unname(thr_part))
}

#' Reshape long -> wide on DIFgroup, add MaxDiff, attach SE columns
#'
#' @keywords internal
#' @noRd
build_wide_table <- function(long_df, groups, level) {
  id_cols <- if (level == "item") "Item" else c("Item", "Threshold")

  loc_wide <- stats::reshape(
    long_df[, c(id_cols, "DIFgroup", "Location"), drop = FALSE],
    idvar = id_cols,
    timevar = "DIFgroup",
    direction = "wide"
  )
  names(loc_wide) <- sub("^Location\\.", "", names(loc_wide))

  se_wide <- stats::reshape(
    long_df[, c(id_cols, "DIFgroup", "SE"), drop = FALSE],
    idvar = id_cols,
    timevar = "DIFgroup",
    direction = "wide"
  )
  names(se_wide) <- sub("^SE\\.", "SE_", names(se_wide))

  out <- merge(loc_wide, se_wide, by = id_cols, sort = FALSE)

  # MaxDiff across the per-group columns (excludes "All")
  group_mat <- as.matrix(out[, groups, drop = FALSE])
  out$MaxDiff <- apply(group_mat, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) NA_real_ else max(x) - min(x)
  })

  ord <- if (level == "threshold") {
    order(out$Item, out$Threshold)
  } else {
    order(out$Item)
  }
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Render the wide table as a kable
#'
#' @keywords internal
#' @noRd
render_lr_kable <- function(
  table_df,
  lr_summary,
  level,
  cutoff,
  sort,
  is_polytomous,
  n_clause = NULL
) {
  if (sort) {
    table_df <- table_df[order(-table_df$MaxDiff), , drop = FALSE]
    rownames(table_df) <- NULL
  }

  display <- table_df

  # Bold flagged rows by wrapping each visible cell in **...**
  if ("Flagged" %in% names(display)) {
    flag <- isTRUE_vec(display$Flagged)
    display$Flagged <- ifelse(flag, "yes", "no")
    if (any(flag)) {
      data_cols <- setdiff(
        names(display),
        if (level == "item") "Item" else c("Item", "Threshold")
      )
      for (cc in data_cols) {
        col <- display[[cc]]
        if (is.numeric(col)) {
          col <- format(col, nsmall = 3, trim = TRUE)
        }
        display[[cc]] <- ifelse(
          flag,
          paste0("**", col, "**"),
          as.character(col)
        )
      }
    }
  }

  caption <- paste0(
    if (is_polytomous) "Partial Credit Model" else "Rasch Model",
    " split by DIF variable (",
    lr_summary$n_groups,
    " groups). ",
    "Andersen LR \u03c7\u00b2 = ",
    round(lr_summary$LR, 3),
    ", df = ",
    lr_summary$df,
    ", ",
    .format_p(lr_summary$p_value),
    ". ",
    n_clause,
    ", ",
    lr_summary$n_items,
    " items.",
    if (!is.null(cutoff)) {
      paste0(" **Bold** = MaxDiff > ", cutoff, " logits.")
    } else {
      ""
    }
  )

  knitr::kable(display, format = "pipe", row.names = FALSE, caption = caption)
}

#' Coerce maybe-character flag column back to logical for rendering
#'
#' @keywords internal
#' @noRd
isTRUE_vec <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  vapply(x, isTRUE, logical(1L))
}

#' Render a ggplot panel of per-group locations with confidence bars
#'
#' @keywords internal
#' @noRd
render_lr_ggplot <- function(
  plot_df,
  lr_summary,
  level,
  groups,
  conf,
  is_polytomous
) {
  z <- stats::qnorm(1 - (1 - conf) / 2)

  plot_df$DIFgroup <- factor(plot_df$DIFgroup, levels = c(groups, "All"))
  group_only <- plot_df[plot_df$DIFgroup %in% groups, , drop = FALSE]
  all_only <- plot_df[plot_df$DIFgroup == "All", , drop = FALSE]

  group_only$DIFgroup <- droplevels(group_only$DIFgroup)
  diamond_x <- (length(groups) + 1) / 2 + 0.15

  caption <- er2_caption(paste0(
    "Andersen LR \u03c7\u00b2 = ",
    round(lr_summary$LR, 3),
    ", df = ",
    lr_summary$df,
    ", ",
    .format_p(lr_summary$p_value),
    ". Error bars: ",
    round(conf * 100),
    "% CI. ",
    "Black diamond = whole-sample location."
  ))

  if (level == "item") {
    p <- ggplot2::ggplot(
      group_only,
      ggplot2::aes(
        x = .data$DIFgroup,
        y = .data$Location,
        group = .data$Item,
        colour = .data$Item
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data$Location - z * .data$SE,
          ymax = .data$Location + z * .data$SE
        ),
        width = 0.1
      ) +
      ggplot2::geom_point(
        data = all_only,
        ggplot2::aes(x = diamond_x, y = .data$Location),
        shape = 18,
        colour = "black",
        size = 3,
        alpha = 0.6,
        inherit.aes = FALSE
      ) +
      ggplot2::facet_wrap(~Item) +
      ggplot2::labs(
        title = "DIF: item locations by group",
        subtitle = "Item locations are means of threshold locations",
        x = "DIF group",
        y = "Item location (logits)",
        caption = caption
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        panel.spacing = ggplot2::unit(0.7, "cm")
      ) +
      er2_axis_margins() +
      er2_plot_caption()
    return(p)
  }

  # level == "threshold"
  p <- ggplot2::ggplot(
    group_only,
    ggplot2::aes(
      x = .data$DIFgroup,
      y = .data$Location,
      group = .data$Threshold,
      colour = .data$Threshold
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$Location - z * .data$SE,
        ymax = .data$Location + z * .data$SE
      ),
      width = 0.1
    ) +
    ggplot2::geom_errorbar(
      data = all_only,
      ggplot2::aes(
        x = diamond_x,
        ymin = .data$Location - z * .data$SE,
        ymax = .data$Location + z * .data$SE
      ),
      width = 0.1,
      colour = "darkgrey",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = all_only,
      ggplot2::aes(x = diamond_x, y = .data$Location),
      shape = 18,
      colour = "black",
      size = 3,
      alpha = 0.6,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~Item) +
    ggplot2::labs(
      title = "DIF: threshold locations by group",
      x = "DIF group",
      y = "Threshold location (logits)",
      caption = caption
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.spacing = ggplot2::unit(0.7, "cm")
    ) +
    er2_axis_margins() +
    er2_plot_caption()

  p
}
