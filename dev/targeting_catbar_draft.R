## ---------------------------------------------------------------------------
## Draft: category-band ("stacked bar") alternative to the bottom panel of
## RMtargeting(). Nothing here modifies easyRasch2.
##
## Each item is one horizontal bar spanning the logit scale, partitioned into
## response-category bands. Boundaries are placed at the MODAL-CATEGORY
## transitions (method "modal"), which reduce to the Andrich thresholds when
## those are ordered and collapse never-modal categories to a tick when they
## are not. "andrich" and "thurstonian" are kept for comparison only.
##
## ci_style controls how the threshold uncertainty is drawn:
##   "errorbar" capped error bars below the bar, coloured by the category the
##              threshold gives entry to
##   "twotone"  as above, but each half takes the colour of the category on
##              that side of the threshold
##   "stagger"  capped error bars, offset vertically by threshold index so
##              overlapping intervals never collide
##   "inbar"    no separate row: a bracket drawn on the bar itself, spanning
##              the CI, anchored at the join it belongs to
##   "blur"     no separate row: the band join itself becomes a gradient whose
##              width is the confidence interval
##   "none"
## ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

OUT <- Sys.getenv("CATBAR_OUT", unset = "figs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else grDevices::png

## --- threshold extraction (mirrors RMtargeting's CML path) -----------------

andrich_thresholds <- function(data) {
  fit <- psychotools::pcmodel(data)
  tp  <- psychotools::threshpar(fit, vcov = TRUE)
  se  <- sqrt(diag(attr(tp, "vcov")))
  loc <- unlist(lapply(tp, as.numeric), use.names = FALSE)
  loc <- loc - mean(loc, na.rm = TRUE)
  data.frame(
    Item      = factor(rep(names(tp), lengths(tp)), levels = names(tp)),
    Threshold = unlist(lapply(tp, seq_along), use.names = FALSE),
    Location  = loc,
    SE        = as.numeric(se),
    row.names = NULL
  )
}

## --- boundary machinery ----------------------------------------------------

## Modal-category boundaries.
## For a PCM item with Andrich thresholds t_1..t_m, log P(X = k | theta) is
## k*theta - S_k + const, with S_k = sum_{j <= k} t_j and S_0 = 0. Category k
## is the most likely response where the line k*theta - S_k is the upper
## envelope of the m+1 lines, i.e. where (k, S_k) lies on the LOWER CONVEX
## HULL of {(k, S_k)}. The boundary between two consecutive hull vertices
## k < l falls at (S_l - S_k)/(l - k), the mean of the thresholds they span.
## So a disordered run is pooled by averaging, and the categories skipped over
## are those that are never the most likely response. Closed form, no grid.
modal_boundaries <- function(tau) {
  S <- c(0, cumsum(tau))
  k <- seq_along(S) - 1L
  hull <- 1L
  for (i in 2:length(S)) {
    while (length(hull) >= 2L) {
      a <- hull[length(hull) - 1L]; b <- hull[length(hull)]
      if ((S[b] - S[a]) * (k[i] - k[a]) >= (S[i] - S[a]) * (k[b] - k[a])) {
        hull <- hull[-length(hull)]
      } else break
    }
    hull <- c(hull, i)
  }
  vert <- k[hull]
  list(cats    = vert,
       dropped = setdiff(k, vert),
       bounds  = diff(S[hull]) / diff(k[hull]))
}

## Thurstonian thresholds: theta where P(X >= k) = 0.5.
thurstonian <- function(tau) {
  m <- length(tau)
  pk <- function(th) {
    lin <- c(0, cumsum(th - tau)); p <- exp(lin - max(lin)); p / sum(p)
  }
  vapply(seq_len(m), function(k) {
    stats::uniroot(function(th) sum(pk(th)[(k + 1L):(m + 1L)]) - 0.5,
                   interval = c(-30, 30), tol = 1e-8)$root
  }, numeric(1))
}

item_bands <- function(tau, xlim, method) {
  m <- length(tau)
  if (method == "andrich") {
    b <- tau; cats <- 0:m; dropped <- integer(0)
  } else if (method == "thurstonian") {
    b <- thurstonian(tau); cats <- 0:m; dropped <- integer(0)
  } else {
    mb <- modal_boundaries(tau)
    b <- mb$bounds; cats <- mb$cats; dropped <- mb$dropped
  }
  list(bands = data.frame(cat = cats, xmin = c(xlim[1], b), xmax = c(b, xlim[2])),
       dropped = dropped)
}

band_data <- function(thr, xlim, method) {
  sp <- split(thr$Location, thr$Item)
  bands <- list(); drops <- list()
  for (it in names(sp)) {
    r <- item_bands(sp[[it]], xlim, method)
    r$bands$Item <- it
    bands[[it]] <- r$bands
    if (length(r$dropped)) {
      cum <- c(0, cumsum(sp[[it]]))
      drops[[it]] <- data.frame(
        Item = it, cat = r$dropped,
        x = vapply(r$dropped, function(cc) {
          keep <- r$bands$cat
          lo <- max(keep[keep < cc]); hi <- min(keep[keep > cc])
          (cum[hi + 1L] - cum[lo + 1L]) / (hi - lo)
        }, numeric(1)))
    }
  }
  list(bands = do.call(rbind, bands),
       drops = if (length(drops)) do.call(rbind, drops) else NULL)
}

disorder_table <- function(thr) {
  sp <- split(thr$Location, thr$Item)
  do.call(rbind, lapply(names(sp), function(it) {
    tau <- sp[[it]]; d <- which(diff(tau) < 0)
    if (!length(d)) return(NULL)
    data.frame(Item = it, k = d, lo = tau[d + 1L], hi = tau[d],
               gap = tau[d] - tau[d + 1L], row.names = NULL)
  }))
}

## --- palettes --------------------------------------------------------------
pal_mako    <- function(n) viridisLite::mako(n, begin = 0.90, end = 0.20)
pal_viridis <- function(n) viridisLite::viridis(n, end = 0.92, direction = -1)

## Pale category fills vanish as a thin line on white, so darken them until
## they clear a luminance floor. Colours that already pass are untouched.
ci_pal <- function(fills, floor = 130) {
  vapply(fills, function(cl) {
    v <- as.numeric(grDevices::col2rgb(cl)[, 1])
    while (sum(v * c(.299, .587, .114)) > floor) v <- v * 0.78
    grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
  }, character(1), USE.NAMES = FALSE)
}

## colour shared by the collapsed-category tick and the reversal spans
rev_col <- "#B2182B"

## "Note." caption in the package's style (er2_caption() does this in easyRasch2)
er2_note <- function(txt, width = 115) {
  paste0("Note. ", paste(strwrap(paste(txt, collapse = " "), width = width),
                         collapse = "\n"))
}

lab_col <- function(fills) ifelse(
  apply(grDevices::col2rgb(fills), 2, function(v) sum(v * c(.299, .587, .114))) > 145,
  "grey15", "white")

## --- the panel -------------------------------------------------------------

catbar_panel <- function(thr,
                         xlim       = c(-4, 4),
                         method     = c("modal", "andrich", "thurstonian"),
                         cat_labels = NULL,
                         labels_in  = c("number", "text", "none"),
                         palette    = pal_mako,
                         ci_style   = c("errorbar", "twotone", "stagger",
                                        "inbar", "blur", "none"),
                         ci_level   = 0.95,
                         fade_edges = TRUE,
                         sort_items = c("data", "location"),
                         bar_h      = 0.52,
                         row_gap    = NULL,
                         title      = NULL,
                         caption    = NULL) {
  method     <- match.arg(method)
  labels_in  <- match.arg(labels_in)
  ci_style   <- match.arg(ci_style)
  sort_items <- match.arg(sort_items)

  ## The collapsed-category label needs headroom above the bar. Leave row_gap
  ## at NULL and it opens up only when there is something to label.
  if (is.null(row_gap)) {
    row_gap <- if (any(vapply(split(thr$Location, thr$Item),
                              function(tau) length(modal_boundaries(tau)$dropped) > 0L,
                              logical(1)))) 1.18 else 1
  }

  ncat  <- max(thr$Threshold) + 1L
  if (is.null(cat_labels)) cat_labels <- as.character(0:(ncat - 1L))
  fills <- palette(ncat)

  lev <- if (sort_items == "location") {
    rev(names(sort(tapply(thr$Location, thr$Item, mean))))
  } else rev(levels(thr$Item))
  yof <- function(x) as.numeric(factor(as.character(x), levels = lev)) * row_gap

  thr$Item <- factor(as.character(thr$Item), levels = lev)
  bd <- band_data(thr, xlim, method)
  b  <- bd$bands
  b$y <- yof(b$Item); b$fill <- fills[b$cat + 1L]

  z  <- stats::qnorm(1 - (1 - ci_level) / 2)
  t2 <- thr
  t2$y    <- yof(t2$Item)
  t2$lo   <- t2$Location - z * t2$SE
  t2$hi   <- t2$Location + z * t2$SE
  cfills    <- ci_pal(fills)
  t2$col_hi <- cfills[t2$Threshold + 1L]  # category entered at this threshold
  t2$col_lo <- cfills[t2$Threshold]       # category left behind

  hh <- bar_h / 2
  p  <- ggplot()

  ## open-ended outer bands fade out
  if (fade_edges) {
    is_out <- b$cat == 0 | b$cat == max(b$cat)
    ns <- 40
    fade <- do.call(rbind, lapply(which(is_out), function(i) {
      r <- b[i, ]; br <- seq(r$xmin, r$xmax, length.out = ns + 1L)
      data.frame(y = r$y, fill = r$fill,
                 alpha = if (r$cat == 0) seq(0.22, 1, length.out = ns)
                         else seq(1, 0.22, length.out = ns),
                 xmin = br[-length(br)], xmax = br[-1])
    }))
    p <- p +
      geom_rect(data = fade, aes(xmin = xmin, xmax = xmax, ymin = y - hh,
                                 ymax = y + hh, fill = I(fill), alpha = I(alpha))) +
      geom_rect(data = b[!is_out, ], aes(xmin = xmin, xmax = xmax, ymin = y - hh,
                                         ymax = y + hh, fill = I(fill)))
  } else {
    p <- p + geom_rect(data = b, aes(xmin = xmin, xmax = xmax, ymin = y - hh,
                                     ymax = y + hh, fill = I(fill)))
  }

  ## band joins: hard white hairline, or a CI-wide gradient under "blur"
  jn <- b[b$xmin > xlim[1], ]
  if (ci_style == "blur") {
    ## match each modal boundary to the threshold(s) it came from and use the
    ## widest of their CIs as the transition width
    jn$key <- paste(jn$Item, round(jn$xmin, 8))
    t2$key <- paste(t2$Item, round(t2$Location, 8))
    w <- vapply(seq_len(nrow(jn)), function(i) {
      cand <- t2[t2$Item == jn$Item[i], ]
      k <- which.min(abs(cand$Location - jn$xmin[i]))
      z * cand$SE[k]
    }, numeric(1))
    ns <- 60
    grad <- do.call(rbind, lapply(seq_len(nrow(jn)), function(i) {
      r <- jn[i, ]
      lo <- r$xmin - w[i]; hi <- r$xmin + w[i]
      br <- seq(lo, hi, length.out = ns + 1L)
      c_lo <- fills[r$cat]        # band to the left
      c_hi <- fills[r$cat + 1L]   # this band
      ramp <- grDevices::colorRampPalette(c(c_lo, c_hi))(ns)
      data.frame(y = r$y, xmin = br[-length(br)], xmax = br[-1], fill = ramp)
    }))
    p <- p + geom_rect(data = grad, show.legend = FALSE,
                       aes(xmin = xmin, xmax = xmax, ymin = y - hh,
                           ymax = y + hh, fill = I(fill)))
  } else {
    p <- p + geom_segment(data = jn, show.legend = FALSE,
                          aes(x = xmin, xend = xmin, y = y - hh, yend = y + hh),
                          colour = "white", linewidth = 0.6)
  }

  ## in-band labels where the band is wide enough
  if (labels_in != "none") {
    lb <- b
    lb$txt <- if (labels_in == "text") cat_labels[lb$cat + 1L] else as.character(lb$cat)
    x0 <- pmax(lb$xmin, xlim[1]); x1 <- pmin(lb$xmax, xlim[2])
    lb$xc <- (x0 + x1) / 2
    need <- if (labels_in == "text") 0.088 * nchar(lb$txt) + 0.15 else 0.30
    lb <- lb[(x1 - x0) > need, ]
    p <- p + geom_text(data = lb, aes(x = xc, y = y, label = txt),
                       colour = lab_col(lb$fill), size = 3.1, show.legend = FALSE)
  }

  ## collapsed (never-modal) categories: tick in that category's own colour
  if (!is.null(bd$drops)) {
    dr <- bd$drops; dr$y <- yof(dr$Item)
    dr$fill <- fills[dr$cat + 1L]; dr$tcol <- cfills[dr$cat + 1L]
    ## several collapsed categories can land on the same boundary, so nudge
    ## their labels apart
    dr <- dr[order(dr$Item, dr$cat), ]
    key <- paste(dr$Item, round(dr$x, 6))
    dr$lab <- stats::ave(dr$cat, key, FUN = function(v) NA)  # placeholder
    dr$lab <- vapply(key, function(kk) paste(dr$cat[key == kk], collapse = ","),
                     character(1), USE.NAMES = FALSE)
    dr <- dr[!duplicated(key), ]
    ## red, to tie the tick to the reversal span that explains it. A white
    ## halo keeps it visible against the dark end of the palette.
    p <- p +
      geom_segment(data = dr, show.legend = FALSE,
                   aes(x = x, xend = x, y = y - hh - 0.06, yend = y + hh + 0.06),
                   colour = "white", linewidth = 2.8, lineend = "butt") +
      geom_segment(data = dr, show.legend = FALSE,
                   aes(x = x, xend = x, y = y - hh - 0.04, yend = y + hh + 0.04),
                   colour = rev_col, linewidth = 1.7, lineend = "butt") +
      geom_text(data = dr, show.legend = FALSE,
                aes(x = x, y = y + hh + 0.14, label = lab),
                colour = rev_col, size = 2.7, fontface = "bold",
                hjust = 0.5, vjust = 0)
  }

  ## --- threshold uncertainty ------------------------------------------------
  if (ci_style %in% c("errorbar", "twotone", "stagger")) {
    cap <- 0.085
    if (ci_style == "stagger") {
      t2$yy <- t2$y - hh - 0.14 - (t2$Threshold - 1L) * 0.105
      cap <- 0.052
    } else {
      t2$yy <- t2$y - hh - 0.17
    }
    ## white halo separates intervals that overlap each other
    p <- p + geom_errorbar(
      data = t2, show.legend = FALSE,
      aes(y = yy, xmin = lo, xmax = hi),
      width = cap * 2, linewidth = 2.1, colour = "white", orientation = "y")

    if (ci_style == "twotone") {
      p <- p +
        geom_errorbar(data = t2, show.legend = FALSE,
                      aes(y = yy, xmin = lo, xmax = Location, colour = I(col_lo)),
                      width = cap * 2, linewidth = 0.8, orientation = "y") +
        geom_errorbar(data = t2, show.legend = FALSE,
                      aes(y = yy, xmin = Location, xmax = hi, colour = I(col_hi)),
                      width = cap * 2, linewidth = 0.8, orientation = "y")
    } else {
      p <- p + geom_errorbar(
        data = t2, show.legend = FALSE,
        aes(y = yy, xmin = lo, xmax = hi, colour = I(col_hi)),
        width = cap * 2, linewidth = 0.8, orientation = "y")
    }
    p <- p + geom_point(data = t2, aes(x = Location, y = yy, colour = I(col_hi)),
                        show.legend = FALSE, size = 1.15)
  }

  ## bracket drawn on the bar itself: no extra row, and each interval stays
  ## anchored to the join it belongs to
  if (ci_style == "inbar") {
    p <- p +
      geom_errorbar(data = t2, show.legend = FALSE,
                    aes(y = y, xmin = lo, xmax = hi),
                    width = bar_h * 0.62, linewidth = 2.0, colour = "grey15",
                    orientation = "y") +
      geom_errorbar(data = t2, show.legend = FALSE,
                    aes(y = y, xmin = lo, xmax = hi),
                    width = bar_h * 0.62, linewidth = 0.8, colour = "white",
                    orientation = "y")
  }

  ## reversal span, above the bar
  dis <- disorder_table(thr)
  if (!is.null(dis) && nrow(dis)) {
    dis <- dis[order(dis$Item, dis$k), ]
    dis$y   <- yof(dis$Item)
    dis$row <- stats::ave(dis$k, dis$Item, FUN = seq_along)   # stack per item
    ## sits under the CI row, next to the two intervals that are out of order
    dis$yy  <- dis$y - hh - 0.32 - (dis$row - 1L) * 0.15
    p <- p +
      geom_segment(data = dis, show.legend = FALSE,
                   aes(x = lo, xend = hi, y = yy, yend = yy),
                   colour = rev_col, linewidth = 0.8,
                   arrow = arrow(ends = "both", length = unit(3, "pt"),
                                 type = "closed")) +
      geom_text(data = dis, show.legend = FALSE,
                aes(x = hi, y = yy, label = sprintf("%.2f", gap)),
                colour = rev_col, size = 2.5, hjust = -0.3, vjust = 0.42)
  }

  ## default caption, describing only the elements actually drawn
  if (is.null(caption)) {
    bits <- c(switch(method,
      modal = "Bands span the locations at which each response category is the most likely response.",
      andrich = "Band boundaries are the Andrich thresholds.",
      thurstonian = "Band boundaries are the Thurstonian thresholds, the locations at which P(X >= k) = 0.5."))
    bits <- c(bits, "The outer bands fade because they are open-ended.")
    if (ci_style %in% c("errorbar", "twotone", "stagger")) {
      bits <- c(bits, sprintf(paste("Points and intervals below each bar are the Andrich",
                                    "thresholds with %d%% confidence intervals, coloured by",
                                    "the category the threshold gives entry to."),
                              round(ci_level * 100)))
    }
    if (!is.null(bd$drops)) {
      bits <- c(bits, paste("A red tick marks a category that is never the most likely",
                            "response at any location, labelled above with the category",
                            "number."))
    }
    if (!is.null(dis) && nrow(dis)) {
      bits <- c(bits, paste("Red arrows give the size in logits of each reversal between",
                            "adjacent thresholds."))
    }
    caption <- er2_note(bits)
  }

  ## legend keyed off cat_labels
  leg <- data.frame(lab = factor(cat_labels, levels = cat_labels))
  p <- p +
    geom_point(data = leg, aes(x = xlim[1] - 100, y = 1, colour = lab), size = 0) +
    scale_colour_manual(values = stats::setNames(fills, levels(leg$lab)),
                        name = NULL, drop = FALSE,
                        guide = guide_legend(nrow = 1, override.aes = list(
                          size = 4.6, shape = 22, fill = fills, stroke = 0)))

  nrev <- if (!is.null(dis) && nrow(dis)) max(table(dis$Item)) else 0L
  ybot <- if (ci_style == "stagger") 0.55 + 0.105 * (max(t2$Threshold) - 1L)
          else max(0.55, 0.32 + 0.15 * nrev + 0.15)
  ytop <- max(0.55, hh + 0.34)

  p +
    scale_y_continuous(breaks = seq_along(lev) * row_gap, labels = lev,
                       expand = expansion(add = c(ybot, ytop))) +
    scale_x_continuous(breaks = seq(xlim[1], xlim[2], 1), expand = c(0, 0)) +
    coord_cartesian(xlim = xlim) +
    labs(x = "Location (logit scale)", y = NULL, title = title, caption = caption) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 11),
          plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey25"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(margin = margin(t = 10)))
}

sv <- function(file, plot, w = 9, h = 5.2) {
  ggsave(file.path(OUT, file), plot, width = w, height = h, dpi = 160,
         bg = "white", device = png_dev)
}

## ===========================================================================
## Demos
## ===========================================================================

phq_labs <- c("Not at all", "Several days",
              "More than half the days", "Nearly every day")

data(phq9, package = "easyRasch2")
thr_phq <- andrich_thresholds(phq9[, 1:9])
pl <- easyRasch2::RMtargeting(phq9[, 1:9], output = "list")

full <- function(ci, ttl = NULL) {
  pl$p1 / pl$p2 /
    catbar_panel(thr_phq, cat_labels = phq_labs, ci_style = ci, title = ttl) +
    patchwork::plot_layout(heights = c(3, 2, 5))
}

## 1. as requested: capped error bars, coloured by the category entered
sv("01_full_errorbar.png", full("errorbar"), h = 8.6)

## 2. two-tone: each half takes the colour of the category on that side
sv("02_full_twotone.png", full("twotone"), h = 8.6)

## 3. staggered by threshold index, so intervals cannot overlap at all
sv("03_full_stagger.png", full("stagger"), h = 9.2)

## 4. no separate row: bracket on the bar at each join
sv("04_full_inbar.png", full("inbar"), h = 8.6)

## 5. no separate row: the join itself becomes a CI-wide gradient
sv("05_full_blur.png", full("blur"), h = 8.6)

## 5. the four side by side, bands panel only
sv("06_ci_styles.png",
   (catbar_panel(thr_phq, cat_labels = phq_labs, ci_style = "errorbar",
                 title = "A. Capped error bars, coloured by the category entered") /
    catbar_panel(thr_phq, cat_labels = phq_labs, ci_style = "twotone",
                 title = "B. Two-tone: colour of the category on each side of the threshold") /
    catbar_panel(thr_phq, cat_labels = phq_labs, ci_style = "stagger",
                 title = "C. Staggered by threshold index: overlap impossible") /
    catbar_panel(thr_phq, cat_labels = phq_labs, ci_style = "inbar",
                 title = "D. Bracket on the bar: one row per item, CI anchored at its own join") /
    catbar_panel(thr_phq, cat_labels = phq_labs, ci_style = "blur",
                 title = "E. No CI row: the band join is a gradient as wide as the CI")) +
     patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom"),
   h = 23)

## --- disordered example ----------------------------------------------------

sim_pcm <- function(n, taus, seed = 1) {
  set.seed(seed)
  th <- rnorm(n)
  d <- as.data.frame(sapply(taus, function(tau) {
    vapply(th, function(x) {
      lin <- c(0, cumsum(x - tau)); p <- exp(lin - max(lin)); p <- p / sum(p)
      sample.int(length(p), 1L, prob = p) - 1L
    }, integer(1))
  }))
  names(d) <- names(taus); d
}

taus_dis <- list(
  i1 = c(-2.0, -1.0,  0.6), i2 = c(-1.4, -0.5,  0.9), i3 = c(-0.3,  0.2,  1.1),
  i4 = c(-1.2, -1.6,  0.4),   # mild reversal of t1/t2
  i5 = c( 0.9, -0.4,  1.4),   # strong reversal of t1/t2
  i6 = c(-0.8,  0.3,  1.6)
)
thr_dis <- andrich_thresholds(sim_pcm(1200, taus_dis, seed = 7))

sv("07_disordered_errorbar.png",
   catbar_panel(thr_dis, ci_style = "errorbar",
                caption = paste("i4 and i5 have reversed thresholds 1 and 2.",
                                "The error bars sit out of order, and category 1",
                                "is never the most likely response,\nso its band",
                                "collapses to a tick.")),
   h = 5.4)
sv("08_disordered_stagger.png", catbar_panel(thr_dis, ci_style = "stagger"), h = 5.8)
sv("09_disordered_inbar.png", catbar_panel(thr_dis, ci_style = "inbar"), h = 5.2)

## --- 7-category stress test ------------------------------------------------
## i5: one reversal  -> category 3 never modal
## i6: two separate reversals -> categories 1 and 4 never modal
## i7: a run of three out of order -> categories 2 and 3 never modal (adjacent)

taus7 <- list(
  i1 = c(-2.4, -1.5, -0.7,  0.1,  0.9,  1.9),
  i2 = c(-2.0, -1.2, -0.4,  0.3,  1.1,  2.1),
  i3 = c(-2.6, -1.7, -0.9,  0.0,  0.8,  1.7),
  i4 = c(-1.8, -1.0, -0.2,  0.5,  1.3,  2.3),
  i5 = c(-2.2, -1.4,  0.2, -0.4,  1.0,  2.0),
  i6 = c(-1.4, -2.0, -0.6,  0.9,  0.2,  1.8),
  i7 = c(-2.2,  0.6,  0.0, -0.6,  1.2,  2.2)
)

sim_pcm7 <- function(n, taus, sd = 1.5, seed = 11) {
  set.seed(seed)
  th <- rnorm(n, sd = sd)
  d <- as.data.frame(sapply(taus, function(tau) {
    vapply(th, function(x) {
      lin <- c(0, cumsum(x - tau)); p <- exp(lin - max(lin)); p <- p / sum(p)
      sample.int(length(p), 1L, prob = p) - 1L
    }, integer(1))
  }))
  names(d) <- names(taus); d
}

dat7 <- sim_pcm7(2000, taus7)
cat("\n7-category example, per-item category counts:\n")
print(sapply(dat7, function(x) tabulate(x + 1L, 7L)))

thr7 <- andrich_thresholds(dat7)
cat("\nReversals:\n"); print(disorder_table(thr7))
cat("\nNever-modal categories:\n")
for (it in names(taus7)) {
  d <- modal_boundaries(thr7$Location[thr7$Item == it])$dropped
  cat(sprintf("  %s: %s\n", it, if (length(d)) paste(d, collapse = ", ") else "-"))
}

labs7 <- c("Never", "Almost never", "Rarely", "Sometimes",
           "Often", "Almost always", "Always")
pl7 <- easyRasch2::RMtargeting(dat7, output = "list")

sv("10_seven_cat_full.png",
   pl7$p1 / pl7$p2 /
     catbar_panel(thr7, cat_labels = labs7, ci_style = "errorbar",
                  row_gap = 1.25) +
     patchwork::plot_layout(heights = c(3, 2, 6)),
   w = 10, h = 9.4)
