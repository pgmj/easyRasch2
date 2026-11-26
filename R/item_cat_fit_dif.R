# Internal: RMitemCatFitPlot() with `dif_var`.
#
# Each group's curve is drawn against predictions from one CML fit to all
# respondents. The band comes from permutations of the group labels within
# total scores. Under the Rasch model without DIF, an item is independent of
# the group given the total score, so the permutation null is exact whatever
# the group ability distributions (Kreiner, 1987; Kreiner & Christensen,
# 2011). A permutation keeps every group's count at every score, so each
# group's curve has the same points in every permutation and nothing is
# refitted. The band is simultaneous for each curve. The study behind it is in
# dev/catfit_dif/ (PLAN.md, "Validating the permutation bands").

#' Okabe-Ito colours, line types and shapes for up to four DIF groups
#' @keywords internal
#' @noRd
.catdif_styles <- function() {
  list(colour = c("#E69F00", "#0072B2", "#009E73", "#CC79A7"),
       linetype = c("solid", "dashed", "dotted", "dotdash"),
       shape = c(21L, 24L, 22L, 23L))
}

#' Data for the DIF path: event indicators and predictions per score
#'
#' @param prep Output of `.catfit_prepare()`.
#' @param grp Integer group (1..G) per row of `prep$X`.
#' @return List used by the layout and band functions.
#' @keywords internal
#' @noRd
.catdif_prepare <- function(prep, grp) {
  X <- prep$X
  r <- prep$r
  m <- lengths(prep$thr)
  ev_item <- rep(seq_along(m), m)
  ev_k <- sequence(m)
  Y <- vapply(seq_along(ev_item), function(e) {
    as.integer(X[, ev_item[e]] >= ev_k[e])
  }, integer(nrow(X)))
  Y <- matrix(Y, nrow = nrow(X))
  K <- nrow(prep$P[[1L]]) - 1L
  Pev <- vapply(seq_along(ev_item), function(e) {
    Pi <- prep$P[[ev_item[e]]]
    rowSums(Pi[, (ev_k[e] + 1L):ncol(Pi), drop = FALSE])
  }, numeric(K + 1L))
  Pev <- matrix(Pev, nrow = K + 1L)
  sc <- sort(unique(r))
  list(X = X, Y = Y, r = r, grp = grp, G = max(grp), sc = sc,
       s_idx = match(r, sc), Pev = Pev[sc + 1L, , drop = FALSE],
       ev_item = ev_item, ev_k = ev_k, m = m, P = prep$P)
}

#' Respondent counts per group x score (G x S) and column sums of `M` per
#' group x score ((G * S) x ncol(M), row (g - 1) * S + s)
#' @keywords internal
#' @noRd
.catdif_counts <- function(M, s_idx, grp, G, S) {
  key <- (grp - 1L) * S + s_idx
  sums <- matrix(0, G * S, ncol(M))
  rs <- rowsum(M, key, reorder = TRUE)
  sums[as.integer(rownames(rs)), ] <- rs
  list(n = matrix(tabulate(key, G * S), G, S, byrow = TRUE), sums = sums)
}

#' Weighted pool-adjacent-violators on points ordered by the predictor
#'
#' Agrees with `.pav()` on the individual data when each point is a score
#' (checked in the tests).
#' @param v Values (shares) in predictor order.
#' @param w Positive weights (respondents).
#' @return The fitted value at each point.
#' @keywords internal
#' @noRd
.pav_w <- function(v, w) {
  J <- length(v)
  bv <- numeric(J)
  bw <- numeric(J)
  bl <- integer(J)
  t <- 0L
  for (j in seq_len(J)) {
    t <- t + 1L
    bv[t] <- v[j]
    bw[t] <- w[j]
    bl[t] <- 1L
    while (t > 1L && bv[t - 1L] > bv[t]) {
      bv[t - 1L] <- (bv[t - 1L] * bw[t - 1L] + bv[t] * bw[t]) /
        (bw[t - 1L] + bw[t])
      bw[t - 1L] <- bw[t - 1L] + bw[t]
      bl[t - 1L] <- bl[t - 1L] + bl[t]
      t <- t - 1L
    }
  }
  rep(bv[seq_len(t)], bl[seq_len(t)])
}

#' Permute group labels within total scores
#' @keywords internal
#' @noRd
.catdif_permute <- function(grp, s_idx) {
  ord <- order(s_idx)
  out <- grp
  out[ord] <- grp[ord][order(s_idx[ord], stats::runif(length(grp)))]
  out
}

# -- Layouts: fixed points of each curve and a value function of the labels --

#' Calibration: item x event x group, one point per score with respondents
#' @keywords internal
#' @noRd
.catdif_calibration <- function(pr) {
  S <- length(pr$sc)
  n0 <- .catdif_counts(pr$Y, pr$s_idx, pr$grp, pr$G, S)$n
  lay <- do.call(rbind, lapply(seq_along(pr$ev_item), function(e) {
    do.call(rbind, lapply(seq_len(pr$G), function(g) {
      s <- which(n0[g, ] > 0)
      if (length(s) == 0L) return(NULL)
      data.frame(item = pr$ev_item[e], k = pr$ev_k[e], event = e, group = g,
                 s = s, score = pr$sc[s], n = n0[g, s],
                 predicted = pr$Pev[s, e])
    }))
  }))
  lay$curve <- (lay$event - 1L) * pr$G + lay$group
  fun <- function(grp) {
    a <- .catdif_counts(pr$Y, pr$s_idx, grp, pr$G, S)
    out <- numeric(nrow(lay))
    for (i in split(seq_len(nrow(lay)), lay$curve)) {
      g <- lay$group[i[1L]]
      w <- a$n[g, lay$s[i]]
      out[i] <- .pav_w(a$sums[(g - 1L) * S + lay$s[i], lay$event[i[1L]]] / w,
                       w)
    }
    out
  }
  list(lay = lay, fun = fun)
}

#' Adjacent categories: item x pair x group, one point per score where the
#' group has respondents in the pair
#' @keywords internal
#' @noRd
.catdif_adjacent <- function(pr) {
  S <- length(pr$sc)
  m <- pr$m
  cat_item <- rep(seq_along(m), m + 1L)
  cat_c <- sequence(m + 1L) - 1L
  Cm <- vapply(seq_along(cat_item), function(j) {
    as.integer(pr$X[, cat_item[j]] == cat_c[j])
  }, integer(nrow(pr$X)))
  Cm <- matrix(Cm, nrow = nrow(pr$X))
  col_of <- function(i, c) which(cat_item == i & cat_c == c)
  n0 <- .catdif_counts(Cm, pr$s_idx, pr$grp, pr$G, S)$n
  lay <- do.call(rbind, lapply(seq_along(m), function(i) {
    Pi <- pr$P[[i]][pr$sc + 1L, , drop = FALSE]
    do.call(rbind, lapply(seq_len(m[i]), function(j) {
      p_exp <- Pi[, j + 1L] / (Pi[, j] + Pi[, j + 1L])
      do.call(rbind, lapply(seq_len(pr$G), function(g) {
        s <- which(n0[g, ] > 0)
        if (length(s) == 0L) return(NULL)
        data.frame(item = i, pair_j = j, group = g, s = s,
                   score = pr$sc[s], expected = p_exp[s],
                   lo_col = col_of(i, j - 1L), hi_col = col_of(i, j))
      }))
    }))
  }))
  value <- function(grp, lay) {
    Cg <- .catdif_counts(Cm, pr$s_idx, grp, pr$G, S)$sums
    row <- (lay$group - 1L) * S + lay$s
    hi <- Cg[cbind(row, lay$hi_col)]
    tot <- hi + Cg[cbind(row, lay$lo_col)]
    list(share = ifelse(tot > 0, hi / tot, NA_real_), n_pair = tot)
  }
  # Points with no respondents in the pair in the observed data are dropped.
  v0 <- value(pr$grp, lay)
  lay$n_pair <- v0$n_pair
  lay <- lay[!is.na(v0$share), , drop = FALSE]
  key <- paste(lay$item, lay$pair_j, lay$group)
  lay$curve <- match(key, unique(key))
  list(lay = lay, fun = function(grp) value(grp, lay)$share)
}

#' Residuals: item x group, one point per score group x category
#'
#' E and V come from the pooled conditional probabilities and the group's
#' counts per score, which a permutation leaves unchanged, so only the
#' observed counts move. Cells with V = 0 are dropped.
#' @param sg Integer score group (1..Q) per respondent.
#' @keywords internal
#' @noRd
.catdif_residuals <- function(pr, sg) {
  S <- length(pr$sc)
  Q <- max(sg)
  m <- pr$m
  sg_of_score <- sg[match(seq_len(S), pr$s_idx)]
  n0 <- .catdif_counts(pr$Y, pr$s_idx, pr$grp, pr$G, S)$n
  lay <- do.call(rbind, lapply(seq_along(m), function(i) {
    Pi <- pr$P[[i]][pr$sc + 1L, , drop = FALSE]
    do.call(rbind, lapply(seq_len(pr$G), function(g) {
      E <- rowsum(n0[g, ] * Pi, sg_of_score, reorder = TRUE)
      V <- rowsum(n0[g, ] * Pi * (1 - Pi), sg_of_score, reorder = TRUE)
      q <- as.integer(rownames(E))
      data.frame(item = i, group = g, q = rep(q, times = m[i] + 1L),
                 category = rep(0:m[i], each = length(q)),
                 expected = as.vector(E), V = as.vector(V))
    }))
  }))
  lay <- lay[lay$V > 0, , drop = FALSE]
  key <- paste(lay$item, lay$group)
  lay$curve <- match(key, unique(key))
  sizes <- pr$G * Q * (m + 1L)
  offset <- c(0L, cumsum(sizes))[seq_along(m)]
  flat <- offset[lay$item] +
    ((lay$group - 1L) * Q + lay$q - 1L) * (m[lay$item] + 1L) +
    lay$category + 1L
  observed <- function(grp) {
    key <- (grp - 1L) * Q + sg
    unlist(lapply(seq_along(m), function(i) {
      tabulate((key - 1L) * (m[i] + 1L) + pr$X[, i] + 1L, nbins = sizes[i])
    }))[flat]
  }
  lay$observed <- observed(pr$grp)
  list(lay = lay,
       fun = function(grp) (observed(grp) - lay$expected) / sqrt(lay$V))
}

# -- Bands --------------------------------------------------------------------

#' Simultaneous permutation band for each curve
#'
#' Each point is standardised by its permutation mean and SD. The curve's
#' statistic is the largest |z| over its points. The band at a point is
#' mean +- c SD, where c is the (m + 1)-th largest of the curve's permuted
#' maxima, m = floor(alpha (B + 1)) - 1. A point is outside when its |z|
#' exceeds c, which happens for some point of a curve exactly when the
#' curve's permutation p, (1 + #{T* >= T}) / (B + 1), is at most alpha. That
#' is the rule validated in dev/catfit_dif/.
#'
#' A point that is empty in a permutation (adjacent pairs) gets z = 0 there,
#' and its mean and SD come from the permutations where it exists. A point
#' that never moves (SD 0, e.g. a certain event) has a band of zero width
#' and is never outside.
#'
#' @param pr Output of `.catdif_prepare()`.
#' @param layout List(lay, fun) from a `.catdif_*()` layout function.
#' @param B Number of permutations.
#' @param conf_level Coverage of the band.
#' @return `layout$lay` with columns `value`, `perm_mean`, `lower`,
#'   `upper`, `outside`.
#' @keywords internal
#' @noRd
.catdif_bands <- function(pr, layout, B, conf_level) {
  lay <- layout$lay
  n_row <- nrow(lay)
  obs <- layout$fun(pr$grp)
  # vapply() returns the n_row x B matrix directly. The moments and maxima
  # are then computed one curve at a time, so the only full-size object is
  # this matrix (memory grows with rows x B, about 180 MB for 22,000 rows at
  # B = 1000).
  Pm <- vapply(seq_len(B), function(b) {
    layout$fun(.catdif_permute(pr$grp, pr$s_idx))
  }, numeric(n_row))
  if (n_row == 1L) Pm <- matrix(Pm, nrow = 1L)
  alpha <- 1 - conf_level
  m_allow <- max(floor(alpha * (B + 1)) - 1, 0)
  mu <- sdv <- crit <- numeric(n_row)
  for (i in split(seq_len(n_row), lay$curve)) {
    P <- Pm[i, , drop = FALSE]
    nb <- pmax(rowSums(!is.na(P)), 1)
    mu[i] <- rowSums(P, na.rm = TRUE) / nb
    D <- (P - mu[i])^2
    D[is.na(D)] <- 0
    # Centred SD: the mean-of-squares form leaves rounding noise where every
    # permutation gives the same value.
    sdv[i] <- sqrt(rowSums(D) / nb)
    Z <- abs(P - mu[i]) / sdv[i]
    Z[sdv[i] <= 1e-12 | is.na(Z)] <- 0
    Tm <- apply(Z, 2L, max)
    crit[i] <- sort(Tm, decreasing = TRUE)[m_allow + 1L]
  }
  moves <- sdv > 1e-12
  z_obs <- abs(obs - mu) / sdv
  z_obs[!moves | is.na(z_obs)] <- 0
  lay$value <- obs
  lay$perm_mean <- mu
  # The drawn band is clipped to the values the permutations reached (and
  # the observed value). For a cell the model expects to be almost empty,
  # mean +- c * SD reaches values the statistic cannot take. The clip never
  # changes `outside`: an observed value beyond mu + c * SD stays beyond the
  # clipped edge, and one inside stays inside.
  pmin_obs <- pmin(suppressWarnings(apply(Pm, 1L, min, na.rm = TRUE)), obs,
                   na.rm = TRUE)
  pmax_obs <- pmax(suppressWarnings(apply(Pm, 1L, max, na.rm = TRUE)), obs,
                   na.rm = TRUE)
  pmin_obs[!is.finite(pmin_obs)] <- -Inf
  pmax_obs[!is.finite(pmax_obs)] <- Inf
  lay$lower <- pmax(mu - crit * sdv, pmin_obs)
  lay$upper <- pmin(mu + crit * sdv, pmax_obs)
  lay$outside <- moves & z_obs > crit + 1e-12
  lay
}

# -- Plots --------------------------------------------------------------------

#' Group scales shared by the three DIF plots
#' @keywords internal
#' @noRd
.catdif_scales <- function(group_labels) {
  st <- .catdif_styles()
  G <- length(group_labels)
  list(
    ggplot2::scale_colour_manual(NULL, values = st$colour[seq_len(G)],
                                 breaks = group_labels, drop = FALSE),
    ggplot2::scale_fill_manual(NULL, values = st$colour[seq_len(G)],
                               breaks = group_labels, drop = FALSE),
    ggplot2::scale_linetype_manual(NULL, values = st$linetype[seq_len(G)],
                                   breaks = group_labels, drop = FALSE),
    ggplot2::scale_shape_manual(NULL, values = st$shape[seq_len(G)],
                                breaks = group_labels, drop = FALSE)
  )
}

#' @keywords internal
#' @noRd
.catdif_calibration_gg <- function(tab, strip, group_labels) {
  out <- tab[tab$outside, , drop = FALSE]
  ggplot2::ggplot(tab, ggplot2::aes(.data$predicted, group = .data$group)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40",
                         linetype = "dotted") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                                      fill = .data$group), alpha = 0.2) +
    ggplot2::geom_step(ggplot2::aes(y = .data$cep, colour = .data$group,
                                    linetype = .data$group),
                       linewidth = 0.6, direction = "mid") +
    ggplot2::geom_point(ggplot2::aes(y = .data$cep, size = .data$n,
                                     fill = .data$group, shape = .data$group),
                        colour = "grey20", stroke = 0.3) +
    ggplot2::geom_point(data = out, ggplot2::aes(y = .data$cep),
                        shape = 1, size = 4, colour = "black", stroke = 0.7) +
    .catdif_scales(group_labels) +
    ggplot2::scale_size_area("Respondents", max_size = 2.5) +
    ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1),
                                labels = c("0", ".5", "1")) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1),
                                labels = c("0", ".5", "1")) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$event)) +
    ggplot2::labs(x = "Predicted probability given the total score",
                  y = "Conditional event probability")
}

#' @keywords internal
#' @noRd
.catdif_adjacent_gg <- function(tab, strip, group_labels) {
  tab$pair <- paste("Categories", tab$pair)
  out <- tab[tab$outside, , drop = FALSE]
  expd <- unique(tab[, c("item", "pair", "score", "expected")])
  ggplot2::ggplot(tab, ggplot2::aes(.data$score, group = .data$group)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                                      fill = .data$group), alpha = 0.2) +
    ggplot2::geom_line(data = expd, ggplot2::aes(.data$score, .data$expected),
                       inherit.aes = FALSE, colour = "grey30",
                       linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(y = .data$observed, size = .data$n_pair,
                                     fill = .data$group, shape = .data$group),
                        colour = "grey20", stroke = 0.3) +
    ggplot2::geom_point(data = out, ggplot2::aes(y = .data$observed),
                        shape = 1, size = 4, colour = "black", stroke = 0.7) +
    .catdif_scales(group_labels) +
    ggplot2::scale_size_area("Respondents in the pair", max_size = 3) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$pair)) +
    ggplot2::labs(x = "Total score", y = "Share in the upper category")
}

#' @keywords internal
#' @noRd
.catdif_residuals_gg <- function(tab, strip, group_labels) {
  G <- length(group_labels)
  # Bars of one category stand side by side with a small gap. The tile
  # behind each bar (its band) is slightly wider than the bar, and
  # neighbouring tiles do not overlap.
  bar_w <- min(0.3, 0.7 / G)
  slot <- bar_w + 0.03
  dodge <- ggplot2::position_dodge(width = slot * G)
  # The y-axis covers the bars and at least +-4, which holds the usual bands.
  # A band that reaches further (a cell the model expects to be almost
  # empty) is cut at the panel edge.
  ylim <- range(c(tab$z, -4, 4), na.rm = TRUE)
  ggplot2::ggplot(tab, ggplot2::aes(factor(.data$category), .data$z,
                                    fill = .data$group,
                                    group = .data$group)) +
    ggplot2::geom_tile(
      ggplot2::aes(y = (.data$lower + .data$upper) / 2,
                   height = .data$upper - .data$lower),
      fill = "grey85", width = bar_w + 0.02, position = dodge, na.rm = TRUE
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_col(ggplot2::aes(colour = .data$outside,
                                   linewidth = .data$outside),
                      width = bar_w, position = dodge, na.rm = TRUE) +
    # A black point at the tip of a bar outside its band, so the bar keeps
    # its group fill. No legend: the caption explains it.
    ggplot2::geom_point(
      ggplot2::aes(y = ifelse(.data$outside, .data$z, NA_real_)),
      shape = 21, fill = "black", colour = "white", size = 1.8, stroke = 0.4,
      position = dodge, na.rm = TRUE, show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(NULL, values = .catdif_styles()$colour[
      seq_len(G)], breaks = group_labels, drop = FALSE) +
    ggplot2::scale_colour_manual(NULL, values = c(`FALSE` = "grey40",
                                                  `TRUE` = "black"),
                                 guide = "none") +
    ggplot2::scale_linewidth_manual(NULL, values = c(`FALSE` = 0.2,
                                                     `TRUE` = 0.4),
                                    guide = "none") +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$score_group)) +
    ggplot2::labs(x = "Category", y = "Standardised residual")
}

# -- Driver ---------------------------------------------------------------------

#' The DIF path of RMitemCatFitPlot()
#'
#' @param prep Output of `.catfit_prepare()` on the rows with a DIF group.
#' @param dif_factor Factor of DIF groups, one per row of `prep$X`, unused
#'   levels dropped.
#' @param n_total Rows supplied to RMitemCatFitPlot().
#' @param n_dif_missing Rows dropped for a missing DIF group. Rows with no
#'   responses are covered by "complete cases" in `prep$quals`.
#' @return A data.frame or a ggplot.
#' @keywords internal
#' @noRd
.catdif_run <- function(prep, dif_factor, type, grouping, class_intervals,
                        score_breaks, iterations, conf_level, item_labels,
                        label_wrap, output, n_total, n_dif_missing) {
  G <- nlevels(dif_factor)
  if (G < 2L) {
    stop("`dif_var` must have at least 2 groups among the respondents used ",
         "(complete cases with a non-extreme total score).", call. = FALSE)
  }
  if (G > 4L) {
    stop("`dif_var` has ", G, " groups. RMitemCatFitPlot() shows at most 4. ",
         "Combine groups, or compare them in pairs.", call. = FALSE)
  }
  grp <- as.integer(dif_factor)
  n_g <- tabulate(grp, G)
  group_labels <- paste0(levels(dif_factor), " (n = ", n_g, ")")
  pr <- .catdif_prepare(prep, grp)

  if (type == "residuals") {
    K <- nrow(prep$P[[1L]]) - 1L
    if (grouping == "manual" && any(score_breaks > K)) {
      stop("`score_breaks` contains values above the maximum possible total ",
           "score (", K, ").", call. = FALSE)
    }
    bin <- .cicc_class_bins(prep$r, grouping, class_intervals, score_breaks)
    sg <- .score_group_labels(prep$r, bin)
    layout <- .catdif_residuals(pr, as.integer(sg))
  } else if (type == "adjacent") {
    layout <- .catdif_adjacent(pr)
  } else {
    layout <- .catdif_calibration(pr)
  }
  b <- .catdif_bands(pr, layout, iterations, conf_level)

  grp_f <- factor(group_labels[b$group], levels = group_labels)
  tab <- switch(
    type,
    calibration = {
      d <- data.frame(item = b$item, event = paste0("X \u2265 ", b$k),
                      group = grp_f, score = b$score, predicted = b$predicted,
                      cep = b$value, n = b$n, perm_mean = b$perm_mean,
                      lower = pmin(pmax(b$lower, 0), 1),
                      upper = pmin(pmax(b$upper, 0), 1),
                      outside = b$outside, curve = b$curve)
      # Scores with the same prediction (certain events, P = 0 or 1, at
      # extreme scores) have fixed responses and identical values: merged.
      key <- paste(d$curve, format(d$predicted, digits = 15))
      n_sum <- tapply(d$n, key, sum)
      d <- d[!duplicated(key), , drop = FALSE]
      d$n <- as.vector(n_sum[paste(d$curve, format(d$predicted,
                                                   digits = 15))])
      d$curve <- NULL
      d$score <- NULL
      d
    },
    adjacent = data.frame(
      item = b$item, pair = paste0(b$pair_j, " vs ", b$pair_j - 1L),
      group = grp_f, score = b$score, n_pair = b$n_pair,
      observed = b$value, expected = b$expected, perm_mean = b$perm_mean,
      lower = pmin(pmax(b$lower, 0), 1), upper = pmin(pmax(b$upper, 0), 1),
      outside = b$outside),
    residuals = data.frame(
      item = b$item,
      score_group = factor(levels(sg)[b$q], levels = levels(sg)),
      category = b$category, group = grp_f, observed = b$observed,
      expected = b$expected, z = b$value, perm_mean = b$perm_mean,
      lower = b$lower, upper = b$upper, outside = b$outside)
  )
  tab$iterations_used <- iterations
  strip_lab <- vapply(item_labels, function(l) {
    paste(strwrap(l, width = label_wrap), collapse = "\n")
  }, character(1L), USE.NAMES = FALSE)
  rownames(tab) <- NULL
  if (output == "dataframe") {
    tab$item <- factor(item_labels[tab$item], levels = item_labels)
    return(tab)
  }
  tab$item <- factor(strip_lab[tab$item], levels = strip_lab)

  quals <- c(prep$quals,
             if (n_dif_missing > 0L) "missing DIF group excluded")
  n_clause <- paste0(.n_caption(prep$n_used, n_total, quals), ".")
  pct <- paste0(format(100 * conf_level), "%")
  band <- paste0(
    "the range of each group's ", if (type == "residuals") "residuals" else
      "curve", " expected without DIF, simultaneous ", pct, " for each ",
    if (type == "residuals") "item and group" else "curve", ", from ",
    iterations, " permutations of group labels within total scores. ",
    if (type == "residuals") "Bars with a black outline and a black point at the tip" else "Ringed points",
    " lie outside their band. Bands are not corrected across ",
    if (type == "residuals") "items and groups" else "curves",
    ", so some cross by chance. ")
  perm_note <- if (iterations < 1000L) {
    paste0("With fewer than 1000 permutations, the band edges depend on ",
           "the seed. ")
  } else {
    ""
  }
  caption <- switch(
    type,
    calibration = paste0(
      "PAV-adjusted conditional event probability for each group against ",
      "the predicted probability of the event given the total score, from ",
      "one fit to all respondents, sized by the number of respondents. ",
      "Bands: ", band, "Dotted line: perfect calibration. "),
    adjacent = paste0(
      "For each group, among respondents in either of two adjacent ",
      "categories at a total score, the share in the upper one, sized by ",
      "the number in the pair. Grey line: expected under the model fitted ",
      "to all respondents. Bands: ", band),
    residuals = paste0(
      "For each group, observed minus expected number of respondents in ",
      "each category, divided by its standard deviation, within groups of ",
      "total scores, with expected counts from one fit to all respondents. ",
      "Grey boxes: ", band,
      .cicc_interval_caption(grouping, class_intervals, score_breaks, bin),
      " ")
  )
  caption <- paste0(caption, perm_note,
                    "Test DIF with RMdifAnchor() or RMdifGamma(). ", n_clause)

  strip <- ggplot2::label_value
  p <- switch(
    type,
    calibration = .catdif_calibration_gg(tab, strip, group_labels),
    adjacent = .catdif_adjacent_gg(tab, strip, group_labels),
    residuals = .catdif_residuals_gg(tab, strip, group_labels)
  )
  p +
    ggplot2::labs(caption = er2_caption(caption)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom",
                   strip.text.y = ggplot2::element_text(angle = 0)) +
    er2_plot_caption() +
    er2_axis_margins()
}
