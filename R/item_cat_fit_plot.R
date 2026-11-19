#' Category Functioning Given the Total Score
#'
#' Four diagnostic plots of how each item's response categories are used,
#' all built on the probability of each category given the respondent's
#' total score, \eqn{P(X_i = j \mid r)}. Under the partial credit model
#' this probability depends on the item parameters alone, so no person
#' locations are estimated and the comparisons are not affected by the
#' error in them. Items are shown in rows.
#'
#' @param data A data.frame or matrix of item responses, items scored from
#'   0. Respondents with missing responses are left out (complete cases),
#'   since the total score must cover the same items for everyone.
#' @param type Character. Which plot to draw:
#'   * `"residuals"` (default): standardised residuals of the number of
#'     respondents in each category, within groups of total scores.
#'   * `"adjacent"`: for each pair of adjacent categories, the share of
#'     respondents in the upper one, against the total score.
#'   * `"restscore"`: the mean and SD of the rest score among respondents
#'     in each category, observed minus expected, against simulated
#'     intervals.
#'   * `"calibration"`: PAV-adjusted calibration of the cumulative events
#'     \eqn{X_i \ge k} against the model's predicted probability.
#' @param fit Optional result of [RMitemCatFit()] with `output = "list"`,
#'   made from the same `data`. Its simulations are used for the rest-score and
#'   calibration plots, so `iterations`, `parallel`, `n_cores` and `seed` are
#'   not needed, and each item's corrected p is shown under its label.
#' @param score_groups Integer. Number of total-score groups for
#'   `type = "residuals"`, formed with roughly equal numbers of
#'   respondents. Default `5L`.
#' @param iterations Integer. Number of simulated datasets for
#'   `type = "restscore"` and `type = "calibration"` when `fit` is not
#'   given. Default `200L`.
#' @param conf_level Numeric in (0, 1). Coverage of the pointwise bands and
#'   intervals. Default `0.95`.
#' @param item_labels Optional character vector of item labels for the
#'   row strips, one per column of `data`. Default: column names.
#' @param label_wrap Integer. Characters per line in the row strips.
#'   Default `15L`.
#' @param parallel,n_cores Logical / integer. Run the simulations in
#'   parallel via \pkg{mirai} when available. Default sequential. When
#'   `n_cores` is `NULL`, `getOption("mc.cores")` is used, and two workers
#'   when that is unset too.
#' @param seed Optional integer for reproducible simulations. See
#'   [easyRasch2-reproducibility].
#' @param output Character. `"ggplot"` (default) or `"dataframe"` for the
#'   data behind the plot.
#'
#' @details
#' Respondents with the lowest or highest possible total score are left
#' out, since their responses are fixed by the score. The item parameters
#' are estimated by conditional maximum likelihood (`psychotools`).
#'
#' A plain comparison of observed and expected category counts would show
#' nothing. The category counts are sufficient statistics for the item
#' parameters, so at the CML estimates the expected counts equal the
#' observed ones. Differences show only once respondents are split by
#' total score, which every plot here does.
#'
#' **Residuals.** For each item, category and score group, the observed
#' number of respondents minus the expected number, divided by its
#' standard deviation under the model. Dotted lines at \eqn{\pm 1.96} mark
#' values that would be unusual for a single cell. The plot has many
#' cells, so some values beyond them are expected even when the model
#' holds. The expected counts treat the item parameters as known.
#'
#' **Adjacent categories.** Among respondents in category \eqn{j - 1} or
#' \eqn{j} at total score \eqn{r}, the share in \eqn{j}. Under the model
#' this share rises with \eqn{r}. A flat observed curve suggests that
#' respondents do not distinguish the two categories. The band is the
#' pointwise binomial range given the number of respondents in the pair.
#'
#' **Rest score.** The mean and SD of the rest score (total score
#' minus the item's score) among respondents in each category, observed
#' minus expected. A category used by respondents with too low or too
#' high rest scores, or with too wide a spread, stands out. Intervals come
#' from datasets simulated from the model given each respondent's total
#' score, with the item parameters re-estimated in each, so the
#' uncertainty in the estimates is included.
#'
#' **Calibration.** For each cumulative event \eqn{X_i \ge k}, the model
#' predicts \eqn{p = P(X_i \ge k \mid r)}. The observed events are
#' regressed on \eqn{p} by pool-adjacent-violators (PAV), which gives the
#' conditional event probability, the best non-decreasing fit to the
#' observed events (Dimitriadis et al., 2021). This is the PAV-adjusted
#' reliability diagram recommended for ordinal outcomes by Säilynoja et
#' al. (2025), with the prediction taken from the total score instead of
#' a person estimate. Agreement puts the curve on the diagonal. The band
#' is the pointwise range of the same curve over datasets simulated given
#' the total scores at the estimated item parameters. A curve that is too
#' flat means the event separates respondents less well than the model
#' predicts, and too steep means better. Because CML reproduces each
#' category's count, the curve cannot lie above or below the diagonal
#' overall, only cross it.
#'
#' All bands and intervals are pointwise. They show where to look, not
#' whether an item fits, and with many points some fall outside even when
#' the model holds.
#'
#' @return
#' * `output = "ggplot"`: a [ggplot2::ggplot] object, one row of panels
#'   per item.
#' * `output = "dataframe"`: the data behind the plot, with columns
#'   `item` and, by `type`:
#'   * `"residuals"`: `score_group`, `category`, `observed`, `expected`,
#'     `z`.
#'   * `"adjacent"`: `pair`, `score`, `n_pair`, `observed`, `expected`,
#'     `lower`, `upper`.
#'   * `"restscore"`: `statistic`, `category`, `observed`, `lower`, `upper`,
#'     `iterations_used`.
#'   * `"calibration"`: `event`, `predicted`, `cep`, `n`, `lower`,
#'     `upper`.
#'
#' @references
#' Dimitriadis, T., Gneiting, T., & Jordan, A. I. (2021). Stable reliability
#' diagrams for probabilistic classifiers. *Proceedings of the National
#' Academy of Sciences, 118*(8), e2016191118.
#' \doi{10.1073/pnas.2016191118}
#'
#' Säilynoja, T., Johnson, A. R., Martin, O. A., & Vehtari, A. (2025).
#' Recommendations for visual predictive checks in Bayesian workflow.
#' *arXiv*. \doi{10.48550/arXiv.2503.01509}
#'
#' Glas, C. A. W. (1988). The derivation of some tests for the Rasch model
#' from the multinomial distribution. *Psychometrika, 53*(4), 525--546.
#' \doi{10.1007/BF02294405}
#'
#' @seealso [RMitemCatFit()] for the test, [RMitemCatProb()] for the model-implied category curves,
#'   [RMitemICCPlot()] for conditional item characteristic curves.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("eRm", quietly = TRUE)) {
#'   data(pcmdat2, package = "eRm")
#'   RMitemCatFitPlot(pcmdat2)
#'   RMitemCatFitPlot(pcmdat2, type = "adjacent")
#'   RMitemCatFitPlot(pcmdat2, type = "calibration", iterations = 100,
#'                      seed = 1)
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
RMitemCatFitPlot <- function(
  data,
  type = c("residuals", "adjacent", "restscore", "calibration"),
  fit = NULL,
  score_groups = 5L,
  iterations = 200L,
  conf_level = 0.95,
  item_labels = NULL,
  label_wrap = 15L,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  output = c("ggplot", "dataframe")
) {
  type <- match.arg(type)
  output <- match.arg(output)
  if (output == "ggplot" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for output = \"ggplot\".",
         call. = FALSE)
  }
  if (!is.numeric(score_groups) || length(score_groups) != 1L ||
        !is.finite(score_groups) || score_groups != round(score_groups) ||
        score_groups < 2) {
    stop("`score_groups` must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(iterations) || length(iterations) != 1L ||
        !is.finite(iterations) || iterations < 20) {
    stop("`iterations` must be a single integer >= 20.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
        conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.null(fit) && !inherits(fit, "RMitemCatFit")) {
    stop("`fit` must be the result of RMitemCatFit(..., output = \"list\").",
         call. = FALSE)
  }
  prep <- .catfit_prepare(data, who = "RMitemCatFitPlot()")
  if (!is.null(fit) && !identical(prep$X, fit$prep$X)) {
    stop("`fit` was made from different data than `data`.", call. = FALSE)
  }
  if (is.null(item_labels) && !is.null(fit)) item_labels <- fit$item_labels
  item_labels <- .catfit_item_labels(item_labels, prep$item_names)
  X <- prep$X
  r <- prep$r
  P <- prep$P
  n_clause <- paste0(.n_caption(prep$n_used, prep$n_total, prep$quals), ".")
  alpha <- 1 - conf_level

  # Simulations: taken from `fit` when given, otherwise run here for the two
  # types that need them.
  S <- NULL
  n_sims <- NA_integer_
  if (type %in% c("restscore", "calibration")) {
    what <- if (type == "restscore") "rest" else "cep"
    if (!is.null(fit)) {
      S <- fit$sims[[what]]
      iterations <- fit$iterations
    } else {
      iterations <- as.integer(iterations)
      if (!is.null(seed)) set.seed(seed)
      S <- .catfit_simulate(prep, iterations, what, parallel, n_cores)[[what]]
    }
    if (is.null(S)) {
      stop("None of the simulated datasets could be refitted.", call. = FALSE)
    }
    n_sims <- ncol(S)
  }

  tab <- switch(
    type,
    residuals = .catplot_residuals(X, r, P, score_groups),
    adjacent = .catplot_adjacent(X, r, P, alpha),
    restscore = .catplot_restscore(X, r, P, S, alpha),
    calibration = .catplot_calibration(X, r, P, S, alpha)
  )
  # Row strips: the item label, and the corrected item p when `fit` is given.
  # Labels are wrapped here, so the p stays on a line of its own.
  strip_lab <- vapply(item_labels, function(l) {
    paste(strwrap(l, width = label_wrap), collapse = "\n")
  }, character(1L), USE.NAMES = FALSE)
  if (!is.null(fit)) {
    strip_lab <- paste0(strip_lab, "\n",
                        vapply(fit$table$padj, .format_p, character(1L)))
  }
  tab$item <- factor(strip_lab[tab$item], levels = strip_lab)
  tab <- tab[, c("item", setdiff(names(tab), "item"))]
  rownames(tab) <- NULL
  if (output == "dataframe") {
    tab$item <- factor(item_labels[as.integer(tab$item)], levels = item_labels)
    return(tab)
  }

  pct <- paste0(format(100 * conf_level), "%")
  fit_clause <- if (is.null(fit)) "" else paste0(
    " Under each item: its p from RMitemCatFit(), corrected across items ",
    "with the ", .correction_label(fit$correction), ".")
  caption <- switch(
    type,
    residuals = paste(
      "Observed minus expected number of respondents in each category,",
      "divided by its standard deviation, within groups of total scores.",
      "Expected counts are from the probability of each category given the",
      "total score. Dotted lines at \u00b11.96. With many cells, some values",
      "beyond them are expected even when the model holds.", n_clause),
    adjacent = paste(
      "Points: among respondents in either of two adjacent categories at a",
      "total score, the share in the upper one, sized by the number of",
      "respondents in the pair. Line: expected under the model. Band:",
      paste0("pointwise ", pct, " binomial range."), n_clause),
    restscore = paste0(
      "Mean and SD of the rest score among respondents in each category, ",
      "observed minus expected. Bars: pointwise ", pct, " range in ",
      n_sims, " datasets simulated given each respondent's total ",
      "score, with item parameters re-estimated in each. ", n_clause,
      .attrition_clause(n_sims, iterations)),
    calibration = paste0(
      "PAV-adjusted conditional event probability against the predicted ",
      "probability of the event given the total score, sized by the number ",
      "of respondents. Band: pointwise ", pct, " range in ", n_sims,
      " datasets simulated given the total scores. Dotted line: perfect ",
      "calibration. ", n_clause)
  )
  caption <- paste0(caption, fit_clause)

  strip <- ggplot2::label_value
  p <- switch(
    type,
    residuals = .catplot_residuals_gg(tab, strip),
    adjacent = .catplot_adjacent_gg(tab, strip),
    restscore = .catplot_restscore_gg(tab, strip),
    calibration = .catplot_calibration_gg(tab, strip)
  )
  p +
    ggplot2::labs(caption = er2_caption(caption)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom",
                   strip.text.y = ggplot2::element_text(angle = 0)) +
    er2_plot_caption() +
    er2_axis_margins()
}

# ---------------------------------------------------------------------
# Internal: conditional probabilities and sampling
# ---------------------------------------------------------------------

#' P(X_i = j | r) for every item
#'
#' Ratios of elementary symmetric functions: the item's category weight
#' times the ESF of the other items at score r - j, over the ESF of all
#' items at r. Uses `.esf_convolve()` (direct polynomial multiplication,
#' all terms positive).
#' @param thr_list Threshold vectors, one per item.
#' @return List of (K + 1) x (m_i + 1) matrices, rows r = 0..K.
#' @keywords internal
#' @noRd
.cond_cat_probs <- function(thr_list) {
  w <- lapply(thr_list, function(t) exp(-c(0, cumsum(t))))
  g_all <- Reduce(.esf_convolve, w, 1)
  K <- length(g_all) - 1L
  lapply(seq_along(w), function(i) {
    g_wo <- Reduce(.esf_convolve, w[-i], 1)
    m <- length(w[[i]]) - 1L
    P <- matrix(0, K + 1L, m + 1L)
    for (j in 0:m) {
      s <- 0:K - j
      ok <- s >= 0 & s < length(g_wo)
      P[ok, j + 1L] <- w[[i]][j + 1L] * g_wo[s[ok] + 1L] / g_all[ok]
    }
    P
  })
}

#' Sample response patterns given each respondent's total score
#'
#' Vectorised over respondents: items are drawn in turn from their
#' distribution given the remaining score, using suffix ESFs (the same
#' scheme as `.sim_conditional()`, which draws one pattern at a time).
#' @param thr_list Threshold vectors, one per item.
#' @param totals Integer total scores, all strictly between 0 and the
#'   maximum.
#' @return Integer matrix, one row per total.
#' @keywords internal
#' @noRd
.sim_given_totals <- function(thr_list, totals) {
  w <- lapply(thr_list, function(t) exp(-c(0, cumsum(t))))
  k <- length(w)
  suffix <- vector("list", k + 1L)
  suffix[[k + 1L]] <- 1
  for (i in k:1) suffix[[i]] <- .esf_convolve(suffix[[i + 1L]], w[[i]])
  n <- length(totals)
  out <- matrix(0L, n, k)
  s <- totals
  for (i in seq_len(k)) {
    m <- length(w[[i]]) - 1L
    g <- suffix[[i + 1L]]
    cp <- matrix(0, n, m + 1L)
    for (j in 0:m) {
      idx <- s - j
      ok <- idx >= 0 & idx < length(g)
      cp[ok, j + 1L] <- w[[i]][j + 1L] * g[idx[ok] + 1L]
    }
    for (j in seq_len(m)) cp[, j + 1L] <- cp[, j + 1L] + cp[, j]
    u <- stats::runif(n) * cp[, m + 1L]
    x <- as.integer(rowSums(u > cp[, seq_len(m), drop = FALSE]))
    out[, i] <- x
    s <- s - x
  }
  out
}

# ---------------------------------------------------------------------
# Internal: tables behind the four plots
# ---------------------------------------------------------------------

#' Observed and expected counts per item, score and category
#' @keywords internal
#' @noRd
.cat_score_table <- function(X, r, P) {
  K <- nrow(P[[1L]]) - 1L
  n_r <- tabulate(r + 1L, nbins = K + 1L)
  do.call(rbind, lapply(seq_along(P), function(i) {
    m <- ncol(P[[i]]) - 1L
    O <- matrix(table(factor(r, levels = 0:K), factor(X[, i], levels = 0:m)),
                K + 1L)
    d <- expand.grid(score = 0:K, category = 0:m)
    d$item <- i
    d$n_r <- n_r[d$score + 1L]
    d$observed <- O[cbind(d$score + 1L, d$category + 1L)]
    d$p <- P[[i]][cbind(d$score + 1L, d$category + 1L)]
    d$expected <- d$n_r * d$p
    d[d$n_r > 0, ]
  }))
}

#' Total-score groups with roughly equal numbers of respondents
#' @return Factor of group labels, one per element of `r`.
#' @keywords internal
#' @noRd
.score_group_factor <- function(r, groups) {
  sc <- sort(unique(r))
  n_s <- tabulate(match(r, sc), length(sc))
  cum <- cumsum(n_s) / sum(n_s)
  g <- pmin(groups, findInterval(cum - 1e-9,
                                 seq(0, 1, length.out = groups + 1L)[-1L]) + 1L)
  g <- match(g, unique(g))
  # Labels built without names: tapply() names would be re-encoded, which
  # warns for the en dash in a non-UTF-8 locale.
  lab <- vapply(seq_len(max(g)), function(j) {
    s <- sc[g == j]
    if (min(s) == max(s)) paste("Score", s[1L])
    else paste0("Scores ", min(s), "\u2013", max(s))
  }, character(1L))
  factor(lab[g[match(r, sc)]], levels = lab)
}

#' @keywords internal
#' @noRd
.catplot_residuals <- function(X, r, P, groups) {
  tab <- .cat_score_table(X, r, P)
  sg <- .score_group_factor(r, groups)
  tab$score_group <- sg[match(tab$score, r)]
  tab$v <- tab$n_r * tab$p * (1 - tab$p)
  agg <- stats::aggregate(cbind(observed, expected, v) ~ item + score_group +
                            category, data = tab, FUN = sum)
  agg$z <- ifelse(agg$v > 0, (agg$observed - agg$expected) / sqrt(agg$v),
                  NA_real_)
  agg$v <- NULL
  agg[order(agg$item, agg$score_group, agg$category), ]
}

#' @keywords internal
#' @noRd
.catplot_adjacent <- function(X, r, P, alpha) {
  tab <- .cat_score_table(X, r, P)
  do.call(rbind, lapply(split(tab, tab$item), function(d) {
    m <- max(d$category)
    do.call(rbind, lapply(seq_len(m), function(j) {
      lo <- d[d$category == j - 1L, c("score", "observed", "expected")]
      hi <- d[d$category == j, c("score", "observed", "expected")]
      a <- merge(lo, hi, by = "score", suffixes = c("_lo", "_hi"))
      a <- a[a$expected_lo + a$expected_hi > 0, ]
      a$n_pair <- a$observed_lo + a$observed_hi
      a <- a[a$n_pair > 0, ]
      p_exp <- a$expected_hi / (a$expected_lo + a$expected_hi)
      data.frame(
        item = d$item[1L],
        pair = paste0(j, " vs ", j - 1L),
        score = a$score,
        n_pair = a$n_pair,
        observed = a$observed_hi / a$n_pair,
        expected = p_exp,
        lower = stats::qbinom(alpha / 2, a$n_pair, p_exp) / a$n_pair,
        upper = stats::qbinom(1 - alpha / 2, a$n_pair, p_exp) / a$n_pair
      )
    }))
  }))
}

#' Rest-score mean and SD per item and category, observed minus expected
#'
#' Expected values weight each rest score by the expected count at that
#' total score. A category no respondent used gives NA.
#' @return data.frame with item, statistic, category, diff.
#' @keywords internal
#' @noRd
.restscore_stats <- function(X, r, P) {
  tab <- .cat_score_table(X, r, P)
  tab$rest <- tab$score - tab$category
  tab <- tab[tab$rest >= 0, ]
  wm <- function(w, x) sum(w * x) / sum(w)
  ws <- function(w, x) sqrt(sum(w * (x - wm(w, x))^2) / sum(w))
  out <- do.call(rbind, lapply(split(tab, list(tab$item, tab$category),
                                     drop = TRUE), function(d) {
    o_ok <- sum(d$observed) > 0
    data.frame(
      item = d$item[1L], category = d$category[1L],
      statistic = c("mean", "sd"),
      diff = c(if (o_ok) wm(d$observed, d$rest) - wm(d$expected, d$rest)
               else NA_real_,
               if (o_ok) ws(d$observed, d$rest) - ws(d$expected, d$rest)
               else NA_real_)
    )
  }))
  out[order(out$item, out$statistic, out$category), ]
}

#' @keywords internal
#' @noRd
.catplot_restscore <- function(X, r, P, S, alpha) {
  obs <- .restscore_stats(X, r, P)
  data.frame(
    item = obs$item,
    statistic = obs$statistic,
    category = obs$category,
    observed = obs$diff,
    lower = apply(S, 1L, stats::quantile, alpha / 2, na.rm = TRUE),
    upper = apply(S, 1L, stats::quantile, 1 - alpha / 2, na.rm = TRUE),
    iterations_used = ncol(S)
  )
}

#' Predicted P(X_i >= k | r) for each respondent, item and event
#' @return List over items of n x m_i matrices.
#' @keywords internal
#' @noRd
.event_forecasts <- function(P, r) {
  lapply(P, function(Pi) {
    m <- ncol(Pi) - 1L
    E <- vapply(seq_len(m), function(k) {
      rowSums(Pi[, (k + 1L):(m + 1L), drop = FALSE])
    }, numeric(nrow(Pi)))
    matrix(E, ncol = m)[r + 1L, , drop = FALSE]
  })
}

#' Weighted PAV of y on p, ties in p pooled first
#'
#' Returns the conditional event probability (CEP) per distinct value of
#' p. Agrees with `stats::isoreg()` on the tie-pooled data.
#' @return data.frame with p, cep, n.
#' @keywords internal
#' @noRd
.pav <- function(p, y) {
  u <- sort(unique(p))
  idx <- match(p, u)
  n <- tabulate(idx, length(u))
  v <- as.numeric(rowsum(y, idx, reorder = TRUE)[, 1L] / n)
  # Stack of blocks: value, weight, number of distinct p values.
  bv <- bw <- numeric(0)
  bl <- integer(0)
  for (j in seq_along(u)) {
    bv <- c(bv, v[j]); bw <- c(bw, n[j]); bl <- c(bl, 1L)
    while (length(bv) > 1L && bv[length(bv) - 1L] > bv[length(bv)]) {
      t <- length(bv)
      bv[t - 1L] <- (bv[t - 1L] * bw[t - 1L] + bv[t] * bw[t]) /
        (bw[t - 1L] + bw[t])
      bw[t - 1L] <- bw[t - 1L] + bw[t]
      bl[t - 1L] <- bl[t - 1L] + bl[t]
      bv <- bv[-t]; bw <- bw[-t]; bl <- bl[-t]
    }
  }
  data.frame(p = u, cep = rep(bv, bl), n = n)
}

#' CEPs of one simulated dataset at the observed forecasts
#' @return Numeric vector, CEPs stacked over items and events in the order
#'   of `.catplot_calibration()`.
#' @keywords internal
#' @noRd
.pav_replicate <- function(Xb, ev) {
  unlist(lapply(seq_along(ev), function(i) {
    lapply(seq_len(ncol(ev[[i]])), function(k) {
      .pav(ev[[i]][, k], as.numeric(Xb[, i] >= k))$cep
    })
  }))
}

#' @keywords internal
#' @noRd
.catplot_calibration <- function(X, r, P, S, alpha) {
  ev <- .event_forecasts(P, r)
  obs <- do.call(rbind, lapply(seq_along(ev), function(i) {
    do.call(rbind, lapply(seq_len(ncol(ev[[i]])), function(k) {
      d <- .pav(ev[[i]][, k], as.numeric(X[, i] >= k))
      data.frame(item = i, event = paste0("X \u2265 ", k),
                 predicted = d$p, cep = d$cep, n = d$n)
    }))
  }))
  obs$lower <- apply(S, 1L, stats::quantile, alpha / 2)
  obs$upper <- apply(S, 1L, stats::quantile, 1 - alpha / 2)
  obs
}

# ---------------------------------------------------------------------
# Internal: plots. Facet labeller and theme are added by the caller.
# ---------------------------------------------------------------------

#' @keywords internal
#' @noRd
.catplot_residuals_gg <- function(tab, strip) {
  tab$sign <- factor(ifelse(tab$z >= 0, "More than expected",
                            "Fewer than expected"),
                     levels = c("More than expected", "Fewer than expected"))
  ggplot2::ggplot(tab, ggplot2::aes(factor(.data$category), .data$z,
                                    fill = .data$sign)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_hline(yintercept = c(-1.96, 1.96), colour = "grey40",
                        linetype = "dotted") +
    ggplot2::geom_col(width = 0.7, colour = "grey30", linewidth = 0.2,
                      na.rm = TRUE) +
    ggplot2::scale_fill_manual(NULL, values = c(
      "More than expected" = "#e34948", "Fewer than expected" = "#0d366b"
    ), drop = FALSE) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$score_group)) +
    ggplot2::labs(x = "Category", y = "Standardised residual")
}

#' @keywords internal
#' @noRd
.catplot_adjacent_gg <- function(tab, strip) {
  tab$pair <- paste("Categories", tab$pair)
  ggplot2::ggplot(tab, ggplot2::aes(.data$score)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "grey85") +
    ggplot2::geom_line(ggplot2::aes(y = .data$expected,
                                    linetype = "Expected"),
                       colour = "#1c5cab", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(y = .data$observed, size = .data$n_pair,
                                     shape = "Observed")) +
    ggplot2::scale_size_area("Respondents in the pair", max_size = 3) +
    ggplot2::scale_shape_manual(NULL, values = 1) +
    ggplot2::scale_linetype_manual(NULL, values = "solid") +
    ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$pair)) +
    ggplot2::labs(x = "Total score", y = "Share in the upper category")
}

#' @keywords internal
#' @noRd
.catplot_restscore_gg <- function(tab, strip) {
  tab$statistic <- factor(c(mean = "Rest-score mean",
                            sd = "Rest-score SD")[tab$statistic],
                          levels = c("Rest-score mean", "Rest-score SD"))
  tab$outside <- factor(
    ifelse(tab$observed < tab$lower | tab$observed > tab$upper,
           "Outside the range", "Inside the range"),
    levels = c("Inside the range", "Outside the range")
  )
  ggplot2::ggplot(tab, ggplot2::aes(factor(.data$category))) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40",
                        linetype = "dotted") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lower,
                                         ymax = .data$upper),
                            colour = "grey60", linewidth = 2.5,
                            na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(y = .data$observed,
                                     shape = .data$outside),
                        size = 2.2, na.rm = TRUE) +
    ggplot2::scale_shape_manual(NULL, values = c(1, 16), drop = FALSE) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$statistic),
                        scales = "free_y") +
    ggplot2::labs(x = "Category", y = "Observed minus expected")
}

#' @keywords internal
#' @noRd
.catplot_calibration_gg <- function(tab, strip) {
  ggplot2::ggplot(tab, ggplot2::aes(.data$predicted)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey40",
                         linetype = "dotted") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "grey85") +
    ggplot2::geom_step(ggplot2::aes(y = .data$cep, linetype = "PAV-adjusted"),
                       colour = "#1c5cab", linewidth = 0.6,
                       direction = "mid") +
    ggplot2::geom_point(ggplot2::aes(y = .data$cep, size = .data$n,
                                     shape = "PAV-adjusted"),
                        colour = "#1c5cab") +
    ggplot2::scale_size_area("Respondents", max_size = 2.5) +
    ggplot2::scale_shape_manual(NULL, values = 16) +
    ggplot2::scale_linetype_manual(NULL, values = "solid") +
    ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1),
                                labels = c("0", ".5", "1")) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1),
                                labels = c("0", ".5", "1")) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::facet_grid(labeller = strip, rows = ggplot2::vars(.data$item),
                        cols = ggplot2::vars(.data$event)) +
    ggplot2::labs(x = "Predicted probability given the total score",
                  y = "Conditional event probability")
}
