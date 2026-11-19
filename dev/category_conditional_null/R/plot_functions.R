# Diagnostic plots for the category-by-score conditional table (prototype,
# 2026-11-18). Built on R/functions.R. Every expected value is the CML
# conditional P(X_i = j | r), so no theta enters any plot.
#
# Note: a plain bar plot of observed vs expected category counts is
# uninformative under CML. The category counts are the sufficient statistics
# for the item's parameters, so at the CML estimate the marginal expected
# counts equal the observed ones exactly. Discrepancies only show once the
# table is split by total score.

library(ggplot2)
if (!exists(".pav")) source("R/calibration_functions.R")

# Condition names wrapped so row strips are not cut off.
.wrap_cond <- labeller(cond = label_wrap_gen(18))

# Facets: columns from `col`, and one row per condition if the data hold a
# `cond` column (used to stack example datasets).
.facet <- function(data, col) {
  if ("cond" %in% names(data)) facet_grid(rows = vars(cond), cols = vars(.data[[col]]),
                                          labeller = .wrap_cond)
  else facet_grid(cols = vars(.data[[col]]))
}

# Long table: one row per score level r and category j of item i.
cond_table <- function(X, delta_list, i) {
  K <- sum(lengths(delta_list))
  r <- rowSums(X)
  keep <- r > 0 & r < K
  X <- X[keep, , drop = FALSE]; r <- r[keep]
  P <- cond_cat_prob(delta_list, i)
  m <- ncol(P) - 1L
  n_r <- tabulate(r + 1L, nbins = K + 1L)
  O <- matrix(table(factor(r, levels = 0:K), factor(X[, i], levels = 0:m)), K + 1L)
  out <- expand.grid(score = 0:K, category = 0:m)
  out$n_r <- n_r[out$score + 1L]
  out$obs <- O[cbind(out$score + 1L, out$category + 1L)]
  out$p_exp <- P[cbind(out$score + 1L, out$category + 1L)]
  out <- out[out$n_r > 0 & out$score - out$category >= 0 &
               out$score - out$category <= K - m, ]
  out$exp <- out$n_r * out$p_exp
  out$dev <- with(out, ifelse(p_exp > 0 & p_exp < 1,
                              (obs - exp) / sqrt(n_r * p_exp * (1 - p_exp)), 0))
  out
}

# Score groups with roughly equal numbers of persons, cut on the total score.
score_groups <- function(tab, groups = 5) {
  nr <- unique(tab[, c("score", "n_r")])
  nr <- nr[order(nr$score), ]
  cum <- cumsum(nr$n_r) / sum(nr$n_r)
  g <- pmin(groups, findInterval(cum - 1e-9, seq(0, 1, length.out = groups + 1)[-1]) + 1)
  lab <- tapply(nr$score, g, function(s) paste0(min(s), "-", max(s)))
  stats::setNames(factor(lab[g], levels = lab), nr$score)
}

# 1. Conditional category probability plot: observed share of each category
#    at each total score against the CML expectation, with the pointwise
#    95% binomial band given n_r (exact under the model, ignoring the
#    uncertainty in the item parameters).
plot_cond_prob <- function(tab, min_n = 1) {
  tab <- tab[tab$n_r >= min_n, ]
  tab$lo <- stats::qbinom(.025, tab$n_r, tab$p_exp) / tab$n_r
  tab$hi <- stats::qbinom(.975, tab$n_r, tab$p_exp) / tab$n_r
  tab$p_obs <- tab$obs / tab$n_r
  tab$cat_lab <- paste("Category", tab$category)
  tab$outside <- factor(ifelse(tab$p_obs < tab$lo | tab$p_obs > tab$hi,
                               "Outside 95% band", "Inside 95% band"),
                        levels = c("Inside 95% band", "Outside 95% band"))
  ggplot(tab, aes(score)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
    geom_line(aes(y = p_exp, linetype = "Expected (CML)"), colour = "#1c5cab",
              linewidth = .6) +
    geom_point(aes(y = p_obs, shape = outside), size = 1.6) +
    scale_shape_manual(NULL, values = c(1, 16), drop = FALSE) +
    scale_linetype_manual(NULL, values = "solid") +
    .facet(tab, "cat_lab") +
    labs(x = "Total score", y = "Proportion in category") +
    theme_bw() + theme(legend.position = "bottom")
}

# 2. Bar plot by score group: observed count per category (bars) against the
#    expected count and its 95% interval (point and bar). Within each group
#    the expected counts sum to the observed group size.
group_table <- function(tab, groups = 5) {
  sg <- score_groups(tab, groups)
  tab$group <- sg[as.character(tab$score)]
  agg <- do.call(rbind, lapply(split(tab, list(tab$group, tab$category), drop = TRUE),
    function(d) data.frame(group = d$group[1], category = d$category[1],
                           obs = sum(d$obs), exp = sum(d$exp),
                           v = sum(d$n_r * d$p_exp * (1 - d$p_exp)))))
  agg$lo <- pmax(0, agg$exp - 1.96 * sqrt(agg$v))
  agg$hi <- agg$exp + 1.96 * sqrt(agg$v)
  agg$z <- (agg$obs - agg$exp) / sqrt(agg$v)
  # Group index, kept as an integer so datasets with different score ranges
  # can be stacked; the score range is in `group`.
  agg$gi <- as.integer(agg$group)
  agg$group <- as.character(agg$group)
  agg
}

.group_lab <- function(gi) {
  g <- max(gi)
  labs <- c("Lowest scores", paste("Score group", 2:(g - 1)), "Highest scores")
  factor(labs[gi], levels = labs)
}

plot_group_bars <- function(agg) {
  agg$group_lab <- .group_lab(agg$gi)
  ggplot(agg, aes(factor(category))) +
    geom_col(aes(y = obs, fill = "Observed"), width = .7, colour = "grey30",
             linewidth = .2) +
    geom_errorbar(aes(ymin = lo, ymax = hi, linetype = "Expected, 95% interval"),
                  width = .3, colour = "#0d366b") +
    geom_point(aes(y = exp, shape = "Expected, 95% interval"), colour = "#0d366b",
               size = 1.8) +
    scale_fill_manual(NULL, values = "grey80") +
    scale_shape_manual(NULL, values = 18) +
    scale_linetype_manual(NULL, values = "solid") +
    .facet(agg, "group_lab") +
    labs(x = "Category", y = "Number of persons") +
    theme_bw() + theme(legend.position = "bottom")
}

# Variant: the same table as standardised residuals (observed - expected) /
#    SD, so every cell is on one scale. Dotted lines at +/- 1.96.
plot_group_resid <- function(agg) {
  agg$group_lab <- .group_lab(agg$gi)
  agg$sign <- factor(ifelse(agg$z >= 0, "More than expected", "Fewer than expected"),
                     levels = c("More than expected", "Fewer than expected"))
  ggplot(agg, aes(factor(category), z, fill = sign)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_hline(yintercept = c(-1.96, 1.96), colour = "grey40", linetype = "dotted") +
    geom_col(width = .7, colour = "grey30", linewidth = .2) +
    scale_fill_manual(NULL, values = c("More than expected" = "#e34948",
                                       "Fewer than expected" = "#0d366b")) +
    .facet(agg, "group_lab") +
    labs(x = "Category", y = "Standardised residual") +
    theme_bw() + theme(legend.position = "bottom")
}

# 3. Heatmap of the Glas deviates per score level and category. Diverging
#    blue-grey-red, clipped at +/- 3.
plot_deviate_map <- function(tab) {
  ggplot(tab, aes(score, factor(category), fill = pmax(-3, pmin(3, dev)))) +
    geom_tile(colour = "white", linewidth = .3) +
    scale_fill_gradient2("Deviate", low = "#2a78d6", mid = "#f0efec",
                         high = "#e34948", midpoint = 0, limits = c(-3, 3)) +
    (if ("cond" %in% names(tab)) facet_grid(rows = vars(cond), labeller = .wrap_cond)) +
    labs(x = "Total score", y = "Category") +
    theme_bw() + theme(panel.grid = element_blank())
}

# 4. Adjacent-category plot: among persons in category j - 1 or j at total
#    score r, the share in j, against its CML expectation. Under the PCM this
#    rises with r. Interchangeable categories give a flat observed curve.
adjacent_table <- function(tab) {
  m <- max(tab$category)
  do.call(rbind, lapply(1:m, function(j) {
    a <- tab[tab$category == j - 1, c("score", "obs", "exp", "n_r")]
    b <- tab[tab$category == j, c("score", "obs", "exp")]
    d <- merge(a, b, by = "score", suffixes = c("_lo", "_hi"))
    d$n_pair <- d$obs_lo + d$obs_hi
    d$p_exp <- d$exp_hi / (d$exp_lo + d$exp_hi)
    d <- d[d$n_pair > 0, ]
    d$p_obs <- d$obs_hi / d$n_pair
    d$lo <- stats::qbinom(.025, d$n_pair, d$p_exp) / d$n_pair
    d$hi <- stats::qbinom(.975, d$n_pair, d$p_exp) / d$n_pair
    d$step <- paste0(j, " vs ", j - 1)
    d
  }))
}

plot_adjacent <- function(adj) {
  adj$step_lab <- paste("Categories", adj$step)
  ggplot(adj, aes(score)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
    geom_line(aes(y = p_exp, linetype = "Expected (CML)"), colour = "#1c5cab",
              linewidth = .6) +
    geom_point(aes(y = p_obs, size = n_pair, shape = "Observed")) +
    scale_size_area("Persons in the pair", max_size = 3, breaks = c(5, 20, 50)) +
    scale_shape_manual(NULL, values = 1) +
    scale_linetype_manual(NULL, values = "solid") +
    .facet(adj, "step_lab") +
    labs(x = "Total score", y = "Share in the upper category") +
    theme_bw() + theme(legend.position = "bottom")
}

# 5. Focused statistics with their null intervals: observed minus expected
#    rest-score mean and SD per category, against the 2.5-97.5% range of the
#    exact conditional bootstrap (patterns sampled given totals, refitted).
boot_stats <- function(X, i, B = 199) {
  delta <- fit_thresholds(X)
  K <- sum(lengths(delta))
  r <- rowSums(X)
  totals <- r[r > 0 & r < K]
  obs <- .boot_vec(item_stats(X, delta, i))
  bt <- t(replicate(B, {
    Xb <- sample_given_totals(delta, totals)
    db <- tryCatch(fit_thresholds(Xb), error = function(e) NULL)
    if (is.null(db) || any(lengths(db) != lengths(delta))) rep(NA_real_, length(obs))
    else .boot_vec(item_stats(Xb, db, i))
  }))
  bt <- bt[stats::complete.cases(bt), , drop = FALSE]
  data.frame(name = names(obs), obs = obs,
             lo = apply(bt, 2, stats::quantile, .025),
             hi = apply(bt, 2, stats::quantile, .975),
             B_used = nrow(bt), row.names = NULL)
}

plot_focused <- function(bs, which = c("mean", "sd")) {
  bs$stat <- sub("_\\d$", "", bs$name)
  bs$category <- as.integer(sub(".*_", "", bs$name))
  bs <- bs[bs$stat %in% which, ]
  bs$stat <- factor(c(mean = "Rest-score mean", sd = "Rest-score SD",
                      tail = "Count in score tails")[bs$stat])
  bs$outside <- factor(ifelse(bs$obs < bs$lo | bs$obs > bs$hi,
                              "Outside null interval", "Inside null interval"),
                       levels = c("Inside null interval", "Outside null interval"))
  ggplot(bs, aes(factor(category))) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "dotted") +
    geom_linerange(aes(ymin = lo, ymax = hi), colour = "grey60", linewidth = 2.5) +
    geom_point(aes(y = obs, shape = outside), size = 2.4) +
    scale_shape_manual(NULL, values = c(1, 16), drop = FALSE) +
    (if ("cond" %in% names(bs)) facet_grid(rows = vars(cond), cols = vars(stat), scales = "free_y",
                                             labeller = .wrap_cond)
     else facet_grid(cols = vars(stat), scales = "free_y")) +
    labs(x = "Category", y = "Observed minus expected") +
    theme_bw() + theme(legend.position = "bottom")
}

# 6. PAV-adjusted calibration of the cumulative events X_i >= k (Dimitriadis,
#    Gneiting & Jordan 2021; Säilynoja et al. 2025, sec. 4.3 and 6), with the
#    forecast taken from the CML conditional P(X_i >= k | r) instead of a
#    theta-based prediction. Consistency band: replicate responses sampled
#    exactly given each person's total score at the fitted item parameters,
#    PAV of each replicate in the same ordering, pointwise 2.5-97.5%.

# .pav() lives in R/calibration_functions.R.

pav_table <- function(X, delta_list, i, B = 199, prob = .95) {
  K <- sum(lengths(delta_list))
  r <- rowSums(X)
  keep <- r > 0 & r < K
  X <- X[keep, , drop = FALSE]; r <- r[keep]
  P <- cond_cat_prob(delta_list, i)
  m <- ncol(P) - 1L
  Xrep <- replicate(B, sample_given_totals(delta_list, r)[, i])  # n x B
  do.call(rbind, lapply(1:m, function(k) {
    pk <- rowSums(P[, (k + 1):(m + 1), drop = FALSE])[r + 1L]
    obs <- .pav(pk, as.numeric(X[, i] >= k))
    rep_cep <- sapply(seq_len(B), function(b) .pav(pk, as.numeric(Xrep[, b] >= k))$cep)
    obs$lo <- apply(rep_cep, 1, stats::quantile, (1 - prob) / 2)
    obs$hi <- apply(rep_cep, 1, stats::quantile, 1 - (1 - prob) / 2)
    obs$event <- paste0("X ≥ ", k)
    obs
  }))
}

plot_pav <- function(pt) {
  ggplot(pt, aes(p)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = "dotted") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
    geom_step(aes(y = cep, linetype = "PAV-adjusted CEP"), colour = "#1c5cab",
              linewidth = .6, direction = "mid") +
    geom_point(aes(y = cep, size = w, shape = "PAV-adjusted CEP"), colour = "#1c5cab") +
    scale_size_area("Persons at this prediction", max_size = 2.5, breaks = c(5, 20, 50)) +
    scale_shape_manual(NULL, values = 16) +
    scale_linetype_manual(NULL, values = "solid") +
    .facet(pt, "event") +
    coord_equal() +
    scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", ".5", "1")) +
    scale_y_continuous(breaks = c(0, .5, 1), labels = c("0", ".5", "1")) +
    labs(x = "Predicted P(event | total score), CML", y = "Conditional event probability") +
    theme_bw() + theme(legend.position = "bottom")
}
