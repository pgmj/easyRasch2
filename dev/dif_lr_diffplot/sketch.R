suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(ggplot2)
out <- commandArgs(TRUE)[1]

# 10 items, 4 categories, n = 900, three groups. Real DIF: I3 is 0.8 logits
# harder in group "b", I7 0.7 logits easier in group "c".
set.seed(42)
k <- 10
base <- lapply(seq(-1.5, 1.5, length.out = k), function(l) l + c(-1, 0, 1))
grp <- rep(c("a", "b", "c"), each = 300)
th <- rnorm(900)
sim_group <- function(shift) {
  d <- base
  for (nm in names(shift)) d[[as.integer(nm)]] <- d[[as.integer(nm)]] + shift[[nm]]
  d
}
dat <- rbind(
  sim_partial_score(base, th[grp == "a"]),
  sim_partial_score(sim_group(list(`3` = 0.8)), th[grp == "b"]),
  sim_partial_score(sim_group(list(`7` = -0.7)), th[grp == "c"])
)
dat <- as.data.frame(dat); names(dat) <- paste0("I", 1:k)

# Differences against a reference group, SE = sqrt(SE_g^2 + SE_ref^2)
# (independent samples). One row per item (or threshold) and comparison.
dif_diff <- function(res, ref, conf = 0.95, cutoff = 0.5) {
  grps <- setdiff(sub("^SE_", "", grep("^SE_", names(res), value = TRUE)), "All")
  others <- setdiff(grps, ref)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  thr <- "Threshold" %in% names(res)
  lab <- if (thr) res$Threshold else res$Item
  do.call(rbind, lapply(others, function(g) {
    d <- res[[g]] - res[[ref]]
    se <- sqrt(res[[paste0("SE_", g)]]^2 + res[[paste0("SE_", ref)]]^2)
    data.frame(Item = res$Item, Label = lab, Comparison = paste0(g, " - ", ref), Diff = d,
               SE = se, lower = d - z * se, upper = d + z * se,
               Beyond = abs(d) > cutoff)
  }))
}

plot_diff <- function(df, cutoff = 0.5, conf = 0.95, title) {
  by_item <- any(df$Label != df$Item)
  df$Label <- factor(df$Label, levels = rev(unique(df$Label)))
  df$Item <- factor(df$Item, levels = unique(df$Item))
  df$Beyond <- factor(ifelse(df$Beyond, paste0("|difference| > ", cutoff), "within cutoff"),
                      levels = c("within cutoff", paste0("|difference| > ", cutoff)))
  ggplot(df, aes(x = Diff, y = Label)) +
    geom_vline(xintercept = 0, colour = "grey40") +
    geom_vline(xintercept = c(-cutoff, cutoff), linetype = "dashed", colour = "grey55") +
    geom_errorbar(aes(xmin = lower, xmax = upper, colour = Beyond, linetype = Beyond), width = 0.25, orientation = "y") +
    geom_point(aes(colour = Beyond, shape = Beyond), size = 2.6) +
    scale_colour_manual(values = c("black", "#D55E00"), drop = FALSE) +
    scale_shape_manual(values = c(16, 17), drop = FALSE) +
    scale_linetype_manual(values = c("solid", "longdash"), drop = FALSE) +
    # Threshold level: one row of panels per item, so each item's thresholds
    # sit together in their own box, with the item name as a left strip.
    (if (by_item) {
      facet_grid(Item ~ Comparison, scales = "free_y", space = "free_y",
                 switch = "y")
    } else {
      facet_wrap(~Comparison)
    }) +
    labs(x = "Difference in location (logits)", y = NULL, colour = NULL, shape = NULL, linetype = NULL,
         title = title,
         caption = paste0("Bars: ", 100 * conf, "% CI of the difference, SE = sqrt(SE_g^2 + SE_ref^2).\n",
                          "Dashed lines: the ", cutoff, "-logit cutoff.")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", strip.text = element_text(size = 12),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold"),
          panel.spacing.y = unit(0.35, "lines"))
}

ab <- dat[grp != "c", ]; g_ab <- grp[grp != "c"]
r_item2 <- RMdifLR(ab, dif_var = g_ab, output = "dataframe")
r_thr2 <- RMdifLR(ab, dif_var = g_ab, level = "threshold", output = "dataframe")
r_item3 <- RMdifLR(dat, dif_var = grp, output = "dataframe")

ggsave(file.path(out, "diff_item_2groups.png"),
       plot_diff(dif_diff(r_item2, "a"), title = "Item level, two groups (real DIF in I3)"),
       width = 6.5, height = 5, dpi = 150)
ggsave(file.path(out, "diff_threshold_2groups.png"),
       plot_diff(dif_diff(r_thr2, "a"), title = "Threshold level, two groups (real DIF in I3)"),
       width = 6.5, height = 8, dpi = 150)
ggsave(file.path(out, "diff_item_3groups.png"),
       plot_diff(dif_diff(r_item3, "a"), title = "Item level, three groups, reference a\n(real DIF in I3 for b, in I7 for c)"),
       width = 9, height = 5, dpi = 150)
print(dif_diff(r_item3, "a"), digits = 2)
