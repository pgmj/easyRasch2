# Example diagnostic plots, one simulated dataset per condition (2026-11-18).
# Run from category_conditional_null/:
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript R/plot_examples.R
# Writes results/plots/*.png and results/plot_examples.rds (data behind the
# plots, including the bootstrap intervals, so the qmd never recomputes).

source("R/study_functions.R")
source("R/plot_functions.R")

N <- 500; K_ITEMS <- 10; B <- 199
conds <- c(A1_pcm_disordered = "PCM, disordered (null)",
           SWAP_1_2 = "Categories 1, 2 swapped",
           S_inner_a0 = "Categories 1, 2 interchangeable",
           S_inner_a2 = "Step 2 discrimination 2",
           DK_narrow_25 = "Don't know 25%, narrow middle")

out_rds <- "results/plot_examples.rds"
if (file.exists(out_rds)) {
  ex <- readRDS(out_rds)
} else {
  ex <- lapply(seq_along(conds), function(c) {
    set.seed(20261118 + c)
    X <- gen_data(names(conds)[c], N, K_ITEMS)
    delta <- fit_thresholds(X)
    tab <- cond_table(X, delta, 1)
    list(tab = tab, agg = group_table(tab), adj = adjacent_table(tab),
         boot = boot_stats(X, 1, B), thresholds = delta[[1]])
  })
  names(ex) <- names(conds)
  saveRDS(ex, out_rds)
  saveRDS(ex, sub("\\.rds$", "_2026-11-18.rds", out_rds))
}

stack <- function(part) {
  do.call(rbind, lapply(names(ex), function(nm) {
    d <- ex[[nm]][[part]]; d$cond <- factor(conds[nm], levels = conds); d
  }))
}

dir.create("results/plots", showWarnings = FALSE)
th <- theme(strip.text.y = element_text(size = 7))
save <- function(p, f, h = 9) ggsave(file.path("results/plots", f), p + th,
                                     width = 9, height = h, dpi = 200)

save(plot_cond_prob(stack("tab")), "1_cond_prob.png")
save(plot_group_bars(stack("agg")), "2a_group_bars.png")
save(plot_group_resid(stack("agg")), "2b_group_resid.png")
save(plot_deviate_map(stack("tab")), "3_deviate_map.png", h = 7)
save(plot_adjacent(stack("adj")), "4_adjacent.png")
save(plot_focused(stack("boot"), c("mean", "sd")), "5_focused.png", h = 8)

# Thresholds of the focus item in each example, for the captions.
print(t(sapply(ex, `[[`, "thresholds")), digits = 2)
print(sapply(ex, function(e) e$boot$B_used[1]))
