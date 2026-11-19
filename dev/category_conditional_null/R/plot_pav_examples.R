# PAV calibration demo on the same five datasets as R/plot_examples.R
# (same seeds, data regenerated). Writes results/plots/6_pav.png and
# results/plot_pav_examples.rds.
source("R/study_functions.R")
source("R/plot_functions.R")
N <- 500; K_ITEMS <- 10
conds <- c(A1_pcm_disordered = "PCM, disordered (null)",
           SWAP_1_2 = "Categories 1, 2 swapped",
           S_inner_a0 = "Categories 1, 2 interchangeable",
           S_inner_a2 = "Step 2 discrimination 2",
           DK_narrow_25 = "Don't know 25%, narrow middle")
out_rds <- "results/plot_pav_examples.rds"
if (file.exists(out_rds)) pv <- readRDS(out_rds) else {
  pv <- do.call(rbind, lapply(seq_along(conds), function(c) {
    set.seed(20261118 + c)
    X <- gen_data(names(conds)[c], N, K_ITEMS)
    delta <- fit_thresholds(X)
    d <- pav_table(X, delta, 1)
    d$cond <- factor(conds[c], levels = conds); d
  }))
  saveRDS(pv, out_rds); saveRDS(pv, sub("\\.rds$", "_2026-11-18.rds", out_rds))
}
ggsave("results/plots/6_pav.png", plot_pav(pv) + theme(strip.text.y = element_text(size = 6)),
       width = 8, height = 10, dpi = 200)
# share of points outside the band, per condition and event (descriptive)
print(aggregate(cbind(out = cep < lo | cep > hi) ~ cond + event, pv, mean), digits = 2)
