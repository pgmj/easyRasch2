# Collect results/cells/*.rds (undated files) into one tidy rep-level table
# and a cell-level rejection table. Saves results/summary.rds (+ dated copy).
#   Rscript R/summarise.R
# SMOKE=1 reads results/smoke_cells and writes results/smoke_summary.rds.

smoke <- identical(Sys.getenv("SMOKE"), "1")
indir <- if (smoke) "results/smoke_cells" else "results/cells"
files <- list.files(indir, pattern = "_k(5|10)\\.rds$", full.names = TRUE)
alpha <- 0.05

reps <- dplyr::bind_rows(lapply(files, function(f) {
  x <- readRDS(f)
  cl <- x$cell
  dplyr::bind_rows(lapply(x$results, function(r) {
    base <- data.frame(cond = cl$cond, null = cl$null, n = cl$n, k = cl$k,
                       rep = r$rep, B = x$B,
                       skipped = if (is.null(r$skipped)) "error" else r$skipped)
    if (!is.na(base$skipped)) return(base)
    pc <- r$p_cat
    base$B_used <- r$B_used
    base$p_index1 <- r$p_item[1]
    base$p_index_other_min <- min(r$p_item[-1])
    base$any_other_holm <- any(r$p_item_holm[-1] < alpha)
    base$item1_holm <- r$p_item_holm[1] < alpha
    for (nm in names(pc)) base[[paste0("p_", nm)]] <- pc[[nm]]
    for (j in 1:3) base[[paste0("z", j, j + 1)]] <- r$z1[j]
    base$disorder_any <- any(r$z1 < -stats::qnorm(1 - alpha))
    base$disorder_pt <- any(diff(r$thresholds1) < 0)
    for (j in 1:4) base[[paste0("p_auc", j - 1, j)]] <- r$auc1[j, "p"]
    base
  }))
}))
reps$cond <- factor(reps$cond, levels = unique(c(
  "A0_pcm_ordered", "A1_pcm_disordered", "S_outer_a0", "S_inner_a0",
  "S_inner_a05", "S_inner_a2", "SWAP_1_2", "DK_wide_05", "DK_wide_10",
  "DK_wide_25", "DK_narrow_05", "DK_narrow_10", "DK_narrow_25", "ERS_30")))

# Guard against the 2026-11-16 RNG bug: every replicate should be a distinct
# dataset, so the focus-item threshold estimates should all differ.
dup <- sapply(files, function(f) {
  x <- readRDS(f)
  th <- vapply(x$results, function(r) if (is.null(r$thresholds1)) paste("skipped", r$rep)
                else paste(signif(r$thresholds1, 6), collapse = " "), "")
  length(unique(th)) / length(th)
})
if (any(dup < 0.99)) warning("cells with repeated datasets:\n",
                             paste(basename(files)[dup < 0.99], collapse = "\n"))

out <- list(reps = reps, alpha = alpha, date = Sys.Date())
stem <- if (smoke) "results/smoke_summary" else "results/summary"
saveRDS(out, paste0(stem, ".rds"))
saveRDS(out, sprintf("%s_%s.rds", stem, Sys.Date()))
message("reps: ", nrow(reps), ", skipped: ", sum(!is.na(reps$skipped)))
