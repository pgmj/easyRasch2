# Collect results/cal_cells/*.rds (undated) into one rep-level table.
# Saves results/calibration_summary.rds (+ dated copy).
#   Rscript R/summarise_calibration.R      (SMOKE=1 for the smoke cells)
smoke <- identical(Sys.getenv("SMOKE"), "1")
indir <- if (smoke) "results/cal_smoke_cells" else "results/cal_cells"
files <- list.files(indir, pattern = "_k(5|10)\\.rds$", full.names = TRUE)
alpha <- 0.05

reps <- dplyr::bind_rows(lapply(files, function(f) {
  x <- readRDS(f); cl <- x$cell
  # A NULL result is a forked worker that was killed (e.g. an interrupted run).
  n_null <- sum(vapply(x$results, is.null, TRUE))
  if (n_null > 0) warning(cl$cell, ": ", n_null, " empty results dropped. Rerun the cell.")
  x$results <- Filter(Negate(is.null), x$results)
  dplyr::bind_rows(lapply(x$results, function(r) {
    base <- data.frame(cond = cl$cond, null = cl$null, n = cl$n, k = cl$k,
                       rep = r$rep, B = x$B, data_id = r$data_id,
                       skipped = if (is.null(r$skipped)) "error" else r$skipped)
    if (!is.na(base$skipped)) return(base)
    base$B_used <- r$B_used
    if (r$B_used == 0) return(base)
    for (e in 1:4) {
      base[[paste0("mcb_", e)]] <- r$mcb[1, e]
      base[[paste0("p_raw_", e)]] <- r$p_raw[1, e]
    }
    base$p_item1 <- r$p_item[1]
    base$item1_holm <- r$p_item_holm[1] < alpha
    base$any_other_holm <- any(r$p_item_holm[-1] < alpha)
    base$any_other_raw <- any(r$p_item[-1] < alpha)  # per-item min-p, no Holm
    # Item p of every item, for null calibration over all items.
    for (i in seq_along(r$p_item)) base[[paste0("p_item_all_", i)]] <- r$p_item[i]
    base
  }))
}))
reps$cond <- factor(reps$cond, levels = unique(c(
  "A0_pcm_ordered", "A1_pcm_disordered", "S_outer_a0", "S_inner_a0",
  "S_inner_a05", "S_inner_a2", "SWAP_1_2", "DK_wide_05", "DK_wide_10",
  "DK_wide_25", "DK_narrow_05", "DK_narrow_10", "DK_narrow_25", "ERS_30")))

# Guard against repeated datasets within a cell (2026-11-16 RNG bug).
dup <- tapply(reps$data_id, paste(reps$cond, reps$n, reps$k),
              function(d) length(unique(d)) / length(d))
if (any(dup < 0.99)) warning("cells with repeated datasets: ",
                             paste(names(dup)[dup < 0.99], collapse = ", "))

out <- list(reps = reps, alpha = alpha, date = Sys.Date())
stem <- if (smoke) "results/cal_smoke_summary" else "results/calibration_summary"
saveRDS(out, paste0(stem, ".rds"))
saveRDS(out, sprintf("%s_%s.rds", stem, Sys.Date()))
message("reps: ", nrow(reps), ", skipped: ", sum(!is.na(reps$skipped)),
        ", min distinct share per cell: ", round(min(dup), 3))
