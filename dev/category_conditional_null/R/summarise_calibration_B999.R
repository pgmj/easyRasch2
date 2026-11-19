# Summary of the large-B null check (R/run_calibration_B999.R).
#   Rscript R/summarise_calibration_B999.R   (SMOKE=1 for the smoke cells)
# Writes results/calibration_b999_summary.rds (+ dated copy).

smoke <- identical(Sys.getenv("SMOKE"), "1")
dir <- if (smoke) "results/cal_b999_smoke_cells" else "results/cal_b999_cells"
files <- list.files(dir, pattern = "_k5\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)

items <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f); cl <- x$cell
  nulls <- vapply(x$results, is.null, TRUE)
  if (any(nulls)) warning(cl$cell, ": ", sum(nulls), " NULL results dropped (worker killed?)")
  do.call(rbind, lapply(x$results[!nulls], function(r) {
    if (!is.na(r$skipped) || is.null(r$p_item)) return(NULL)
    data.frame(cond = cl$cond, n = cl$n, k = cl$k, B = x$B, rep = r$rep,
               B_used = r$B_used, item = seq_along(r$p_item),
               p_minp = r$p_item_minp, p_maxT = r$p_item, padj = r$padj)
  }))
}))

skipped <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  data.frame(cond = x$cell$cond, n = x$cell$n, reps = x$reps,
             skipped = sum(vapply(x$results, function(r) is.null(r) || !is.na(r$skipped) ||
                                    is.null(r$p_item), TRUE)))
}))

key <- c("cond", "n")
cells <- do.call(rbind, lapply(split(items, items[key], drop = TRUE), function(d) {
  ds <- split(d, d$rep)
  fw <- vapply(ds, function(z) any(z$padj <= .05), TRUE)
  data.frame(cond = d$cond[1], n = d$n[1], B = d$B[1], datasets = length(ds),
             min_B_used = min(d$B_used),
             maxT_01 = mean(d$p_maxT <= .01), maxT_05 = mean(d$p_maxT <= .05),
             maxT_10 = mean(d$p_maxT <= .10),
             minp_01 = mean(d$p_minp <= .01), minp_05 = mean(d$p_minp <= .05),
             fwer_05 = mean(fw), fwer_se = sqrt(mean(fw) * (1 - mean(fw)) / length(ds)))
}))
rownames(cells) <- NULL

# Strict rejection (p < alpha), the easyRasch2 convention, and the same
# datasets with only the first 199 replicates (the B of the rerun). The
# subsample keeps the studentisation by the full-pool moments, so it is an
# approximation of a B = 199 run.
source("R/calibration_functions.R")
strict <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  R <- Filter(function(r) !is.null(r) && is.na(r$skipped) && !is.null(r$Tm), x$results)
  sub <- t(vapply(R, function(r) {
    T <- r$Tm[1:min(200L, nrow(r$Tm)), , drop = FALSE]
    pi <- colMeans(sweep(T, 2L, T[1L, ], `>=`))
    pj <- pmax(.stepdown_maxT(T[1L, ], T[-1L, , drop = FALSE]), pi)
    c(any(pj < .05), mean(pi < .05), mean(pi < .01))
  }, numeric(3)))
  p <- unlist(lapply(R, `[[`, "p_item"))
  pm <- unlist(lapply(R, `[[`, "p_item_minp"))
  data.frame(cond = x$cell$cond, n = x$cell$n, k = x$cell$k, B = x$B, datasets = length(R),
             item_01 = mean(p < .01), item_05 = mean(p < .05),
             minp_05 = mean(pm < .05),
             fwer = mean(vapply(R, function(r) any(r$padj < .05), TRUE)),
             sub199_fwer = mean(sub[, 1]), sub199_item_05 = mean(sub[, 2]),
             sub199_item_01 = mean(sub[, 3]))
}))
strict$se <- sqrt(.05 * .95 / strict$datasets)
print(strict, digits = 3)
print(skipped)
print(cells, digits = 3)

out <- list(cells = cells, strict = strict, items = items, skipped = skipped, date = Sys.Date())
f <- if (smoke) "results/calibration_b999_smoke_summary.rds" else "results/calibration_b999_summary.rds"
saveRDS(out, f)
saveRDS(out, sub("\\.rds$", paste0("_", Sys.Date(), ".rds"), f))
