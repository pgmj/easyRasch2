# Summary of the specificity study (run_specificity_study.R).
# Per dataset: infit p and corrected p for every item from the cell's oracle
# reference set, and RMitemCatFit()'s item p and corrected p. Per condition
# and cell: the share of datasets in which the focus item is flagged
# (corrected p < .05) by each test, by both, by one only, and in which any
# other item is flagged. Infit flags are also split by direction. `events`
# is the number of an item's events with uncorrected p < .05.
#
#   CORES=10 Rscript R/summarise_specificity.R           (SMOKE=1 for the smoke cells)
# Writes results/specificity_summary.rds (+ dated copy).

source("R/specificity_functions.R")
library(parallel)
smoke <- identical(Sys.getenv("SMOKE"), "1")
cores <- as.integer(Sys.getenv("CORES", "4"))
alpha <- 0.05
outdir <- if (smoke) "results/spec_smoke_cells" else "results/spec_cells"
refdir <- if (smoke) "results/spec_smoke_cells/ref" else "results/spec_infit_ref"

refs <- list()
for (f in list.files(refdir, pattern = "^ref_.*\\.rds$", full.names = TRUE)) {
  r <- readRDS(f)
  refs[[sprintf("%s_n%d_k%d", r$base, r$n, r$k)]] <- r$ref
}

files <- list.files(outdir, pattern = "_k(5|10)\\.rds$", full.names = TRUE)
items <- do.call(rbind, mclapply(files, function(f) {
  x <- readRDS(f); cl <- x$cell
  ref <- refs[[sprintf("%s_n%d_k%d", x$base, cl$n, cl$k)]]
  do.call(rbind, lapply(x$results, function(r) {
    if (is.null(r$infit) || is.null(r$catfit_padj)) return(NULL)
    ip <- infit_p(r$infit, ref)
    data.frame(cond = cl$cond, type = cl$type, null = cl$null, n = cl$n, k = cl$k,
               rep = r$rep, item = seq_along(r$infit), infit = r$infit,
               infit_p = ip$p, infit_padj = ip$padj,
               catfit_p = r$catfit_p, catfit_padj = r$catfit_padj,
               events = rowSums(r$p_raw < alpha, na.rm = TRUE),
               ref_size = nrow(ref), B_used = r$B_used)
  }))
}, mc.cores = cores))
stopifnot(nrow(items) > 0)

items$inf_flag <- items$infit_padj < alpha
items$cat_flag <- items$catfit_padj < alpha
foc <- items[items$item == 1, ]
oth <- items[items$item != 1, ]
any_oth <- stats::aggregate(cbind(inf_any = inf_flag, cat_any = cat_flag) ~
                              cond + n + k + rep, oth, any)
foc <- merge(foc, any_oth, by = c("cond", "n", "k", "rep"))

rates <- do.call(rbind, lapply(split(foc, list(foc$cond, foc$n, foc$k), drop = TRUE),
  function(d) data.frame(
    cond = d$cond[1], type = d$type[1], null = d$null[1], n = d$n[1], k = d$k[1],
    datasets = nrow(d),
    infit = mean(d$inf_flag),
    infit_under = mean(d$inf_flag & d$infit > 1),
    infit_over = mean(d$inf_flag & d$infit < 1),
    catfit = mean(d$cat_flag),
    both = mean(d$inf_flag & d$cat_flag),
    catfit_only = mean(d$cat_flag & !d$inf_flag),
    infit_only = mean(d$inf_flag & !d$cat_flag),
    infit_any_other = mean(d$inf_any),
    catfit_any_other = mean(d$cat_any),
    mean_infit = mean(d$infit))))
rates <- rates[order(rates$type != "pcm", rates$cond, rates$n, rates$k), ]

files_all <- list.files(outdir, pattern = "_k(5|10)\\.rds$", full.names = TRUE)
skipped <- do.call(rbind, lapply(files_all, function(f) {
  x <- readRDS(f)
  data.frame(cell = x$cell$cell, reps = x$reps,
             used = sum(vapply(x$results, function(r)
               !is.null(r$infit) && !is.null(r$catfit_padj), TRUE)))
}))

out <- list(rates = rates, skipped = skipped, alpha = alpha, date = Sys.Date(),
            focus_items = foc)
fn <- if (smoke) "results/specificity_smoke_summary.rds" else "results/specificity_summary.rds"
saveRDS(out, fn)
saveRDS(out, sub("\\.rds$", paste0("_", Sys.Date(), ".rds"), fn))
message("datasets used: ", sum(skipped$used), " of ", sum(skipped$reps))
