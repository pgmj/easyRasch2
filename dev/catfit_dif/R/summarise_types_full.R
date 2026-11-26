# Summarise run_types_full.R cells (adjacent and residuals).
#   LANG=en_US.UTF-8 Rscript dev/catfit_dif/R/summarise_types_full.R
`%||%` <- function(a, b) if (is.null(a)) b else a
files <- list.files("dev/catfit_dif/results/permband_types_cells",
                    pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)
summ <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  ok <- Filter(function(r) !isTRUE(r$failed), x$reps)
  hashes <- vapply(x$reps, function(r) r$hash %||% NA_character_, "")
  a <- do.call(rbind, lapply(ok, `[[`, "res"))
  s <- stats::aggregate(cbind(points, points_out, curves, any_point, simult) ~
                          type + alpha, data = a, FUN = sum)
  ppc <- stats::aggregate(points_per_curve ~ type + alpha, data = a,
                          FUN = mean)
  s <- merge(s, ppc)
  data.frame(cell = x$cell$cell, impact = x$cell$impact,
             groups = x$cell$groups, n = x$cell$n, k = x$cell$k, reps_ok = length(ok),
             distinct = length(unique(stats::na.omit(hashes))),
             type = s$type, alpha = s$alpha,
             pointwise = s$points_out / s$points,
             any_point = s$any_point / s$curves,
             simultaneous = s$simult / s$curves,
             curves = s$curves, points_per_curve = s$points_per_curve)
}))
out <- "dev/catfit_dif/results/permband_types_summary.rds"
saveRDS(summ, out)
saveRDS(summ, sub("\\.rds$", paste0("_", format(Sys.Date(), "%Y%m%d"), ".rds"),
                  out))
if (any(summ$distinct < summ$reps_ok)) warning("repeated datasets in a cell")
print(format(summ, digits = 3), row.names = FALSE)
