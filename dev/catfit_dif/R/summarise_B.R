# Summarise study B cells into one slim table per cell and alpha.
#   SMOKE=1 LANG=en_US.UTF-8 Rscript dev/catfit_dif/R/summarise_B.R
smoke <- identical(Sys.getenv("SMOKE"), "1")
in_dir <- if (smoke) "dev/catfit_dif/results/permband_B_smoke_cells" else
  "dev/catfit_dif/results/permband_B_cells"
files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)

summ <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  ok <- Filter(function(r) !isTRUE(r$failed), x$reps)
  hashes <- vapply(x$reps, function(r) r$hash %||% NA_character_, "")
  n_distinct <- length(unique(stats::na.omit(hashes)))
  if (n_distinct < length(x$reps)) {
    warning(x$cell$cell, ": ", n_distinct, " distinct datasets of ",
            length(x$reps))
  }
  pts <- do.call(rbind, lapply(ok, `[[`, "points"))
  cvs <- do.call(rbind, lapply(seq_along(ok), function(i) {
    cbind(rep = i, ok[[i]]$curves)
  }))
  do.call(rbind, lapply(c(0.05, 0.10), function(a) {
    p <- pts[pts$alpha == a, ]
    by_bin <- tapply(p$outside, p$nbin, sum) / tapply(p$points, p$nbin, sum)
    nc <- nrow(cvs)
    data.frame(
      cell = x$cell$cell, impact = x$cell$impact, groups = x$cell$groups,
      n = x$cell$n, k = x$cell$k, alpha = a,
      reps_ok = length(ok), reps_failed = length(x$reps) - length(ok),
      distinct = n_distinct, B = x$B,
      point_rate = sum(p$outside) / sum(p$points),
      point_rate_1_4 = by_bin["1-4"], point_rate_5_9 = by_bin["5-9"],
      point_rate_10_19 = by_bin["10-19"], point_rate_20_49 = by_bin["20-49"],
      point_rate_50p = by_bin["50+"],
      curve_any_point = mean(cvs$min_p_point <= a),
      curve_simult = mean(cvs$p_sim <= a),
      curve_simult_se = sqrt(mean(cvs$p_sim <= a) *
                               (1 - mean(cvs$p_sim <= a)) / nc),
      curves = nc,
      points_per_curve = mean(cvs$n_points_var),
      row.names = NULL
    )
  }))
}))

out <- if (smoke) "dev/catfit_dif/results/permband_B_smoke_summary.rds" else
  "dev/catfit_dif/results/permband_B_summary.rds"
saveRDS(summ, out)
if (!smoke) {
  saveRDS(summ, sub("\\.rds$", paste0("_", format(Sys.Date(), "%Y%m%d"),
                                      ".rds"), out))
}
print(format(summ, digits = 3), row.names = FALSE)
