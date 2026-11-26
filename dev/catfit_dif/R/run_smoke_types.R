# Band smoke for all three plot types (calibration as regression check
# against the B smoke, adjacent, residuals). Same cells and seeds as the B
# smoke. Run from the package root (container, about 10 min on 4 cores):
#   NCORES=4 LANG=en_US.UTF-8 Rscript dev/catfit_dif/R/run_smoke_types.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/catfit_dif/R/perm_functions.R")

cores <- as.integer(Sys.getenv("NCORES", "4"))
B <- 999L
R <- 200L
alphas <- c(0.05, 0.10)
groups <- list(two_80_20 = c(0.8, 0.2), three_equal = rep(1 / 3, 3))
# Cell ids as in run_B.R, so the datasets are the B smoke datasets.
cells <- data.frame(id = c(6L, 44L), impact = c(1, 0.5),
                    groups = c("two_80_20", "three_equal"),
                    n = c(250L, 500L), k = c(5L, 10L),
                    stringsAsFactors = FALSE)
cells$cell <- sprintf("c%02d_imp%.1f_%s_n%d_k%d", cells$id, cells$impact,
                      cells$groups, cells$n, cells$k)

one_rep <- function(row, rep) {
  set.seed(20261125L + 10000L * row$id + rep)
  d <- gen_null(row$n, row$k, row$impact, groups[[row$groups]])
  hash <- digest::digest(d$X)
  pr <- prep_dif(d$X, d$grp)
  if (is.null(pr)) return(list(failed = TRUE, hash = hash))
  setups <- list(calibration = calib_setup(pr), adjacent = adjacent_setup(pr),
                 residuals = residual_setup(pr, 5L))
  res <- lapply(names(setups), function(nm) {
    st <- setups[[nm]]
    b <- perm_bands_generic(pr, st$lay, st$fun, B = B)
    pv <- b$points$p_point[b$points$varies]
    data.frame(type = nm, alpha = alphas,
               points = length(pv),
               points_out = vapply(alphas, function(a) sum(pv <= a), 0),
               curves = nrow(b$curves),
               any_point = vapply(alphas, function(a) {
                 sum(b$curves$min_p_point <= a)
               }, 0),
               simult = vapply(alphas, function(a) {
                 sum(b$curves$p_sim <= a)
               }, 0),
               points_per_curve = mean(b$curves$n_points_var))
  })
  list(failed = FALSE, hash = hash, res = do.call(rbind, res))
}

out <- list()
for (j in seq_len(nrow(cells))) {
  row <- cells[j, ]
  tc <- Sys.time()
  reps <- parallel::mclapply(seq_len(R), function(rep) {
    tryCatch(one_rep(row, rep), error = function(e) {
      list(failed = TRUE, error = conditionMessage(e))
    })
  }, mc.cores = cores)
  out[[row$cell]] <- list(cell = row, reps = reps)
  cat(sprintf("%s done in %.1f min (failed: %d)\n", row$cell,
              as.numeric(difftime(Sys.time(), tc, units = "mins")),
              sum(vapply(reps, `[[`, logical(1), "failed"))))
}
saveRDS(out, "dev/catfit_dif/results/permband_types_smoke.rds")

summ <- do.call(rbind, lapply(out, function(x) {
  ok <- Filter(function(r) !isTRUE(r$failed), x$reps)
  hashes <- vapply(x$reps, function(r) r$hash %||% NA_character_, "")
  a <- do.call(rbind, lapply(ok, `[[`, "res"))
  s <- stats::aggregate(cbind(points, points_out, curves, any_point, simult) ~
                          type + alpha, data = a, FUN = sum)
  ppc <- stats::aggregate(points_per_curve ~ type + alpha, data = a,
                          FUN = mean)
  s <- merge(s, ppc)
  data.frame(cell = x$cell$cell, reps_ok = length(ok),
             distinct = length(unique(stats::na.omit(hashes))),
             type = s$type, alpha = s$alpha,
             pointwise = s$points_out / s$points,
             any_point = s$any_point / s$curves,
             simultaneous = s$simult / s$curves,
             curves = s$curves, points_per_curve = s$points_per_curve)
}))
saveRDS(summ, "dev/catfit_dif/results/permband_types_smoke_summary.rds")
print(format(summ, digits = 3), row.names = FALSE)
