# Null coverage of the simulated residual ranges in RMitemCatFitPlot(type =
# "residuals", fit = RMitemCatFit(...)). Under the PCM a pointwise 95% range
# should leave about 5% of cells outside. Compared with the share beyond
# +-1.96 in the same cells.
#
# Run from easyRasch2/:  Rscript dev/category_conditional_null/R/residual_band_check.R [datasets] [cores]
# Container: 100 datasets, 4 cores (about 3 min).

args <- commandArgs(trailingOnly = TRUE)
NDATA <- if (length(args) >= 1) as.integer(args[1]) else 100L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 4L
N <- 400L; K <- 8L; M <- 3L; ITER <- 200L

suppressMessages(devtools::load_all(".", quiet = TRUE))
out_dir <- "dev/category_conditional_null/results"
# Same generator and thresholds as residual_null_check.R.
src <- readLines("dev/category_conditional_null/R/residual_null_check.R")
eval(parse(text = src[grep("^true_thr <- ", src):(grep("^one_rep <- ", src) - 1L)]))
thr <- true_thr(K, M)

t0 <- Sys.time()
res <- parallel::mclapply(seq_len(NDATA), function(i) {
  set.seed(20261121 + i)
  X <- as.data.frame(sim_pcm(stats::rnorm(N), thr))
  names(X) <- paste0("I", seq_len(K))
  fit <- tryCatch(suppressMessages(suppressWarnings(
    RMitemCatFit(X, iterations = ITER, seed = i, output = "list"))),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  d <- RMitemCatFitPlot(X, type = "residuals", fit = fit, output = "dataframe")
  d <- d[!is.na(d$z) & !is.na(d$lower), ]
  data.frame(dataset = i, expected = d$expected, z = d$z,
             lower = d$lower, upper = d$upper)
}, mc.cores = CORES)
skipped <- sum(vapply(res, is.null, TRUE))
cells <- do.call(rbind, res)
cells$e_bin <- cut(cells$expected, c(0, 5, 20, Inf),
                   labels = c("E < 5", "5-20", "> 20"), right = FALSE)
summ <- function(d) data.frame(
  cells = nrow(d),
  below_range = mean(d$z < d$lower), above_range = mean(d$z > d$upper),
  below_196 = mean(d$z < -1.96), above_196 = mean(d$z > 1.96))
out <- list(
  settings = list(datasets = NDATA, skipped = skipped, n = N, k = K, m = M,
                  iterations = ITER, score_groups = 5L, conf_level = 0.95),
  overall = summ(cells),
  by_bin = cbind(e_bin = levels(cells$e_bin),
                 do.call(rbind, lapply(split(cells, cells$e_bin), summ))),
  per_dataset = do.call(rbind, lapply(split(cells, cells$dataset), function(d)
    data.frame(dataset = d$dataset[1], outside = mean(d$z < d$lower | d$z > d$upper),
               beyond_196 = mean(abs(d$z) > 1.96))))
)
saveRDS(out, file.path(out_dir, "residual_band_check.rds"))
saveRDS(out, file.path(out_dir, paste0("residual_band_check_", Sys.Date(), ".rds")))
message("done, skipped ", skipped, ", ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
print(out$overall); print(out$by_bin)
