## Coverage of fractional weighted bootstrap (FWB) intervals for conditional
## infit MSQ
##
## Question: does an FWB percentile interval for an item's observed infit MSQ
## cover the item's population MSQ, for misfitting items, for fitting items
## distorted by the cascade, and under the null? And do simultaneous intervals
## (one max-type constant over items, the infit analogue of Q3max*) cover all
## items at once? This is the interval an effect-size floor on infit needs.
##
## Datasets are regenerated from conditional_rmsd_partA.qmd with its seeds, and
## each is checked against Part A's cached infit (must match). The population
## MSQ of each item is its Part A mean at n = 20,000 (300 replications).
##
## Run from easyRasch2/dev:  source("fwb_msq_coverage.R")
## Results are cached per cell in fwb_msq_coverage/ (undated and dated copies),
## so an interrupted run resumes at the next unfinished cell. Report:
## fwb_msq_coverage.qmd, which only reads the cache.

suppressMessages({
  library(parallel)
})

NCORES  <- 10
B       <- 200          # FWB draws per dataset
R_COV   <- 200          # datasets per cell
COV_NS  <- c(500, 2000, 5000)
COV_CONDS <- c("null", "xdim50_k1", "xdim15_k3", "xdim75_k3")
OUT_DIR <- "fwb_msq_coverage"
dir.create(OUT_DIR, showWarnings = FALSE)

## ---- Part A definitions (generator, conditions, seeds, LPE calibration) ----
## Evaluates conditional_rmsd_partA.qmd up to its `run` chunk. Everything it
## needs from its own calibration chunks is read from conditional_rmsd/.
local({
  tmp <- file.path(tempdir(), "partA_defs.R")
  knitr::purl("conditional_rmsd_partA.qmd", output = tmp, quiet = TRUE)
  code <- readLines(tmp)
  stop_at <- grep("^## ----run", code)[1]
  invisible(capture.output(
    for (e in parse(text = code[seq_len(stop_at - 1)])) eval(e, globalenv())
  ))
})
stopifnot(exists("gen_data"), exists("conds"), exists("misfit_idx"))

## ---- population MSQ per item: Part A mean at n = 20,000 --------------------
pop_msq <- do.call(rbind, lapply(COV_CONDS, function(cc) {
  d <- readRDS(file.path("conditional_rmsd", sprintf("%s_n20000.rds", cc)))
  agg <- aggregate(infit ~ item, data = d, FUN = mean)
  data.frame(cond = cc, item = agg$item, pop = agg$infit)
}))

## ---- FWB draws of conditional infit ---------------------------------------
## Weighted CML refit, then infit from the conditional moments with the same
## weights. With unit weights this equals iarm::out_infit() (checked
## 2026-10-12 to three decimals).
infit_weighted <- function(Xm, w) {
  m <- psychotools::pcmodel(as.data.frame(Xm), weights = w)
  thr <- lapply(psychotools::threshpar(m), as.numeric)
  mom <- .cond_moments(thr)
  s <- rowSums(Xm)
  keep <- s > 0 & s < mom$max_score
  E <- mom$E[s[keep] + 1, , drop = FALSE]
  V <- mom$V[s[keep] + 1, , drop = FALSE]
  wk <- w[keep]
  colSums(wk * (Xm[keep, , drop = FALSE] - E)^2) / colSums(wk * V)
}

fwb_draws <- function(Xm, B) {
  n <- nrow(Xm)
  out <- matrix(NA_real_, B, ncol(Xm))
  for (b in seq_len(B)) {
    w <- stats::rexp(n)
    w <- w / mean(w)
    out[b, ] <- tryCatch(infit_weighted(Xm, w),
                         error = function(e) rep(NA_real_, ncol(Xm)))
  }
  out
}

one_dataset <- function(ci, n, r, cached_infit) {
  set.seed(1e7 * ci + 100 * (n %/% 100) + r)   # Part A's seed for this dataset
  cc <- conds[ci, ]
  X <- gen_data(cc$type, cc$k, cc$par, n)
  Xm <- as.matrix(X)
  obs <- infit_weighted(Xm, rep(1, n))
  match_err <- max(abs(obs - cached_infit))
  set.seed(5e8 + 1e6 * ci + 100 * (n %/% 100) + r)
  draws <- fwb_draws(Xm, B)
  ok <- stats::complete.cases(draws)
  draws <- draws[ok, , drop = FALSE]
  lo <- apply(draws, 2, stats::quantile, 0.025)
  hi <- apply(draws, 2, stats::quantile, 0.975)
  sdv <- apply(draws, 2, stats::sd)
  # simultaneous: one constant over items on the standardised scale
  tmax <- apply(abs(sweep(draws, 2, obs)) / rep(sdv, each = nrow(draws)), 1, max)
  cmax <- stats::quantile(tmax, 0.95)
  data.frame(rep = r, item = seq_len(ncol(Xm)), obs = obs, fwb_mean = colMeans(draws),
             fwb_sd = sdv, lo = lo, hi = hi,
             slo = obs - cmax * sdv, shi = obs + cmax * sdv, cmax = cmax,
             draws_ok = sum(ok), match_err = match_err)
}

run_cell <- function(cc_name, n) {
  f <- file.path(OUT_DIR, sprintf("%s_n%d.rds", cc_name, n))
  if (file.exists(f)) return(invisible(readRDS(f)))
  ci <- match(cc_name, conds$cond)
  partA <- readRDS(file.path("conditional_rmsd", sprintf("%s_n%d.rds", cc_name, n)))
  t0 <- Sys.time()
  res <- mclapply(seq_len(R_COV), function(r) {
    cached <- partA$infit[partA$rep == r][order(partA$item[partA$rep == r])]
    tryCatch(one_dataset(ci, n, r, cached), error = function(e) conditionMessage(e))
  }, mc.cores = NCORES)
  ok <- vapply(res, is.data.frame, logical(1))
  d <- do.call(rbind, res[ok])
  d$cond <- cc_name; d$n <- n
  d$misfit <- d$item %in% misfit_idx[seq_len(conds$k[ci])]
  d <- merge(d, pop_msq, by = c("cond", "item"))
  attr(d, "failed") <- unlist(res[!ok])
  attr(d, "minutes") <- as.numeric(Sys.time() - t0, units = "mins")
  attr(d, "B") <- B
  saveRDS(d, f)
  saveRDS(d, file.path(OUT_DIR, sprintf("%s_n%d_%s.rds", cc_name, n, Sys.Date())))
  cat(sprintf("%s n=%d: %d datasets, %.1f min, max match error %.2g\n",
              cc_name, n, sum(ok), attr(d, "minutes"), max(d$match_err)))
  invisible(d)
}

## Smoke test: Sys.setenv(FWB_SMOKE = "1") before sourcing defines everything
## but runs nothing.
if (!identical(Sys.getenv("FWB_SMOKE"), "1")) {
  for (n in COV_NS) for (cc in COV_CONDS) run_cell(cc, n)
  cat("done\n")
}
