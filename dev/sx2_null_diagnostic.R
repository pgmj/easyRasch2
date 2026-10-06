## S-X2 null diagnostic
##
## Design: dev/sx2-null-diagnostic-design.md. Null-only check of why
## mirt::itemfit(fit_stats = "S_X2") rejects too often under a true Rasch
## model in the July 2026 S-X2 study (removed, see git history). Variants:
##   V0  Rasch, latent variance estimated, default itemfit (mincell = 1)
##   V1  V0 with p recomputed on df + 1 (H1: df one too small)
##   V2  same fit, mincell = 5 (H2: sparse cells)
##   V3  Rasch, latent variance fixed at its true value 2.25 (H1 vs H3)
##
## Run: source() from easyRasch2/dev. About 15 minutes on 10 cores. Results
## are cached in sx2_null_diagnostic/results.rds; delete it to re-run.

suppressMessages({
  library(mirt)
  library(parallel)
})

R_REPS   <- 2000
NS       <- c(150, 500, 2000)
THETA_SD <- 1.5
NCORES   <- 10
OUT_DIR  <- "sx2_null_diagnostic"
RES_FILE <- file.path(OUT_DIR, "results.rds")

## ---- item sets ---------------------------------------------------------
set.seed(20250727)                     # as make_master() in the July 2026 S-X2 study (removed, see git history)
orig_locs <- runif(20, -2, 2)

item_sets <- list(
  dich_even = as.list(seq(-2, 2, length.out = 20)),
  dich_orig = as.list(orig_locs),
  poly      = lapply(seq(-1.5, 1.5, length.out = 9),
                     function(l) l + c(-1.2, 0, 1.2))
)

## ---- PCM generator, as restscore_asymptotic_null.qmd -------------------
sim_pcm <- function(theta, thr) {
  n <- length(theta)
  sapply(thr, function(d) {
    m <- length(d)
    eta <- outer(theta, 0:m) -
      matrix(c(0, cumsum(d)), n, m + 1, byrow = TRUE)
    p <- exp(eta - apply(eta, 1, max))
    p <- p / rowSums(p)
    cp <- p %*% upper.tri(diag(m + 1), diag = TRUE)
    rowSums(runif(n) > cp)
  })
}

## ---- one replication ---------------------------------------------------
sx2_cols <- function(f) cbind(S = f$S_X2, df = f$df.S_X2, p = f$p.S_X2)

one_rep <- function(rep, set, n) {
  set.seed(1e7 * match(set, names(item_sets)) + 10 * n + rep)
  thr <- item_sets[[set]]
  d <- as.data.frame(sim_pcm(rnorm(n, 0, THETA_SD), thr))
  names(d) <- paste0("I", seq_along(thr))

  # every category observed in every item, as in the restscore study
  full <- all(mapply(function(x, m) all(tabulate(x + 1, m + 1) > 0),
                     d, lengths(thr)))
  if (!full) return(NULL)

  tryCatch({
    m0 <- mirt(d, 1, itemtype = "Rasch", verbose = FALSE)
    v0 <- sx2_cols(itemfit(m0, fit_stats = "S_X2"))
    v2 <- sx2_cols(itemfit(m0, fit_stats = "S_X2", mincell = 5))

    sv <- mirt(d, 1, itemtype = "Rasch", pars = "values")
    sv$value[sv$name == "COV_11"] <- THETA_SD^2
    sv$est[sv$name == "COV_11"]   <- FALSE
    m3 <- mirt(d, 1, itemtype = "Rasch", pars = sv, verbose = FALSE)
    v3 <- sx2_cols(itemfit(m3, fit_stats = "S_X2"))

    out <- cbind(v0, v2, v3)
    colnames(out) <- paste0(rep(c("S", "df", "p"), 3), "_",
                            rep(c("v0", "v2", "v3"), each = 3))
    cbind(out, var_hat = coef(m0, simplify = TRUE)$cov[1, 1])
  }, error = function(e) NULL)
}

## ---- run -------------------------------------------------------------
grid <- expand.grid(rep = seq_len(R_REPS), n = NS, set = names(item_sets),
                    stringsAsFactors = FALSE)

if (file.exists(RES_FILE)) {
  sim <- readRDS(RES_FILE)
} else {
  dir.create(OUT_DIR, showWarnings = FALSE)
  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(grid)),
                  function(i) one_rep(grid$rep[i], grid$set[i], grid$n[i]),
                  mc.cores = NCORES, mc.preschedule = TRUE)
  sim <- list(res = res, grid = grid, item_sets = item_sets,
              mirt_version = as.character(packageVersion("mirt")),
              elapsed_min = as.numeric(Sys.time() - t0, units = "mins"))
  saveRDS(sim, RES_FILE)
}
cat("elapsed:", round(sim$elapsed_min, 1), "min\n")
