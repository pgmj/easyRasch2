# How often do empty response categories occur in observed data, in the
# package's parametric simulations, in row-bootstrap resamples, and within
# DIF subgroups? Code review 2026-11-10, item_parameters.R unit, point 2 of
# the shared category check decision.
#
# Run from the package root:
#   Rscript dev/code-review/category_frequency/category_frequency.R [R] [B] [cores]
# Container smoke size: R = 60, B = 40, cores = 4.

args  <- commandArgs(TRUE)
R_obs <- if (length(args) >= 1) as.integer(args[1]) else 60L
B     <- if (length(args) >= 2) as.integer(args[2]) else 40L
cores <- if (length(args) >= 3) as.integer(args[3]) else 4L

suppressMessages(devtools::load_all(".", quiet = TRUE))
ns <- asNamespace("easyRasch2")
out_dir <- "dev/code-review/category_frequency"

# Shape of each item: "ok", "top" (only the top category empty, allowed),
# "fail" (constant item, or an empty lowest or interior category).
item_shape <- function(x, m) {
  x <- x[!is.na(x)]
  if (length(unique(x)) <= 1L) return("fail")
  tab <- tabulate(x + 1L, nbins = m + 1L)
  if (all(tab > 0L)) return("ok")
  empty <- which(tab == 0L) - 1L
  top_run <- rev(cumprod(rev(tab == 0L)))  # 1 for categories in an empty top run
  if (all(top_run[empty + 1L] == 1L)) "top" else "fail"
}
data_shape <- function(d, m) {
  s <- vapply(seq_len(ncol(d)), function(j) item_shape(d[, j], m), "")
  c(fail = any(s == "fail"), any_empty = any(s != "ok"),
    n_items_empty = sum(s != "ok"))
}

gen_thresholds <- function(k, m) {
  lapply(seq_len(k), function(i) {
    loc <- stats::runif(1, -1, 1)
    spread <- stats::runif(1, 0.8, 1.4)
    sort(loc + spread * seq(-1.5, 1.5, length.out = m) + stats::rnorm(m, 0, 0.2))
  })
}

one_dataset <- function(n, k, m, mu) {
  thr <- gen_thresholds(k, m)
  d <- ns$sim_partial_score(thr, stats::rnorm(n, mu, 1))
  colnames(d) <- paste0("I", seq_len(k))
  obs <- data_shape(d, m)
  res <- c(obs_fail = obs[["fail"]], obs_any_empty = obs[["any_empty"]],
           par_any_empty = NA, par_fail = NA, boot_fail = NA,
           split50_any_empty = NA, split50_items = NA,
           split70_any_empty = NA, split70_items = NA)
  if (obs[["fail"]]) return(res)

  # Parametric DGP used by the *Cutoff() functions: observed CML thresholds and
  # resampled finite WLE locations (.wle_theta_pool), responses from
  # sim_partial_score().
  pool <- ns$.wle_theta_pool(d)
  par <- t(replicate(B, {
    s <- ns$sim_partial_score(pool$thr_list, sample(pool$thetas, n, replace = TRUE))
    data_shape(s, m)[c("any_empty", "fail")]
  }))
  boot <- replicate(B, data_shape(d[sample.int(n, n, replace = TRUE), , drop = FALSE], m)[["fail"]])
  split <- function(p) {
    sh <- t(replicate(B, {
      g <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(p, 1 - p))
      a <- data_shape(d[g, , drop = FALSE], m); b <- data_shape(d[!g, , drop = FALSE], m)
      # items excluded by eRm::LRtest: any category empty in either group
      ex <- vapply(seq_len(k), function(j)
        item_shape(d[g, j], m) != "ok" || item_shape(d[!g, j], m) != "ok", TRUE)
      c(any = a[["any_empty"]] || b[["any_empty"]], items = sum(ex))
    }))
    colMeans(sh)
  }
  s50 <- split(0.5); s70 <- split(0.7)
  res[c("par_any_empty", "par_fail", "boot_fail")] <-
    c(mean(par[, "any_empty"]), mean(par[, "fail"]), mean(boot))
  res[c("split50_any_empty", "split50_items", "split70_any_empty", "split70_items")] <-
    c(s50[["any"]], s50[["items"]], s70[["any"]], s70[["items"]])
  res
}

design <- expand.grid(n = c(100L, 200L, 400L), k = c(6L, 10L), m = c(3L, 4L),
                      mu = c(0, -1), KEEP.OUT.ATTRS = FALSE)

RNGkind("L'Ecuyer-CMRG"); set.seed(20261110)
t0 <- Sys.time()
cells <- lapply(seq_len(nrow(design)), function(r) {
  p <- design[r, ]
  reps <- parallel::mclapply(seq_len(R_obs), function(i)
    suppressWarnings(one_dataset(p$n, p$k, p$m, p$mu)), mc.cores = cores, mc.set.seed = TRUE)
  M <- do.call(rbind, reps)
  cat(sprintf("cell %2d/%d n=%d k=%d cats=%d mu=%g done (%.1f min)\n", r, nrow(design),
              p$n, p$k, p$m + 1, p$mu, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  cbind(p[rep(1, nrow(M)), ], M)
})
raw <- do.call(rbind, cells)

summ <- stats::aggregate(
  cbind(obs_fail, obs_any_empty, par_any_empty, par_fail, boot_fail,
        split50_any_empty, split50_items, split70_any_empty, split70_items) ~ n + k + m + mu,
  data = transform(raw, obs_fail = as.numeric(obs_fail), obs_any_empty = as.numeric(obs_any_empty)),
  FUN = function(z) mean(z, na.rm = TRUE), na.action = stats::na.pass)
summ$categories <- summ$m + 1L

stamp <- format(Sys.Date(), "%Y%m%d")
saveRDS(list(raw = raw, summary = summ, R = R_obs, B = B), file.path(out_dir, "category_frequency.rds"))
saveRDS(list(raw = raw, summary = summ, R = R_obs, B = B), file.path(out_dir, paste0("category_frequency_", stamp, ".rds")))
print(summ[order(summ$mu, summ$categories, summ$k, summ$n),
           c("n", "k", "categories", "mu", "obs_fail", "obs_any_empty", "par_any_empty",
             "par_fail", "boot_fail", "split50_any_empty", "split50_items",
             "split70_any_empty", "split70_items")], digits = 2, row.names = FALSE)
