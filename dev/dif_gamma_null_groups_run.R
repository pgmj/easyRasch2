# Runner for dif_gamma_null_groups.qmd (same chunks). Resumable: skips finished conditions.
#| label: setup
suppressMessages(devtools::load_all("..", quiet = TRUE))
source("dif_partgam_fast.R")
library(parallel)

out_dir <- "dif_gamma_null_groups"
dir.create(out_dir, showWarnings = FALSE)
R_reps <- 300L
B_iter <- 400L
n_true <- 2000L
cores <- 10L

item_sets <- list(
  poly9 = lapply(seq(-1.5, 1.5, length.out = 9), function(l) l + c(-0.8, 0, 0.8)),
  dich15 = as.list(seq(-2, 2, length.out = 15))
)
for (nm in names(item_sets)) {
  names(item_sets[[nm]]) <- paste0("I", seq_along(item_sets[[nm]]))
}
grid <- expand.grid(
  items = names(item_sets),
  delta = c(0, 0.5, 1),
  sizes = c("150/150", "480/120"),
  stringsAsFactors = FALSE
)
grid$id <- sprintf("%s_d%s_%s", grid$items, grid$delta, sub("/", "-", grid$sizes))

#| label: generators
# Vectorised PCM simulation (distributionally identical to the package's
# sim_partial_score() and psychotools::rrm(), much faster).
sim_pcm_fast <- function(theta, thr_list) {
  n <- length(theta)
  vapply(thr_list, function(thr) {
    m <- length(thr)
    cum <- outer(theta, 0:m) - matrix(c(0, cumsum(thr)), n, m + 1L, byrow = TRUE)
    p <- exp(cum - apply(cum, 1, max))
    cp <- t(apply(p / rowSums(p), 1, cumsum))
    as.integer(rowSums(stats::runif(n) > cp[, -(m + 1L), drop = FALSE]))
  }, integer(n))
}

valid_sim <- function(X, n_cats) {
  if (all(n_cats == 2L)) return(all(colSums(X) >= 8L))
  all(vapply(seq_len(ncol(X)), function(j) all(tabulate(X[, j] + 1L, nbins = n_cats[j]) > 0L), logical(1)))
}

gen_observed <- function(cond) {
  n <- as.integer(strsplit(cond$sizes, "/")[[1]])
  thr <- item_sets[[cond$items]]
  theta <- c(stats::rnorm(n[1], 0, 1.5), stats::rnorm(n[2], cond$delta * 1.5, 1.5))
  g <- factor(rep(c("g1", "g2"), n))
  X <- sim_pcm_fast(theta, thr)
  colnames(X) <- names(thr)
  list(X = X, g = g)
}

gam <- function(X, g) {
  out <- partgam_dif_fast(X, g)
  stats::setNames(out$gamma, out$Item)
}

#| label: nulls
# One dataset: observed gamma plus B draws from each null.
null_draws <- function(obs, B) {
  X <- obs$X
  g <- obs$g
  n <- nrow(X)
  n_cats <- apply(X, 2, max) + 1L
  pool <- .wle_theta_pool(X)
  thr <- pool$thr_list
  wle <- .estimate_thetas(X, thr, method = "WLE")$theta
  props <- as.numeric(table(g)) / n
  lev <- levels(g)
  wle_g <- lapply(lev, function(l) { w <- wle[g == l]; w[is.finite(w)] })
  n_g <- as.integer(table(g))
  cond_dl <- list(thr_list = thr, cond_groups = .cond_groups(X, thr),
                  sample_n = n, n_items = ncol(X))
  score <- rowSums(X)
  strata <- split(seq_len(n), score)

  draw <- function(method) {
    repeat {
      if (method == "random") {
        th <- sample(pool$thetas, n, replace = TRUE)
        Xs <- sim_pcm_fast(th, thr)
        gs <- factor(sample(lev, n, replace = TRUE, prob = props), levels = lev)
        if (valid_sim(Xs, n_cats)) break
      } else if (method == "groupwise") {
        th <- unlist(lapply(seq_along(lev), function(k) sample(wle_g[[k]], n_g[k], replace = TRUE)))
        Xs <- sim_pcm_fast(th, thr)
        gs <- factor(rep(lev, n_g), levels = lev)
        if (valid_sim(Xs, n_cats)) break
      } else if (method == "conditional") {
        Xs <- as.matrix(.sim_cond_dataset(cond_dl))
        gs <- g
        break
      } else {
        Xs <- X
        gs <- g
        for (s in strata) if (length(s) > 1L) gs[s] <- gs[s[sample.int(length(s))]]
        break
      }
    }
    colnames(Xs) <- colnames(X)
    gam(Xs, gs)
  }
  methods <- c("random", "groupwise", "conditional", "permutation")
  sims <- lapply(stats::setNames(methods, methods), function(mt) {
    do.call(rbind, lapply(seq_len(B), function(b) draw(mt)))
  })
  list(observed = gam(X, g), sims = sims)
}

summarise_dataset <- function(nd) {
  do.call(rbind, lapply(names(nd$sims), function(mt) {
    pv <- .bootstrap_pvalues(nd$observed, nd$sims[[mt]], correction = "fwer", tail = "two.sided")
    data.frame(
      method = mt,
      item = pv$name,
      p = pv$p,
      padj = pv$padj,
      null_sd = apply(nd$sims[[mt]], 2, stats::sd)[pv$name],
      stringsAsFactors = FALSE
    )
  }))
}

#| label: run
RNGkind("L'Ecuyer-CMRG")
for (i in seq_len(nrow(grid))) {
  cond <- grid[i, ]
  f <- file.path(out_dir, paste0(cond$id, ".rds"))
  if (file.exists(f)) next
  t0 <- Sys.time()
  set.seed(1000L + i)
  true_g <- do.call(rbind, mclapply(seq_len(n_true), function(r) {
    o <- gen_observed(cond)
    gam(o$X, o$g)
  }, mc.cores = cores))
  true_sd <- apply(true_g, 2, stats::sd)
  res <- mclapply(seq_len(R_reps), function(r) {
    o <- gen_observed(cond)
    out <- tryCatch(summarise_dataset(null_draws(o, B_iter)), error = function(e) NULL)
    if (!is.null(out)) out$rep <- r
    out
  }, mc.cores = cores)
  res <- do.call(rbind, res)
  res$true_sd <- true_sd[res$item]
  res <- cbind(cond[rep(1, nrow(res)), c("items", "delta", "sizes")], res)
  attr(res, "runtime") <- Sys.time() - t0
  saveRDS(res, f)
}
