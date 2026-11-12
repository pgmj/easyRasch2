# Does the *Cutoff() rule "discard a simulated dataset unless every response
# category is used" change the item-fit null, and what should replace it?
# Code review 2026-11-11 (open finding from item_parameters.md, frequency
# study in ../category_frequency/).
#
# Background. All six *Cutoff() functions (infit, restscore, Q3, gamma LD,
# residual PCA, CFA) simulate datasets from the observed CML thresholds and
# discard any dataset in which an item misses one of the categories the
# observed fit has. Off target that discards 56-93% of iterations at n = 100.
# The discarded share is lost: RMitemInfitCutoff(iterations = 400) can return
# 30 iterations, or stop with "All simulation iterations failed".
#
# The observed data always use every category their own fit has (that is how
# the fit gets its threshold count), so the current rule conditions the null
# on an event the observed data satisfy. That suggests the rule is right in
# principle and only wasteful. This study checks that, for the two item-fit
# functions, with three arms on the same simulated datasets:
#
#   A  current rule: of the requested iterations, keep those with every
#      category used. Fewer iterations than requested.
#   B  allow empty top categories: keep every iteration that .cml_fit()
#      accepts (it still stops on an empty lowest or interior category). The
#      refit then has fewer thresholds for that item than the observed fit.
#   C  current rule, topped up: keep simulating until the requested number of
#      iterations use every category (at most cap x requested draws).
#
# Prediction: A and C have the same calibration, C with more stable cutoffs
# and no total failures; B can drift, since part of its null comes from
# datasets unlike the observed one.
#
# Outcomes per cell and arm: iterations used, share of runs with no usable
# iteration, familywise false-flag rate under the default bootstrap p-value
# flagging (correction = "fwer", alpha = .05) and under the interval (HDCI)
# flagging, and the cutoffs themselves. Data are simulated from the PCM, so
# every flag is false.
#
# Each arm calls the package's own simulation worker (run_single_infit_sim(),
# run_single_restscore_sim()), patched only so that the category rule can be
# switched, and the package's own RMitemInfit() / RMitemRestscore() for the
# flags. A fidelity check at the start stops the run unless arm A reproduces
# RMitemInfitCutoff() and RMitemRestscoreCutoff() exactly.
#
# Run from the package root, on the host (not the container):
#   nohup Rscript dev/code-review/cutoff_empty_categories/cutoff_empty_categories.R \
#     300 10 400 > dev/code-review/cutoff_empty_categories/run.log 2>&1 &
# Arguments: R (datasets per cell), cores, iterations (per cutoff call),
# cap (C's draw limit as a multiple of iterations), output folder.
# Container smoke size: R = 4, cores = 4, iterations = 40.

args       <- commandArgs(TRUE)
R_obs      <- if (length(args) >= 1) as.integer(args[1]) else 4L
cores      <- if (length(args) >= 2) as.integer(args[2]) else 4L
iterations <- if (length(args) >= 3) as.integer(args[3]) else 40L
cap        <- if (length(args) >= 4) as.integer(args[4]) else 10L
out_dir    <- if (length(args) >= 5) args[5] else "dev/code-review/cutoff_empty_categories"

suppressMessages(devtools::load_all(".", quiet = TRUE))
ns <- asNamespace("easyRasch2")
dir.create(out_dir, showWarnings = FALSE)
q <- function(expr) suppressWarnings(suppressMessages(expr))

# --- Patched workers ----------------------------------------------------------
# The only change: the discard line becomes
#   if (.note_empty()) return("validation_failed: not all categories represented")
# .note_empty() records that the dataset misses a category and returns
# flag$skip: TRUE discards as the package does (no fit), FALSE goes on to the
# refit, where .cml_fit() still stops on an empty lowest or interior category.
flag <- new.env()
flag$empty <- FALSE
flag$skip <- TRUE
patch_worker <- function(fn) {
  target <- 'return("validation_failed: not all categories represented")'
  txt <- deparse(fn, width.cutoff = 500L)
  hits <- lengths(regmatches(txt, gregexpr(target, txt, fixed = TRUE)))
  if (sum(hits) != 1L) stop("expected one discard line in the worker, found ", sum(hits))
  txt <- sub(target, paste0("if (.note_empty()) ", target), txt, fixed = TRUE)
  out <- eval(parse(text = txt))
  env <- new.env(parent = ns)
  env$.note_empty <- function() {
    flag$empty <- TRUE
    flag$skip
  }
  environment(out) <- env
  out
}
worker <- list(
  infit = patch_worker(ns$run_single_infit_sim),
  restscore = patch_worker(ns$run_single_restscore_sim)
)
# each function's default DGP, as users get it
dgp_default <- c(infit = "resample", restscore = "conditional")

# --- Mirrors of the *Cutoff() set-up and summary -------------------------------
build_data_list <- function(d, dgp) {
  dm <- as.matrix(d)
  pool <- ns$.wle_theta_pool(dm)
  dl <- list(dgp = dgp, type = "polytomous", thr_list = pool$thr_list,
             n_items = ncol(dm), sample_n = nrow(dm), item_names = colnames(dm))
  if (dgp == "resample") {
    dl$thetas <- pool$thetas
    dl$deltaslist <- pool$thr_list
  } else {
    dl$cond_groups <- ns$.cond_groups(dm, pool$thr_list)
  }
  list(dl = dl, thetas = pool$thetas)
}

make_cutoff <- function(stat, ok_results, d, thetas, dgp) {
  iter_dfs <- lapply(seq_along(ok_results), function(i) {
    df <- ok_results[[i]]
    df$iteration <- i
    df
  })
  res <- do.call(rbind, iter_dfs)
  items <- colnames(d)
  if (stat == "restscore") {
    res <- res[, c("iteration", "Item", "Observed", "Expected", "Difference")]
    rownames(res) <- NULL
    ic <- do.call(rbind, lapply(items, function(it) {
      x <- res$Difference[res$Item == it]
      x <- x[is.finite(x)]
      h <- ggdist::hdci(x, .width = 0.95)
      data.frame(Item = it, diff_low = h[1L, 1L], diff_high = h[1L, 2L],
                 stringsAsFactors = FALSE, row.names = NULL)
    }))
  } else {
    rownames(res) <- NULL
    ic <- do.call(rbind, lapply(unique(res$Item), function(it) {
      s <- res[res$Item == it, ]
      hi <- ggdist::hdci(s$InfitMSQ, .width = 0.95)
      ho <- ggdist::hdci(s$OutfitMSQ, .width = 0.95)
      data.frame(Item = it, infit_low = hi[1L, 1L], infit_high = hi[1L, 2L],
                 outfit_low = ho[1L, 1L], outfit_high = ho[1L, 2L],
                 stringsAsFactors = FALSE, row.names = NULL)
    }))
  }
  rownames(ic) <- NULL
  list(results = res, item_cutoffs = ic, actual_iterations = length(ok_results),
       requested_iterations = iterations, sample_n = nrow(d),
       sample_n_total = nrow(d), sample_has_na = FALSE,
       sample_summary = summary(thetas), item_names = items,
       cutoff_method = "hdci", hdci_width = 0.95, dgp = dgp)
}

# --- The three arms from one stream of simulated datasets --------------------
run_arms <- function(stat, d, cut_seed) {
  dgp <- dgp_default[[stat]]
  b <- build_data_list(d, dgp)
  set.seed(cut_seed)
  seeds_main <- sample.int(.Machine$integer.max, iterations)
  set.seed(cut_seed + 1L)
  seeds_extra <- sample.int(.Machine$integer.max, (cap - 1L) * iterations)
  call_worker <- function(s, skip) {
    flag$empty <- FALSE
    flag$skip <- skip
    r <- q(worker[[stat]](s, b$dl))
    list(res = r, ok = is.list(r), empty = flag$empty)
  }
  # requested iterations, every fittable one refitted (arms A and B)
  main <- lapply(seeds_main, call_worker, skip = FALSE)
  ok <- vapply(main, `[[`, TRUE, "ok")
  empty <- vapply(main, `[[`, TRUE, "empty")
  keep_a <- ok & !empty
  # arm C tops A up from further draws, discarding as the package does
  c_res <- lapply(main[keep_a], `[[`, "res")
  draws <- iterations
  for (s in seeds_extra) {
    if (length(c_res) >= iterations) break
    r <- call_worker(s, skip = TRUE)
    draws <- draws + 1L
    if (r$ok) c_res[[length(c_res) + 1L]] <- r$res
  }
  arms <- list(
    A = lapply(main[keep_a], `[[`, "res"),
    B = lapply(main[ok], `[[`, "res"),
    C = c_res
  )
  meta <- c(n_empty = sum(empty), n_unfit = sum(!ok), c_draws = draws)
  list(arms = lapply(arms, function(x) {
    if (length(x) == 0L) NULL else make_cutoff(stat, x, d, b$thetas, dgp)
  }), meta = meta)
}

flag_table <- function(stat, d, cutoff) {
  f <- if (stat == "infit") RMitemInfit else RMitemRestscore
  p <- q(f(d, cutoff = cutoff, output = "dataframe"))
  iv <- q(f(d, cutoff = cutoff, p_value = FALSE, output = "dataframe"))
  if (stat == "infit") {
    data.frame(Item = p$Item, value = p$Infit_MSQ, low = p$Infit_low,
               high = p$Infit_high, p = p$p_infit, padj = p$padj_infit,
               flag_p = nzchar(p$Flagged), flag_int = nzchar(iv$Flagged))
  } else {
    data.frame(Item = p$Item, value = p$Difference, low = p$Diff_low,
               high = p$Diff_high, p = p$p_restscore, padj = p$padj_restscore,
               flag_p = nzchar(p$Flagged), flag_int = nzchar(iv$Flagged))
  }
}

# --- Observed data, as in the frequency study ---------------------------------
gen_thresholds <- function(k, m) {
  lapply(seq_len(k), function(i) {
    loc <- stats::runif(1, -1, 1)
    spread <- stats::runif(1, 0.8, 1.4)
    sort(loc + spread * seq(-1.5, 1.5, length.out = m) + stats::rnorm(m, 0, 0.2))
  })
}
# Redraw until the data pass the package's category check, as a user's data
# must before any *Cutoff() call. Empty top categories are kept.
observed_data <- function(n, k, m, mu) {
  for (tries in 1:200) {
    d <- as.data.frame(ns$sim_partial_score(gen_thresholds(k, m), stats::rnorm(n, mu, 1)))
    colnames(d) <- paste0("I", seq_len(k))
    if (!inherits(try(ns$.check_response_categories(d), silent = TRUE), "try-error")) {
      return(list(d = d, tries = tries))
    }
  }
  stop("no passing dataset in 200 draws")
}

one_rep <- function(p, rep_seed) {
  set.seed(rep_seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  od <- observed_data(p$n, p$k, p$m, p$mu)
  d <- od$d
  cut_seed <- sample.int(1e8, 1L)
  top_empty <- vapply(d, function(x) max(x) < p$m, TRUE)
  rows <- list()
  for (stat in c("infit", "restscore")) {
    ra <- run_arms(stat, d, cut_seed)
    for (arm in names(ra$arms)) {
      co <- ra$arms[[arm]]
      base <- data.frame(n = p$n, k = p$k, cats = p$m + 1L, mu = p$mu,
                         rep = rep_seed, stat = stat, arm = arm,
                         iters = if (is.null(co)) 0L else co$actual_iterations,
                         n_empty = ra$meta[["n_empty"]], n_unfit = ra$meta[["n_unfit"]],
                         c_draws = ra$meta[["c_draws"]], obs_tries = od$tries,
                         n_top_empty = sum(top_empty))
      ft <- if (is.null(co)) {
        data.frame(Item = colnames(d), value = NA_real_, low = NA_real_, high = NA_real_,
                   p = NA_real_, padj = NA_real_, flag_p = NA, flag_int = NA)
      } else {
        flag_table(stat, d, co)
      }
      ft$item_top_empty <- top_empty[ft$Item]
      rows[[length(rows) + 1L]] <- cbind(base[rep(1L, nrow(ft)), ], ft, row.names = NULL)
    }
  }
  do.call(rbind, rows)
}

# --- Fidelity check: arm A must equal the package -----------------------------
fidelity <- function() {
  set.seed(7, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  d <- observed_data(150, 6, 3, -1)$d # off target, so some draws are discarded
  for (stat in c("infit", "restscore")) {
    ra <- run_arms(stat, d, 99L)
    f <- if (stat == "infit") RMitemInfitCutoff else RMitemRestscoreCutoff
    pk <- q(f(d, iterations = iterations, parallel = FALSE, seed = 99L))
    a <- ra$arms$A
    chk <- c(
      iterations = isTRUE(all.equal(a$actual_iterations, pk$actual_iterations)),
      results = isTRUE(all.equal(a$results, pk$results, check.attributes = FALSE)),
      cutoffs = isTRUE(all.equal(a$item_cutoffs, pk$item_cutoffs, check.attributes = FALSE)),
      flags = isTRUE(all.equal(flag_table(stat, d, a), flag_table(stat, d, pk)))
    )
    cat(sprintf("fidelity %-9s A = %d of %d iterations (B %d, C %d, %d draws): %s\n",
                stat, a$actual_iterations, iterations, ra$arms$B$actual_iterations,
                ra$arms$C$actual_iterations, ra$meta[["c_draws"]],
                paste(names(chk), ifelse(chk, "ok", "DIFFERS"), collapse = ", ")))
    if (!all(chk)) stop("arm A does not reproduce the package; fix the mirror before running")
  }
}
fidelity()

# --- Run ------------------------------------------------------------------------
design <- expand.grid(n = c(100L, 200L, 400L), k = c(6L, 10L), mu = c(0, -1),
                      KEEP.OUT.ATTRS = FALSE)
design$m <- ifelse(design$k == 6L, 3L, 4L) # 6 items x 4 categories, 10 x 5

t0 <- Sys.time()
cells <- lapply(seq_len(nrow(design)), function(r) {
  p <- design[r, ]
  set.seed(20261111L + r)
  rep_seeds <- sample.int(.Machine$integer.max, R_obs)
  reps <- parallel::mclapply(rep_seeds, function(s) {
    tryCatch(one_rep(p, s), error = function(e) {
      data.frame(n = p$n, k = p$k, cats = p$m + 1L, mu = p$mu, rep = s,
                 error = conditionMessage(e))
    })
  }, mc.cores = cores, mc.preschedule = FALSE)
  errs <- vapply(reps, function(x) "error" %in% names(x), TRUE)
  if (any(errs)) {
    cat(sprintf("  %d rep(s) failed: %s\n", sum(errs), reps[[which(errs)[1L]]]$error))
  }
  out <- do.call(rbind, reps[!errs])
  cat(sprintf("cell %2d/%d n=%d items=%d cats=%d mu=%g done (%.1f min)\n", r,
              nrow(design), p$n, p$k, p$m + 1L, p$mu,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(out, file.path(out_dir, sprintf("cell_%02d.rds", r))) # partial results
  out
})
raw <- do.call(rbind, cells)

# --- Summary --------------------------------------------------------------------
# one row per run (dataset x statistic x arm); a familywise false flag is any
# flagged item. NA when the arm had no usable iteration.
run_key <- interaction(raw[c("n", "k", "mu", "stat", "arm", "rep")], drop = TRUE)
per_run <- do.call(rbind, lapply(split(raw, run_key), function(s) {
  data.frame(s[1L, c("n", "k", "cats", "mu", "stat", "arm", "rep", "iters")],
             fwer_p = any(s$flag_p), fwer_int = any(s$flag_int))
}))
summ <- do.call(rbind, lapply(
  split(per_run, per_run[c("n", "k", "mu", "stat", "arm")], drop = TRUE),
  function(s) data.frame(
    n = s$n[1], items = s$k[1], cats = s$cats[1], mu = s$mu[1], stat = s$stat[1],
    arm = s$arm[1], runs = nrow(s),
    iters_mean = mean(s$iters), iters_min = min(s$iters),
    none_pct = 100 * mean(s$iters == 0),
    fwer_p_pct = 100 * mean(s$fwer_p, na.rm = TRUE),
    fwer_int_pct = 100 * mean(s$fwer_int, na.rm = TRUE)
  )
))
summ <- summ[order(summ$stat, summ$mu, summ$items, summ$n, summ$arm), ]
rownames(summ) <- NULL

res <- list(raw = raw, summary = summ, design = design, R = R_obs,
            iterations = iterations, cap = cap,
            minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
            session = utils::sessionInfo())
saveRDS(res, file.path(out_dir, "cutoff_empty_categories.rds"))
saveRDS(res, file.path(out_dir, paste0("cutoff_empty_categories_", Sys.Date(), ".rds")))
print(summ, digits = 3, row.names = FALSE)
cat(sprintf("Total %.1f min\n", res$minutes))
