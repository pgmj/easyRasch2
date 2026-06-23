# Short / low-reliability scale: where resampling theta-hat should over-disperse
suppressMessages(pkgload::load_all(".", quiet = TRUE)); set.seed(7)
N <- 300L; M_truth <- 2000L; K <- 25L; B <- 300L
diff_i   <- seq(-1.6, 1.6, length.out = 5L)          # 5 dichotomous items
thr_list <- as.list(diff_i); sigma <- 1.5
sim_null <- function() as.data.frame(sim_partial_score(thr_list, rnorm(N, 0, sigma)))
stat_mm  <- function(D){ m <- .q3_residual_matrix(D,"CML"); max(m,na.rm=TRUE)-mean(m,na.rm=TRUE) }
# reliability of this scale (for context)
cat(sprintf("Short scale: 5 dich items, latent SD=%.2f, N=%d\n", sigma, N))

truth <- vapply(seq_len(M_truth), function(i) stat_mm(sim_null()), numeric(1))
P <- quantile(truth, c(.90,.95,.99))
cat(sprintf("True null max-mean: mean=%.3f P95=%.3f P99=%.3f\n\n", mean(truth), P[2], P[3]))

res <- data.frame(stat=numeric(K), cut_r=numeric(K), cut_c=numeric(K))
for (k in seq_len(K)) { Dk <- sim_null(); res$stat[k] <- stat_mm(Dk)
  res$cut_r[k] <- RMlocdepQ3Cutoff(Dk, iterations=B, parallel=FALSE, seed=k, dgp="resample")$suggested_cutoff
  res$cut_c[k] <- RMlocdepQ3Cutoff(Dk, iterations=B, parallel=FALSE, seed=k, dgp="conditional")$suggested_cutoff }
p99 <- P[3]
cat(sprintf("True P99=%.3f (nominal FPR 0.01):\n", p99))
print(data.frame(DGP=c("resample","conditional"),
  mean_cutoff=round(c(mean(res$cut_r),mean(res$cut_c)),3),
  bias_vs_P99=round(c(mean(res$cut_r)-p99, mean(res$cut_c)-p99),3),
  realized_fpr=round(c(mean(res$stat>res$cut_r), mean(res$stat>res$cut_c)),3)))
