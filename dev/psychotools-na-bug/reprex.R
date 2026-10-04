# Minimal reproducible example: pcmodel() and rsmodel() return wrong estimates
# when no response pattern is complete and all NA patterns have the same shape.
#
# psychotools 0.7.6, R 4.x. Needs only psychotools; eRm is used at the end as an
# independent CML reference and can be skipped.

library(psychotools)
packageVersion("psychotools")

# Partial credit data: 8 items, 3 categories, n = 1000 -------------------------
set.seed(1)
n <- 1000
k <- 8
theta <- rnorm(n)
loc <- seq(-1.5, 1.5, length.out = k)
rpcm <- function(theta, thr) {
  eta <- cbind(0, outer(theta, thr, "-"))
  eta <- t(apply(eta, 1, cumsum))
  p <- exp(eta - apply(eta, 1, max))
  p <- p / rowSums(p)
  apply(p, 1, function(pr) sample(seq_along(pr) - 1L, 1L, prob = pr))
}
y <- sapply(loc, function(l) rpcm(theta, l + c(-0.5, 0.5)))
colnames(y) <- paste0("I", 1:k)

# Every respondent skips exactly one item, so no row is complete -------------
y_na <- y
y_na[cbind(1:n, sample(1:k, n, replace = TRUE))] <- NA
sum(complete.cases(y_na)) # 0

# The complete-data fit is fine, the incomplete one is not --------------------
m_full <- pcmodel(y)
m_na <- pcmodel(y_na) # warning: could not invert Hessian
logLik(m_full) # -4581
logLik(m_na) # positive: +1188555
round(itempar(m_full), 2)
round(itempar(m_na), 2) # up to 0.6 logits off the complete-data estimates

# The same data with five complete rows added back works ----------------------
y_na5 <- y_na
y_na5[1:5, ] <- y[1:5, ]
round(itempar(pcmodel(y_na5)), 2)

# rsmodel() has the same problem ----------------------------------------------
logLik(rsmodel(y_na)) # positive: +828067

# Cause -----------------------------------------------------------------------
# In the missing-data branch the per-pattern parameter lists are built as
#   mapply(split, esf_par, mapply(function(x, y) rep(1:x, y), m_i, oj_max_i))
# (pcmodel: in cloglik(), agrad() and the `full` ESF computation; rsmodel: the
# corresponding three calls with rep.int and oj_i). When every NA pattern has
# the same number of items and categories, both mapply() calls simplify their
# results from lists to matrices, so elementary_symmetric_functions() receives
# the wrong structure. A complete-case pattern has a different length, which
# keeps the results as lists and hides the problem. Dichotomous data are not
# affected in practice because each item has a single parameter.
#
# Fix: SIMPLIFY = FALSE on both mapply() calls, in all three places.

patch <- function(fun) {
  src <- deparse(fun)
  h <- grep("mapply(function(x, ", src, fixed = TRUE)
  h <- h[grepl("\\)\\)$", src[h + 1])]
  src[h + 1] <- sub("\\)\\)$", ", SIMPLIFY = FALSE), SIMPLIFY = FALSE)", src[h + 1])
  f <- eval(parse(text = src))
  environment(f) <- asNamespace("psychotools")
  f
}
pcmodel_fixed <- patch(pcmodel)
rsmodel_fixed <- patch(rsmodel)

m_fix <- pcmodel_fixed(y_na)
logLik(m_fix) # -3813, fewer observations than the complete data
round(itempar(m_fix), 2) # within 0.04 of the complete-data estimates

# Unchanged where the original already worked
all.equal(coef(pcmodel(y)), coef(pcmodel_fixed(y)))
y_mcar <- y
y_mcar[matrix(runif(n * k) < 0.1, n)] <- NA
all.equal(coef(pcmodel(y_mcar)), coef(pcmodel_fixed(y_mcar)))
all.equal(coef(rsmodel(y_mcar)), coef(rsmodel_fixed(y_mcar)))

# Independent check against eRm (optional) ------------------------------------
if (requireNamespace("eRm", quietly = TRUE)) {
  tp <- unlist(threshpar(m_fix))
  te <- as.vector(t(eRm::thresholds(eRm::PCM(y_na))$threshtable[[1]][, -1]))
  print(max(abs((tp - mean(tp)) - (te - mean(te))))) # about 2e-5
}
