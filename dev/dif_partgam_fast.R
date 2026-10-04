# Prototype: vectorised partial gamma for DIF. Development code for
# dev/dif-impact-design.md (Part A, question 4), not package code.
#
# Generalises .partgam_one() in R/ld_partgam.R from an m x m table (item by
# item) to an mx x my table (item category by group), conditioned on the full
# score. For one stratum with count matrix N (rows: item categories, columns:
# groups in factor-level order), Gx[a, a'] = 1 iff a' > a and Gy likewise:
#   A = Gx N            A[a, b] = sum over a' > a of N[a', b]
#   C = sum(N * (A t(Gy)))   pairs where the other observation is higher on both
#   D = sum(N * (A Gy))      higher on the item, lower on the group code
# Each unordered pair is counted once. gamma = (C - D) / (C + D) pooled over
# strata, Davis (1967).
#
# Coding follows iarm::partgam_DIF(): complete cases on items and grouping
# variable, score = sum over all items (the item included), group order = the
# factor levels (character and numeric variables are converted with factor(),
# so alphabetical or numeric order). Positive gamma: higher-coded groups score
# higher on the item at the same total score.

.partgam_dif_one <- function(x, y, z, mx, my, Gx, Gy, tGy) {
  zc <- z - min(z) + 1L
  nz <- max(zc)
  counts <- tabulate((zc - 1L) * mx * my + y * mx + x + 1L, nbins = mx * my * nz)
  dim(counts) <- c(mx, my, nz)
  conc <- 0
  disc <- 0
  for (s in seq_len(nz)) {
    N <- counts[, , s]
    if (sum(N) < 2L) next
    A <- Gx %*% N
    conc <- conc + sum(N * (A %*% tGy))
    disc <- disc + sum(N * (A %*% Gy))
  }
  total <- conc + disc
  if (total == 0) NA_real_ else (conc - disc) / total
}

partgam_dif_fast <- function(data, dif_var) {
  X <- as.matrix(data)
  storage.mode(X) <- "integer"
  items <- colnames(X)
  if (is.null(items)) items <- paste0("I", seq_len(ncol(X)))
  f <- if (is.factor(dif_var)) dif_var else factor(dif_var)
  ok <- stats::complete.cases(X) & !is.na(f)
  X <- X[ok, , drop = FALSE]
  f <- f[ok]
  y <- as.integer(f) - 1L
  my <- nlevels(f)
  mx <- max(X) + 1L
  score <- rowSums(X)
  upper <- function(m) outer(seq_len(m), seq_len(m), function(a, b) as.numeric(b > a))
  Gx <- upper(mx)
  Gy <- upper(my)
  tGy <- t(Gy)
  data.frame(
    Item = items,
    gamma = vapply(seq_len(ncol(X)), function(i) {
      .partgam_dif_one(X[, i], y, score, mx, my, Gx, Gy, tGy)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

# Brute-force reference straight from the definition: within each stratum of
# the total score, count every unordered pair of respondents as concordant or
# discordant on (item, group code). Used where iarm errors (constant item).
partgam_dif_brute <- function(data, dif_var) {
  X <- as.matrix(data)
  f <- if (is.factor(dif_var)) dif_var else factor(dif_var)
  ok <- stats::complete.cases(X) & !is.na(f)
  X <- X[ok, , drop = FALSE]
  y <- as.integer(f[ok])
  score <- rowSums(X)
  vapply(seq_len(ncol(X)), function(i) {
    C <- 0
    D <- 0
    for (s in unique(score)) {
      idx <- which(score == s)
      if (length(idx) < 2L) next
      sx <- sign(outer(X[idx, i], X[idx, i], "-"))
      sy <- sign(outer(y[idx], y[idx], "-"))
      prod <- (sx * sy)[upper.tri(sx)]
      C <- C + sum(prod > 0)
      D <- D + sum(prod < 0)
    }
    if (C + D == 0) NA_real_ else (C - D) / (C + D)
  }, numeric(1))
}

# iarm reference with its console output suppressed. Returns NULL on error.
partgam_dif_iarm <- function(data, dif_var) {
  out <- NULL
  utils::capture.output(
    out <- tryCatch(iarm::partgam_DIF(as.data.frame(data), dif_var),
                    error = function(e) NULL)
  )
  if (is.null(out)) return(NULL)
  as.numeric(out$gamma)
}
