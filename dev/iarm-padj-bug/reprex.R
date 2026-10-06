# iarm 0.4.3: partgam_DIF() and partgam_LD() return a Bonferroni correction
# in the adjusted p-value column, whatever `p.adj` method is requested.
#
# Inside the item (or item-pair) loop, each p-value is adjusted on its own:
#   pkorr <- p.adjust(pvalue, method = padj, n = l * k)          # partgam_DIF
#   pkorr <- p.adjust(pvalue, method = padj, n = (k * (k - 1)))  # partgam_LD
# With a single p-value, BH, Holm, Hochberg and Bonferroni all return
# min(1, n * p). item_restscore(), out_infit() and boot_fit() adjust the
# whole vector at once and are not affected.

library(iarm)
packageVersion("iarm")

set.seed(1)
n <- 400
k <- 6
theta <- rnorm(n)
grp <- factor(rep(c("a", "b"), each = n / 2))
X <- as.data.frame(sapply(seq(-1.5, 1.5, length.out = k), function(b) {
  rbinom(n, 1, plogis(theta - b - 0.6 * (grp == "b") * (b == -1.5)))
}))
names(X) <- paste0("I", 1:k)

res <- partgam_DIF(X, grp) # default p.adj = "BH"
p <- as.numeric(res$pvalue)
padj_iarm <- as.numeric(res[[6]])

data.frame(
  item = res$Item,
  p = signif(p, 3),
  iarm_padj_BH = signif(padj_iarm, 3),
  bonferroni = signif(pmin(1, k * p), 3),
  BH = signif(p.adjust(p, method = "BH"), 3)
)

# The "BH" column is the Bonferroni correction
all.equal(padj_iarm, pmin(1, k * p)) # TRUE

# and the requested method makes no difference
res_holm <- partgam_DIF(X, grp, p.adj = "holm")
all.equal(as.numeric(res_holm[[6]]), padj_iarm) # TRUE

# partgam_LD() behaves the same way over its k * (k - 1) tests
ld <- partgam_LD(X)
p_ld <- as.numeric(ld[[1]][[5]])
all.equal(as.numeric(ld[[1]][[6]]), pmin(1, k * (k - 1) * p_ld)) # TRUE

# Possible fix: collect the p-values in the loop and adjust once afterwards,
# e.g. for partgam_DIF()
#   result$pkorr <- p.adjust(result$pvalue, method = padj)
# and for partgam_LD() over all k * (k - 1) tests in both tables together,
# then derive the significance symbols from the adjusted values.
