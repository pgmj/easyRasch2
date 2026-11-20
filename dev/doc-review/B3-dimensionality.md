# B3 Dimensionality: help-page review

2026-11-20. RMdimCFACutoff, RMdimCFA, RMdimCFAPlot, RMdimResidualPCA,
RMdimResidualPCACutoff, RMdimMartinLof, RMdimMartinLofResiduals. Roxygen
only. Rd regenerated, checkRd clean (UTF-8 note for "Löf" only), Rd2txt
scanned, test-documentation passes, all seven examples run
(`run.donttest = TRUE`).

## Changes

| Page | Changes |
|---|---|
| all | Sentence-case titles (the CFA cutoff title was 90 characters, now "Simulated reference for one-factor CFA fit and loadings"), `[fn()]` links, dashes and prose semicolons removed (including example comments, also in RMitemRestscoreBoot from B2), `n_cores` wording aligned with B1/B2, "Monte-Carlo"/"MC" spelled out. |
| RMdimCFACutoff | `data`: complete cases. `iterations`: 1000 for p-values in RMdimCFA. Rosseel (2012) was listed but not cited, now cited for lavaan. Datasets-that-fail paragraph now says what happens when all fail (extreme-only category named). See also adds RMdimResidualPCA and RMdimMartinLof. |
| RMdimCFA | `data`: complete cases. `alpha` and `correction` "used only with p_value = TRUE". Westfall & Young (1993) listed but not cited, now cited. "marginal" now "uncorrected". |
| RMdimCFAPlot | Value now says what the segments are: fixed 66% and 99.8% quantile intervals that do not follow `percentile` (F13). |
| RMdimResidualPCA | Description and rule-of-thumb paragraph tightened, em-dash removed. `output`: "loadings" alias explained. History "(the previous MLE partition dropped extreme-score cases)" removed. See also adds RMdimCFA and RMdimMartinLof. |
| RMdimResidualPCACutoff | `data`: complete cases (the code uses `na.omit`, the page did not say). "Silently dropped" kept in substance ("left out without a message"), since this function does not report attrition (F16). |
| RMdimMartinLof | Title sentence case. Description without "--" dashes and "e.g.". Partition list items no longer joined by semicolons. p-value element rewritten. **Example bug fixed**: it claimed to build two-dimensional data, but the generator ignored theta (`sample(0:4)`, "toy simulation here"). It now simulates PCM data on two correlated latent variables. The test now rejects (p at the floor, subscale WLE r = .41). Tjur (1982) was cited but missing from the references, added (please check the entry). See also adds RMdimMartinLofResiduals and RMdimCFA. |
| RMdimMartinLofResiduals | Title says what it is. "section7" typo. Arguments rewritten without dashes. `flag_threshold` now says these are single-cell flags and some are expected under unidimensionality (F18). |

## Findings in the code (not changed)

- **F13** RMdimCFAPlot draws the loadings segments at fixed 0.1/99.9% and
  16.7/83.3% quantiles, whatever `percentile` is. The flag colour follows
  `percentile`, so the two can disagree. The outer segment could use
  `loading_cutoffs`.
- **F14** RMdimCFAPlot marks flagged loadings by colour alone, red against
  the usual orange (`sienna2`). Red and orange are hard to tell apart with
  colour vision deficiency. Add a shape or fill difference.
- **F15** RMdimResidualPCA returns the variance partition only as
  `attr(result, "variance_partition")` in the data.frame output (attributes
  rule).
- **F16** RMdimResidualPCACutoff and RMdimMartinLof drop failed simulations
  without a message, while the infit and restscore functions report
  attrition in the caption.
- **F17** RMdimCFACutoff lists Yuan & Bentler (2000), robust ML for
  non-normal missing data, which the page never cites and which does not
  seem to concern WLSMV. Drop it, or say what it supports?
- **F18** RMdimMartinLofResiduals flags cells at a fixed |residual| > 2
  (single-cell, no correction). Documented. A note on how many to expect,
  or a simulated or corrected flag, would fit the package's line on fixed
  cutoffs.
- **F19** RMdimCFA's `p_value` defaults to `FALSE`, while RMitemInfit and
  RMitemRestscore default to `NULL` (use the p-values when the full cutoff
  object is given). For B10.
