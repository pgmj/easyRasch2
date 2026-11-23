# /code-review medium R/ (2026-11-23)

Target: uncommitted diff of `R/` on top of a5de817 (12 files; help-page
review group E, F21 `Flagged` rename).
Level: medium. Six findings; each checked by hand afterwards.

| # | File | Finding | Verdict |
|---|------|---------|---------|
| M1 | R/martin_lof.R:1099 | RMdimMartinLofResiduals() example still runs `df[df$flagged, ]`; column is now `Flagged` and only exists with `flag_threshold` (default NULL), so the example silently prints 0 rows | Real. Open |
| M2 | jamovi personfit.b.R:155,200; locdepgamma.b.R:305; partgamdif.b.R:191 | read lowercase `$flagged` | Real, but already logged under F21 for the next pin bump (jamovi pinned to b45f3a3). F21 note missed **partgamdif.b.R:191** (RMdifGamma). Not an easyRasch2 fix |
| M3 | R/cfa_cutoff.R:887 | loadings `Flagged` character, fit `Flagged` logical | Deliberate, F21 rule (direction string for two-sided item-level tests). Dismissed |
| M4 | R/cfa_cutoff.R:827 / jamovi cfacutoff.b.R:202 | NA "not tested" shows blank in jamovi | Not a regression: old "not tested" string also mapped to "" via isTRUE(as.logical()). Pre-existing jamovi nit |
| M5 | R/item_cat_fit_plot.R:611 | fallback warning twice | False positive: lines 255 (restscore/calibration) and 278 (residuals with fit) are exclusive branches |
| M6 | R/residual_pca.R:32 | `p_value` doc says default flags "against the cutoff (the 99th percentile)", but a single number is also accepted | Real, small doc nit. Open |

## Fixes (2026-11-23, approved by user)

- M1: example now passes `flag_threshold = 2` and uses `df[df$Flagged, ]`;
  run with run_examples, prints two flagged cells.
- M6: `p_value` doc says "the 99th percentile when the full
  RMdimResidualPCACutoff() result is given".
- M2: partgamdif.b.R:191 added to the F21 note in dev/doc-review/FINDINGS.md.
- Rd regenerated (RMdimMartinLofResiduals.Rd, RMdimResidualPCA.Rd only).
