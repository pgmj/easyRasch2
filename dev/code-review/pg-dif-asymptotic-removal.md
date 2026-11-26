# Code review: uncommitted changes after f44db3e (v9009)

Date: 2026-11-24. Level: medium. Scope: the uncommitted working tree. The main
change takes the asymptotic p-value and CI out of RMdifGamma(). The tree also
has edits to RMitemCatFit and the CICC plot.

## Findings

| ID | Severity | File | Finding | Status |
|----|----------|------|---------|--------|
| PG1 | important | DESCRIPTION:3 | The version changed from `1.4.0.9009` to `1.4.9010`, which drops the patch component. `1.4.9010` sorts above any future `1.4.x` release. It should be `1.4.0.9010`. Confirmed with git diff. | fixed by user |
| PG2 | important | R/dif_partgam.R:326 | RMdifGamma() no longer returns `lower`, `upper`, `padj_bh` or `Significance`. easyRasch2jmv `R/partgamdif.b.R` reads these columns at lines 178-193 and 294-298 (confirmed). The jamovi module breaks when its SHA pin moves. This belongs to jmv plan G5 and needs a decision before crossing the package boundary. | deferred to jmv plan G5 |
| PG3 | nit | R/cicc_plot.R:56 | `@return` says RMdifGamma() gives p-values corrected across items. A plain call gives no p-values now; the user has to pass an RMdifGammaCutoff() result as `cutoff`. | fixed |
| PG4 | nit | R/cicc_plot.R:582, NEWS.md:92 | The caption and the NEWS entry say "Test DIF with RMdifGamma()". They should mention the cutoff, as the RMdifGamma caption already does. | fixed |
| PG5 | nit | R/cicc_plot.R:122 | The text describes "`p_value = TRUE` or a cutoff" as two routes. `p_value = TRUE` without a cutoff object is an error, so there is only one route. | fixed |
| PG6 | nit | tests/testthat/testthat-problems.rds | A leftover file from a failed test run on 2026-11-24 09:20. It is not in .gitignore. Delete it and add it to .gitignore and .Rbuildignore. | fixed |
