## Submission summary

This is the first CRAN submission of **easyRasch2**, an R package that
streamlines Rasch measurement theory workflows (local dependence,
conditional infit, item-restscore, DIF via LR test / partial gamma /
Rasch trees, dimensionality via Martin-Löf / residual PCA / CFA,
targeting, reliability, and score-to-WLE transformations). It is the
successor of the GitHub-only `easyRasch` package and rebuilds the API
around a consistent `RM*()` naming scheme with full test coverage and
simulation-based cutoffs for several diagnostics.

## Test environments

* Local: macOS 14.6, R 4.5.3 — `R CMD check --as-cran` clean
  (0 errors / 0 warnings / 0 notes).
* R-hub v2 (`rhub::rhub_check()`), Linux container — clean
  (0 errors / 0 warnings / 0 notes).
* `devtools::check_win_*()` only 1 NOTE for New submission

## R CMD check results

```
0 errors | 0 warnings | 0 notes
```

## CRAN policy notes

* **License**: GPL (>= 3). Contributed code from the `raschtreeMH` and
  `effecttree` packages (Mirka Henninger, Jan Radek) is adapted under
  its original MIT licence; contributors are credited in `Authors@R`
  and in the relevant source files (`R/dif_tree.R`). The MIT terms
  are compatible with GPL (>= 3).
* **Examples**:
  * Fast-running examples (well under 5 s) are unwrapped.
  * Examples that exceed CRAN's 5 s budget but run correctly are wrapped
    in `\donttest{}` so they are still executed by
    `R CMD check --run-donttest` and by reverse-dependency checks.
  * `\dontrun{}` is reserved for examples that genuinely cannot run on
    CRAN: parallel back-ends via `mirai`, multiple-imputation workflows
    via `mice`, and the `mirt`-based residual-PCA simulations whose
    runtime exceeds reasonable check budgets.
* **Suggested packages** are all guarded with `requireNamespace(...,
  quietly = TRUE)` before use, with a clear `stop()` message instructing
  the user how to install them. The package's hard dependencies
  (`eRm`, `psychotools`, `mirt`, `knitr`, `parallel`, `stats`, `utils`,
  `rlang`) are all on CRAN.
* **Writing to user filespace**: no examples, tests, or vignettes write
  to the user's home directory, working directory, or any location
  outside `tempdir()`.
* **No `T`/`F` shorthand**, no `<<-` into the global environment, no
  `library()`/`require()` calls inside functions, no `cat()`/`print()`
  to stdout outside of intended user-facing output (which uses
  `message()` or `knitr::kable()`).
* **Parallel processing**: all functions that support parallel
  evaluation (via `mirai`) default to **sequential** execution and only
  use multiple workers when the user explicitly opts in. The CRAN
  worker limit of 2 cores is respected; examples and tests never
  request more than one worker.
* **Reproducibility**: every simulation-based function accepts a `seed`
  argument; tests set seeds explicitly.

## Downstream dependencies

There are currently no reverse dependencies on CRAN (first submission).
A companion jamovi module, `easyRasch2jmv`, calls into this package and
is maintained in parallel; it is **not** a CRAN package and therefore
imposes no reverse-dependency burden.

## Additional notes for the reviewer

* The package targets researchers working with ordinal/dichotomous
  questionnaire data. It deliberately re-exports a small layer of
  opinionated wrappers around `eRm`, `psychotools`, `iarm`, `mirt`,
  and `lavaan` rather than duplicating their estimation routines.
* Several diagnostics provide parametric-bootstrap cutoff alternatives
  to rule-of-thumb thresholds (e.g. `RMitemInfitCutoff()`,
  `RMdifGammaCutoff()`, `RMlocdepGammaCutoff()`, `RMlocdepQ3Cutoff()`,
  `RMdimResidualPCACutoff()`, `RMdimCFACutoff()`). These are the most computationally
  expensive entry points and account for the majority of `\donttest{}`
  / `\dontrun{}` examples.
* Spelling: the package uses British spellings in some user-facing
  prose ("behaviour", "centred") consistent with the underlying
  statistical literature. No misspellings are expected from
  `urlchecker::url_check()` or `devtools::spell_check()`.
* On win-builder (R-devel), `R CMD check` flags the following words in
  DESCRIPTION as "possibly misspelled". All are legitimate technical
  terms or abbreviations and are spelled correctly:
    - **CFA**  -- standard abbreviation for *confirmatory factor
      analysis*.
    - **Rasch** -- proper noun, the surname of Georg Rasch, after
      whom the family of measurement models is named.
    - **infit** -- a standard Rasch-model fit statistic (the
      *information-weighted* mean-square residual; Smith, 1986).
    - **unidimensionality** -- the standard psychometric term for the
      assumption that a single latent trait underlies item responses.
    - **et / al** -- parts of the inline citation
      *"Christensen et al. (2021)"*.


Thank you for your time reviewing this submission.
