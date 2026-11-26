# Specificity study: host run (2026-11-23)

**Question.** Does `RMitemCatFit()` flag items whose categories work as
scored but which misfit in other ways, and does conditional infit tell the
two apart? This tests the workflow the help pages, vignette and book now
recommend: check item fit first, then read a flag as a category problem.

## Design

- Stylised design of the main study: focus item 1 at location 0, steps
  (-1.5, -0.5, 0.5, 1.5), k - 1 PCM items over [-1, 1], theta ~ N(0, 1),
  5 categories. n = 250, 500. k = 5, 10.
- Conditions: all 14 of the main study, plus four calibrated to the real
  phq9 misfit (`R/infit_calibration_stylised.R`):

  | Condition | Mechanism | Focus infit |
  |---|---|---|
  | `U_multi_r75` | second dimension, rho = .75 | 1.22-1.31 |
  | `U_slope_055` | slope 0.55 | 1.23-1.33 |
  | `O_ld_d1` | response dependence on the item nearest 0, d = 1 | 0.80-0.83 |
  | `O_slope_16` | slope 1.6 | 0.80-0.84 |

- Each dataset gets both tests. `RMitemCatFit()`: main-study implementation
  (`analyse_calibration()`), max-T over events, Westfall-Young across items.
  Infit: `easyRasch2::RMitemInfit()` values, p-values from an oracle
  reference set per n, k and base thresholds (1000 PCM datasets at the true
  parameters), studentised two-sided, Westfall-Young across items.
- Replication as the main calibration run: PCM 500 datasets (B = 199),
  others 200 (B = 99).

## Run

From `easyRasch2/dev/category_conditional_null/`, with the dev easyRasch2
(1.4.0.9007 or later) and iarm installed:

```sh
nohup env CORES=10 Rscript R/run_specificity_study.R > results/run_specificity_study.log 2>&1 &
# when "all cells done" is in the log:
CORES=10 Rscript R/summarise_specificity.R
```

`CORES` sets the number of workers (default 4). The run resumes: finished
cells in `results/spec_cells/` and reference sets in
`results/spec_infit_ref/` are skipped. Both folders are gitignored. The
summary is `results/specificity_summary.rds` (+ dated copy), to be tracked.

**Expected runtime: about 2 h on 10 cores.** The main calibration run took
3.3 h on 4 cores for these conditions at n = 250 and 500. The four new
conditions add about 2.5 core-hours, infit about 1 core-hour (0.2 s per
dataset) and the reference sets about 30 core-minutes.

Seeds are fixed per cell (from the cell name) and per reference set. As in
the main study, the random-number streams depend on the number of cores.

## Smoke test (container, 2026-11-23)

`SMOKE=1 CORES=4`: n = 250, k = 5, 4-8 datasets, B = 19, 50 reference
datasets, 5 s. All 28 datasets used. Mean focus infit matched the
calibration (O_ld_d1 0.82, O_slope_16 0.81, U_multi_r75 1.23, U_slope_055
1.29). `RMitemCatFit()` flagged nothing, as expected with B = 19 (smallest p
= .05).
