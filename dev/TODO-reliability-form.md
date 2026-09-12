# Decided: which reliability formula `RMreliability()` reports

Status: **decided 2026-09-11, option 3 (switch), implemented.** Raised when
`RMreliabilityCurve()` was added. This note keeps the evidence and the
reasoning; see the closing section for what was done.

## The two formulas

With `sigma` the MML estimate of the latent SD and `SEM(theta) = 1/sqrt(I(theta))`:

| | formula | range |
|---|---|---|
| subtractive (Green/Lord) | `1 - SEM^2 / sigma^2` | unbounded below |
| ratio | `sigma^2 / (sigma^2 + SEM^2)` | (0, 1) |

`.marginal_rxx()` uses the subtractive form, so it is what `RMreliability()`
prints in its "Marginal" row (floored at 0). `RMreliabilityCurve()` uses the
ratio form for its reliability axis.

## Why the curve went with the ratio form

1. **It is bounded.** The subtractive form is negative whenever
   `SEM(theta) > sigma`, which is common at the tails of a short scale and at
   the floor of a skewed one. A curve that dives below zero at both ends is
   unreadable, and flooring it hides the part of the scale the figure exists
   to show.
2. **It matches what sigma is.** The MML estimate is the SD of the latent
   trait, i.e. true-score SD, so `sigma^2 + SEM^2` is the implied observed
   variance and the ratio is the true/observed proportion. The subtractive
   form treats sigma^2 as the observed variance instead.
3. **It tracks observed-score reliability far better.** Milanzi et al. (2015,
   eq. 12) criticise exactly the subtractive coefficient. Their design was
   replicated in `dev/milanzi_check.R` (24 binary items, b ~ U(-4,4),
   2PL a ~ N(2, 0.8^2), N = 400, no models fitted, 1000 replications per cell),
   and the comparison extended to polytomous Rasch data.

### Part A, binary, replicating their Table 1

Mean absolute error against the exact expected-sum-score reliability
`Var(mu) / (Var(mu) + Var(eps))`:

| Model | sigma^2 | exact | subtractive | ratio | ratio (curve avg) | Taylor | % subtractive < 0 |
|---|---|---|---|---|---|---|---|
| 1PL | 0.25 | .416 | -.427 | .415 | .416 | .417 | 99.7 |
| 1PL | 1 | .737 | .636 | .734 | .735 | .739 | 0 |
| 1PL | 4 | .912 | .889 | .900 | .904 | .918 | 0 |
| 2PL | 0.25 | .533 | .227 | .575 | .579 | .527 | 16.8 |
| 2PL | 1 | .821 | .805 | .839 | .841 | .809 | 0 |
| 2PL | 4 | .947 | .913 | .922 | .942 | .943 | 0.2 |

Overall MAE: **subtractive 0.222, ratio 0.018**, ratio-of-the-curve 0.015,
Milanzi's own Taylor approximation 0.012.

### Part B, polytomous Rasch, the package's territory

PCM, 4 categories, 500 replications per cell, using the package's
`.test_information()` (a vectorised equivalent was checked against it to
3.4e-15).

| items | sigma | exact | subtractive | ratio | MAE subtractive | MAE ratio | % < 0 |
|---|---|---|---|---|---|---|---|
| 5 | 0.5 | .423 | -.414 | .418 | .837 | .006 | 99.8 |
| 10 | 0.5 | .593 | .303 | .591 | .290 | .003 | 0.2 |
| 20 | 0.5 | .745 | .655 | .744 | .090 | .002 | 0 |
| 5 | 1.0 | .733 | .591 | .711 | .142 | .022 | 0 |
| 10 | 1.0 | .846 | .803 | .836 | .042 | .010 | 0 |
| 20 | 1.0 | .916 | .903 | .912 | .013 | .005 | 0 |
| 5 | 1.5 | .850 | .742 | .796 | .108 | .054 | 0 |
| 10 | 1.5 | .918 | .880 | .893 | .038 | .025 | 0 |
| 20 | 1.5 | .957 | .942 | .945 | .016 | .012 | 0 |

Overall MAE: **subtractive 0.175, ratio 0.015**, ratio-of-the-curve 0.010.

### What the replication settles, and what it does not

**Settled.** The earlier recomputation by inverting their published rho_f was
sound. The identity `ratio = 1 / (2 - rho_f)` holds to 2.2e-16 in every
replication, and the reading of their eq. 12 as the subtractive coefficient is
confirmed by reproducing its negative value in the low-variance 1PL cell. The
choice of sample `Var(theta)` over the true sigma^2 changes nothing (.415 vs
.416 and so on). The latent-correlation column reproduces their published
values to three decimals (1PL: .645/.879/.967 against .646/.879/.967), which
anchors the replication, and their headline claim that the Taylor
approximation is closest to exact is reproduced.

**Not a clean sweep.** Holding one item draw fixed across cells, as they did,
the ratio form is closer than the subtractive in 100% of draws for all three
1PL cells and 99% / 47% / 96.5% for the 2PL cells at sigma^2 = 0.25 / 1 / 4.
The 2PL sigma^2 = 1 cell is a coin flip. The overall gap is driven entirely by
the low-information cases. Where both forms are well behaved, choosing between
them on accuracy grounds is not possible, and the argument rests on
boundedness instead.

**Level differences from their Table 1.** Their table is one realization and
several of its values sit above the 97.5th percentile of the replication
spread, most visibly for the 2PL. This does not bear on the comparison, which
is paired within replication, but it does mean the simulated means should not
be read as failed reproductions of their published numbers.

**One finding that touched the code, now acted on.** Averaging the ratio form
over the curve beats the ratio-of-averages in both parts (0.015 vs 0.018
binary, 0.010 vs 0.015 polytomous). `RMreliabilityCurve()`'s `marginal_ratio`
was switched to the density-weighted mean of the curve on 2026-09-11, and the
dashed reference line on the reliability axis follows it. `sigma^2 / (sigma^2 +
x)` is convex in x, so the curve mean is the larger of the two, and the exact
value sits above both. `marginal_green` still uses the ratio of averages,
because it has to keep matching `.marginal_rxx()`.

## The decision

**Option 3: switch to the ratio form, and rename the row so nothing is silently
redefined.** `RMreliability()` now prints "Marginal (curve mean)".

What settled it was not the accuracy margin on its own but the coherence of the
table. PSI is `1 - mean(SEM^2)/Var(theta_hat)`, and since
`Var(theta_hat) ~ sigma^2 + mean(SEM^2)`, PSI is already the ratio form
computed from the observed spread of the WLE estimates. The Marginal row was
the other form. So the gap between the two rows, which the help page told users
to read as evidence of an off-target or non-normal sample, was substantially an
artefact of the differing formulas:

| data | items | PSI | ratio form | Green | PSI - Green |
|---|---|---|---|---|---|
| phq9 | 9 | .839 | .878 | .862 | -.023 |
| raschdat1 | 20 | .725 | .766 | .695 | +.029 |
| phq9[, 1:5] | 5 | .661 | .814 | .771 | -.111 |
| simulated, theta ~ N(0, .45^2) | 6 | .452 | .506 | .025 | **+.427** |

The last row is data generated to be well targeted and normal. A user following
the old help text would have read that .427 gap as a serious targeting problem.
With both rows on the same coefficient, the remaining gap is what the help text
always claimed it was: one coefficient computed from the observed spread, the
other from the fitted latent density.

**Why not option 1 (leave it).** The documented interpretation of the
PSI-to-marginal gap was unsafe, and the simulation above shows the estimator
itself is poor and unbounded.

**Why not option 2 (add a row).** A ratio-form row beside PSI would largely
duplicate it, and shipping a coefficient that has just been shown to leave
(0, 1) invites people to report it.

**On the naming.** "Marginal reliability" names a specific formula in the
literature and in mirt and TAM, so changing the formula under the same label
would have been an interoperability trap. The row label carries "(curve mean)"
so the difference is visible in the output itself, not only in the docs.

**What was implemented.**

- `.marginal_rxx()` is the density-weighted mean of the ratio curve. The floor
  at 0 is gone, being unnecessary for a bounded quantity.
- The `RMreliability()` row is "Marginal (curve mean)". Details rewritten and
  the pre-1.2.0 behaviour named.
- `RMreliabilityCurve()`'s `marginal_ratio` now agrees with it; the superseded
  subtractive value stays as `marginal_green`, labelled "Green/Lord,
  superseded" in the summary table, for users comparing against easyRasch2
  <= 1.2.0 or against other software.
- NEWS carries a breaking-change entry saying results move and by how much.
- Regression tests in both test files pin the formula and the boundedness.

## References

Green, B. F., Bock, R. D., Humphreys, L. G., Linn, R. L., & Reckase, M. D.
(1984). Technical Guidelines for Assessing Computerized Adaptive Tests.
*Journal of Educational Measurement, 21*(4), 347-360.

Milanzi, E., Molenberghs, G., Alonso, A., Verbeke, G., & De Boeck, P. (2015).
Reliability measures in item response theory: Manifest versus latent
correlation functions. *British Journal of Mathematical and Statistical
Psychology, 68*(1), 43-64. doi:10.1111/bmsp.12033
