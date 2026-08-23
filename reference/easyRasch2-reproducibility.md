# Random number generation and reproducibility in easyRasch2

The functions that simulate, bootstrap or resample take a `seed`
argument, and most also take `parallel`. This topic describes what
`seed` guarantees, what the default `seed = NULL` inherits from the
calling session, what a call does to that session's random number
generator, and the one situation in which the package overrides a
deliberate choice of generator.

## How seeding works

A single call proceeds in three steps.

1.  `set.seed(seed)` is applied once in the calling session.

2.  A vector of per-iteration seeds is drawn from that stream, one per
    Monte Carlo iteration, bootstrap replicate or respondent.

3.  Each iteration seeds itself from its own value before drawing
    anything.

Because each iteration carries its own seed, iterations are independent
of the order in which they finish. This is what makes `parallel = TRUE`
and `parallel = FALSE` return **identical** results for the same `seed`,
rather than merely equivalent ones, and what makes a result independent
of `n_cores`.

## The default `seed = NULL` inherits the session's stream

`seed` defaults to `NULL` everywhere, and a `NULL` seed sets nothing.
Step 2 above then draws the per-iteration seeds from whatever stream the
session is already on, so a
[`set.seed()`](https://rdrr.io/r/base/Random.html) earlier in the script
reproduces the call just as passing `seed` does. One
[`set.seed()`](https://rdrr.io/r/base/Random.html) at the top of an
analysis therefore covers every function listed below, and there is no
need to give each call a seed of its own.

The two routes differ in what they are robust to. An explicit `seed`
pins one call whatever runs before it. A session-level
[`set.seed()`](https://rdrr.io/r/base/Random.html) pins the script as a
sequence, so inserting, removing or reordering an earlier call that
draws random numbers changes every result after it. Use an explicit
`seed` where a single result has to be reproducible on its own, such as
a published cutoff.

## The generator is pinned inside iterations

Parallel iterations run in `mirai` daemons, and a daemon starts under
`RNGkind("L'Ecuyer-CMRG")` while an ordinary R session uses the
`"Mersenne-Twister"` default. Seeding alone is therefore not enough: the
same integer seed produces a different stream in a daemon than in the
calling session. Each iteration consequently pins the generator as well
as the seed, with the equivalent of

    set.seed(seed, kind = "Mersenne-Twister",
             normal.kind = "Inversion", sample.kind = "Rejection")

which are R's defaults since 3.6.0.

Two consequences worth knowing:

- If you have deliberately selected a non-default generator with
  [`RNGkind()`](https://rdrr.io/r/base/Random.html), the simulation
  inside these functions still uses the defaults above. Your choice is
  not honoured there.

- After such a call the session is left on the default generator. Wrap
  the call in
  [`withr::with_preserve_seed()`](https://withr.r-lib.org/reference/with_seed.html)
  or save and restore [`RNGkind()`](https://rdrr.io/r/base/Random.html)
  yourself if you need your setting back.

This is a deliberate trade. Pinning the generator is what allows the
parallel and sequential paths to agree exactly, and a non-default
generator is an unusual thing to want inside a bootstrap whose purpose
is to be reproducible.

## The session's random state is advanced

These functions consume random numbers from the calling session and do
not restore `.Random.seed` afterwards, so code that runs after a call
will not see the stream it would have seen without one. Use
[`withr::with_preserve_seed()`](https://withr.r-lib.org/reference/with_seed.html)
if that matters.

## Functions this applies to

[`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/reference/RMdifGammaCutoff.md),
[`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/reference/RMdimCFACutoff.md),
[`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/reference/RMdimMartinLof.md),
[`RMdimResidualPCACutoff()`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCACutoff.md),
[`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitCutoff.md),
[`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitCutoffMI.md),
[`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/reference/RMitemRestscoreBoot.md),
[`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/reference/RMlocdepGammaCutoff.md),
[`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3cutoff.md),
[`RMpersonFit()`](https://pgmj.github.io/easyRasch2/reference/RMpersonFit.md)
and
[`RMreliability()`](https://pgmj.github.io/easyRasch2/reference/RMreliability.md).

Functions that involve no simulation take no `seed` and are unaffected.

[`RMUreliability()`](https://pgmj.github.io/easyRasch2/reference/RMUreliability.md)
is the one resampling function without a `seed`. It splits the draw
columns at random, and reproduces from the session's stream in the way
described above.
[`RMreliability()`](https://pgmj.github.io/easyRasch2/reference/RMreliability.md)
calls it unseeded on purpose, so that each of the `rmu_iter` repetitions
it averages uses a different split.

## A note on other sources of non-determinism

A reproducible `seed` does not make every downstream number identical
across machines. Model fitting can converge to marginally different
optima under different BLAS implementations, and `mirt`-based routines
advance their own internal state. Where a function needs to defend
against that it re-seeds explicitly; see for example the plausible-value
block in
[`RMreliability()`](https://pgmj.github.io/easyRasch2/reference/RMreliability.md).

## See also

[easyRasch2-renaming](https://pgmj.github.io/easyRasch2/reference/easyRasch2-renaming.md)
