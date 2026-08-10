#' Random number generation and reproducibility in easyRasch2
#'
#' Every function in the package that simulates, bootstraps or resamples
#' takes a `seed` argument, and most also take `parallel`. This topic
#' describes what `seed` guarantees, what it does to the calling session's
#' random number generator, and the one situation in which the package
#' overrides a deliberate choice of generator.
#'
#' @section How seeding works:
#' A single call proceeds in three steps.
#'
#' \enumerate{
#'   \item `set.seed(seed)` is applied once in the calling session.
#'   \item A vector of per-iteration seeds is drawn from that stream, one
#'     per Monte Carlo iteration, bootstrap replicate or respondent.
#'   \item Each iteration seeds itself from its own value before drawing
#'     anything.
#' }
#'
#' Because each iteration carries its own seed, iterations are independent
#' of the order in which they finish. This is what makes
#' `parallel = TRUE` and `parallel = FALSE` return **identical** results
#' for the same `seed`, rather than merely equivalent ones, and what makes
#' a result independent of `n_cores`.
#'
#' @section The generator is pinned inside iterations:
#' Parallel iterations run in `mirai` daemons, and a daemon starts under
#' `RNGkind("L'Ecuyer-CMRG")` while an ordinary R session uses the
#' `"Mersenne-Twister"` default. Seeding alone is therefore not enough: the
#' same integer seed produces a different stream in a daemon than in the
#' calling session. Each iteration consequently pins the generator as well
#' as the seed, with the equivalent of
#'
#' ```
#' set.seed(seed, kind = "Mersenne-Twister",
#'          normal.kind = "Inversion", sample.kind = "Rejection")
#' ```
#'
#' which are R's defaults since 3.6.0.
#'
#' Two consequences worth knowing:
#'
#' \itemize{
#'   \item If you have deliberately selected a non-default generator with
#'     [RNGkind()], the simulation inside these functions still uses the
#'     defaults above. Your choice is not honoured there.
#'   \item After such a call the session is left on the default generator.
#'     Wrap the call in [withr::with_preserve_seed()] or save and restore
#'     [RNGkind()] yourself if you need your setting back.
#' }
#'
#' This is a deliberate trade. Pinning the generator is what allows the
#' parallel and sequential paths to agree exactly, and a non-default
#' generator is an unusual thing to want inside a bootstrap whose purpose
#' is to be reproducible.
#'
#' @section The session's random state is advanced:
#' These functions consume random numbers from the calling session and do
#' not restore `.Random.seed` afterwards, so code that runs after a call
#' will not see the stream it would have seen without one. Use
#' [withr::with_preserve_seed()] if that matters.
#'
#' @section Functions this applies to:
#' [RMdifGammaCutoff()], [RMdimCFACutoff()], [RMdimMartinLof()],
#' [RMdimResidualPCACutoff()], [RMitemInfitCutoff()],
#' [RMitemInfitCutoffMI()], [RMitemRestscoreBoot()],
#' [RMlocdepGammaCutoff()], [RMlocdepQ3Cutoff()], [RMpersonFit()] and
#' [RMreliability()].
#'
#' Functions that involve no simulation take no `seed` and are unaffected.
#'
#' @section A note on other sources of non-determinism:
#' A reproducible `seed` does not make every downstream number identical
#' across machines. Model fitting can converge to marginally different
#' optima under different BLAS implementations, and `mirt`-based routines
#' advance their own internal state. Where a function needs to defend
#' against that it re-seeds explicitly; see for example the plausible-value
#' block in `RMreliability()`.
#'
#' @seealso [easyRasch2-renaming]
#'
#' @name easyRasch2-reproducibility
#' @aliases easyRasch2-reproducibility
#' @keywords internal
NULL
