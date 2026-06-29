# Contributing to easyRasch2

Contributions, bug reports, and questions are all welcome.

## Reporting issues

Please report bugs and unexpected behaviour on the
[issue tracker](https://github.com/pgmj/easyRasch2/issues). A good report
includes:

- a minimal reproducible example (e.g. using a built-in dataset such as the
  bundled `phq9` data or `eRm::pcmdat2`),
- the output of `sessionInfo()`, and
- the installed version (`packageVersion("easyRasch2")`).

Searching existing issues first helps avoid duplicates.

## Contributing code

1. Open an issue describing the change before starting substantial work, so
   scope and approach can be discussed.
2. Fork the repository and create a feature branch.
3. Follow the existing code style and the package's `RM`-prefixed naming
   conventions.
4. Add or update tests under `tests/testthat/` and run `devtools::test()` and
   `devtools::check()` until clean.
5. Update roxygen2 documentation and, where relevant, `NEWS.md`.
6. Open a pull request referencing the issue.

## Seeking support

For usage questions, open a
[GitHub issue](https://github.com/pgmj/easyRasch2/issues) or contact the
maintainer at <pgmj@pm.me>. Broader Rasch-analysis tutorials are linked from
the package [website](https://pgmj.github.io/easyRasch2/).

## Code of conduct

Please keep all project interactions respectful and constructive.
