# Contributing to biomed

1. Create a focused branch from `main`.
2. Add or update roxygen2 documentation and tests with each behavior change.
3. Run `devtools::document()`, `devtools::test()`, and `devtools::check()`.
4. Confirm that no patient data, secrets, local paths, or large generated files
   are included in the diff.
5. Open a pull request describing the scientific purpose, assumptions, and
   validation performed.

Changes to clinical endpoint definitions or statistical defaults must include
tests and a clear rationale. Breaking changes should be documented in
`NEWS.md`.
