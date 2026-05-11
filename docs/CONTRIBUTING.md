# Contributing to BetaDanish

Contributions to BetaDanish are welcome.

## Reporting bugs

When reporting a bug, please include:

- a minimal reproducible example,
- your R version,
- your operating system,
- the BetaDanish version,
- the full error message or warning.

## Pull requests

Before submitting a pull request, please ensure that:

- [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  passes without errors or warnings,
- new functions include roxygen2 documentation,
- new features include tests,
- examples are small enough to run quickly.

## Development workflow

Recommended checks:

``` r

devtools::document()
devtools::test()
devtools::check()
```
