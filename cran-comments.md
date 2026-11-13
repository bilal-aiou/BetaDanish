## R CMD check results

0 errors | 0 warnings | 0 notes

* Local: Windows, R release, `devtools::check(args = "--as-cran")`.
* win-builder: R-devel (2026-11-12 r90881 ucrt), 1 NOTE (see below).

The NOTE lists possibly misspelled words in DESCRIPTION. All are spelled
correctly: 'Allama' is part of the university name in the cited thesis
(Allama Iqbal Open University), 'Cramer', 'von' and 'Mises' name the
Cramer-von Mises statistic, and 'Kullback' and 'Leibler' name the
Kullback-Leibler divergence.

## This submission

This is a major update, 0.3.0 -> 1.0.0. It brings the estimation,
goodness-of-fit, interval and Bayesian routines into line with the final
analysis code of the doctoral thesis on which the package is based. All
changes, including the behaviour changes, are listed in NEWS.md.

The package now cites a single reference, the doctoral thesis
(Ahmad, 2026). The thesis has no DOI, so it is cited by author, year and
title.

## Reverse dependencies

There are currently no reverse dependencies on CRAN.
