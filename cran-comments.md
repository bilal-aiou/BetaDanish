## Test environments

* local Ubuntu 24.04, R 4.4.x
* GitHub Actions: ubuntu-latest (R-release, R-devel), macOS-latest (R-release), windows-latest (R-release)
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0-1 notes

All notes (if any) relate to the resubmission status or new-maintainer
declaration, which we have addressed.

## Reverse dependencies

There are no reverse dependencies for this package.

## Notes for CRAN

This is the second release of the BetaDanish package. The 0.2.0 release
contains a substantial methodological expansion (Bayesian estimation,
bound-constrained competing-risks estimation with Aalen-Johansen
comparators and Gray's test, moments-with-existence checks, Shannon
entropy, order statistics, mean residual life, stress-strength
reliability, hazard-shape classifier, bootstrap confidence intervals,
and a simulation-study runner) along with bug fixes to the AFT and cure
summary methods (delta-method back-transform).

All new heavyweight dependencies (MCMCpack, coda, cmprsk, MASS,
flexsurv) are declared in Suggests with requireNamespace() guards at
the call sites.
