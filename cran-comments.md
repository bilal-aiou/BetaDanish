## Test environments

* local: Windows 11 x64, R 4.5.2
* win-builder: R-devel and R-release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the reviewer

This is a feature release of an existing package (0.2.0 -> 0.3.0). Two changes are user-visible and worth flagging.

**A dataset has been removed.** The `brain_cancer` dataset that shipped in
0.1.0 and 0.2.0 is no longer included. This was done at the request of the
maintainer's doctoral supervisor. No functionality depends on it; the
affected example now uses the `melanoma` dataset. A new dataset,
`guinea_pig` (Bjerkedal 1960), has been added in its place.

**One signature has changed.** `bd_entropy_shannon()` previously computed the
entropy by quadrature with arguments `(a, b, c, k, subdivisions, rel.tol)`.
It now uses a closed-form expression with arguments
`(a, b, c, k, terms, method, rel.tol, subdivisions)`. Named calls are
unaffected, and the previous behaviour remains available as
`method = "quadrature"`. The change is recorded in NEWS.md.

The remaining changes are additive: structural properties (moments,
entropies, mean residual life, stress-strength reliability, order
statistics), penalized and grouped-likelihood estimation, profile and Wald
intervals, and a file-driven analysis entry point.

## Reverse dependencies

There are no reverse dependencies on CRAN.
