# Compare Model-Based CIF to the Aalen-Johansen Estimator

Computes the nonparametric Aalen-Johansen CIF (via cmprsk) for each
competing-risks cause, overlays it on the fitted Beta-Danish CIF, and
returns the Aalen-Johansen times/estimates, the fitted CIF values on a
common grid, and Gray's CIF-equality test where applicable.

## Usage

``` r
cif_compare(fit, tmax = NULL, n_grid = 160, plot = TRUE)
```

## Arguments

- fit:

  A fitted object of class `"bd_competing"`.

- tmax:

  Optional upper time for evaluation; default the 95th percentile of
  observed times.

- n_grid:

  Number of time points on the evaluation grid (default 160).

- plot:

  Logical; if `TRUE` (default) draws a panel of overlays.

## Value

A list with elements `tgrid`, `cif_fit` (data frame long format),
`cif_aj` (data frame long format) and optionally `gray_test`.

## Details

Requires cmprsk (Suggests).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
T1 <- rbetadanish(200, 1.2, 1.5, 1.0, 0.4)
T2 <- rbetadanish(200, 1.0, 2.0, 1.0, 0.2)
C  <- stats::rexp(200, 0.05)
time  <- pmin(T1, T2, C)
cause <- ifelse(time == C, 0L, ifelse(T1 <= T2, 1L, 2L))
fit <- fit_bd_competing(time = time, cause = cause)
cif_compare(fit)   # requires cmprsk to be installed
} # }
```
