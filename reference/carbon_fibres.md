# Breaking Stress of Carbon Fibres

Breaking stress (in Gba) of 100 carbon fibre specimens. This dataset
exhibits a unimodal (increasing-then-decreasing) hazard pattern that
classical distributions like the Weibull cannot adequately capture.

## Usage

``` r
carbon_fibres
```

## Format

A data frame with 100 rows and 2 columns:

- time:

  Breaking stress in Gba

- status:

  Event indicator (1 = event occurred)

## Source

Nichols, M. D., & Padgett, W. J. (2006). A bootstrap control chart for
Weibull percentiles. Quality and Reliability Engineering International,
22(2), 141-151.

## Examples

``` r
data(carbon_fibres)
# \donttest{
fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = carbon_fibres)
#> Warning: Parameters a, b, c, and k must be strictly positive.
#> Warning: Parameters a, b, c, and k must be strictly positive.
#> Warning: pbeta(*, log.p=TRUE) -> bpser(a=374115, b=13.5121, x=0.940558,...) underflow to -Inf
#> Warning: pbeta(*, log.p=TRUE) -> bpser(a=374115, b=13.5121, x=0.973071,...) underflow to -Inf
#> Warning: Parameters a, b, c, and k must be strictly positive.
# }
```
