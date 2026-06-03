# Density of the r-th Order Statistic

Evaluates the probability density function of the r-th order statistic
from a sample of size n drawn from the Beta-Danish distribution.

## Usage

``` r
bd_order_stat_pdf(x, r, n, a, b, c, k, log = FALSE)
```

## Arguments

- x:

  Numeric vector of time points.

- r:

  Integer order (1 = minimum, n = maximum).

- n:

  Integer sample size.

- a, b, c, k:

  Positive parameters of the Beta-Danish distribution.

- log:

  Logical; if `TRUE` return the log-density.

## Value

Numeric vector (or its log).

## Examples

``` r
tgrid <- seq(0.01, 5, length.out = 50)
bd_order_stat_pdf(tgrid, r = 5, n = 20,
                  a = 1.5, b = 2.5, c = 2, k = 1)
#>  [1] 8.944348e-21 8.310447e-07 1.270589e-03 4.737029e-02 3.491252e-01
#>  [6] 1.036766e+00 1.737659e+00 1.984234e+00 1.733529e+00 1.250355e+00
#> [11] 7.844586e-01 4.442210e-01 2.332267e-01 1.158070e-01 5.520174e-02
#> [16] 2.554856e-02 1.158160e-02 5.177243e-03 2.294292e-03 1.012081e-03
#> [21] 4.458734e-04 1.966756e-04 8.703806e-05 3.870603e-05 1.731811e-05
#> [26] 7.803650e-06 3.544031e-06 1.623116e-06 7.499654e-07 3.497109e-07
#> [31] 1.646072e-07 7.822078e-08 3.752872e-08 1.817968e-08 8.891680e-09
#> [36] 4.390730e-09 2.188835e-09 1.101467e-09 5.594529e-10 2.867703e-10
#> [41] 1.483284e-10 7.740530e-11 4.074819e-11 2.163564e-11 1.158474e-11
#> [46] 6.254444e-12 3.404139e-12 1.867551e-12 1.032560e-12 5.752628e-13
```
