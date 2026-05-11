# The Beta-Danish Distribution

Density, distribution function, quantile function, hazard function, and
random generation for the four-parameter Beta-Danish distribution.

## Usage

``` r
dbetadanish(x, a, b, c, k, log = FALSE)

pbetadanish(q, a, b, c, k, lower.tail = TRUE, log.p = FALSE)

qbetadanish(p, a, b, c, k, lower.tail = TRUE, log.p = FALSE)

rbetadanish(n, a, b, c, k)

hbetadanish(x, a, b, c, k, log = FALSE)
```

## Arguments

- x, q:

  Vector of quantiles (time points).

- a:

  Shape parameter (beta generator). Set \`a = 1\` for the 3-parameter
  submodel.

- b:

  Shape parameter (beta generator / tail weight).

- c:

  Shape parameter (baseline shape).

- k:

  Scale parameter (baseline scale).

- log, log.p:

  Logical; if TRUE, probabilities/densities are given as log.

- lower.tail:

  Logical; if TRUE (default), probabilities are P\[X \<= x\], otherwise
  P\[X \> x\].

- p:

  Vector of probabilities.

- n:

  Number of observations to generate.

## Value

\`dbetadanish\` gives the density, \`pbetadanish\` gives the
distribution function, \`qbetadanish\` gives the quantile function,
\`hbetadanish\` gives the hazard function, and \`rbetadanish\` generates
random deviates.

## Details

The Beta-Danish distribution is a highly flexible lifetime distribution
capable of modeling decreasing, increasing, unimodal, and bathtub-shaped
hazard rates.

## References

Ahmad, B., & Danish, M. Y. (2026). Development and Characterization of a
Flexible Three-Parameter Lifetime Distribution.

## Examples

``` r
# Density
dbetadanish(x = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#> [1] 0.1087591

# CDF
pbetadanish(q = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#> [1] 0.102199

# Hazard
hbetadanish(x = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#> [1] 0.1211394

# Random generation
rbetadanish(n = 10, a = 1.5, b = 2, c = 3, k = 0.5)
#>  [1]  1.942771 11.545763 11.700365 80.185448 42.730097  4.566257  5.361773
#>  [8]  3.848609  2.630686  6.275763
```
