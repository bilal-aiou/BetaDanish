# Simulate Data from the Beta-Danish Distribution

Generates synthetic survival data from the Beta-Danish distribution,
with optional right-censoring.

## Usage

``` r
simulate_bd_data(n, a, b, c, k, censor_rate = 0, seed = NULL)
```

## Arguments

- n:

  Integer; number of observations to simulate.

- a, b, c, k:

  Numeric; parameters of the Beta-Danish distribution.

- censor_rate:

  Numeric; rate parameter for the exponential censoring distribution. If
  \`0\` (default), no censoring is applied.

- seed:

  Integer; optional seed for reproducibility.

## Value

A data frame with columns \`time\` and \`status\`.

## Examples

``` r
# Simulate complete data
dat <- simulate_bd_data(n = 100, a = 1.5, b = 2, c = 3, k = 0.5)

# Simulate censored data
dat_cens <- simulate_bd_data(n = 100, a = 1.5, b = 2, c = 3, k = 0.5, censor_rate = 0.1)
```
