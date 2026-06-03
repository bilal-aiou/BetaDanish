# Simulate Beta-Danish Cure Data

Generates synthetic survival data from a Beta-Danish mixture or
promotion-time cure model, incorporating covariates.

## Usage

``` r
simulate_bd_cure_data(
  n,
  type = c("mixture", "promotion"),
  a = 1,
  b = 2,
  c = 1.5,
  delta,
  gamma,
  X,
  Z,
  target_censor = 0.3,
  seed = NULL
)
```

## Arguments

- n:

  Integer; number of observations.

- type:

  Character; \`"mixture"\` or \`"promotion"\`.

- a, b, c:

  Numeric; baseline shape parameters.

- delta:

  Numeric vector; coefficients for the latency scale \`k\`.

- gamma:

  Numeric vector; coefficients for the incidence/cure component.

- X:

  Matrix; design matrix for latency (must match length of \`delta\`).

- Z:

  Matrix; design matrix for incidence (must match length of \`gamma\`).

- target_censor:

  Numeric; target proportion of censoring to calibrate the exponential
  censoring rate. Default is 0.3.

- seed:

  Integer; optional seed.

## Value

A list containing the simulated \`data\` (time, status), the \`cured\`
indicator, and the true parameters.
