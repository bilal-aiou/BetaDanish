# Fit Beta-Danish Cure Models

Fits mixture and promotion-time (non-mixture) cure models using the
Beta-Danish AFT baseline.

## Usage

``` r
fit_bd_cure(
  formula_aft,
  formula_cure,
  data,
  type = c("mixture", "promotion"),
  n_starts = 10,
  method = "BFGS"
)
```

## Arguments

- formula_aft:

  A formula for the latency component (e.g., \`Surv(time, status) ~
  age\`).

- formula_cure:

  A one-sided formula for the incidence/cure component (e.g., \`~
  treatment\`).

- data:

  A data frame containing the variables.

- type:

  Character; either \`"mixture"\` or \`"promotion"\`.

- n_starts:

  Integer; number of random starts for optimization.

- method:

  Optimization method passed to \`maxLik\`.

## Value

An object of class \`bd_cure\`.

## Details

In the \*\*mixture\*\* model, the population is split into susceptible
and cured fractions. The susceptible probability is modeled via logistic
regression: \`pi = exp(Z

In the \*\*promotion-time\*\* (non-mixture) model, the cure fraction is
derived from a latent Poisson process of clonogenic cells: \`theta =
exp(Z The cure fraction is \`exp(-theta)\`.
