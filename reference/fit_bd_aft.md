# Fit Beta-Danish AFT Regression Model

Fits an Accelerated Failure Time (AFT) regression model using the
Complementary Exponentiated Danish (CED) baseline (Beta-Danish with
a=1).

## Usage

``` r
fit_bd_aft(formula, data, n_starts = 10, method = "BFGS")
```

## Arguments

- formula:

  A survival formula (e.g., \`Surv(time, status) ~ age + treatment\`).

- data:

  A data frame containing the variables.

- n_starts:

  Integer; number of random starts for optimization.

- method:

  Optimization method passed to \`maxLik\`.

## Value

An object of class \`bd_aft\`.

## Details

To ensure identifiability, the shape parameter \`a\` is fixed to 1. The
scale parameter \`k\` is linked to covariates via \`k_i = exp(X_i
Positive coefficients in \`delta\` indicate a larger \`k\`, which
corresponds to shorter survival times (accelerated failure).
