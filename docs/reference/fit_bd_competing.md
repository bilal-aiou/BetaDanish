# Fit Beta-Danish Competing Risks Model

Fits a parametric competing risks model assuming independent latent
failure times, where each cause-specific baseline follows a 4-parameter
Beta-Danish distribution.

## Usage

``` r
fit_bd_competing(time, cause, n_starts = 5, method = "BFGS")
```

## Arguments

- time:

  Numeric vector of observed times.

- cause:

  Integer vector of event causes. \`0\` indicates right-censored, and
  \`1, 2, ..., m\` indicate specific event causes.

- n_starts:

  Integer; number of random starts for the joint optimization.

- method:

  Optimization method passed to \`maxLik\`.

## Value

An object of class \`bd_competing\`.

## Details

Under the assumption of independent latent failure times, the joint
likelihood factorizes. The function first fits independent Beta-Danish
models for each cause (treating other causes as censored) to find robust
starting values, then optimizes the joint likelihood.
