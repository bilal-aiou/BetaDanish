# Bayesian Estimation for the Beta-Danish Distribution

Samples from the posterior of the Beta-Danish or Exponentiated Danish
parameters using a random-walk Metropolis sampler with vague
\\\Gamma(0.01, 0.01)\\ priors on the positive parameters.

## Usage

``` r
bayes_betadanish(
  time,
  status = NULL,
  submodel = TRUE,
  burnin = 5000,
  mcmc = 15000,
  tune = 0.5,
  theta_init = NULL,
  seed = NULL,
  verbose = 0
)
```

## Arguments

- time:

  Numeric vector of observed times.

- status:

  Numeric vector of event indicators (1 = event, 0 = right-censored).

- submodel:

  Logical; `TRUE` for the 3-parameter ED submodel, `FALSE` for the
  4-parameter full model.

- burnin:

  Burn-in iterations.

- mcmc:

  Post-burnin iterations.

- tune:

  Random-walk tuning parameter.

- theta_init:

  Optional starting values on the log scale.

- seed:

  Optional integer seed.

- verbose:

  Integer; passed to MCMCmetrop1R (0 = silent).

## Value

An object of class `"bd_bayes"` with components `draws` (mcmc object),
`summary`, `HPD`, `submodel`, `call`.

## Details

Requires MCMCpack and coda (Suggests).

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
dat <- rbetadanish(100, a = 1.5, b = 2, c = 3, k = 0.5)
fit <- bayes_betadanish(time = dat, submodel = TRUE,
                        burnin = 500, mcmc = 1500)
fit$summary
} # }
```
