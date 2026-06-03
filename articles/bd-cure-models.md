# Cure Models with the Beta-Danish Distribution

## Two cure formulations

The package implements two cure formulations on the Exponentiated Danish
kernel:

- **Mixture cure** splits the population into a susceptible fraction
  `pi(Z)` modelled by logistic regression on the incidence covariates,
  and a cured fraction `1 - pi(Z)`.
- **Promotion-time cure** derives the cure fraction from a latent
  Poisson process of clonogenic cells with intensity
  `theta(Z) = exp(Z' gamma)`; the cure fraction is `exp(-theta(Z))`.

Both are fitted via
[`fit_bd_cure()`](https://bilal-aiou.github.io/BetaDanish/reference/fit_bd_cure.md).

## Quick example

Both fits below use simulated data with a known cure structure.

``` r

library(BetaDanish)
set.seed(2026)
n <- 250
z <- stats::rbinom(n, 1, 0.5)
pi_susc <- stats::plogis(0.3 + 0.7 * z)
cured   <- stats::rbinom(n, 1, 1 - pi_susc) == 1
T_true  <- ifelse(cured, Inf,
                  rbetadanish(n, a = 1, b = 2, c = 1.5, k = 0.4))
C   <- stats::rexp(n, 0.04)
time   <- pmin(T_true, C)
status <- ifelse(T_true <= C, 1, 0)
dat <- data.frame(time = time, status = status, z = z)
cat("Sample size:", n, " Censoring rate:", round(mean(status == 0), 2), "\n")
#> Sample size: 250  Censoring rate: 0.45
```

## Mixture cure model

``` r

fit_mix <- fit_bd_cure(
  formula_aft  = survival::Surv(time, status) ~ 1,
  formula_cure = ~ z,
  data         = dat,
  type         = "mixture",
  n_starts     = 3
)
summary(fit_mix)
#> 
#> Call:
#> fit_bd_cure(formula_aft = survival::Surv(time, status) ~ 1, formula_cure = ~z, 
#>     data = dat, type = "mixture", n_starts = 3)
#> 
#> Beta-Danish Cure Model (mixture)
#> 
#>                   Estimate Std. Error z value Pr(>|z|)  
#> log_b              1.26173    0.52936  2.3835  0.01715 *
#> log_c              0.43147    0.19780  2.1814  0.02916 *
#> delta_(Intercept) -1.39644    0.76201 -1.8326  0.06687 .
#> gamma_(Intercept)  0.22464    0.19410  1.1573  0.24713  
#> gamma_z            0.56245    0.29202  1.9260  0.05410 .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-Likelihood: -416.2507
```

## Promotion-time cure model

``` r

fit_prom <- fit_bd_cure(
  formula_aft  = survival::Surv(time, status) ~ 1,
  formula_cure = ~ z,
  data         = dat,
  type         = "promotion",
  n_starts     = 3
)
summary(fit_prom)
#> 
#> Call:
#> fit_bd_cure(formula_aft = survival::Surv(time, status) ~ 1, formula_cure = ~z, 
#>     data = dat, type = "promotion", n_starts = 3)
#> 
#> Beta-Danish Cure Model (promotion)
#> 
#>                   Estimate Std. Error z value Pr(>|z|)  
#> log_b              1.56319    0.79236  1.9728  0.04851 *
#> log_c              0.39322    0.17026  2.3095  0.02091 *
#> delta_(Intercept) -2.03331    0.94620 -2.1489  0.03164 *
#> gamma_(Intercept) -0.19425    0.12821 -1.5151  0.12975  
#> gamma_z            0.33526    0.17070  1.9640  0.04953 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-Likelihood: -416.1921
```

## See also

- [`?fit_bd_cure`](https://bilal-aiou.github.io/BetaDanish/reference/fit_bd_cure.md)
  for full documentation
- `?bd_bootstrap_ci` for bootstrap confidence intervals
- [`?plot.bd_cure`](https://bilal-aiou.github.io/BetaDanish/reference/plot.bd_aft.md)
  for Cox-Snell residual plots
