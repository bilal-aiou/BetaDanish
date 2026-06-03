# BetaDanish: The Beta-Danish Distribution for Lifetime Data Analysis

<!-- badges: start -->
[![R-CMD-check](https://github.com/bilal-aiou/BetaDanish/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bilal-aiou/BetaDanish/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.r-project.org/Licenses/GPL-3)
<!-- badges: end -->

The **BetaDanish** R package provides a comprehensive suite of tools for survival and reliability analysis using the four-parameter Beta-Danish distribution and its three-parameter Exponentiated Danish (ED) submodel.

Developed at the Department of Statistics, Allama Iqbal Open University (AIOU), Islamabad, the package addresses the limitations of classical lifetime models by accommodating decreasing, increasing, unimodal, and bathtub-shaped hazard rates in a single unified family.

## Why use BetaDanish?

| Model | Typical hazard shape | Flexibility |
|---|---|---|
| Exponential | Constant | Low |
| Weibull | Increasing or decreasing | Moderate |
| Gamma | Flexible but limited | Moderate |
| Log-normal | Non-monotone | Moderate |
| Log-logistic | Non-monotone | Moderate |
| Beta-Danish | Increasing, decreasing, unimodal, bathtub | High |

## Features (v0.2.0)

* **Core distribution**: numerically stable d/p/q/r/h/s/logS for the Beta-Danish and the Danish baseline.
* **MLE inference**: complete and right-censored data with delta-method standard errors, multi-start optimization, and AIC/BIC stored at fit time.
* **Bayesian inference**: `bayes_betadanish()` random-walk Metropolis sampler (requires `MCMCpack`).
* **AFT regression** and **mixture / promotion-time cure models**.
* **Competing risks**: bound-constrained multi-start L-BFGS-B with formula interface, Aalen-Johansen overlay, and Gray's test (`cif_compare()`).
* **Structural properties**: moments with existence checks, probability-weighted moments, mean residual life, Shannon entropy, order statistics, stress-strength reliability, and a Glaser-type hazard-shape classifier.
* **Diagnostics**: survival, hazard, density, P-P, Q-Q, and Cox-Snell residual plots for `betadanish`, `bd_aft`, and `bd_cure` fits.
* **Goodness-of-fit**: KS, Cramer-von Mises, Anderson-Darling, plus a 7-distribution `flexsurv` comparator.
* **Tooling**: bootstrap confidence intervals (`bd_bootstrap_ci`) and a finite-sample simulation-study runner (`bd_mle_study`).

## Installation

```r
# Development version
# install.packages("devtools")
devtools::install_github("bilal-aiou/BetaDanish")

# Optional Suggests for full functionality:
install.packages(c("MCMCpack", "coda", "cmprsk", "flexsurv", "MASS"))
```

## Quick Start

```r
library(BetaDanish)
data("remission")

# Fit the four-parameter Beta-Danish model
fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = remission)
summary(fit)
plot(fit, type = "all")        # 6 diagnostic panels including Cox-Snell

# Test the submodel hypothesis a = 1 against the full model
fit_sub <- fit_betadanish(survival::Surv(time, status) ~ 1,
                          data = remission, submodel = TRUE)
compare_models(fit, fit_sub)

# Benchmark against seven standard distributions
compare_distributions(fit)
```

## Bayesian estimation

```r
fit_bayes <- bayes_betadanish(
  time = remission$time, status = remission$status,
  submodel = TRUE, burnin = 2000, mcmc = 5000, seed = 1)
fit_bayes$summary
```

## Cure models

```r
data("transplant")
fit_cure <- fit_bd_cure(
  formula_aft  = survival::Surv(time, status) ~ 1,
  formula_cure = ~ group,
  data         = transplant,
  type         = "mixture")
summary(fit_cure)
plot(fit_cure)                 # Cox-Snell residuals
```

## Competing risks

```r
set.seed(1)
T1 <- rbetadanish(300, 1.2, 1.5, 1.0, 0.4)
T2 <- rbetadanish(300, 1.0, 2.0, 1.0, 0.2)
C  <- stats::rexp(300, 0.05)
time <- pmin(T1, T2, C)
cause <- ifelse(time == C, 0L, ifelse(T1 <= T2, 1L, 2L))
fit_cr <- fit_bd_competing(time = time, cause = cause, n_restarts = 5)
cif_compare(fit_cr)            # Aalen-Johansen overlay + Gray's test
```

## Theoretical properties

```r
bd_moments(a = 1.5, b = 2.5, c = 2, k = 1)
bd_entropy_shannon(1.5, 2.5, 2, 1)
bd_mrl(c(0.5, 1, 2), 1.5, 2.5, 2, 1)
bd_stress_strength(c(1.5, 2.5, 2, 1), c(1.5, 2.5, 2, 0.5))
bd_hazard_shape(a = 0.5, b = 4, c = 1.2, k = 1)$shape   # "bathtub"
```

## Citation

Ahmad, B., & Danish, M. Y. (2025). The Beta-Danish distribution for lifetime data analysis. *Journal of Applied Mathematics, Statistics and Informatics*, 21(1). <https://doi.org/10.2478/jamsi-2025-0010>

```r
citation("BetaDanish")
```

## License

GPL-3
