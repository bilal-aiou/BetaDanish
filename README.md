# BetaDanish: The Beta-Danish Distribution for Lifetime Data Analysis

<!-- badges: start -->
[![R-CMD-check](https://github.com/bilal-aiou/BetaDanish/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bilal-aiou/BetaDanish/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.r-project.org/Licenses/GPL-3)
<!-- badges: end -->

The **BetaDanish** R package provides a comprehensive suite of tools for survival and reliability analysis using the four-parameter Beta-Danish distribution and its three-parameter Exponentiated Danish (ED) submodel.

## Features (v0.2.0)

* **Core distribution**: numerically stable d/p/q/r/h/s/logS for the Beta-Danish and the Danish baseline.
* **MLE inference**: complete and right-censored data with delta-method standard errors and multi-start optimization.
* **Bayesian inference**: `bayes_betadanish()` random-walk Metropolis sampler (requires `MCMCpack`).
* **AFT regression** and **mixture / promotion-time cure models**.
* **Competing risks**: bound-constrained multi-start L-BFGS-B with Aalen-Johansen overlay and Gray's test (`cif_compare()`).
* **Structural properties**: moments with existence checks, probability-weighted moments, mean residual life, Shannon entropy, order statistics, stress-strength reliability, and a hazard-shape classifier.
* **Diagnostics**: survival, hazard, density, P-P, Q-Q, and Cox-Snell residual plots.
* **Tooling**: bootstrap confidence intervals and a finite-sample simulation-study runner.

## Installation

```r
# install.packages("devtools")
devtools::install_github("bilal-aiou/BetaDanish")
install.packages(c("MCMCpack", "coda", "cmprsk", "flexsurv", "MASS"))
```

## Quick Start

```r
library(BetaDanish)
data("remission")
fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = remission)
summary(fit)
plot(fit, type = "all")
```

## Citation

Ahmad, B., & Danish, M. Y. (2025). The Beta-Danish distribution for lifetime data analysis. *Journal of Applied Mathematics, Statistics and Informatics*, 21(1). <https://doi.org/10.2478/jamsi-2025-0010>

## License

GPL-3
