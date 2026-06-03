# BetaDanish <a href="https://bilal-aiou.github.io/BetaDanish/"><img src="man/figures/logo.png" align="right" height="139" alt="BetaDanish website" /></a>

<!-- badges: start -->
[![R-CMD-check](https://github.com/bilal-aiou/BetaDanish/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bilal-aiou/BetaDanish/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bilal-aiou/BetaDanish/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/bilal-aiou/BetaDanish/actions/workflows/pkgdown.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.r-project.org/Licenses/GPL-3)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

## Overview

The **BetaDanish** package implements the four-parameter Beta-Danish
distribution and its three-parameter Exponentiated Danish (ED)
submodel for survival, reliability, and lifetime-data analysis. The
distribution was introduced by Ahmad and Danish (2025) and offers a
flexible alternative to classical lifetime models such as the
Weibull, gamma, and log-normal, while accommodating monotonic,
unimodal, and bathtub-shaped hazards within a single parametric
family.

Beyond the core distribution, the package provides a comprehensive
toolkit for modern survival modelling:

- Maximum-likelihood and **Bayesian** inference for complete and
  right-censored data
- **Accelerated failure time (AFT)** regression
- Mixture and **promotion-time cure** models
- **Competing risks** analysis with Aalen-Johansen comparisons and
  Gray's test
- Closed-form moments, **Shannon entropy**, order statistics,
  mean residual life, hazard-shape classification, and
  stress-strength reliability
- A complete diagnostic toolkit: survival, hazard, density, P-P, Q-Q
  and **Cox-Snell residual** plots

## Installation

You can install the development version of BetaDanish from
[GitHub](https://github.com/bilal-aiou/BetaDanish) with:

```r
# install.packages("devtools")
devtools::install_github("bilal-aiou/BetaDanish")
```

Optional packages used by advanced modules (declared in `Suggests`):

```r
install.packages(c("MCMCpack", "coda", "cmprsk", "flexsurv", "MASS"))
```

## The Beta-Danish Distribution

### Density function

For $t > 0$, the Beta-Danish pdf is

$$f(t \mid a, b, c, k) = \frac{1}{B(a,b)}\, k\, c\, t^{c-1} \left[1 - \exp(-k\, t^{c})\right]^{ac - 1} \exp(-k\, b\, t^{c}),$$

with positive parameters $a, b, c, k > 0$. The case $a = 1$ yields
the three-parameter Exponentiated Danish (ED) submodel used
throughout the case studies.

### Survival, hazard, and CDF

The package provides numerically stable d/p/q/r/h functions analogous
to the base R distribution interface:

| Function | Purpose |
|---|---|
| `dbetadanish()` | Density (PDF) |
| `pbetadanish()` | Cumulative distribution (CDF) and survival |
| `qbetadanish()` | Quantile function |
| `rbetadanish()` | Random number generation |
| `hbetadanish()` | Hazard rate |

## Quick Start

```r
library(BetaDanish)

# 1. Evaluate the density, CDF, hazard
dbetadanish(x = 2, a = 1.5, b = 2, c = 3, k = 0.5)
pbetadanish(q = 2, a = 1.5, b = 2, c = 3, k = 0.5)
hbetadanish(x = 2, a = 1.5, b = 2, c = 3, k = 0.5)

# 2. Simulate data and fit a model
set.seed(123)
dat <- simulate_bd_data(n = 200, a = 1.5, b = 2, c = 3, k = 0.5,
                        censor_rate = 0.2)
fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat)
summary(fit)

# 3. Diagnostic plots
plot(fit, type = "all")

# 4. Compare against standard distributions
compare_distributions(fit)
```

## Featured Examples

### Bayesian estimation

```r
data("remission")
fit_bayes <- bayes_betadanish(
  time     = remission$time,
  status   = remission$status,
  submodel = TRUE,
  burnin   = 2000,
  mcmc     = 5000,
  seed     = 1
)
print(fit_bayes)
```

See the **Bayesian Estimation** vignette for full details.

### AFT regression

```r
data("brain_cancer")
fit_aft <- fit_bd_aft(
  survival::Surv(Survtime, Survstatus) ~ Age + Grade + Surgery,
  data = brain_cancer
)
summary(fit_aft)
plot(fit_aft)   # Cox-Snell residual diagnostic
```

### Cure models

Both **mixture** and **promotion-time** cure models on the
Exponentiated Danish kernel:

```r
fit_cure <- fit_bd_cure(
  formula_aft  = survival::Surv(time, status) ~ 1,
  formula_cure = ~ group,
  data         = mydata,
  type         = "mixture"
)
summary(fit_cure)
```

### Competing risks

```r
fit_cr <- fit_bd_competing(time = time, cause = cause)
res <- cif_compare(fit_cr, plot = TRUE)
res$gray_test
```

The fitted Beta-Danish cumulative incidence functions are overlaid
against the nonparametric Aalen-Johansen estimator, with Gray's test
reported for cause equality.

## Structural Properties

| Function | Property |
|---|---|
| `bd_entropy_shannon()` | Shannon (differential) entropy |
| `bd_order_stat_pdf()` | r-th order statistic density |

## Built-in Datasets

| Dataset | n | Description |
|---|---|---|
| `remission` | 128 | Bladder cancer remission times |
| `carbon_fibres` | 100 | Breaking stress of carbon fibres (Gba) |
| `transplant` | 91 | Bone marrow transplant survival |
| `aarset` | 50 | Aarset device failure times (bathtub hazard) |
| `leukemia` | 23 | Acute myelogenous leukemia survival |
| `melanoma` | 205 | Malignant melanoma post-surgery |
| `brain_cancer` | varies | Brain cancer survival with comorbidities |

## Vignettes

The package ships with five comprehensive tutorials:

1. **Introduction to BetaDanish** — Overview, parameterization,
   first-time-user walkthrough
2. **Brain Cancer Case Study** — Complete real-world analysis
3. **Bayesian Estimation** — Posterior inference with
   `bayes_betadanish()`
4. **Competing Risks** — Joint cause-specific modelling with CIF
   comparison
5. **Cure Models** — Mixture and promotion-time formulations

Browse them with:

```r
browseVignettes("BetaDanish")
```

## Citation

If you use BetaDanish in published work, please cite:

> Ahmad, B., & Danish, M. Y. (2025). The Beta-Danish distribution
> for lifetime data analysis. *Journal of Applied Mathematics,
> Statistics and Informatics*, 21(1).
> <https://doi.org/10.2478/jamsi-2025-0010>

A BibTeX entry is available via:

```r
citation("BetaDanish")
```

## Getting Help

- **Bug reports and feature requests:** please open an issue at
  <https://github.com/bilal-aiou/BetaDanish/issues>
- **Function reference:** see the Reference tab on this site, or
  type `?fit_betadanish` in R
- **Vignettes:** see the Articles tab on this site, or
  `vignette(package = "BetaDanish")`

## Authors

**Bilal Ahmad** *(maintainer)* — Allama Iqbal Open University,
Islamabad, Pakistan. <bilalahmad.imcbh9@gmail.com>

**Muhammad Yameen Danish** — Allama Iqbal Open University,
Islamabad, Pakistan. <yameen.danish@aiou.edu.pk>

## License

GPL-3
