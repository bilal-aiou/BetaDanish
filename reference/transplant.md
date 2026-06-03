# Bone Marrow Transplant Survival

Survival times for 91 patients with refractory acute lymphoblastic
leukemia who received either an allogeneic or autologous bone marrow
transplant. This dataset includes right-censoring and a treatment
covariate, making it ideal for demonstrating cure-rate models and AFT
regression.

## Usage

``` r
transplant
```

## Format

A data frame with 91 rows and 3 columns:

- time:

  Survival time in days

- status:

  Event indicator (1 = death/relapse, 0 = censored)

- group:

  Treatment group (0 = Allogeneic, 1 = Autologous)

## Source

Klein, J. P., & Moeschberger, M. L. (2003). Survival Analysis:
Techniques for Censored and Truncated Data (2nd ed.). Springer.

## Examples

``` r
data(transplant)
# \donttest{
# Fit a model with a covariate
fit <- fit_bd_aft(survival::Surv(time, status) ~ group, data = transplant)
# }
```
