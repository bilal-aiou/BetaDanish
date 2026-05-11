# Malignant Melanoma Survival After Surgery

Survival times for 205 patients with malignant melanoma after surgery.
This rich clinical dataset includes multiple covariates and heavy
censoring.

## Usage

``` r
melanoma
```

## Format

A data frame with 205 rows and 6 columns:

- time:

  Survival time in days

- status:

  Event indicator (1 = died from melanoma, 0 = alive, 2 = other death)

- thickness:

  Tumor thickness in mm

- sex:

  Patient sex (1 = male, 0 = female)

- age:

  Patient age in years

- ulcer:

  Ulceration indicator (1 = present, 0 = absent)

- year:

  Year of operation

## Source

Andersen, P. K., Borgan, O., Gill, R. D., & Keiding, N. (1993).
Statistical Models Based on Counting Processes. Springer.

## Examples

``` r
data(melanoma)
# \donttest{
# Treat status 1 as event, others as censored
melanoma$event <- ifelse(melanoma$status == 1, 1, 0)
fit <- fit_betadanish(survival::Surv(time, event) ~ age + thickness, data = melanoma)
#> Warning: bgrat(a=21.0971, b=0.000245906, x=1, y=9.88131e-323): z=2.03555e-321, b*z == 0 underflow, hence inaccurate pbeta()
#> Warning: bgrat(a=21.0971, b=0.000245906, x=1, y=9.88131e-323): z=2.03555e-321, b*z == 0 underflow, hence inaccurate pbeta()
# }
```
