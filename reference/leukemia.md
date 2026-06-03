# Acute Myelogenous Leukemia Survival

Survival times (in weeks) for 23 patients with acute myelogenous
leukemia. A classic, small dataset perfect for fast testing of censored
data workflows.

## Usage

``` r
leukemia
```

## Format

A data frame with 23 rows and 3 columns:

- time:

  Survival time in weeks

- status:

  Event indicator (1 = event, 0 = censored)

- group:

  Treatment group (Maintained vs Non-maintained)

## Source

Miller, R. G. (1997). Survival Analysis. Wiley.

## Examples

``` r
data(leukemia)
# \donttest{
fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = leukemia)
# }
```
