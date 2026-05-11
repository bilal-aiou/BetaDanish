# Brain Cancer Survival Data

A comprehensive dataset of 500 brain cancer patients, including survival
times, censoring status, and multiple clinical covariates. This dataset
was used to demonstrate Accelerated Failure Time (AFT) regression and
Cure-Rate models using the Beta-Danish distribution.

## Usage

``` r
brain_cancer
```

## Format

A data frame with 500 rows and 16 columns:

- ID:

  Patient identifier

- Gender:

  Patient gender (1 = Male, 0 = Female)

- Age:

  Age group (1 = Young, 2 = Middle, 3 = Old)

- Area:

  Geographic area (1 = Urban, 0 = Rural)

- FH:

  Family history of cancer (1 = Yes, 0 = No)

- CMH:

  Comorbid history (1 = Yes, 0 = No)

- Grade:

  Tumor grade (1 = I/II, 2 = III, 3 = IV)

- Surgery:

  Surgical intervention (1 = Yes, 0 = No)

- Radiotherapy:

  Radiotherapy treatment (1 = Yes, 0 = No)

- Chemotherapy:

  Chemotherapy treatment (1 = Yes, 0 = No)

- Treatment:

  Treatment type

- Morphology:

  Tumor morphology

- Survstatus:

  Survival status (1 = Event/Death, 0 = Censored)

- Survtime:

  Survival time in months

- Types:

  Tumor types classification

- Morphology1:

  Alternative morphology classification

## Source

Atomic Energy Cancer Hospital (NORI), Islamabad, Pakistan.

## Examples

``` r
data(brain_cancer)
# \donttest{
# Fit an AFT model using the brain cancer data
fit <- fit_bd_aft(survival::Surv(Survtime, Survstatus) ~ Age + Grade + Surgery,
                  data = brain_cancer, n_starts = 2)
#> Warning: Parameters a, b, c, and k must be strictly positive.
#> Warning: Parameters a, b, c, and k must be strictly positive.
summary(fit)
#> 
#> Call:
#> fit_bd_aft(formula = survival::Surv(Survtime, Survstatus) ~ Age + 
#>     Grade + Surgery, data = brain_cancer, n_starts = 2)
#> 
#> Beta-Danish AFT Model
#> 
#>                    Estimate Std. Error z value  Pr(>|z|)    
#> log_b              1.324576   0.441683  2.9989 0.0027093 ** 
#> log_c              0.485017   0.129954  3.7322 0.0001898 ***
#> delta_(Intercept) -4.965820   0.633670 -7.8366 4.629e-15 ***
#> delta_Age          0.515220   0.083159  6.1956 5.806e-10 ***
#> delta_Grade        0.621435   0.082143  7.5653 3.869e-14 ***
#> delta_Surgery     -1.403098   0.240828 -5.8261 5.672e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-Likelihood: -937.1851 
# }
```
