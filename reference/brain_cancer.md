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
summary(fit)
#> 
#> Call:
#> fit_bd_aft(formula = survival::Surv(Survtime, Survstatus) ~ Age + 
#>     Grade + Surgery, data = brain_cancer, n_starts = 2)
#> 
#> Beta-Danish AFT Model
#> 
#>                    Estimate Std. Error z value  Pr(>|z|)    
#> log_b              1.324748   0.425678  3.1121 0.0018577 ** 
#> log_c              0.485262   0.128174  3.7860 0.0001531 ***
#> delta_(Intercept) -4.964901   0.617402 -8.0416 8.867e-16 ***
#> delta_Age          0.515316   0.083177  6.1955 5.812e-10 ***
#> delta_Grade        0.621136   0.082252  7.5516 4.298e-14 ***
#> delta_Surgery     -1.403160   0.241036 -5.8214 5.837e-09 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Log-Likelihood: -937.1852 
# }
```
