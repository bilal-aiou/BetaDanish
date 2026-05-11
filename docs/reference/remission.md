# Bladder Cancer Remission Times

Remission times (in months) for 128 bladder cancer patients. This is a
complete (uncensored) sample widely used in lifetime distribution
literature to demonstrate decreasing or right-skewed hazard rates.

## Usage

``` r
remission
```

## Format

A data frame with 128 rows and 2 columns:

- time:

  Remission time in months

- status:

  Event indicator (1 = event occurred)

## Source

Lee, E. T., & Wang, J. W. (2003). Statistical Methods for Survival Data
Analysis (3rd ed.). Wiley.

## Examples

``` r
data(remission)
# \donttest{
fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = remission)
#> Warning: Parameters a, b, c, and k must be strictly positive.
#> Warning: Parameters a, b, c, and k must be strictly positive.
#> Warning: Parameters a, b, c, and k must be strictly positive.
#> Warning: Parameters a, b, c, and k must be strictly positive.
#> Warning: Parameters a, b, c, and k must be strictly positive.
summary(fit)
#> 
#> Call:
#> fit_betadanish(formula = survival::Surv(time, status) ~ 1, data = remission)
#> 
#> Beta-Danish Distribution Fit
#> Model: Full 4-Parameter Model 
#> 
#>    Estimate Std. Error Lower 95% Upper 95% z value Pr(>|z|)   
#> a  0.686572   0.790045 -0.861916  2.235059  0.8690 0.384831   
#> b  4.078300   1.530117  1.079271  7.077328  2.6654 0.007691 **
#> c  2.196511   2.561668 -2.824359  7.217381  0.8575 0.391194   
#> k  0.082971   0.077796 -0.069508  0.235451  1.0665 0.286186   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-Likelihood: -409.9137 
#> AIC: 827.8274  | BIC: 839.2356 
# }
```
