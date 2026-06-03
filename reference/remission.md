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
#> a  0.686568   0.828451 -0.937197  2.310333  0.8287  0.40725   
#> b  4.078195   1.497987  1.142139  7.014250  2.7224  0.00648 **
#> c  2.196543   2.685695 -3.067419  7.460505  0.8179  0.41343   
#> k  0.082975   0.079826 -0.073483  0.239433  1.0394  0.29860   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Log-Likelihood: -409.9137 
#> AIC: 827.8274  | BIC: 839.2356 
# }
```
