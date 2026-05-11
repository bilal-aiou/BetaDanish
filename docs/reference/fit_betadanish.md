# Fit the Beta-Danish Distribution to Survival Data

Fits the Beta-Danish distribution using Maximum Likelihood Estimation
(MLE). Supports both complete and right-censored data via
\`survival::Surv\` objects.

## Usage

``` r
fit_betadanish(formula, data, submodel = FALSE, n_starts = 10, method = "BFGS")
```

## Arguments

- formula:

  A formula object, with the response on the left of a \`~\` operator,
  and the terms on the right. The response must be a survival object as
  returned by the \`Surv\` function. Use \`~ 1\` for models without
  covariates.

- data:

  A data frame containing the variables named in the formula.

- submodel:

  Logical; if \`TRUE\`, fits the 3-parameter submodel by fixing \`a =
  1\`.

- n_starts:

  Integer; the number of random starting points to use for the
  optimization to ensure global convergence. Default is 10.

- method:

  Character; the optimization method passed to \`maxLik\`. Default is
  "BFGS".

## Value

An object of S3 class \`"betadanish"\`, containing the parameter
estimates, log-likelihood, variance-covariance matrix, and convergence
diagnostics.

## Details

The optimization is performed on the log-transformed parameters to
strictly enforce positivity constraints. The returned coefficients and
variance-covariance matrix are transformed back to the natural scale
using the Delta method.

## Examples

``` r
# \donttest{
# Simulate some data
set.seed(123)
sim_time <- rbetadanish(100, a = 1.5, b = 2, c = 3, k = 0.5)
sim_status <- sample(c(0, 1), 100, replace = TRUE, prob = c(0.2, 0.8))
dat <- data.frame(time = sim_time, status = sim_status)

# Fit the 4-parameter model
fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat)
#> Warning: Parameters a, b, c, and k must be strictly positive.
summary(fit)
#> 
#> Call:
#> fit_betadanish(formula = survival::Surv(time, status) ~ 1, data = dat)
#> 
#> Beta-Danish Distribution Fit
#> Model: Full 4-Parameter Model 
#> 
#>    Estimate Std. Error Lower 95% Upper 95% z value  Pr(>|z|)    
#> a   1.53758    3.90797  -6.12205   9.19721  0.3934 0.6939895    
#> b   1.94070    0.49901   0.96263   2.91876  3.8891 0.0001006 ***
#> c   2.59131    6.77172 -10.68126  15.86388  0.3827 0.7019671    
#> k   0.36910    0.33163  -0.28089   1.01910  1.1130 0.2657089    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-Likelihood: -269.7122 
#> AIC: 547.4243  | BIC: 557.845 

# Fit the 3-parameter submodel
fit_sub <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat, submodel = TRUE)
summary(fit_sub)
#> 
#> Call:
#> fit_betadanish(formula = survival::Surv(time, status) ~ 1, data = dat, 
#>     submodel = TRUE)
#> 
#> Beta-Danish Distribution Fit
#> Model: 3-Parameter Submodel (a=1) 
#> 
#>   Estimate Std. Error Lower 95% Upper 95% z value  Pr(>|z|)    
#> b  1.94466    0.49437   0.97570   2.91362  3.9336 8.367e-05 ***
#> c  4.04255    1.86684   0.38353   7.70156  2.1654   0.03035 *  
#> k  0.40928    0.28778  -0.15477   0.97333  1.4222   0.15497    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-Likelihood: -269.7296 
#> AIC: 545.4592  | BIC: 553.2747 
# }
```
