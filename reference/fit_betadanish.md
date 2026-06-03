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
#>   Estimate Std. Error Lower 95% Upper 95% z value  Pr(>|z|)    
#> a  1.53758    3.33535  -4.99971   8.07487  0.4610 0.6448023    
#> b  1.94070    0.50087   0.95900   2.92240  3.8747 0.0001068 ***
#> c  2.59131    5.76955  -8.71700  13.89962  0.4491 0.6533337    
#> k  0.36910    0.31223  -0.24287   0.98108  1.1822 0.2371461    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
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
#>   Estimate Std. Error Lower 95% Upper 95% z value Pr(>|z|)    
#> b  1.94466    0.49477   0.97491   2.91442  3.9304 8.48e-05 ***
#> c  4.04255    1.89669   0.32504   7.76005  2.1314  0.03306 *  
#> k  0.40928    0.29166  -0.16237   0.98093  1.4033  0.16053    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Log-Likelihood: -269.7296 
#> AIC: 545.4592  | BIC: 553.2747 
# }
```
