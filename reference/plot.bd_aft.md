# Cox-Snell Residual Plot for AFT and Cure Fits

Diagnostic Cox-Snell residual plot for a fitted AFT or cure model.

## Usage

``` r
# S3 method for class 'bd_aft'
plot(x, ...)

# S3 method for class 'bd_cure'
plot(x, ...)
```

## Arguments

- x:

  A fitted `"bd_aft"` or `"bd_cure"` object.

- ...:

  Further graphical parameters.

## Value

Invisibly returns `x`.

## Examples

``` r
# \donttest{
set.seed(42)
n <- 300
x <- stats::rnorm(n)
k <- exp(-0.5 - 0.3 * x)
t_sim <- rbetadanish(n, a = 1, b = 2, c = 1.5, k = k)
status <- stats::rbinom(n, 1, 0.85)  # ~15% censoring
dat <- data.frame(time = t_sim, status = status, x = x)
fit <- fit_bd_aft(survival::Surv(time, status) ~ x, data = dat,
                  n_starts = 5)
plot(fit)
#> Warning: Cox-Snell plot: fitted coefficients contain NA/NaN; cannot compute residuals.
# }
```
