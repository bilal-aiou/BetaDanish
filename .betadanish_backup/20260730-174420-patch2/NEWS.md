# BetaDanish 0.2.0.9000 (development version)

## Correctness fixes

* `plot.bd_aft()` and `plot.bd_cure()` now produce the Cox-Snell residual
  plot. Both looked up the shape parameters as `coefficients["b"]` and
  `coefficients["c"]`, but `fit_bd_aft()` and `fit_bd_cure()` store them as
  `log_b` and `log_c`. The lookup returned `NA`, the internal guard caught it,
  and the functions returned without drawing anything. The lookup is fixed and
  the values are exponentiated back to the natural scale.

* `dbetadanish()` is now accurate in the far right tail. The term
  \eqn{\log\{1 - G(t)\}} was formed as `log1p(-exp(log_G))`, which pins at
  about `log(1e-16)` once `G` rounds to one, flooring the log-density near
  -36.8 regardless of its true value. It is now formed as
  `log(-expm1(-c * log1p(1/(k*t))))`, which holds full relative precision for
  every `k*t > 0`.

* `pbetadanish()` computes the survival function through the beta mirror
  identity \eqn{1 - I_y(a,b) = I_{1-y}(b,a)}, so a probability near one is
  never subtracted from one.

* `qbetadanish()` obtains \eqn{1 - u} directly from
  \eqn{1 - q\beta(p; a, b) = q\beta(p; b, a)} with the tail flag reversed. The
  previous route computed `y^(1/c) / (k * (1 - y^(1/c)))`, which loses all
  significant digits as `u` approaches one. Round-trip accuracy now holds to
  `p = 1 - 1e-13`.

* `hbetadanish()` no longer substitutes `-700` for an exhausted log-survival.
  A positive density with zero survival gives `Inf`, an honest divergent
  hazard, and the indeterminate case gives `NaN`.

* The density, distribution, quantile, survival and hazard functions now
  recycle their parameters element-wise against `x`, `q` or `p`, following the
  usual convention for R distribution functions.

* `fit_betadanish()` records `npar`, `nobs`, `nevent`, `AIC` and `BIC` on the
  fitted object, and `report_betadanish()` computes information criteria
  through the `logLik` method. Both previously omitted AIC and BIC silently.

* `summary.bd_aft()` and `summary.bd_cure()` report shape parameters on the
  natural scale with delta-method standard errors and an exponentiated
  log-scale confidence interval, separately from the regression coefficients.
  No Wald test is reported for a shape parameter, since the implied null lies
  outside the parameter space.

* `read_survival_data()` selects columns by name rather than by position, so a
  covariate named `time` or `status` no longer collides with the response.
  Covariates keep their original type. The returned data frame carries a
  `bd_data_report` attribute.

* `extract_surv_data()` validates times and the event indicator up front and
  reports a clear error rather than failing inside the optimiser.

## New

* `sbetadanish()`, an explicit survival function.

* `fit_betadanish()` gains `check_identifiability`. It warns when the optimiser
  reports a poor code, when the information matrix is singular, when
  \eqn{\hat b} is within two standard errors of the \eqn{b = 1}
  non-identifiability ridge, and when the fitted correlation between
  \eqn{\hat a} and \eqn{\hat c} exceeds 0.95 in absolute value. The
  Identifiability section of `?fit_betadanish` explains each case.

## Documentation

* `?fit_betadanish` documents the \eqn{b = 1} ridge and the lower-tail
  \eqn{(a, c)} confounding.

* `?fit_bd_competing` documents the Tsiatis (1975) non-identifiability result
  and the direction of bias under latent dependence.

* `?BetaDanish` records that \eqn{S(t) \propto t^{-b}} in the upper tail and
  hence that \eqn{E(X^r)} is finite if and only if \eqn{b > r}.

* The `a = 1` submodel is called the Exponentiated Danish (ED) throughout, in
  line with the underlying thesis. "Complementary Exponentiated Danish (CED)"
  has been removed.

* Carbon fibre breaking stress is given in GPa. It was previously written
  "Gba", and the typo was whitelisted in `inst/WORDLIST`, which is why the
  spell check never caught it.

* The `melanoma` help page reported six columns and documented seven.

## Tests

* Regression tests pin the tail behaviour of the distribution functions. The
  survival function is regularly varying with index \eqn{-b}, so the slope of
  \eqn{\log S} against \eqn{\log t} must approach \eqn{-b}. The previous
  implementation cannot pass this test, which is why the tail defects went
  unnoticed.

* Smoke tests added for `fit_bd_aft()`, `fit_bd_cure()`, `plot.bd_aft()` and
  `plot.bd_cure()`, none of which were previously covered.

# BetaDanish 0.2.0

> **Changelog correction.** The 0.2.0 entry below has been rewritten to
> describe only what that release actually contained. As first published it
> also listed mean residual life, hazard-shape classification, stress-strength
> reliability, bootstrap confidence intervals for AFT and cure models, and a
> finite-sample simulation-study runner, none of which were implemented, and it
> reported four bug fixes that had not been applied. Those features are being
> added in the 0.3.0 development series and the fixes are recorded above.

## Major new functionality

* **Bayesian inference**: `bayes_betadanish()` provides random-walk Metropolis
  sampling for the Exponentiated Danish submodel and the full four-parameter
  Beta-Danish model with vague Gamma priors.
* **Competing risks rewrite**: `fit_bd_competing()` uses bound-constrained
  multi-start L-BFGS-B optimisation. `cif_compare()` overlays fitted cumulative
  incidence functions against the Aalen-Johansen estimator and reports Gray's
  test.
* **Structural properties**: Shannon entropy (`bd_entropy_shannon()`, by
  adaptive quadrature) and order-statistic densities (`bd_order_stat_pdf()`).
* **Diagnostics**: Cox-Snell residual plot methods for AFT and cure fits. These
  were shipped but non-functional; see the fix above.

## Vignettes

Three new vignettes were added:

* "Bayesian Estimation with BetaDanish"
* "Competing Risks with the Beta-Danish Distribution"
* "Cure Models with the Beta-Danish Distribution"

## Infrastructure

* Continuous integration via GitHub Actions on four OS/R configurations.
* Test coverage reporting via Codecov.
* Online package website built with pkgdown.
* All `Suggests` packages guarded with `requireNamespace()` at the call sites.

# BetaDanish 0.1.0

* First public release.
* Implements the four-parameter Beta-Danish distribution and its
  three-parameter Exponentiated Danish submodel for survival and reliability
  analysis.
* Maximum-likelihood estimation, goodness-of-fit, model comparison, and
  visualization.
* Built-in datasets: remission, carbon_fibres, transplant, aarset, leukemia,
  melanoma, brain_cancer.
