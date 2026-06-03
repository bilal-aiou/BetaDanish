# Package index

## Distribution functions

Core d/p/q/r/h functions for the Beta-Danish distribution

- [`dbetadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/BetaDanish.md)
  [`pbetadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/BetaDanish.md)
  [`qbetadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/BetaDanish.md)
  [`rbetadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/BetaDanish.md)
  [`hbetadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/BetaDanish.md)
  : The Beta-Danish Distribution

## Maximum likelihood estimation

- [`fit_betadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/fit_betadanish.md)
  : Fit the Beta-Danish Distribution to Survival Data
- [`gof_betadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/gof_betadanish.md)
  : Goodness-of-Fit Statistics for Beta-Danish Models
- [`compare_models()`](https://bilal-aiou.github.io/BetaDanish/reference/compare_models.md)
  : Compare Nested Beta-Danish Models
- [`compare_distributions()`](https://bilal-aiou.github.io/BetaDanish/reference/compare_distributions.md)
  : Compare Beta-Danish with Standard Distributions

## Bayesian estimation

- [`bayes_betadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/bayes_betadanish.md)
  : Bayesian Estimation for the Beta-Danish Distribution

## AFT regression

- [`fit_bd_aft()`](https://bilal-aiou.github.io/BetaDanish/reference/fit_bd_aft.md)
  : Fit Beta-Danish AFT Regression Model

## Cure models

- [`fit_bd_cure()`](https://bilal-aiou.github.io/BetaDanish/reference/fit_bd_cure.md)
  : Fit Beta-Danish Cure Models
- [`simulate_bd_cure_data()`](https://bilal-aiou.github.io/BetaDanish/reference/simulate_bd_cure_data.md)
  : Simulate Beta-Danish Cure Data

## Competing risks

- [`fit_bd_competing()`](https://bilal-aiou.github.io/BetaDanish/reference/fit_bd_competing.md)
  : Fit Beta-Danish Competing Risks Model
- [`cif_betadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/cif_betadanish.md)
  : Compute Cumulative Incidence Function (CIF)
- [`cif_compare()`](https://bilal-aiou.github.io/BetaDanish/reference/cif_compare.md)
  : Compare Model-Based CIF to the Aalen-Johansen Estimator

## Structural properties

- [`bd_entropy_shannon()`](https://bilal-aiou.github.io/BetaDanish/reference/bd_entropy_shannon.md)
  : Shannon Entropy of the Beta-Danish Distribution
- [`bd_order_stat_pdf()`](https://bilal-aiou.github.io/BetaDanish/reference/bd_order_stat_pdf.md)
  : Density of the r-th Order Statistic

## Diagnostics and plots

- [`plot(`*`<betadanish>`*`)`](https://bilal-aiou.github.io/BetaDanish/reference/plot.betadanish.md)
  : Plot Diagnostics for Beta-Danish Fit
- [`plot(`*`<bd_aft>`*`)`](https://bilal-aiou.github.io/BetaDanish/reference/plot.bd_aft.md)
  [`plot(`*`<bd_cure>`*`)`](https://bilal-aiou.github.io/BetaDanish/reference/plot.bd_aft.md)
  : Cox-Snell Residual Plot for AFT and Cure Fits

## S3 methods

- [`summary(`*`<betadanish>`*`)`](https://bilal-aiou.github.io/BetaDanish/reference/summary.betadanish.md)
  : Summary Method for Beta-Danish Fit
- [`print(`*`<betadanish>`*`)`](https://bilal-aiou.github.io/BetaDanish/reference/print.betadanish.md)
  : Print Method for Beta-Danish Fit
- [`print(`*`<summary.betadanish>`*`)`](https://bilal-aiou.github.io/BetaDanish/reference/print.summary.betadanish.md)
  : Print Summary Method for Beta-Danish Fit
- [`coef(`*`<betadanish>`*`)`](https://bilal-aiou.github.io/BetaDanish/reference/coef.betadanish.md)
  : Extract Coefficients
- [`vcov(`*`<betadanish>`*`)`](https://bilal-aiou.github.io/BetaDanish/reference/vcov.betadanish.md)
  : Extract Variance-Covariance Matrix
- [`logLik(`*`<betadanish>`*`)`](https://bilal-aiou.github.io/BetaDanish/reference/logLik.betadanish.md)
  : Extract Log-Likelihood

## Reports and pipelines

- [`analyze_betadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/analyze_betadanish.md)
  : Comprehensive Beta-Danish Analysis Pipeline
- [`report_betadanish()`](https://bilal-aiou.github.io/BetaDanish/reference/report_betadanish.md)
  : Create a compact report from a BetaDanish model fit
- [`read_survival_data()`](https://bilal-aiou.github.io/BetaDanish/reference/read_survival_data.md)
  : Read and Prepare Survival Data
- [`simulate_bd_data()`](https://bilal-aiou.github.io/BetaDanish/reference/simulate_bd_data.md)
  : Simulate Data from the Beta-Danish Distribution

## Datasets

- [`remission`](https://bilal-aiou.github.io/BetaDanish/reference/remission.md)
  : Bladder Cancer Remission Times
- [`carbon_fibres`](https://bilal-aiou.github.io/BetaDanish/reference/carbon_fibres.md)
  : Breaking Stress of Carbon Fibres
- [`transplant`](https://bilal-aiou.github.io/BetaDanish/reference/transplant.md)
  : Bone Marrow Transplant Survival
- [`aarset`](https://bilal-aiou.github.io/BetaDanish/reference/aarset.md)
  : Aarset Device Failure Times
- [`leukemia`](https://bilal-aiou.github.io/BetaDanish/reference/leukemia.md)
  : Acute Myelogenous Leukemia Survival
- [`melanoma`](https://bilal-aiou.github.io/BetaDanish/reference/melanoma.md)
  : Malignant Melanoma Survival After Surgery
- [`brain_cancer`](https://bilal-aiou.github.io/BetaDanish/reference/brain_cancer.md)
  : Brain Cancer Survival Data
