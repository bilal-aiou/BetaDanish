---
title: 'BetaDanish: An R package for heavy-tailed survival analysis with the Beta-Danish distribution'
tags:
  - R
  - survival analysis
  - lifetime distributions
  - heavy tails
  - maximum likelihood
  - cure models
  - competing risks
authors:
  - name: Bilal Ahmad
    corresponding: true
    affiliation: 1
  - name: Muhammad Yameen Danish
    affiliation: 1
affiliations:
  - name: Department of Statistics, Allama Iqbal Open University, Islamabad, Pakistan
    index: 1
date: 13 November 2026
bibliography: paper.bib
---

# Summary

Survival and reliability studies record how long it takes until an event
happens: a patient relapses, a component breaks, a claim is settled. In many
such data sets most times are moderate but a few are very long, so the
distribution has a heavy upper tail. Common lifetime models such as the
Weibull, gamma and log-normal distributions have upper tails that decay
quickly, and they tend to understate the chance of very long survival.

The Beta-Danish distribution [@Ahmad2026] is a four-parameter lifetime model
for such data. It applies the beta-generated construction [@Eugene2002] to a
two-parameter "Danish" baseline, a special case of the Burr III and Dagum
distributions [@Burr1942; @Dagum1977] whose survival function decays like a
power of time. One parameter, $b$, sets how heavy the upper tail is, and
fixing $a = 1$ gives a three-parameter submodel, the Exponentiated Danish
distribution, whose distribution and quantile functions have closed forms.

`BetaDanish` is an R [@R] package that makes the distribution usable in
practice. It evaluates the density, distribution, quantile, survival and
hazard functions accurately far into both tails; fits the model to complete,
right-censored and coarsely recorded samples by maximum likelihood or by
Bayesian simulation; and says when the data cannot determine all four
parameters. Regression extensions cover accelerated failure time, cure and
competing-risks models. The package is available from CRAN
[@BetaDanishCRAN] and is developed openly on GitHub.

# Statement of need

A four-parameter lifetime family is useful only if its fits can be trusted,
and the Beta-Danish family is difficult in three ways that generic fitting
code does not handle. First, its parameters are not equally well determined.
On the stratum $b = 1$ they cannot all be identified, and in the lower tail
the shape parameters $a$ and $c$ act almost only through their product, so
the likelihood has long, nearly flat ridges. Second, for some samples the
likelihood increases without bound as the parameters move toward the edge of
the parameter space. A general-purpose optimizer then stops somewhere along
the way and reports a "maximum" whose standard errors and information
criteria have no meaning. Third, the probabilities that matter for
heavy-tailed data lie far in the tails, where the textbook formulas lose all
precision.

`BetaDanish` is written for statisticians and applied researchers in
biomedicine, reliability engineering and actuarial science who want to use
the family, or its submodel, on real data and to know when the fourth
parameter is warranted. It implements the estimation and inference
procedures developed in the first author's doctoral thesis [@Ahmad2026],
so that the analyses there can be reproduced and applied to new data. A
full-length description of the package with worked examples has been
submitted to the *Journal of Statistical Software* and is under review; this
paper summarizes the software and the reasons for its design.

# State of the field

Several R packages fit parametric lifetime models. `flexsurv` [@Jackson2016]
is the most general platform for parametric survival regression and accepts
user-supplied distributions. `flexsurvcure` [@Amdahl2025] and `smcure`
[@Cai2012] fit cure models. `fitdistrplus` [@DelignetteMuller2015] fits
univariate distributions to complete and censored data, `actuar`
[@Dutang2008] provides the Burr and Dagum families, `gamlss`
[@Stasinopoulos2007] models every distribution parameter as a function of
covariates, and `Newdistns` [@Nadarajah2016] evaluates beta-generated and
related families for a user-supplied baseline. None of them includes the
Beta-Danish family.

We built a dedicated package rather than contribute the family to one of
these frameworks because its difficulties lie in estimation and reporting,
not in the density alone. Supplying the density to `flexsurv` or `gamlss`
would give fits, but not a search that guarantees the four-parameter fit is
at least as good as the nested three-parameter fit, not the distinction
between maxima and suprema at the edge of the parameter space, and not the
diagnostics described below. Adding these would mean changing the estimation
core of a general package for the sake of one family. Instead, `BetaDanish`
follows the base R conventions for its density, distribution, quantile,
random-generation and hazard functions, so the distribution can still be
passed to `flexsurv` as a custom distribution, and it uses the formulas and
objects of `survival` [@Therneau2000]. It complements the existing ecosystem
rather than replacing it.

# Software design

The package is written in base R without compiled code. We traded speed for
portability and ease of maintenance: the models have at most four shape and
scale parameters plus regression coefficients, so the computing cost is
small, and accuracy comes from how quantities are computed rather than from
extended precision. The functions work on the log scale where possible,
evaluate upper-tail probabilities through the mirror identity of the
incomplete beta function instead of subtracting from one, and switch to the
leading term of a series expansion where the baseline underflows.

Estimation is built around nesting. The likelihood is maximized in two
stages: the three-parameter submodel is fitted from several starting points,
and the four-parameter model is then started from points along the $(a, c)$
ridge through the submodel optimum, one of which is the submodel optimum
itself. This costs more likelihood evaluations than a single start, but it
guarantees that the likelihood ratio statistic for $a = 1$ is never
negative. Shape parameters are confined to a wide box. A fit that ends at its
edge is reported as a supremum, without standard errors or information
criteria, and interior fits on the weakly identified ridge are flagged. We
chose to return fewer numbers rather than numbers that look valid but are
not.

Inference follows the same principle. Wald intervals are computed on the log
scale; profile likelihood intervals extend their grid until the profile
crosses the cut-off or reaches the box, and an endpoint that never closes is
reported as open. Estimates can be re-expressed in terms of the identified
product $ac$. Goodness-of-fit statistics (Kolmogorov-Smirnov, Cramer-von
Mises and Anderson-Darling) are computed from the Kaplan-Meier estimator
[@Kaplan1958] on the probability-integral-transform scale, so complete and
censored samples are judged on one scale. The Bayesian sampler is an
adaptive multi-chain Metropolis algorithm [@Haario2001; @Roberts1997],
written in base R to avoid a hard dependency; it reports R-hat and effective
sample sizes and declares convergence only when every parameter passes a
fixed threshold.

Model objects are S3 classes with `print`, `summary`, `coef`, `vcov`,
`logLik` and `plot` methods, so they behave like other R model fits. A
file-based workflow reads a delimited text file or spreadsheet, fits the
requested models and returns tidy tables, for users who prefer not to write
R code. The package has 22 `testthat` test files, is checked on Linux, macOS
and Windows by continuous integration, and is documented in help pages, six
vignettes and a `pkgdown` website.

# Research impact statement

`BetaDanish` grew out of the analysis code of the doctoral thesis that
introduces the Beta-Danish distribution [@Ahmad2026]. Version 1.0.0 brings
every estimation and inference routine into line with the thesis's final
analysis code, so the analyses reported there can be reproduced with the
package. The package is distributed through CRAN, where the current release is
1.0.0 and the package has a DOI, and its development history, with tagged
releases 0.2.0, 0.3.0 and 1.0.0, is public on GitHub. The full-length manuscript submitted to
the *Journal of Statistical Software* comes with a standalone replication
script that reproduces every result and figure from data sets shipped with
the package, including benchmark samples from medicine and materials
science. Together with the test suite and documentation, these materials
make the package's results checkable by anyone with R installed.

# AI usage disclosure

The authors used Claude (Anthropic), including Claude Opus 5.5, to generate
and refactor parts of the package code, to debug it, to draft tests and
documentation, and to draft and edit the text of this paper. The statistical
methods and design decisions are the authors' own, and the authors reviewed,
edited and validated all AI-assisted code and text, including through the
package's test suite and `R CMD check` on several platforms.

# Acknowledgements

The package was developed during the first author's doctoral research at the
Department of Statistics, Allama Iqbal Open University, Islamabad, Pakistan.

# References
