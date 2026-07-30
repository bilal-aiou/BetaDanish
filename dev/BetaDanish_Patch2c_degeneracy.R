## =============================================================================
##  BetaDanish  --  PATCH 2c : DEGENERATE-OPTIMUM GUARD
## =============================================================================
##
##  WHY
##    The Patch 2b test log contained this, from the four-parameter fit on the
##    transplant data:
##
##      Model                     LogLik        Chisq Df Pr(>Chisq)
##      Submodel (3-param)  -4.768013e+02           NA NA         NA
##      Full Model (4-param) 5.210644e+77 1.042129e+78  1          0
##
##    A log-likelihood of +5.2e77 on 91 observations is not a maximum. The
##    optimiser slid down a degenerate ridge and maxLik returned the result as
##    converged. optim_multistart() selects on `fit$maximum > best_ll` behind
##    nothing but an is.finite() test, so the runaway beat every honest optimum
##    and won.
##
##  WHAT THIS PATCH DOES
##    It ports the policy already established in Master_R_Code_for_Thesis.R
##    rather than inventing a new one. Constants and thresholds are the ones
##    used there, so the package and the thesis pipeline now agree:
##
##      .ch6_degenerate()    lim_shape = 500, hard reject      -> .bd_make_accept()
##      .flag_and_sanitize() big = 1e8, weak = 1e6, se_ratio = 50
##      default_starts()     the deterministic start grid      -> .bd_default_starts()
##      "best NON-DEGENERATE optimum over a structured multi-start"
##
##    One guard is new, because the master script does not have it and the
##    +5.2e77 case shows it is needed: a plausibility bound on the fitted
##    log-likelihood per observation.
##
##  SCOPE
##    fit_betadanish(), fit_bd_aft() and fit_bd_cure(). Competing risks is left
##    alone because recommendation 37 rewrites that fitter in Patch 3.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch2c_degeneracy.R")
##  IDEMPOTENT   Yes.
## =============================================================================

if (getRversion() < "4.0.0") stop("This patch needs R >= 4.0.")

.step_n <- 0L
.step <- function(m) { .step_n <<- .step_n + 1L; cat(sprintf("\n[%02d] %s\n", .step_n, m)) }
.ok   <- function(m) cat("     OK   ", m, "\n", sep = "")
.info <- function(m) cat("     ..   ", m, "\n", sep = "")
.warn <- function(m) cat("     WARN ", m, "\n", sep = "")
.die  <- function(...) stop("\n\n*** PATCH ABORTED ***\n", ..., "\n", call. = FALSE)

BACKUP_DIR <- NULL
.backup <- function(p) {
  if (!file.exists(p)) return(invisible(FALSE))
  d <- file.path(BACKUP_DIR, p)
  dir.create(dirname(d), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(p, d, overwrite = TRUE)) .die("Could not back up ", p)
  invisible(TRUE)
}
.put <- function(path, content) {
  .backup(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- strsplit(content, "\n", fixed = TRUE)[[1]]
  while (length(lines) && !nzchar(lines[length(lines)])) lines <- lines[-length(lines)]
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(lines, con = con, sep = "\n")
  .ok(paste("wrote", path))
  invisible(TRUE)
}
.write_lines <- function(path, lines) {
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(lines, con = con, sep = "\n")
}
.sub_in <- function(path, from, to, label) {
  if (!file.exists(path)) { .warn(paste(path, "missing")); return(invisible(FALSE)) }
  txt <- readLines(path, warn = FALSE)
  if (!any(grepl(from, txt, fixed = TRUE))) {
    if (any(grepl(to, txt, fixed = TRUE))) .info(paste(label, "-- already applied"))
    else .warn(paste(label, "-- pattern not found; check by hand"))
    return(invisible(FALSE))
  }
  .backup(path)
  .write_lines(path, gsub(from, to, txt, fixed = TRUE))
  .ok(label)
  invisible(TRUE)
}

cat(strrep("=", 78), "\n")
cat("  BetaDanish  --  Patch 2c : degenerate-optimum guard\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("R/analyze_csv.R")) .die("Patch 2b has not been applied.")
.ok("Patch 2b detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch2c"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  R/utils-internal.R  --  guarded multi-start
## =============================================================================

.step("Rewriting R/utils-internal.R (acceptance predicate, start grid)")

.put("R/utils-internal.R", r"---(#' Robust Multi-Start Optimization using maxLik
#'
#' @param ll_fun Log-likelihood function to maximise.
#' @param start_list A list of named numeric vectors of starting values.
#' @param method Optimisation method (default `"BFGS"`).
#' @param accept Optional predicate taking a fitted `maxLik` object and
#'   returning `TRUE` if the optimum is admissible. Optima for which it returns
#'   `FALSE` are discarded rather than competing for the maximum. `NULL`
#'   (default) accepts anything finite, which is the historical behaviour.
#'
#' @details
#' Selecting purely on the highest log-likelihood is unsafe for this family.
#' The Beta-Danish likelihood has degenerate ridges along which the objective
#' can be driven arbitrarily high without a finite maximiser existing, so an
#' unguarded search will prefer a runaway to every honest optimum. When
#' `accept` is supplied, the best *admissible* optimum is returned instead.
#'
#' The returned object carries three attributes for diagnostics:
#' `bd_starts_ok`, `bd_starts_rejected` and `bd_loglik_spread`, the last being
#' the range of the accepted log-likelihoods.
#'
#' @return The best admissible `maxLik` object, or `NULL` if none qualified.
#' @noRd
optim_multistart <- function(ll_fun, start_list, method = "BFGS", accept = NULL) {
  best_fit <- NULL
  best_ll  <- -Inf
  n_ok     <- 0L
  n_rej    <- 0L
  lls      <- numeric(0)

  for (start_vals in start_list) {
    fit <- tryCatch(
      maxLik::maxLik(logLik = ll_fun, start = start_vals, method = method),
      error = function(e) NULL)

    if (is.null(fit) || is.null(fit$maximum) ||
        is.na(fit$maximum) || !is.finite(fit$maximum)) next

    if (!is.null(accept) && !isTRUE(accept(fit))) {
      n_rej <- n_rej + 1L
      next
    }

    n_ok <- n_ok + 1L
    lls  <- c(lls, fit$maximum)
    if (fit$maximum > best_ll) {
      best_ll  <- fit$maximum
      best_fit <- fit
    }
  }

  if (!is.null(best_fit)) {
    attr(best_fit, "bd_starts_ok")       <- n_ok
    attr(best_fit, "bd_starts_rejected") <- n_rej
    attr(best_fit, "bd_loglik_spread")   <-
      if (length(lls) > 1L) diff(range(lls)) else 0
  }
  best_fit
}

#' Admissibility Predicate for a Beta-Danish Optimum
#'
#' Ports the policy of `.ch6_degenerate()` in the thesis master script. That
#' function documents the mechanism: the density contains
#' \eqn{w = kx/(1 + kx)} raised to the power \eqn{c}, and as \eqn{k \to \infty}
#' and \eqn{c \to \infty} with \eqn{c/k} fixed, \eqn{w \to 1} and
#' \eqn{w^{c} \to \exp\{-(c/k)/x\}}: the family collapses to a different,
#' Frechet-type limiting distribution. That limit is a genuine ridge in the
#' likelihood, and an optimiser started in the wrong place slides down it.
#'
#' The master script records a real instance: an ED AFT fit returning
#' `b = 0.82`, `c = 296209.74` with log-likelihood -934.68, against a true
#' optimum of -907.19. Its rule, adopted verbatim here, is that any fit whose
#' shape parameters explode past `lim_shape` is rejected however good its
#' log-likelihood looks.
#'
#' The per-observation log-likelihood bound is additional. The master script
#' does not have it, but a four-parameter fit on the transplant data returned
#' a log-likelihood of `+5.2e77` on 91 observations, which no shape bound alone
#' would have caught.
#'
#' @param n Number of observations, used for the per-observation bound. Pass 0
#'   to disable that check.
#' @param log_shape,log_scale Names of the shape and scale parameters on the
#'   optimisation (log) scale.
#' @param lim_shape Shape parameters above this are degenerate. Default 500,
#'   the master script's value.
#' @param lim_scale Scale parameters outside `[1/lim_scale, lim_scale]` are
#'   degenerate. Default 1e8, the master script's `big`.
#' @param lim_ll_per_obs Mean log-density above this is implausible for real
#'   data and indicates a runaway rather than a fit.
#'
#' @return A function of one argument, suitable as `accept` in
#'   [optim_multistart()].
#' @noRd
.bd_make_accept <- function(n,
                            log_shape = c("log_a", "log_b", "log_c"),
                            log_scale = "log_k",
                            lim_shape = 500,
                            lim_scale = 1e8,
                            lim_ll_per_obs = 50) {
  force(n)
  function(fit) {
    if (is.null(fit) || is.null(fit$maximum) || !is.finite(fit$maximum))
      return(FALSE)
    if (n > 0 && fit$maximum / n > lim_ll_per_obs) return(FALSE)

    est <- fit$estimate
    if (is.null(est) || any(!is.finite(est))) return(FALSE)

    sh <- est[names(est) %in% log_shape]
    if (length(sh) && any(exp(sh) > lim_shape)) return(FALSE)

    sc <- est[names(est) %in% log_scale]
    if (length(sc) && (any(exp(sc) > lim_scale) ||
                       any(exp(sc) < 1 / lim_scale))) return(FALSE)

    TRUE
  }
}

#' Deterministic Starting Grid
#'
#' The shape combinations are those of `default_starts()` in the thesis master
#' script, which were chosen against the datasets in the dissertation. The
#' scale is set from the data rather than fixed, so the grid transfers to
#' arbitrary time units.
#'
#' A deterministic grid matters more than the number of random starts. The
#' master script uses fixed, sanely placed starts and rarely runs away; the
#' package used `runif(0.5, 5)` on every shape and did.
#'
#' @param submodel Logical; `TRUE` for the three-parameter ED submodel.
#' @param k_base Data-driven scale, typically the reciprocal of the mean
#'   observed event time.
#' @return A list of named numeric vectors on the log scale.
#' @noRd
.bd_default_starts <- function(submodel, k_base) {
  shapes <- if (submodel) {
    list(c(b = 5.4, c = 2.3), c(b = 5.0, c = 2.3), c(b = 4.5, c = 2.3),
         c(b = 2.0, c = 2.0), c(b = 4.0, c = 1.5), c(b = 1.5, c = 2.5),
         c(b = 8.0, c = 2.0), c(b = 3.0, c = 1.0))
  } else {
    list(c(a = 3.0, b = 4.0,  c = 2.3), c(a = 3.0, b = 24.0, c = 2.3),
         c(a = 5.7, b = 2.5,  c = 2.3), c(a = 1.0, b = 2.0,  c = 2.0),
         c(a = 2.0, b = 2.0,  c = 1.0), c(a = 0.7, b = 4.0,  c = 1.5),
         c(a = 3.0, b = 5.0,  c = 2.0), c(a = 1.2, b = 10.0, c = 3.0),
         c(a = 0.5, b = 20.0, c = 2.0))
  }

  out <- list()
  for (kf in c(1, 0.25)) {
    for (s in shapes) {
      v <- c(s, k = k_base * kf)
      out[[length(out) + 1L]] <- stats::setNames(log(v), paste0("log_", names(v)))
    }
  }
  out
}

#' Log-Space Addition
#' Computes log(exp(x) + exp(y)) stably
#' @noRd
logspace_add <- function(logx, logy) {
  m <- pmax(logx, logy)
  m + log(exp(logx - m) + exp(logy - m))
}
)---")

## =============================================================================
##  R/fit_models.R  --  use the grid and the guard
## =============================================================================

.step("Patching R/fit_models.R to use the guarded multi-start")

.sub_in("R/fit_models.R",
"  start_list <- vector(\"list\", n_starts)
  for (i in seq_len(n_starts)) {
    core <- c(log_b = log(stats::runif(1, 0.5, 5)),
              log_c = log(stats::runif(1, 0.5, 5)),
              log_k = log(k_base * stats::runif(1, 0.5, 2)))
    start_list[[i]] <- if (submodel) core else
      c(log_a = log(stats::runif(1, 0.5, 5)), core)
  }

  fit <- optim_multistart(ll_fun, start_list, method = method)
  if (is.null(fit))
    stop(\"Optimisation failed to converge from any starting point. Try a \",
         \"larger n_starts, or method = \\\"NM\\\".\", call. = FALSE)",
"  ## A deterministic grid first: it is what keeps the search out of the
  ## degenerate ridge. n_starts adds random starts on top of it.
  start_list <- .bd_default_starts(submodel, k_base)
  for (i in seq_len(n_starts)) {
    core <- c(log_b = log(stats::runif(1, 0.5, 5)),
              log_c = log(stats::runif(1, 0.5, 5)),
              log_k = log(k_base * stats::runif(1, 0.5, 2)))
    start_list[[length(start_list) + 1L]] <- if (submodel) core else
      c(log_a = log(stats::runif(1, 0.5, 5)), core)
  }

  accept <- .bd_make_accept(n = length(time))
  fit <- optim_multistart(ll_fun, start_list, method = method, accept = accept)
  if (is.null(fit))
    stop(\"No admissible optimum was found. Every start either failed or \",
         \"landed on a degenerate ridge, where the shape parameters explode \",
         \"and no finite maximiser exists. Try method = \\\"NM\\\", or fit the \",
         \"ED submodel with submodel = TRUE.\", call. = FALSE)",
        "fit_betadanish: deterministic grid + acceptance predicate")

.sub_in("R/fit_models.R",
"    convergence  = fit$code,
    message      = fit$message,",
"    convergence  = fit$code,
    message      = fit$message,
    starts_ok        = .bd_or(attr(fit, \"bd_starts_ok\"), NA_integer_),
    starts_rejected  = .bd_or(attr(fit, \"bd_starts_rejected\"), NA_integer_),
    loglik_spread    = .bd_or(attr(fit, \"bd_loglik_spread\"), NA_real_),",
        "fit_betadanish: record multi-start diagnostics")

.step("Extending the diagnostics block")

.sub_in("R/fit_models.R",
"  d <- list(
    converged        = isTRUE(fit$convergence %in% c(0L, 1L, 2L)),
    convergence_code = fit$convergence,
    vcov_singular    = anyNA(fit$vcov) || any(!is.finite(diag(fit$vcov))) ||
                       any(diag(fit$vcov) <= 0),
    near_b_ridge     = NA,
    b_distance_se    = NA_real_,
    ac_correlation   = NA_real_
  )",
"  d <- list(
    converged        = isTRUE(fit$convergence %in% c(0L, 1L, 2L)),
    convergence_code = fit$convergence,
    vcov_singular    = anyNA(fit$vcov) || any(!is.finite(diag(fit$vcov))) ||
                       any(diag(fit$vcov) <= 0),
    near_b_ridge     = NA,
    b_distance_se    = NA_real_,
    ac_correlation   = NA_real_,
    starts_ok        = fit$starts_ok,
    starts_rejected  = fit$starts_rejected,
    loglik_spread    = fit$loglik_spread,
    loglik_per_obs   = if (is.null(fit$nobs) || !fit$nobs) NA_real_ else
                         as.numeric(fit$logLik) / fit$nobs
  )",
        "diagnostics: multi-start fields")

.sub_in("R/fit_models.R",
"  if (is.finite(d$ac_correlation) && abs(d$ac_correlation) > 0.95)
    warning(sprintf(paste0(\"The fitted correlation between a-hat and c-hat is \",
                           \"%.3f, so effectively only the product c*a is \",
                           \"identified. Individual estimates of a and c should \",
                           \"not be interpreted.\"),
                    d$ac_correlation), call. = FALSE)

  invisible(NULL)
}",
"  if (is.finite(d$ac_correlation) && abs(d$ac_correlation) > 0.95)
    warning(sprintf(paste0(\"The fitted correlation between a-hat and c-hat is \",
                           \"%.3f, so effectively only the product c*a is \",
                           \"identified. Individual estimates of a and c should \",
                           \"not be interpreted.\"),
                    d$ac_correlation), call. = FALSE)

  if (!is.na(d$starts_rejected) && d$starts_rejected > 0)
    warning(sprintf(paste0(\"%d starting point(s) reached a degenerate ridge \",
                           \"and were discarded. The reported fit is the best \",
                           \"admissible optimum. If this is most of the grid, \",
                           \"the four-parameter model is a poor choice for \",
                           \"these data.\"), d$starts_rejected), call. = FALSE)

  if (!is.na(d$loglik_spread) && d$loglik_spread > 2)
    warning(sprintf(paste0(\"Accepted optima span %.2f log-likelihood units, \",
                           \"so the surface has several local maxima and the \",
                           \"reported fit may not be global. Increase \",
                           \"n_starts.\"), d$loglik_spread), call. = FALSE)

  invisible(NULL)
}",
        "diagnostics: degeneracy and multi-start warnings")

## =============================================================================
##  AFT and cure  --  the failure mode the master script actually recorded
## =============================================================================

.step("Guarding fit_bd_aft and fit_bd_cure")

.sub_in("R/aft_models.R",
        "  fit <- optim_multistart(ll_fun, start_list, method = method)",
        "  ## The master script's recorded degenerate case was an ED AFT fit that\n  ## returned c = 296209.74. Shape explosion is rejected; large regression\n  ## coefficients are not, following .ch6_separated: near-separation is\n  ## information, not an error.\n  fit <- optim_multistart(\n    ll_fun, start_list, method = method,\n    accept = .bd_make_accept(n = length(surv_data$time),\n                             log_shape = c(\"log_b\", \"log_c\"),\n                             log_scale = character(0)))",
        "fit_bd_aft: shape-explosion guard")

.sub_in("R/cure_models.R",
        "  fit <- optim_multistart(ll_fun, start_list, method = method)",
        "  fit <- optim_multistart(\n    ll_fun, start_list, method = method,\n    accept = .bd_make_accept(n = length(time),\n                             log_shape = c(\"log_b\", \"log_c\"),\n                             log_scale = character(0)))",
        "fit_bd_cure: shape-explosion guard")

## =============================================================================
##  report_betadanish  --  surface the new flags
## =============================================================================

.step("Surfacing the new flags in report_betadanish()")

.sub_in("R/report_betadanish.R",
"    if (is.finite(d$ac_correlation) && abs(d$ac_correlation) > 0.95)
      flags <- c(flags, \"(a, c) confounded\")",
"    if (is.finite(d$ac_correlation) && abs(d$ac_correlation) > 0.95)
      flags <- c(flags, \"(a, c) confounded\")
    if (!is.null(d$starts_rejected) && !is.na(d$starts_rejected) &&
        d$starts_rejected > 0)
      flags <- c(flags, sprintf(\"%d degenerate start(s) discarded\",
                                d$starts_rejected))
    if (!is.null(d$loglik_spread) && !is.na(d$loglik_spread) &&
        d$loglik_spread > 2)
      flags <- c(flags, sprintf(\"local optima span %.2f log-lik units\",
                                d$loglik_spread))",
        "report_betadanish: degeneracy flags")

## =============================================================================
##  TESTS
## =============================================================================

.step("Writing tests/testthat/test-degeneracy.R")

.put("tests/testthat/test-degeneracy.R", r"---(## Regression tests for the degenerate-ridge guard.
##
## The motivating failure: a four-parameter fit on the transplant data returned
## a log-likelihood of +5.2e77 on 91 observations, and the unguarded
## multi-start selected it as the maximum.

fake_fit <- function(ll, ...) {
  list(maximum = ll, estimate = c(...))
}

test_that("the acceptance predicate rejects exploded shapes", {
  acc <- BetaDanish:::.bd_make_accept(n = 100)

  ok <- fake_fit(-250, log_a = log(1.5), log_b = log(3),
                 log_c = log(2), log_k = log(0.5))
  expect_true(acc(ok))

  ## The master script's recorded case: c = 296209.74
  blown <- fake_fit(-934.68, log_a = log(1), log_b = log(0.82),
                    log_c = log(296209.74), log_k = log(0.5))
  expect_false(acc(blown))

  ## Exactly at the limit is admissible; past it is not.
  expect_true(acc(fake_fit(-250, log_a = log(1), log_b = log(499),
                           log_c = log(2), log_k = log(0.5))))
  expect_false(acc(fake_fit(-250, log_a = log(1), log_b = log(501),
                            log_c = log(2), log_k = log(0.5))))
})

test_that("the acceptance predicate rejects an implausible log-likelihood", {
  acc <- BetaDanish:::.bd_make_accept(n = 91)
  ## The observed runaway. Shapes alone would not have caught it.
  expect_false(acc(fake_fit(5.210644e77, log_a = log(2), log_b = log(3),
                            log_c = log(2), log_k = log(0.1))))
  expect_false(acc(fake_fit(Inf, log_a = 0, log_b = 0, log_c = 0, log_k = 0)))
  expect_false(acc(fake_fit(NaN, log_a = 0, log_b = 0, log_c = 0, log_k = 0)))
  ## A large but attainable positive log-likelihood stays admissible.
  expect_true(acc(fake_fit(91 * 3, log_a = 0, log_b = log(2),
                           log_c = 0, log_k = 0)))
})

test_that("the acceptance predicate rejects an exploded or vanished scale", {
  acc <- BetaDanish:::.bd_make_accept(n = 100)
  expect_false(acc(fake_fit(-250, log_a = 0, log_b = log(2),
                            log_c = 0, log_k = log(1e9))))
  expect_false(acc(fake_fit(-250, log_a = 0, log_b = log(2),
                            log_c = 0, log_k = log(1e-9))))
})

test_that("non-finite estimates are inadmissible", {
  acc <- BetaDanish:::.bd_make_accept(n = 100)
  expect_false(acc(fake_fit(-250, log_a = NaN, log_b = 0,
                            log_c = 0, log_k = 0)))
  expect_false(acc(NULL))
})

test_that("the deterministic start grid is well formed", {
  for (sub in c(TRUE, FALSE)) {
    g <- BetaDanish:::.bd_default_starts(sub, k_base = 0.02)
    expect_gt(length(g), 10L)
    for (s in g) {
      expect_true(all(is.finite(s)))
      expect_true(all(grepl("^log_", names(s))))
      expect_equal(length(s), if (sub) 3L else 4L)
      ## Every start is itself admissible, so the grid can never seed the ridge.
      expect_true(all(exp(s[names(s) != "log_k"]) < 500))
    }
  }
})

test_that("the four-parameter fit no longer runs away on the transplant data", {
  skip_on_cran()
  data(transplant, package = "BetaDanish", envir = environment())
  set.seed(99)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                   n_starts = 5, check_identifiability = FALSE))

  n <- nrow(transplant)
  ## The bug produced logLik / n = 5.7e75. Anything remotely sane is far below.
  expect_lt(fit$logLik / n, 5)
  expect_gt(fit$logLik, -1e4)
  expect_true(all(fit$coefficients[intersect(names(fit$coefficients),
                                             c("a", "b", "c"))] < 500))
  expect_true(is.finite(fit$AIC))
})

test_that("a likelihood ratio test on real data is no longer absurd", {
  skip_on_cran()
  data(transplant, package = "BetaDanish", envir = environment())
  set.seed(101)
  full <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                   submodel = FALSE, n_starts = 5, check_identifiability = FALSE))
  sub <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                   submodel = TRUE, n_starts = 5, check_identifiability = FALSE))

  ## The full model cannot fit worse than its own submodel by more than
  ## optimiser noise, and one extra parameter cannot buy an unbounded gain.
  chisq <- 2 * (full$logLik - sub$logLik)
  expect_gt(chisq, -1)
  expect_lt(chisq, 100)
  expect_true(is.finite(chisq))
})

test_that("multi-start diagnostics are recorded on the fit", {
  skip_on_cran()
  dat <- simulate_bd_data(120, a = 1, b = 3, c = 2, k = 0.5, seed = 5)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  expect_true(is.numeric(fit$starts_ok))
  expect_gte(fit$starts_ok, 1L)
  expect_true(is.numeric(fit$loglik_spread))
  expect_gte(fit$loglik_spread, 0)
  expect_true(is.numeric(fit$diagnostics$loglik_per_obs))
})
)---")

## =============================================================================
##  NEWS
## =============================================================================

.step("Recording the fix in NEWS.md")

.nw <- readLines("NEWS.md", warn = FALSE)
if (!any(grepl("degenerate ridge", .nw, fixed = TRUE))) {
  .hdr <- grep("^## Correctness fixes\\s*$", .nw)
  if (length(.hdr) == 1L) {
    .backup("NEWS.md")
    .sec <- c(
      "",
      "* **Degenerate optima are no longer reported as fits.** The Beta-Danish",
      "  likelihood has ridges along which the objective can be driven",
      "  arbitrarily high with no finite maximiser: as `k` and `c` grow with",
      "  `c/k` fixed, the family collapses to a Frechet-type limit.",
      "  `optim_multistart()` selected purely on the highest log-likelihood",
      "  behind nothing but a finiteness test, so a runaway beat every honest",
      "  optimum. A four-parameter fit on the `transplant` data returned a",
      "  log-likelihood of `+5.2e77` on 91 observations, and the resulting",
      "  likelihood ratio test reported a chi-squared statistic of `1.0e78`.",
      "",
      "  `fit_betadanish()`, `fit_bd_aft()` and `fit_bd_cure()` now search over",
      "  a deterministic start grid before any random starts, and accept only",
      "  admissible optima: shape parameters bounded by 500, scale within",
      "  `1e-8` to `1e8`, and a plausible per-observation log-likelihood. The",
      "  thresholds and the start grid are those already used by the thesis",
      "  analysis pipeline, so the two now agree.",
      "",
      "* The fitted object records `starts_ok`, `starts_rejected` and",
      "  `loglik_spread`. Warnings fire when starts are discarded as degenerate",
      "  and when the accepted optima span more than two log-likelihood units,",
      "  which is the signal that `n_starts` is too low for the surface.")
    .nw <- append(.nw, .sec, after = .hdr)
    .write_lines("NEWS.md", .nw)
    .ok("NEWS.md updated")
  } else {
    .warn("anchor not found in NEWS.md; add the note by hand")
  }
} else {
  .info("NEWS.md already records it")
}

## =============================================================================
##  VERIFY
## =============================================================================

.step("Parsing all R and test files")
.targets <- c(list.files("R", pattern = "[.]R$", full.names = TRUE),
              list.files("tests", pattern = "[.]R$", recursive = TRUE, full.names = TRUE))
.bad <- character(0)
for (f in .targets) {
  e <- tryCatch({ parse(f); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(e)) .bad <- c(.bad, paste0("  ", f, ": ", e))
}
if (length(.bad)) .die("These files do not parse:\n", paste(.bad, collapse = "\n"),
                       "\n\nBackups: ", BACKUP_DIR)
.ok(sprintf("%d file(s) parse cleanly", length(.targets)))

.step("Confirming the guard is wired into all three fitters")
for (f in c("R/fit_models.R", "R/aft_models.R", "R/cure_models.R")) {
  if (any(grepl(".bd_make_accept", readLines(f, warn = FALSE), fixed = TRUE)))
    .ok(paste(f, "guarded"))
  else
    .warn(paste(f, "-- guard not detected; check by hand"))
}

.step("devtools::document()")
.r <- tryCatch({ devtools::document(); TRUE }, error = function(e) conditionMessage(e))
if (!isTRUE(.r)) .die("document() failed:\n  ", .r, "\n\nBackups: ", BACKUP_DIR)
.ok("documentation regenerated")

.step("devtools::test()")
.t <- tryCatch(devtools::test(), error = function(e) { .warn(conditionMessage(e)); NULL })

.step("devtools::check() -- several minutes, do not interrupt")
.chk <- tryCatch(devtools::check(document = FALSE, args = "--as-cran", error_on = "never"),
                 error = function(e) { .warn(conditionMessage(e)); NULL })

cat("\n", strrep("=", 78), "\n", sep = "")
if (!is.null(.chk)) {
  cat("  CHECK RESULT\n", strrep("=", 78), "\n", sep = "")
  cat(sprintf("  errors=%d  warnings=%d  notes=%d\n",
              length(.chk$errors), length(.chk$warnings), length(.chk$notes)))
  for (nm in c("errors", "warnings", "notes")) {
    if (length(.chk[[nm]])) {
      cat("\n---- ", toupper(nm), " ----\n", sep = "")
      cat(.chk[[nm]], sep = "\n\n")
    }
  }
} else {
  cat("  check() did not complete; run devtools::check() manually.\n")
}

cat("\n", strrep("=", 78), "\n", sep = "")
cat("  PATCH 2c COMPLETE  --  degenerate-optimum guard\n")
cat(strrep("=", 78), "\n\n")
cat("  optim_multistart() gains an acceptance predicate\n")
cat("  .bd_make_accept()  ports .ch6_degenerate (lim_shape = 500) plus a new\n")
cat("                     per-observation log-likelihood bound\n")
cat("  .bd_default_starts() ports default_starts() with a data-driven scale\n")
cat("  fit_betadanish / fit_bd_aft / fit_bd_cure all guarded\n")
cat("  starts_ok, starts_rejected, loglik_spread recorded and warned on\n\n")
cat("  n_starts now means RANDOM starts added on top of the deterministic\n")
cat("  grid, so every fit searches at least 16 sane starting points.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
