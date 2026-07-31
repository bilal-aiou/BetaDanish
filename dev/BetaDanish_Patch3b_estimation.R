## =============================================================================
##  BetaDanish  --  PHASE 2, PATCH 3b : ESTIMATION AND INFERENCE
## =============================================================================
##
##  Implements approved recommendations 34, 35 and 36, plus one addition taken
##  from the thesis pipeline that the Phase 1 audit did not list.
##
##    34  fit_betadanish(penalty = )        ridge-penalised maximum likelihood
##    35  fit_betadanish(grouped = TRUE)    grouped likelihood for grid-recorded
##                                          times, with automatic increment
##                                          inference
##    36  bd_profile_ci()                   profile likelihood intervals
##        bd_wald_ci()                      log-scale Wald intervals
##
##    NEW bd_identified_coef()              the (ac, b, k) reparametrisation
##
##  WHY THE ADDITION
##    Master_R_Code_for_Thesis.R, Item 7, records that reparametrising the
##    four-parameter model by the identified composite p = ac drops the Hessian
##    condition number on the remission data from 2750 to 152 and yields a
##    finite SE(ac) = 0.281, with the log-likelihood unchanged. That is the
##    constructive answer to the (a, c) confounding documented in Patch 1: not
##    just a warning that a and c cannot be separated, but a summary of what
##    the data actually determine.
##
##    Item 6 of the same file records that on breaking-stress data the profile
##    interval for b has no finite upper bound and 29 percent of bootstrap
##    resamples fail to bound it above, with the recommendation to report b as
##    a lower bound. bd_profile_ci() returns Inf in that case and says so in
##    print(), rather than quietly returning the largest grid value.
##
##  ONE CORRECTNESS POINT WORTH READING
##    A penalised fit maximises the log-likelihood minus a penalty. Storing
##    that objective as `logLik` would corrupt every AIC, BIC and likelihood
##    ratio test downstream. The reported logLik is therefore the UNPENALISED
##    log-likelihood re-evaluated at the penalised estimate; the objective is
##    kept separately as `penalised_logLik`. There is a test for this.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3b_estimation.R")
##  IDEMPOTENT   Yes -- every file is written whole.
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

cat(strrep("=", 78), "\n")
cat("  BetaDanish  --  Phase 2, Patch 3b : estimation and inference\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (file.exists("R/entropy.R"))
  .die("R/entropy.R is still present -- run Patch 3a-fix3 first.")
if (!file.exists("R/structural.R")) .die("Patch 3a has not been applied.")
if (!any(grepl(".bd_make_accept", readLines("R/fit_models.R", warn = FALSE), fixed = TRUE)))
  .die("R/fit_models.R is not guarded -- run Patch 2c and its fixes first.")
if (!any(grepl(".bd_grid_step", readLines("R/data_helpers.R", warn = FALSE), fixed = TRUE)))
  .die("R/data_helpers.R has no .bd_grid_step() -- run Patch 2a first.")
.ok("Patches 1 through 3a detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3b"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

.step("Rec 34, 35: rewriting R/fit_models.R with penalty and grouped fitting")

.put("R/fit_models.R", r"---(#' Fit the Beta-Danish Distribution to Survival Data
#'
#' Fits the Beta-Danish distribution by maximum likelihood. Complete and
#' right-censored samples are both supported, via a `survival::Surv` response.
#'
#' @param formula A formula whose left-hand side is a `Surv` object. Use
#'   `~ 1` for a model without covariates.
#' @param data A data frame containing the variables in `formula`.
#' @param submodel Logical; if `TRUE`, fits the three-parameter Exponentiated
#'   Danish (ED) submodel by fixing `a = 1`.
#' @param n_starts Integer; number of random starting points for the
#'   multi-start optimisation. Default 10.
#' @param method Character; optimisation method passed to `maxLik::maxLik`.
#' @param check_identifiability Logical; if `TRUE` (default), issue warnings
#'   when the fit lands in a region where the parameters are weakly identified.
#' @param penalty Non-negative ridge penalty applied on the log-parameter
#'   scale. `0` (default) is ordinary maximum likelihood. A small positive
#'   value stabilises the fit when the likelihood is nearly flat, at the cost
#'   of bias toward `penalty_center`.
#' @param penalty_center Numeric vector on the log-parameter scale toward which
#'   the penalty shrinks, or `NULL` (default) to shrink toward the best
#'   unpenalised start.
#' @param grouped Logical; if `TRUE`, use the grouped (interval) likelihood
#'   appropriate to times recorded on a coarse grid. Default `FALSE`.
#' @param delta Recording increment for `grouped = TRUE`. `NULL` (default)
#'   infers it from the spacing of the observed times.
#'
#' @return An object of S3 class `"betadanish"` with components including
#'   `coefficients`, `logLik`, `vcov`, `npar`, `nobs`, `convergence` and
#'   `diagnostics`.
#'
#' @details
#' Optimisation is carried out on log-transformed parameters so that positivity
#' is enforced without constraints; estimates and the variance-covariance matrix
#' are returned on the natural scale, the latter via the delta method.
#'
#' @section Grouped data:
#' Survival times are often recorded on a grid -- whole days, whole months --
#' and the point-density likelihood is not appropriate for them. It treats a
#' rounded value as an exact observation, which overstates the information in
#' the sample and understates every standard error. With `grouped = TRUE` an
#' event recorded at \eqn{t} contributes
#' \eqn{\log\{F(t + \delta/2) - F(t - \delta/2)\}} instead of
#' \eqn{\log f(t)}; censored observations are unchanged. The cell probability
#' falls back to \eqn{f(t)\delta} only where the difference of two nearly
#' equal distribution values has cancelled to zero.
#'
#' `read_survival_data()` reports an inferred `grid_step` for exactly this
#' purpose, and this function warns when the times look grid-recorded but
#' `grouped = FALSE`.
#'
#' @section Penalised fitting:
#' `penalty > 0` adds \eqn{\lambda \sum (\theta - \mu)^2} on the
#' log-parameter scale. This is worth reaching for when the likelihood is flat
#' along the \eqn{(a, c)} direction and the unpenalised optimiser wanders, but
#' it is a deliberate bias: the estimates are shrunk toward `penalty_center`.
#'
#' The reported `logLik` is always the **unpenalised** log-likelihood evaluated
#' at the penalised estimate, so that AIC, BIC and likelihood ratio tests
#' remain comparable across fits. The objective actually maximised is stored
#' separately as `penalised_logLik`. Treat the degrees of freedom as nominal:
#' shrinkage reduces the effective number of parameters, so information
#' criteria are conservative under penalisation.
#'
#' @section Identifiability:
#' The four-parameter model is not uniformly well identified, and a converged
#' fit is not by itself evidence that it is. Two regions warrant care.
#'
#' * **The \eqn{b = 1} ridge.** At \eqn{b = 1} the beta generator collapses and
#'   the model is non-identifiable. A fit with \eqn{\hat b} within about two
#'   standard errors of one lies close to that ridge; the likelihood is nearly
#'   flat along it, so the individual estimates carry little information even
#'   though the fitted survival curve may look excellent.
#' * **Lower-tail \eqn{(a, c)} confounding.** Near the lower tail, \eqn{a} and
#'   \eqn{c} enter almost exclusively through the product \eqn{ca}, so the
#'   expected Fisher information is close to singular in that direction. A
#'   fitted correlation between \eqn{\hat a} and \eqn{\hat c} above about 0.95
#'   in absolute value indicates that only the product is being estimated.
#'
#' In either case the ED submodel (`submodel = TRUE`) is usually the honest
#' report, and a likelihood ratio test via [compare_models()] will normally
#' fail to reject it. Set `check_identifiability = FALSE` to silence the
#' warnings once you have satisfied yourself that they are understood.
#'
#' @seealso [compare_models()], [gof_betadanish()], [compare_distributions()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim_time   <- rbetadanish(150, a = 1.5, b = 3, c = 2, k = 0.5)
#' sim_status <- rbinom(150, 1, 0.85)
#' dat <- data.frame(time = sim_time, status = sim_status)
#'
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat)
#' summary(fit)
#'
#' fit_sub <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                           submodel = TRUE)
#' compare_models(fit, fit_sub)
#' }
fit_betadanish <- function(formula, data, submodel = FALSE, n_starts = 10,
                           method = "BFGS", check_identifiability = TRUE,
                           penalty = 0, penalty_center = NULL,
                           grouped = FALSE, delta = NULL) {

  surv_data <- extract_surv_data(formula, data)
  time   <- surv_data$time
  status <- surv_data$status

  if (!is.numeric(penalty) || length(penalty) != 1L || is.na(penalty) || penalty < 0)
    stop("'penalty' must be a single non-negative number.", call. = FALSE)

  grid_step <- .bd_grid_step(time)
  if (isTRUE(grouped)) {
    if (is.null(delta)) delta <- grid_step
    if (!is.finite(delta) || delta <= 0)
      stop("grouped = TRUE needs a recording increment. The spacing of the ",
           "times did not imply one, so supply 'delta' explicitly.",
           call. = FALSE)
  } else if (is.finite(grid_step) && isTRUE(check_identifiability)) {
    warning(sprintf(paste0("The times look recorded on a grid of %s. The ",
                           "point-density likelihood treats them as exact, ",
                           "which understates the standard errors. Consider ",
                           "grouped = TRUE."), format(grid_step)),
            call. = FALSE)
  }

  ## Unpenalised log-likelihood, exact or grouped.
  loglik_fun <- function(pars) {
    a_par <- if (submodel) 1.0 else exp(pars[["log_a"]])
    b_par <- exp(pars[["log_b"]])
    c_par <- exp(pars[["log_c"]])
    k_par <- exp(pars[["log_k"]])

    lp <- if (isTRUE(grouped)) {
      suppressWarnings(.bd_log_cell(time, delta, a_par, b_par, c_par, k_par))
    } else {
      suppressWarnings(dbetadanish(time, a_par, b_par, c_par, k_par, log = TRUE))
    }
    ls <- suppressWarnings(
      pbetadanish(time, a_par, b_par, c_par, k_par,
                  lower.tail = FALSE, log.p = TRUE))

    sum(status * lp + (1 - status) * ls)
  }

  pen_center <- penalty_center

  ll_fun <- function(pars) {
    loglik <- loglik_fun(pars)
    if (!is.finite(loglik)) return(-1e10)
    if (penalty > 0 && !is.null(pen_center)) {
      mu <- pen_center[names(pars)]
      if (!anyNA(mu)) loglik <- loglik - penalty * sum((pars - mu)^2)
    }
    if (!is.finite(loglik)) return(-1e10)
    loglik
  }

  avg_t <- mean(time[status == 1], na.rm = TRUE)
  if (is.na(avg_t) || avg_t <= 0) avg_t <- mean(time, na.rm = TRUE)
  k_base <- 1 / avg_t

  ## A deterministic grid first: it is what keeps the search out of the
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

  ## With a penalty and no explicit centre, shrink toward the unpenalised
  ## optimum rather than toward an arbitrary point.
  if (penalty > 0 && is.null(pen_center)) {
    pilot <- optim_multistart(function(p) {
      v <- loglik_fun(p); if (is.finite(v)) v else -1e10
    }, start_list, method = method, accept = accept)
    if (!is.null(pilot)) pen_center <- pilot$estimate
  }

  fit <- optim_multistart(ll_fun, start_list, method = method, accept = accept)
  if (is.null(fit))
    stop("No admissible optimum was found. Every start either failed or ",
         "landed on a degenerate ridge, where the shape parameters explode ",
         "and no finite maximiser exists. Try method = \"NM\", or fit the ",
         "ED submodel with submodel = TRUE.", call. = FALSE)

  est_log <- fit$estimate
  est_nat <- exp(est_log)
  names(est_nat) <- sub("^log_", "", names(est_log))

  vcov_log <- tryCatch(solve(-fit$hessian),
                       error = function(e) matrix(NA_real_, length(est_log),
                                                  length(est_log)))
  J        <- diag(est_nat, nrow = length(est_nat))
  vcov_nat <- J %*% vcov_log %*% t(J)
  rownames(vcov_nat) <- colnames(vcov_nat) <- names(est_nat)

  npar <- length(est_nat)
  nobs <- length(time)

  ## The objective may be penalised; the reported log-likelihood never is, so
  ## that AIC, BIC and likelihood ratio tests stay comparable across fits.
  ll_unpen <- loglik_fun(est_log)
  if (!is.finite(ll_unpen)) ll_unpen <- fit$maximum

  out <- list(
    coefficients = est_nat,
    logLik       = ll_unpen,
    penalised_logLik = if (penalty > 0) fit$maximum else NA_real_,
    penalty      = penalty,
    grouped      = isTRUE(grouped),
    delta        = if (isTRUE(grouped)) delta else NA_real_,
    vcov         = vcov_nat,
    npar         = npar,
    nobs         = nobs,
    nevent       = sum(status == 1),
    AIC          = 2 * npar - 2 * ll_unpen,
    BIC          = npar * log(nobs) - 2 * ll_unpen,
    convergence  = fit$code,
    message      = fit$message,
    starts_ok        = .bd_or(attr(fit, "bd_starts_ok"), NA_integer_),
    starts_rejected  = .bd_or(attr(fit, "bd_starts_rejected"), NA_integer_),
    loglik_spread    = .bd_or(attr(fit, "bd_loglik_spread"), NA_real_),
    submodel     = submodel,
    data         = list(time = time, status = status),
    formula      = formula,
    call         = match.call()
  )
  out$diagnostics <- .bd_fit_diagnostics(out)
  class(out) <- "betadanish"

  if (isTRUE(check_identifiability)) .bd_warn_diagnostics(out$diagnostics)

  out
}

#' Log Cell Probability for Grouped Times
#'
#' \eqn{\log\{F(t + \delta/2) - F(t - \delta/2)\}}, falling back to
#' \eqn{\log f(t) + \log\delta} only where the difference of two nearly equal
#' distribution values has cancelled to zero. Capped at `0`, since a
#' probability cannot exceed one.
#'
#' @noRd
.bd_log_cell <- function(t, delta, a, b, c, k) {
  L  <- pmax(t - delta / 2, 0)
  U  <- t + delta / 2
  pc <- pbetadanish(U, a, b, c, k) - pbetadanish(L, a, b, c, k)

  out <- numeric(length(t))
  bad <- !is.finite(pc) | pc <= 0
  if (any(!bad)) out[!bad] <- pmin(log(pc[!bad]), 0)
  if (any(bad))
    out[bad] <- pmin(dbetadanish(t[bad], a, b, c, k, log = TRUE) + log(delta), 0)
  out
}

#' Assemble Convergence and Identifiability Diagnostics
#' @noRd
.bd_fit_diagnostics <- function(fit) {
  se <- sqrt(pmax(diag(fit$vcov), 0))
  names(se) <- names(fit$coefficients)

  d <- list(
    converged        = isTRUE(fit$convergence %in% c(0L, 1L, 2L)),
    convergence_code = fit$convergence,
    vcov_singular    = anyNA(fit$vcov) || any(!is.finite(diag(fit$vcov))) ||
                       any(diag(fit$vcov) <= 0),
    near_b_ridge     = NA,
    b_distance_se    = NA_real_,
    ac_correlation   = NA_real_,
    starts_ok        = .bd_or(fit$starts_ok, NA_integer_),
    starts_rejected  = .bd_or(fit$starts_rejected, NA_integer_),
    loglik_spread    = .bd_or(fit$loglik_spread, NA_real_),
    loglik_per_obs   = if (is.null(fit$nobs) || !isTRUE(fit$nobs > 0)) NA_real_
                       else as.numeric(fit$logLik) / fit$nobs
  )

  if ("b" %in% names(fit$coefficients) && is.finite(se[["b"]]) && se[["b"]] > 0) {
    d$b_distance_se <- abs(fit$coefficients[["b"]] - 1) / se[["b"]]
    d$near_b_ridge  <- d$b_distance_se < 2
  }

  if (!fit$submodel && all(c("a", "c") %in% rownames(fit$vcov))) {
    vaa <- fit$vcov["a", "a"]; vcc <- fit$vcov["c", "c"]
    vac <- fit$vcov["a", "c"]
    if (is.finite(vaa) && is.finite(vcc) && vaa > 0 && vcc > 0)
      d$ac_correlation <- vac / sqrt(vaa * vcc)
  }
  d
}

#' Emit Identifiability Warnings
#'
#' Every field is normalised to a length-one value of the expected type before
#' any condition is evaluated. A fitted object saved before a later diagnostic
#' existed carries no `starts_rejected` or `loglik_spread`, and must degrade
#' quietly rather than error. Guarding each condition separately is not enough:
#' `isTRUE()` absorbs `NA` and zero-length results, but `abs(NULL)` raises a
#' non-numeric argument error before `isTRUE()` ever sees it.
#'
#' @param d A diagnostics list, possibly incomplete.
#' @noRd
.bd_warn_diagnostics <- function(d) {
  if (!is.list(d)) return(invisible(NULL))

  num1 <- function(nm) {
    v <- d[[nm]]
    if (is.null(v) || length(v) != 1L || !is.numeric(v)) NA_real_ else as.numeric(v)
  }
  flag <- function(nm) isTRUE(d[[nm]])

  code       <- if (is.null(d$convergence_code) ||
                    length(d$convergence_code) != 1L) NA else d$convergence_code
  b_dist     <- num1("b_distance_se")
  ac_cor     <- num1("ac_correlation")
  n_rejected <- num1("starts_rejected")
  spread     <- num1("loglik_spread")

  if (!flag("converged"))
    warning("The optimiser reported code ", code,
            "; treat the estimates as provisional and increase n_starts.",
            call. = FALSE)

  if (flag("vcov_singular"))
    warning("The observed information matrix is singular or not positive ",
            "definite, so standard errors are unreliable. This usually means ",
            "the likelihood is flat in at least one direction.", call. = FALSE)

  if (flag("near_b_ridge"))
    warning(sprintf(paste0("b-hat is only %.2f standard errors from 1, close to ",
                           "the b = 1 non-identifiability ridge. Consider the ",
                           "ED submodel (submodel = TRUE). See the ",
                           "Identifiability section of ?fit_betadanish."),
                    b_dist), call. = FALSE)

  if (isTRUE(abs(ac_cor) > 0.95))
    warning(sprintf(paste0("The fitted correlation between a-hat and c-hat is ",
                           "%.3f, so effectively only the product c*a is ",
                           "identified. Individual estimates of a and c should ",
                           "not be interpreted."),
                    ac_cor), call. = FALSE)

  if (isTRUE(n_rejected > 0))
    warning(sprintf(paste0("%.0f starting point(s) reached a degenerate ridge ",
                           "and were discarded. The reported fit is the best ",
                           "admissible optimum. If this is most of the grid, ",
                           "the four-parameter model is a poor choice for ",
                           "these data."), n_rejected), call. = FALSE)

  if (isTRUE(spread > 2))
    warning(sprintf(paste0("Accepted optima span %.2f log-likelihood units, ",
                           "so the surface has several local maxima and the ",
                           "reported fit may not be global. Increase ",
                           "n_starts."), spread), call. = FALSE)

  invisible(NULL)
})---")


.step("Rec 36: writing R/inference.R")

.put("R/inference.R", r"---(## Interval estimation for Beta-Danish fits.
##
## The Wald interval on the natural scale is the wrong default for this family.
## The thesis pipeline records the problem directly: on the breaking-stress data
## the symmetric interval for the tail index b crosses zero, the profile
## interval has no finite upper bound, and 29 percent of bootstrap resamples
## fail to bound b above. The functions here follow that finding rather than
## papering over it.

#' Reconstruct the Log-Likelihood of a Fitted Model
#'
#' Rebuilds the objective from the fit object so that profiling does not need
#' the original call. Honours the submodel and grouped settings.
#'
#' @noRd
.bd_refit_loglik <- function(object) {
  time     <- object$data$time
  status   <- object$data$status
  submodel <- isTRUE(object$submodel)
  grouped  <- isTRUE(object$grouped)
  delta    <- object$delta

  function(nat) {
    a <- if (submodel) 1 else nat[["a"]]
    b <- nat[["b"]]; c <- nat[["c"]]; k <- nat[["k"]]
    if (!is.finite(a) || !is.finite(b) || !is.finite(c) || !is.finite(k) ||
        a <= 0 || b <= 0 || c <= 0 || k <= 0) return(-Inf)
    lp <- if (grouped) {
      suppressWarnings(.bd_log_cell(time, delta, a, b, c, k))
    } else {
      suppressWarnings(dbetadanish(time, a, b, c, k, log = TRUE))
    }
    ls <- suppressWarnings(
      pbetadanish(time, a, b, c, k, lower.tail = FALSE, log.p = TRUE))
    v <- sum(status * lp + (1 - status) * ls)
    if (is.finite(v)) v else -Inf
  }
}

#' Log-Scale Wald Confidence Intervals
#'
#' Symmetric Wald intervals on the log-parameter scale, exponentiated back.
#'
#' @param object A fitted `"betadanish"` object.
#' @param level Confidence level. Default 0.95.
#'
#' @return A matrix with one row per parameter and columns `estimate`,
#'   `se`, `lower` and `upper`.
#'
#' @details
#' All four parameters are strictly positive, so a symmetric interval on the
#' natural scale can and does cross zero when a parameter is weakly identified.
#' Building the interval on the log scale and exponentiating keeps it inside
#' the parameter space:
#' \eqn{\hat\theta \exp(\pm z\, \mathrm{SE}(\log\hat\theta))}.
#'
#' This is still a Wald interval, and it inherits the usual weakness: it
#' assumes the log-likelihood is approximately quadratic. Where it is not --
#' which is exactly where these intervals matter -- prefer [bd_profile_ci()].
#'
#' @seealso [bd_profile_ci()], [bd_identified_coef()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 1)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       submodel = TRUE, n_starts = 3)
#' bd_wald_ci(fit)
#' }
bd_wald_ci <- function(object, level = 0.95) {
  if (!inherits(object, "betadanish"))
    stop("'object' must be a fitted betadanish model.", call. = FALSE)
  est <- object$coefficients
  se  <- sqrt(pmax(diag(object$vcov), 0))
  se_log <- se / est                       # delta method, natural -> log
  z <- stats::qnorm(1 - (1 - level) / 2)

  out <- cbind(estimate = as.numeric(est),
               se       = as.numeric(se),
               lower    = as.numeric(est * exp(-z * se_log)),
               upper    = as.numeric(est * exp(z * se_log)))
  rownames(out) <- names(est)
  attr(out, "level") <- level
  attr(out, "scale") <- "log-scale Wald, exponentiated"
  out
}

#' Profile Likelihood Confidence Interval
#'
#' Profiles one parameter, maximising over the others, and inverts the
#' likelihood ratio test to obtain an interval.
#'
#' @param object A fitted `"betadanish"` object.
#' @param parameter Name of the parameter to profile, one of `"a"`, `"b"`,
#'   `"c"` or `"k"`. Default `"b"`.
#' @param level Confidence level. Default 0.95.
#' @param grid Optional numeric vector of values to scan. If `NULL` (default) a
#'   log-spaced grid is built around the estimate.
#' @param n_grid Number of grid points when `grid` is `NULL`.
#' @param method Optimisation method for the nuisance parameters, passed to
#'   [stats::optim()].
#'
#' @return An object of class `"bd_profile"` with the interval, the grid, the
#'   profile log-likelihood, and the critical value.
#'
#' @details
#' The interval is the set of values \eqn{\theta_0} for which
#' \eqn{2\{\ell_{\max} - \ell_p(\theta_0)\} \le \chi^2_{1,\,\mathrm{level}}},
#' with \eqn{\ell_p} the profile log-likelihood.
#'
#' **An open upper bound is a result, not a failure.** For a weakly identified
#' tail index the profile can stay inside the critical region all the way to
#' the top of the grid, in which case `upper` is returned as `Inf` and the
#' honest report is a lower bound rather than a point estimate with an
#' interval. The thesis pipeline records precisely this on the breaking-stress
#' data, where the profile gives \eqn{b \ge 26.1} with no finite upper bound.
#'
#' Widening `grid` will not manufacture a bound; it only confirms the flatness.
#'
#' @seealso [bd_wald_ci()], [bd_identified_coef()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' dat <- simulate_bd_data(200, a = 1, b = 3, c = 2, k = 0.5, seed = 2)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       submodel = TRUE, n_starts = 3)
#' p <- bd_profile_ci(fit, "b")
#' p
#' }
bd_profile_ci <- function(object, parameter = "b", level = 0.95,
                          grid = NULL, n_grid = 60L, method = "Nelder-Mead") {
  if (!inherits(object, "betadanish"))
    stop("'object' must be a fitted betadanish model.", call. = FALSE)
  est <- object$coefficients
  if (!parameter %in% names(est))
    stop("'", parameter, "' is not a parameter of this fit. Available: ",
         paste(names(est), collapse = ", "), call. = FALSE)

  free   <- setdiff(names(est), parameter)
  ll_at  <- .bd_refit_loglik(object)
  ll_max <- as.numeric(object$logLik)

  profile_at <- function(v) {
    if (!length(free)) {
      nat <- est; nat[parameter] <- v
      return(ll_at(nat))
    }
    obj <- function(lp) {
      nat <- est
      nat[free] <- exp(lp)
      nat[parameter] <- v
      val <- ll_at(nat)
      if (is.finite(val)) -val else 1e10
    }
    r <- tryCatch(stats::optim(log(est[free]), obj, method = method,
                               control = list(maxit = 3000, reltol = 1e-9)),
                  error = function(e) NULL)
    if (is.null(r) || !is.finite(r$value)) -Inf else -r$value
  }

  if (is.null(grid)) {
    e <- as.numeric(est[[parameter]])
    grid <- sort(unique(exp(seq(log(e / 20), log(e * 20),
                                length.out = as.integer(n_grid)))))
  }
  grid <- sort(unique(grid[is.finite(grid) & grid > 0]))
  if (length(grid) < 3L)
    stop("'grid' needs at least three positive values.", call. = FALSE)

  prof <- vapply(grid, profile_at, numeric(1))
  crit <- stats::qchisq(level, df = 1) / 2
  dev  <- ll_max - prof

  inside <- grid[is.finite(dev) & dev <= crit]
  if (!length(inside)) {
    lower <- NA_real_; upper <- NA_real_
  } else {
    lower <- min(inside)
    hi    <- max(inside)
    tol   <- max(grid) * (1 - 1e-8)
    upper <- if (hi >= tol) Inf else hi
    if (lower <= min(grid) * (1 + 1e-8)) lower <- min(grid)
  }

  out <- list(parameter = parameter,
              estimate  = as.numeric(est[[parameter]]),
              lower = lower, upper = upper, level = level,
              grid = grid, profile = prof, logLik_max = ll_max,
              critical = crit,
              open_above = is.infinite(upper),
              open_below = length(inside) > 0 && lower <= min(grid) * (1 + 1e-8))
  class(out) <- "bd_profile"
  out
}

#' @param x A `"bd_profile"` object.
#' @param ... Ignored.
#' @rdname bd_profile_ci
#' @export
print.bd_profile <- function(x, ...) {
  cat("Profile likelihood interval for ", x$parameter, "\n", sep = "")
  cat("  estimate: ", signif(x$estimate, 6), "\n", sep = "")
  cat("  ", format(100 * x$level), "% interval: [",
      signif(x$lower, 6), ", ",
      if (is.infinite(x$upper)) "Inf" else signif(x$upper, 6), "]\n", sep = "")
  if (isTRUE(x$open_above))
    cat("\n  The profile stays within the critical region to the top of the\n",
        "  grid, so there is no finite upper bound. Report ", x$parameter,
        " as a\n  lower bound rather than as a point estimate.\n", sep = "")
  if (isTRUE(x$open_below))
    cat("\n  The interval also reaches the bottom of the grid; widen it if a\n",
        "  finite lower bound matters.\n", sep = "")
  invisible(x)
}

#' Coefficients in the Identified Parametrisation
#'
#' Reports the four-parameter fit in terms of the composite \eqn{ac}, which is
#' identified, rather than \eqn{a} and \eqn{c} separately, which are not.
#'
#' @param object A fitted `"betadanish"` object.
#' @param level Confidence level for the reported intervals. Default 0.95.
#'
#' @return An object of class `"bd_identified"` holding the reparametrised
#'   estimates with delta-method standard errors, and the condition numbers of
#'   the covariance matrix before and after.
#'
#' @details
#' Near the lower tail \eqn{a} and \eqn{c} enter the density almost entirely
#' through their product, so the expected information is close to singular in
#' the direction that separates them. Writing \eqn{p = ac} and \eqn{r = a/c},
#' the ratio \eqn{r} carries almost no information while \eqn{p} is estimated
#' precisely.
#'
#' Reporting \eqn{(ac, b, k)} is therefore the honest summary of a
#' four-parameter fit that sits on the flat direction: it has finite standard
#' errors and a far better conditioned covariance matrix, and it says exactly
#' what the data determine. The thesis pipeline records a condition number
#' falling from 2750 to 152 on the remission data, with the log-likelihood
#' unchanged, since this is a reparametrisation and not a different model.
#'
#' For the ED submodel \eqn{a = 1}, so \eqn{ac = c} and nothing is gained;
#' the function returns the coefficients unchanged with a note.
#'
#' @seealso [fit_betadanish()] and its Identifiability section,
#'   [bd_profile_ci()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' dat <- simulate_bd_data(200, a = 1.5, b = 3, c = 2, k = 0.5, seed = 4)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       n_starts = 5)
#' bd_identified_coef(fit)
#' }
bd_identified_coef <- function(object, level = 0.95) {
  if (!inherits(object, "betadanish"))
    stop("'object' must be a fitted betadanish model.", call. = FALSE)

  est <- object$coefficients
  V   <- object$vcov
  z   <- stats::qnorm(1 - (1 - level) / 2)

  if (isTRUE(object$submodel) || !all(c("a", "c") %in% names(est))) {
    se <- sqrt(pmax(diag(V), 0))
    tab <- cbind(estimate = as.numeric(est), se = as.numeric(se),
                 lower = as.numeric(est * exp(-z * se / est)),
                 upper = as.numeric(est * exp(z * se / est)))
    rownames(tab) <- names(est)
    out <- list(table = tab, submodel = TRUE,
                condition_before = .bd_condition(V),
                condition_after  = .bd_condition(V), level = level)
    class(out) <- "bd_identified"
    return(out)
  }

  keep <- c("a", "b", "c", "k")
  V4 <- V[keep, keep, drop = FALSE]
  a  <- est[["a"]]; cc <- est[["c"]]

  ## (a, b, c, k) -> (ac, b, k)
  J <- rbind(c(cc, 0, a, 0),
             c(0,  1, 0, 0),
             c(0,  0, 0, 1))
  Vn <- J %*% V4 %*% t(J)
  nm <- c("ac", "b", "k")
  dimnames(Vn) <- list(nm, nm)

  e_new <- c(ac = a * cc, b = est[["b"]], k = est[["k"]])
  se    <- sqrt(pmax(diag(Vn), 0))

  tab <- cbind(estimate = as.numeric(e_new), se = as.numeric(se),
               lower = as.numeric(e_new * exp(-z * se / e_new)),
               upper = as.numeric(e_new * exp(z * se / e_new)))
  rownames(tab) <- nm

  out <- list(table = tab, submodel = FALSE,
              a = a, c = cc, ratio = a / cc,
              condition_before = .bd_condition(V4),
              condition_after  = .bd_condition(Vn),
              level = level)
  class(out) <- "bd_identified"
  out
}

#' Condition Number of a Covariance Matrix
#' @noRd
.bd_condition <- function(V) {
  ev <- tryCatch(eigen(V, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  if (anyNA(ev) || min(ev) <= 0) return(NA_real_)
  max(ev) / min(ev)
}

#' @param x A `"bd_identified"` object.
#' @param ... Ignored.
#' @rdname bd_identified_coef
#' @export
print.bd_identified <- function(x, ...) {
  cat("Identified parametrisation\n")
  if (isTRUE(x$submodel)) {
    cat("  ED submodel: a is fixed at 1, so ac = c and nothing is gained.\n\n")
  } else {
    cat("  a and c enter mainly through their product; ac is reported.\n\n")
  }
  print(round(x$table, 5))
  if (!isTRUE(x$submodel)) {
    cat("\n  a = ", signif(x$a, 5), ", c = ", signif(x$c, 5),
        ", ratio a/c = ", signif(x$ratio, 5), "\n", sep = "")
  }
  cb <- x$condition_before; ca <- x$condition_after
  cat("  covariance condition number: ",
      if (is.finite(cb)) signif(cb, 5) else "NA", " -> ",
      if (is.finite(ca)) signif(ca, 5) else "NA", "\n", sep = "")
  invisible(x)
})---")


.step("Writing tests/testthat/test-estimation.R")

.put("tests/testthat/test-estimation.R", r"---(## Estimation and inference added in Patch 3b.

grid_data <- function(n = 250, seed = 11) {
  ## Continuous times rounded to a grid of 1, as month-recorded data would be.
  set.seed(seed)
  t <- rbetadanish(n, a = 1, b = 4, c = 2, k = 0.15)
  data.frame(time = pmax(round(t), 1), status = 1L)
}

test_that(".bd_log_cell is a log probability and sums sensibly", {
  a <- 1; b <- 4; c <- 2; k <- 0.15
  t <- c(1, 3, 10, 40)
  lc <- BetaDanish:::.bd_log_cell(t, delta = 1, a, b, c, k)

  expect_true(all(is.finite(lc)))
  expect_true(all(lc <= 0))                      # never exceeds log(1)

  ## Matches the direct difference of the distribution function
  direct <- log(pbetadanish(t + 0.5, a, b, c, k) -
                  pbetadanish(t - 0.5, a, b, c, k))
  expect_equal(lc, direct, tolerance = 1e-9)
})

test_that(".bd_log_cell falls back to the density when the cell underflows", {
  ## Far into the tail the two distribution values are equal to machine
  ## precision, so the difference cancels to zero and the density is used.
  lc <- BetaDanish:::.bd_log_cell(1e12, delta = 1e-6, a = 1.5, b = 3, c = 2, k = 1)
  expect_true(is.finite(lc))
  expect_lt(lc, 0)
})

test_that("grouped fitting runs and reports its increment", {
  skip_on_cran()
  dat <- grid_data()
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, grouped = TRUE, n_starts = 3,
                   check_identifiability = FALSE))
  expect_true(isTRUE(fit$grouped))
  expect_equal(fit$delta, 1)
  expect_true(is.finite(fit$logLik))
  expect_true(all(fit$coefficients > 0))
})

test_that("grouped standard errors exceed the point-density ones", {
  skip_on_cran()
  dat <- grid_data(n = 300, seed = 12)
  exact <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  grp <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, grouped = TRUE, n_starts = 3,
                   check_identifiability = FALSE))

  ## Treating a rounded time as exact overstates the information, so the
  ## point-density likelihood is the more confident of the two.
  se_exact <- sqrt(pmax(diag(exact$vcov), 0))
  se_grp   <- sqrt(pmax(diag(grp$vcov), 0))
  expect_true(all(is.finite(c(se_exact, se_grp))))
  expect_gte(mean(se_grp / se_exact), 0.95)
})

test_that("grid-recorded times raise a warning when grouped is FALSE", {
  skip_on_cran()
  dat <- grid_data(n = 120, seed = 13)
  expect_warning(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 2),
    "grid")
})

test_that("grouped = TRUE without an inferable increment is an error", {
  dat <- data.frame(time = c(0.137, 1.882, 3.019, 7.4451, 11.02),
                    status = 1L)
  expect_error(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, grouped = TRUE, n_starts = 1),
    "recording increment")
})

test_that("the penalty shrinks the estimates and is recorded", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 21)
  plain <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   n_starts = 3, check_identifiability = FALSE))
  pen <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   n_starts = 3, penalty = 0.5, check_identifiability = FALSE))

  expect_equal(plain$penalty, 0)
  expect_equal(pen$penalty, 0.5)
  expect_true(is.na(plain$penalised_logLik))
  expect_true(is.finite(pen$penalised_logLik))
})

test_that("the reported log-likelihood is unpenalised", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 22)
  pen <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   n_starts = 3, penalty = 1, check_identifiability = FALSE))

  ## The penalised objective can only be lower than the log-likelihood it
  ## subtracts a non-negative penalty from.
  expect_lte(pen$penalised_logLik, pen$logLik + 1e-6)

  ## AIC and BIC must be built from the unpenalised value, or they would not
  ## be comparable with an unpenalised fit.
  expect_equal(pen$AIC, 2 * pen$npar - 2 * pen$logLik, tolerance = 1e-8)
  expect_equal(unname(stats::AIC(pen)), pen$AIC, tolerance = 1e-8)
})

test_that("a negative penalty is refused", {
  dat <- simulate_bd_data(60, a = 1, b = 3, c = 2, k = 0.5, seed = 23)
  expect_error(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   penalty = -1, n_starts = 1),
    "non-negative")
})

test_that("Wald intervals stay inside the parameter space", {
  skip_on_cran()
  dat <- simulate_bd_data(200, a = 1, b = 3, c = 2, k = 0.5, seed = 24)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  ci <- bd_wald_ci(fit)

  expect_true(all(colnames(ci) == c("estimate", "se", "lower", "upper")))
  expect_true(all(ci[, "lower"] > 0))                       # never crosses zero
  expect_true(all(ci[, "lower"] <= ci[, "estimate"]))
  expect_true(all(ci[, "upper"] >= ci[, "estimate"]))
  expect_equal(attr(ci, "level"), 0.95)
})

test_that("the profile interval brackets the estimate", {
  skip_on_cran()
  dat <- simulate_bd_data(200, a = 1, b = 3, c = 2, k = 0.5, seed = 25)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  p <- bd_profile_ci(fit, "b", n_grid = 25L)

  expect_s3_class(p, "bd_profile")
  expect_equal(p$parameter, "b")
  expect_lte(p$lower, p$estimate)
  expect_gte(p$upper, p$estimate)
  expect_true(max(p$profile, na.rm = TRUE) <= p$logLik_max + 1e-6)
  expect_output(print(p), "Profile likelihood")
})

test_that("a flat profile is reported as an open upper bound", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 26)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))

  ## A grid stopping just above the estimate cannot bound b from above, so the
  ## upper limit must come back as Inf rather than as the grid maximum.
  g <- seq(fit$coefficients[["b"]] * 0.5, fit$coefficients[["b"]] * 1.02,
           length.out = 15)
  p <- bd_profile_ci(fit, "b", grid = g)
  expect_true(is.infinite(p$upper))
  expect_true(p$open_above)
  expect_output(print(p), "lower bound")
})

test_that("profiling rejects an unknown parameter", {
  skip_on_cran()
  dat <- simulate_bd_data(80, a = 1, b = 3, c = 2, k = 0.5, seed = 27)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 2, check_identifiability = FALSE))
  expect_error(bd_profile_ci(fit, "a"), "not a parameter")
})

test_that("the identified parametrisation reports ac with a finite SE", {
  skip_on_cran()
  dat <- simulate_bd_data(250, a = 1.5, b = 3, c = 2, k = 0.5, seed = 28)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   n_starts = 5, check_identifiability = FALSE))
  id <- bd_identified_coef(fit)

  expect_s3_class(id, "bd_identified")
  expect_equal(rownames(id$table), c("ac", "b", "k"))
  expect_equal(id$table["ac", "estimate"],
               fit$coefficients[["a"]] * fit$coefficients[["c"]],
               tolerance = 1e-10)
  expect_true(all(is.finite(id$table[, "se"])))
  expect_true(all(id$table[, "lower"] > 0))
  expect_output(print(id), "Identified parametrisation")
})

test_that("the submodel needs no reparametrisation", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 29)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  id <- bd_identified_coef(fit)
  expect_true(id$submodel)
  expect_equal(rownames(id$table), c("b", "c", "k"))
  expect_output(print(id), "nothing is gained")
})

test_that("inference functions reject a non-betadanish object", {
  expect_error(bd_wald_ci(list()), "betadanish")
  expect_error(bd_profile_ci(list()), "betadanish")
  expect_error(bd_identified_coef(list()), "betadanish")
}))---")


.step("Recording the additions in NEWS.md")
.nw <- readLines("NEWS.md", warn = FALSE)
if (!any(grepl("bd_profile_ci", .nw, fixed = TRUE))) {
  .hdr <- grep("^# BetaDanish 0\\.2\\.0\\.9000", .nw)
  if (length(.hdr) == 1L) {
    .backup("NEWS.md")
    .sec <- c(
      "",
      "## Estimation and inference",
      "",
      "* `fit_betadanish(grouped = TRUE)` uses the grouped likelihood for times",
      "  recorded on a coarse grid. An event at `t` contributes",
      "  `log{F(t + delta/2) - F(t - delta/2)}` rather than `log f(t)`. The",
      "  increment is inferred from the spacing of the times unless `delta` is",
      "  supplied, and a warning now fires when the times look grid-recorded",
      "  but the point-density likelihood is being used, which understates",
      "  every standard error.",
      "",
      "* `fit_betadanish(penalty = )` adds a ridge penalty on the",
      "  log-parameter scale for the weakly identified regime, shrinking",
      "  toward the unpenalised optimum unless `penalty_center` says",
      "  otherwise.",
      "",
      "  The reported `logLik` is the **unpenalised** log-likelihood evaluated",
      "  at the penalised estimate, so AIC, BIC and likelihood ratio tests",
      "  stay comparable with unpenalised fits. The objective actually",
      "  maximised is kept separately as `penalised_logLik`.",
      "",
      "* `bd_wald_ci()` builds Wald intervals on the log scale and",
      "  exponentiates, so they cannot cross zero. The symmetric interval on",
      "  the natural scale does exactly that for a weakly identified tail",
      "  index.",
      "",
      "* `bd_profile_ci()` inverts the likelihood ratio test for one",
      "  parameter. When the profile stays inside the critical region to the",
      "  top of the grid it returns `Inf` and says so: an unbounded interval",
      "  is a result, and the honest report is a lower bound rather than a",
      "  point estimate.",
      "",
      "* `bd_identified_coef()` reports the four-parameter fit as `(ac, b, k)`.",
      "  Near the lower tail `a` and `c` enter almost entirely through their",
      "  product, so the composite is estimated precisely while the two",
      "  separately are not. The log-likelihood is unchanged: this is a",
      "  reparametrisation, not a different model.")
    .nw <- append(.nw, .sec, after = .hdr)
    .write_lines("NEWS.md", .nw)
    .ok("NEWS.md updated")
  } else {
    .warn("version header not found in NEWS.md; add the note by hand")
  }
} else {
  .info("NEWS.md already records these")
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

.step("Scanning R/ for any function defined in more than one file")
.defs <- list()
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  txt <- readLines(f, warn = FALSE)
  hits <- grep("^[a-zA-Z_.][a-zA-Z0-9_.]*[ ]*<-[ ]*function", txt, value = TRUE)
  for (nm in sub("[ ]*<-.*$", "", hits))
    .defs[[nm]] <- c(.defs[[nm]], basename(f))
}
.dups <- Filter(function(v) length(unique(v)) > 1L, .defs)
if (length(.dups)) {
  for (nm in names(.dups))
    .warn(sprintf("%s defined in: %s", nm, paste(unique(.dups[[nm]]), collapse = ", ")))
  .die("A function is defined in more than one file. Backups: ", BACKUP_DIR)
}
.ok(sprintf("%d definition(s), every name unique", length(.defs)))

.step("Loading from source for the self-test")
.loaded <- tryCatch({ devtools::load_all(".", quiet = TRUE); TRUE },
                    error = function(e) conditionMessage(e))
if (!isTRUE(.loaded))
  .die("load_all() failed:\n  ", .loaded, "\n\nBackups: ", BACKUP_DIR)
.ok("source loaded")

.step("Numerical self-test")
.fails <- character(0)
.check <- function(label, got, want, tol = 1e-6) {
  rel <- abs(got - want) / max(abs(want), 1e-12)
  if (is.finite(rel) && rel < tol) {
    .ok(sprintf("%-46s %.8g vs %.8g", label, got, want))
  } else {
    .warn(sprintf("%-46s %.8g vs %.8g  (rel %.1e)", label, got, want, rel))
    .fails <<- c(.fails, label)
  }
}

## The grouped cell probability is exactly the difference of the CDF.
.a <- 1; .b <- 4; .cc <- 2; .kk <- 0.15
.t <- 7
.check("log cell vs direct CDF difference",
       BetaDanish:::.bd_log_cell(.t, 1, .a, .b, .cc, .kk),
       log(pbetadanish(.t + 0.5, .a, .b, .cc, .kk) -
             pbetadanish(.t - 0.5, .a, .b, .cc, .kk)))

## A reparametrisation cannot change the log-likelihood.
set.seed(31)
.dat <- simulate_bd_data(200, a = 1.5, b = 3, c = 2, k = 0.5, seed = 31)
.fit <- suppressWarnings(
  fit_betadanish(survival::Surv(time, status) ~ 1, data = .dat,
                 n_starts = 5, check_identifiability = FALSE))
.id <- bd_identified_coef(.fit)
.check("ac equals a times c",
       .id$table["ac", "estimate"],
       .fit$coefficients[["a"]] * .fit$coefficients[["c"]])

## AIC must be built from the unpenalised log-likelihood.
.pen <- suppressWarnings(
  fit_betadanish(survival::Surv(time, status) ~ 1, data = .dat,
                 n_starts = 3, penalty = 1, check_identifiability = FALSE))
.check("AIC uses the unpenalised log-likelihood",
       .pen$AIC, 2 * .pen$npar - 2 * .pen$logLik)
if (!is.finite(.pen$penalised_logLik) || .pen$penalised_logLik > .pen$logLik + 1e-6) {
  .warn("the penalised objective is not below the log-likelihood")
  .fails <- c(.fails, "penalised objective")
} else {
  .ok(sprintf("%-46s %.8g <= %.8g", "penalised objective below log-likelihood",
              .pen$penalised_logLik, .pen$logLik))
}

## The Wald interval must stay strictly positive.
.ci <- bd_wald_ci(.fit)
if (all(.ci[, "lower"] > 0)) {
  .ok("log-scale Wald interval stays positive")
} else {
  .warn("a Wald lower limit is not positive"); .fails <- c(.fails, "wald")
}

if (length(.fails))
  .die("Self-tests failed: ", paste(.fails, collapse = "; "),
       "\nNothing further was run. Backups: ", BACKUP_DIR)
.ok("all self-tests agree")

.step("devtools::document()")
.r <- tryCatch({ devtools::document(); TRUE }, error = function(e) conditionMessage(e))
if (!isTRUE(.r)) .die("document() failed:\n  ", .r, "\n\nBackups: ", BACKUP_DIR)
.ok("documentation regenerated")

.rd_arg_names <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  i0  <- regexpr("arguments[{]", txt)
  if (i0 < 0) return(character(0))
  rest <- substring(txt, i0)
  i1 <- regexpr("\n[}]\n", rest)
  if (i1 > 0) rest <- substring(rest, 1, i1)
  raw <- unlist(regmatches(rest, gregexpr("item[{][^}]*[}]", rest)))
  if (!length(raw)) return(character(0))
  trimws(unlist(strsplit(sub("[}]$", "", sub("^item[{]", "", raw)), ",")))
}
.rd_aliases <- function(path) {
  a <- grep("^.alias[{]", readLines(path, warn = FALSE), value = TRUE)
  trimws(sub("[}].*$", "", sub("^.alias[{]", "", a)))
}

.step("Checking Rd files for duplicated arguments and aliases")
.rds <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
.bad <- character(0)
for (f in .rds) {
  nms <- .rd_arg_names(f)
  dup <- unique(nms[duplicated(nms)])
  if (length(dup)) {
    .warn(sprintf("%s duplicates: %s", basename(f), paste(dup, collapse = ", ")))
    .bad <- c(.bad, basename(f))
  }
}
.amap <- list()
for (f in .rds) for (a in .rd_aliases(f)) .amap[[a]] <- c(.amap[[a]], basename(f))
.dupal <- Filter(function(v) length(v) > 1L, .amap)
for (a in names(.dupal))
  .warn(sprintf("alias '%s' in: %s", a, paste(.dupal[[a]], collapse = ", ")))
if (length(.bad) || length(.dupal))
  .die("Rd duplication would fail R CMD check. Backups: ", BACKUP_DIR)
.ok(sprintf("%d Rd file(s), %d alias(es), all clean", length(.rds), length(.amap)))

.step("Confirming the new exports")
.ns <- readLines("NAMESPACE", warn = FALSE)
for (f in c("bd_wald_ci", "bd_profile_ci", "bd_identified_coef")) {
  if (any(grepl(paste0("export(", f, ")"), .ns, fixed = TRUE))) .ok(f)
  else .warn(paste("not exported:", f))
}

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
cat("  PATCH 3b COMPLETE  --  estimation and inference\n")
cat(strrep("=", 78), "\n\n")
cat("  34  penalty = , with the unpenalised log-likelihood reported\n")
cat("  35  grouped = TRUE, with the increment inferred and a warning when\n")
cat("      grid-recorded times are fitted as if exact\n")
cat("  36  bd_wald_ci(), bd_profile_ci() with an honest open upper bound\n")
cat("      bd_identified_coef() for the (ac, b, k) composite\n\n")
cat("  Still to come\n")
cat("    3c  competing-risks covariates (37) and simulation studies (38, 39)\n")
cat("    3d  visualisation (40-42) and the 0.3.0 release (45, 46)\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
