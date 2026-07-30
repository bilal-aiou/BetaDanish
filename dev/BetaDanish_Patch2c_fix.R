## =============================================================================
##  BetaDanish  --  PATCH 2c-fix : complete the degeneracy guard
## =============================================================================
##
##  WHAT WENT WRONG IN PATCH 2c
##    Five of its seven substitutions reported "pattern not found". The cause
##    was in my patch helper, not in your package. `.sub_in()` did:
##
##        txt <- readLines(path)                 # a vector of LINES
##        grepl(from, txt, fixed = TRUE)         # matched against each line
##
##    Four of the patterns spanned several lines, so they could never match a
##    single element of that vector. The two single-line ones (aft_models.R,
##    cure_models.R) applied correctly, which is why those two are already
##    done. My pre-flight check tested whole-file matching -- the right
##    semantics, but not what the helper implemented -- so it passed.
##
##  WHAT THIS PATCH DOES
##    No pattern matching at all. R/fit_models.R and R/report_betadanish.R are
##    written out complete, with the Patch 2c changes already applied and
##    verified offline. Wholesale writes are idempotent and cannot half-apply.
##
##  CURRENT STATE IT EXPECTS
##    Applied by Patch 2c:  R/utils-internal.R, R/aft_models.R,
##                          R/cure_models.R, tests/testthat/test-degeneracy.R,
##                          NEWS.md
##    Still to do (this patch):  R/fit_models.R, R/report_betadanish.R
##
##  ALSO
##    The transplant regression test in Patch 2c passed even though the guard
##    was absent, purely because seed 99 happened not to run away. It is
##    replaced here by a sweep over the seeds that did fail, including the ones
##    from the analyze-csv suite where the +5.2e77 result actually appeared.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch2c_fix.R")
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

cat(strrep("=", 78), "\n")
cat("  BetaDanish  --  Patch 2c-fix : complete the degeneracy guard\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
.u <- readLines("R/utils-internal.R", warn = FALSE)
if (!any(grepl(".bd_make_accept", .u, fixed = TRUE)))
  .die("Patch 2c has not been applied -- R/utils-internal.R has no .bd_make_accept().")
if (!any(grepl(".bd_default_starts", .u, fixed = TRUE)))
  .die("R/utils-internal.R has no .bd_default_starts(); re-run Patch 2c first.")
.ok("Patch 2c groundwork detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch2cfix"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

.step("Writing R/fit_models.R complete")

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
                           method = "BFGS", check_identifiability = TRUE) {

  surv_data <- extract_surv_data(formula, data)
  time   <- surv_data$time
  status <- surv_data$status

  ll_fun <- function(pars) {
    a_par <- if (submodel) 1.0 else exp(pars[["log_a"]])
    b_par <- exp(pars[["log_b"]])
    c_par <- exp(pars[["log_c"]])
    k_par <- exp(pars[["log_k"]])

    lp <- suppressWarnings(
      dbetadanish(time, a_par, b_par, c_par, k_par, log = TRUE))
    ls <- suppressWarnings(
      pbetadanish(time, a_par, b_par, c_par, k_par,
                  lower.tail = FALSE, log.p = TRUE))

    loglik <- sum(status * lp + (1 - status) * ls)
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

  out <- list(
    coefficients = est_nat,
    logLik       = fit$maximum,
    vcov         = vcov_nat,
    npar         = npar,
    nobs         = nobs,
    nevent       = sum(status == 1),
    AIC          = 2 * npar - 2 * fit$maximum,
    BIC          = npar * log(nobs) - 2 * fit$maximum,
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
#' @noRd
.bd_warn_diagnostics <- function(d) {
  if (!isTRUE(d$converged))
    warning("The optimiser reported code ", d$convergence_code,
            "; treat the estimates as provisional and increase n_starts.",
            call. = FALSE)

  if (isTRUE(d$vcov_singular))
    warning("The observed information matrix is singular or not positive ",
            "definite, so standard errors are unreliable. This usually means ",
            "the likelihood is flat in at least one direction.", call. = FALSE)

  if (isTRUE(d$near_b_ridge))
    warning(sprintf(paste0("b-hat is only %.2f standard errors from 1, close to ",
                           "the b = 1 non-identifiability ridge. Consider the ",
                           "ED submodel (submodel = TRUE). See the ",
                           "Identifiability section of ?fit_betadanish."),
                    d$b_distance_se), call. = FALSE)

  if (is.finite(d$ac_correlation) && abs(d$ac_correlation) > 0.95)
    warning(sprintf(paste0("The fitted correlation between a-hat and c-hat is ",
                           "%.3f, so effectively only the product c*a is ",
                           "identified. Individual estimates of a and c should ",
                           "not be interpreted."),
                    d$ac_correlation), call. = FALSE)

  if (!is.na(d$starts_rejected) && d$starts_rejected > 0)
    warning(sprintf(paste0("%d starting point(s) reached a degenerate ridge ",
                           "and were discarded. The reported fit is the best ",
                           "admissible optimum. If this is most of the grid, ",
                           "the four-parameter model is a poor choice for ",
                           "these data."), d$starts_rejected), call. = FALSE)

  if (!is.na(d$loglik_spread) && d$loglik_spread > 2)
    warning(sprintf(paste0("Accepted optima span %.2f log-likelihood units, ",
                           "so the surface has several local maxima and the ",
                           "reported fit may not be global. Increase ",
                           "n_starts."), d$loglik_spread), call. = FALSE)

  invisible(NULL)
})---")


.step("Writing R/report_betadanish.R complete")

.put("R/report_betadanish.R", r"---(#' Create a Compact Report from a Beta-Danish Model Fit
#'
#' Collects the headline quantities from a fitted model into a small object with
#' a `print` method. Information criteria are computed from the fitted
#' log-likelihood via the `logLik` method, so they cannot fall out of step with
#' the fit.
#'
#' @param fit A fitted `"betadanish"` object.
#'
#' @return An object of class `"betadanish_report"`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- simulate_bd_data(120, a = 1, b = 3, c = 2, k = 0.5)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       submodel = TRUE)
#' report_betadanish(fit)
#' }
#'
#' @export
report_betadanish <- function(fit) {
  if (is.null(fit) || !inherits(fit, "betadanish"))
    stop("'fit' must be a fitted betadanish object.", call. = FALSE)

  out <- list(
    call         = fit$call,
    coefficients = fit$coefficients,
    submodel     = isTRUE(fit$submodel),
    logLik       = as.numeric(fit$logLik),
    npar         = .bd_or(fit$npar, length(fit$coefficients)),
    nobs         = .bd_or(fit$nobs, length(fit$data$time)),
    AIC          = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
    BIC          = tryCatch(stats::BIC(fit), error = function(e) NA_real_),
    convergence  = fit$convergence,
    diagnostics  = fit$diagnostics
  )

  class(out) <- "betadanish_report"
  out
}

## Not `%||%`: base R gained that operator in 4.4.0, and defining it here would
## mask it for anyone attaching the package.
#' @noRd
.bd_or <- function(x, y) if (is.null(x)) y else x

#' @param x A `"betadanish_report"` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @rdname report_betadanish
#' @export
print.betadanish_report <- function(x, ...) {
  cat("Beta-Danish Model Report\n")
  cat("------------------------\n")
  cat("Model:          ",
      if (x$submodel) "3-parameter ED submodel (a = 1)" else "4-parameter Beta-Danish",
      "\n", sep = "")
  cat("Observations:   ", x$nobs, "  Parameters: ", x$npar, "\n", sep = "")
  cat("Log-likelihood: ", format(round(x$logLik, 4), nsmall = 4), "\n", sep = "")
  cat("AIC:            ", format(round(x$AIC, 4), nsmall = 4), "\n", sep = "")
  cat("BIC:            ", format(round(x$BIC, 4), nsmall = 4), "\n", sep = "")
  cat("Convergence:    ", x$convergence, "\n", sep = "")
  cat("\nEstimates:\n")
  print(round(x$coefficients, 4))

  d <- x$diagnostics
  if (!is.null(d)) {
    flags <- character(0)
    if (isTRUE(d$vcov_singular)) flags <- c(flags, "singular information matrix")
    if (isTRUE(d$near_b_ridge))  flags <- c(flags, "near the b = 1 ridge")
    if (is.finite(d$ac_correlation) && abs(d$ac_correlation) > 0.95)
      flags <- c(flags, "(a, c) confounded")
    if (!is.null(d$starts_rejected) && !is.na(d$starts_rejected) &&
        d$starts_rejected > 0)
      flags <- c(flags, sprintf("%d degenerate start(s) discarded",
                                d$starts_rejected))
    if (!is.null(d$loglik_spread) && !is.na(d$loglik_spread) &&
        d$loglik_spread > 2)
      flags <- c(flags, sprintf("local optima span %.2f log-lik units",
                                d$loglik_spread))
    if (length(flags))
      cat("\nDiagnostic flags: ", paste(flags, collapse = "; "),
          "\n  See the Identifiability section of ?fit_betadanish.\n", sep = "")
  }
  invisible(x)
})---")


.step("Replacing the seed-dependent transplant test")

.put("tests/testthat/test-degeneracy-transplant.R", r"---(
## The Patch 2c version of this test used a single seed and passed even with
## the guard absent, which made it worthless as a regression test. The observed
## +5.2e77 runaway appeared under the analyze-csv seeds, so those are swept
## here alongside the original.

test_that("the four-parameter fit never runs away on the transplant data", {
  skip_on_cran()
  data(transplant, package = "BetaDanish", envir = environment())
  n <- nrow(transplant)

  for (s in c(11, 12, 14, 17, 99, 101, 202)) {
    set.seed(s)
    fit <- suppressWarnings(
      fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                     submodel = FALSE, n_starts = 3,
                     check_identifiability = FALSE))

    ## The bug gave logLik / n = 5.7e75. Any real fit is far below this.
    expect_lt(fit$logLik / n, 5, label = paste("seed", s, "loglik per obs"))
    expect_gt(fit$logLik, -1e4, label = paste("seed", s, "loglik floor"))

    shapes <- fit$coefficients[intersect(names(fit$coefficients),
                                         c("a", "b", "c"))]
    expect_true(all(shapes < 500), label = paste("seed", s, "shapes bounded"))
    expect_true(all(is.finite(fit$coefficients)),
                label = paste("seed", s, "finite estimates"))
  }
})

test_that("the guard is actually in force, not merely unexercised", {
  skip_on_cran()
  data(transplant, package = "BetaDanish", envir = environment())
  set.seed(11)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                   submodel = FALSE, n_starts = 3,
                   check_identifiability = FALSE))

  ## These fields only exist once fit_models.R routes through the guarded
  ## multi-start, so they double as proof the patch applied.
  expect_true(is.numeric(fit$starts_ok))
  expect_gte(fit$starts_ok, 1L)
  expect_true(is.numeric(fit$starts_rejected))
  expect_true(is.numeric(fit$loglik_spread))
  expect_gte(fit$loglik_spread, 0)
  expect_true(is.numeric(fit$diagnostics$loglik_per_obs))
})

test_that("a likelihood ratio test on real data is no longer absurd", {
  skip_on_cran()
  data(transplant, package = "BetaDanish", envir = environment())
  set.seed(11)
  full <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                   submodel = FALSE, n_starts = 3, check_identifiability = FALSE))
  sub <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))

  chisq <- 2 * (full$logLik - sub$logLik)
  expect_true(is.finite(chisq))
  expect_gt(chisq, -1)
  expect_lt(chisq, 100)
})
)---")

.step("Removing the superseded tests from test-degeneracy.R")

.dg <- readLines("tests/testthat/test-degeneracy.R", warn = FALSE)
.txt <- paste(.dg, collapse = "\n")
.drop <- c("the four-parameter fit no longer runs away on the transplant data",
           "a likelihood ratio test on real data is no longer absurd",
           "multi-start diagnostics are recorded on the fit")
.removed <- 0L
for (nm in .drop) {
  pat <- paste0('test_that\\("', nm, '"')
  st <- grep(pat, .dg)
  if (!length(st)) next
  en <- st[1]
  depth <- 0L; started <- FALSE
  for (i in seq(st[1], length(.dg))) {
    depth <- depth + lengths(regmatches(.dg[i], gregexpr("\\{", .dg[i]))) -
                     lengths(regmatches(.dg[i], gregexpr("\\}", .dg[i])))
    if (!started && depth > 0L) started <- TRUE
    if (started && depth == 0L) { en <- i; break }
  }
  while (en < length(.dg) && !nzchar(.dg[en + 1L])) en <- en + 1L
  .dg <- .dg[-(st[1]:en)]
  .removed <- .removed + 1L
}
if (.removed) {
  .backup("tests/testthat/test-degeneracy.R")
  con <- file("tests/testthat/test-degeneracy.R", open = "wb")
  writeLines(.dg, con = con, sep = "\n"); close(con)
  .ok(sprintf("moved %d test block(s) into the new file", .removed))
} else {
  .info("already moved")
}

## =============================================================================
##  VERIFY
## =============================================================================

.step("Confirming the guard is wired into all three fitters")
.all_ok <- TRUE
for (f in c("R/fit_models.R", "R/aft_models.R", "R/cure_models.R")) {
  if (any(grepl(".bd_make_accept", readLines(f, warn = FALSE), fixed = TRUE))) {
    .ok(paste(f, "guarded"))
  } else {
    .warn(paste(f, "-- guard NOT detected")); .all_ok <- FALSE
  }
}
if (!any(grepl(".bd_default_starts", readLines("R/fit_models.R", warn = FALSE), fixed = TRUE))) {
  .warn("R/fit_models.R does not use the deterministic start grid"); .all_ok <- FALSE
} else {
  .ok("R/fit_models.R uses the deterministic start grid")
}
for (fld in c("starts_ok", "starts_rejected", "loglik_spread", "loglik_per_obs")) {
  if (any(grepl(fld, readLines("R/fit_models.R", warn = FALSE), fixed = TRUE)))
    .ok(paste("diagnostic field:", fld))
  else { .warn(paste("missing field:", fld)); .all_ok <- FALSE }
}
if (!.all_ok) .die("The guard did not wire up completely. Backups: ", BACKUP_DIR)

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
cat("  PATCH 2c-fix COMPLETE\n")
cat(strrep("=", 78), "\n\n")
cat("  R/fit_models.R        deterministic grid + acceptance predicate,\n")
cat("                        starts_ok / starts_rejected / loglik_spread,\n")
cat("                        degeneracy and local-optima warnings\n")
cat("  R/report_betadanish.R degeneracy flags surfaced in the report\n")
cat("  tests                 transplant regression swept over seven seeds\n")
cat("                        instead of relying on one lucky one\n\n")
cat("  The analyze-csv suite should no longer print a +5.2e77 log-likelihood.\n")
cat("  Watch that block in the test output.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
