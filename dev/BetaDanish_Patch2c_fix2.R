## =============================================================================
##  BetaDanish  --  PATCH 2c-fix2 : make the diagnostic guards absent-safe
## =============================================================================
##
##  THE GUARD ITSELF IS WORKING. In the Patch 2c-fix run the analyze-csv block
##  changed from
##
##      Full Model (4-param)  5.210644e+77   Chisq 1.042129e+78   p = 0
##  to
##      Full Model (4-param)      -476.8225   Chisq 0             p = 1
##
##  which is the honest result: on the transplant data the fourth parameter
##  buys nothing, exactly as the flat (a, c) direction predicts.
##
##  WHAT BROKE
##    .bd_warn_diagnostics() gained two guards written as
##
##        if (!is.na(d$starts_rejected) && d$starts_rejected > 0)
##
##    That is fine for a real fit, which always carries the field. But
##    test-identifiability.R builds a diagnostics list by hand, without it.
##    is.na(NULL) is logical(0), and `logical(0) && logical(0)` yields NA, so
##    the `if` failed with "missing value where TRUE/FALSE needed". Five tests
##    errored -- none of them about the degeneracy guard, all of them about my
##    guard style.
##
##    The pre-existing ac_correlation guard had the same latent fault:
##    is.finite(NULL) is also logical(0). It never fired only because the test
##    helper happened to supply that field.
##
##  THE FIX
##    Every guard in .bd_warn_diagnostics() now uses isTRUE(), which returns
##    FALSE for NULL, NA and zero-length alike. A fit object saved before the
##    degeneracy guard existed simply reports nothing about it, instead of
##    erroring. A regression test covers exactly that shape.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch2c_fix2.R")
##  IDEMPOTENT   Yes -- both files are written whole.
## =============================================================================

if (getRversion() < "4.0.0") stop("This patch needs R >= 4.0.")

.step_n <- 0L
.step <- function(m) { .step_n <<- .step_n + 1L; cat(sprintf("\n[%02d] %s\n", .step_n, m)) }
.ok   <- function(m) cat("     OK   ", m, "\n", sep = "")
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
cat("  BetaDanish  --  Patch 2c-fix2 : absent-safe diagnostic guards\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!any(grepl(".bd_make_accept", readLines("R/fit_models.R", warn = FALSE), fixed = TRUE)))
  .die("Patch 2c-fix has not been applied -- R/fit_models.R is not guarded.")
.ok("Patch 2c-fix detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch2cfix2"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

.step("Writing R/fit_models.R with isTRUE() guards")

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

  if (isTRUE(abs(d$ac_correlation) > 0.95))
    warning(sprintf(paste0("The fitted correlation between a-hat and c-hat is ",
                           "%.3f, so effectively only the product c*a is ",
                           "identified. Individual estimates of a and c should ",
                           "not be interpreted."),
                    d$ac_correlation), call. = FALSE)

  if (isTRUE(d$starts_rejected > 0))
    warning(sprintf(paste0("%d starting point(s) reached a degenerate ridge ",
                           "and were discarded. The reported fit is the best ",
                           "admissible optimum. If this is most of the grid, ",
                           "the four-parameter model is a poor choice for ",
                           "these data."), d$starts_rejected), call. = FALSE)

  if (isTRUE(d$loglik_spread > 2))
    warning(sprintf(paste0("Accepted optima span %.2f log-likelihood units, ",
                           "so the surface has several local maxima and the ",
                           "reported fit may not be global. Increase ",
                           "n_starts."), d$loglik_spread), call. = FALSE)

  invisible(NULL)
})---")


.step("Updating tests/testthat/test-identifiability.R")

.put("tests/testthat/test-identifiability.R", r"---(## The identifiability diagnostics are tested directly rather than through a
## fitted model, so the test does not depend on where the optimiser lands.

mk_diag <- function(...) {
  base <- list(converged = TRUE, convergence_code = 1L, vcov_singular = FALSE,
               near_b_ridge = FALSE, b_distance_se = 5, ac_correlation = 0.2,
               starts_ok = 12L, starts_rejected = 0L, loglik_spread = 0.1,
               loglik_per_obs = -2.5)
  utils::modifyList(base, list(...))
}

test_that("the b = 1 ridge warning fires when b-hat sits close to one", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(near_b_ridge = TRUE,
                                              b_distance_se = 0.6)),
    "non-identifiability ridge")
})

test_that("a singular information matrix is reported", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(vcov_singular = TRUE)),
    "singular")
})

test_that("(a, c) confounding is reported above the 0.95 threshold", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(ac_correlation = 0.99)),
    "only the product")
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(ac_correlation = -0.98)),
    "only the product")
})

test_that("a poor convergence code is reported", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(converged = FALSE,
                                              convergence_code = 4L)),
    "code 4")
})

test_that("a well-behaved fit produces no diagnostic warnings", {
  expect_silent(BetaDanish:::.bd_warn_diagnostics(mk_diag()))
})

test_that("check_identifiability = FALSE suppresses the warnings", {
  skip_on_cran()
  dat <- simulate_bd_data(120, a = 1, b = 1.1, c = 2, k = 0.5, seed = 99)
  fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                        n_starts = 3, check_identifiability = FALSE)
  expect_s3_class(fit, "betadanish")
  expect_true(is.list(fit$diagnostics))
})

test_that("zero-length input gives zero-length output for every function", {
  z <- numeric(0)
  expect_length(dbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(pbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(qbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(sbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(hbetadanish(z, 1.5, 3, 2, 1), 0L)
  ## A zero-length parameter also collapses the result, as in stats::dnorm.
  expect_length(dbetadanish(1, a = z, b = 3, c = 2, k = 1), 0L)
})

test_that("degenerate starts and local-optima spread are reported", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(starts_rejected = 4L)),
    "degenerate ridge")
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(loglik_spread = 7.5)),
    "several local maxima")
  expect_silent(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(starts_rejected = 0L,
                                              loglik_spread = 0.4)))
})

test_that("a diagnostics list missing the newer fields is tolerated", {
  ## Objects saved before the degeneracy guard existed have no starts_rejected
  ## or loglik_spread. Every guard uses isTRUE(), so an absent field is simply
  ## not reported rather than raising "missing value where TRUE/FALSE needed".
  legacy <- list(converged = TRUE, convergence_code = 1L,
                 vcov_singular = FALSE, near_b_ridge = FALSE,
                 b_distance_se = 5, ac_correlation = 0.2)
  expect_silent(BetaDanish:::.bd_warn_diagnostics(legacy))

  bare <- list(converged = TRUE)
  expect_silent(BetaDanish:::.bd_warn_diagnostics(bare))
}))---")


.step("Confirming no fragile guards remain")
.fm <- readLines("R/fit_models.R", warn = FALSE)
.blk <- .fm[seq(grep("^.bd_warn_diagnostics <- function", .fm)[1], length(.fm))]
.guards <- grep("if \\(", .blk, value = TRUE)
.fragile <- grep("isTRUE", .guards, value = TRUE, invert = TRUE)
for (g in .guards) cat("       ", trimws(g), "\n", sep = "")
if (length(.fragile)) {
  .warn("guard(s) not using isTRUE():")
  for (g in .fragile) cat("        ", trimws(g), "\n", sep = "")
} else {
  .ok(sprintf("all %d guard(s) are absent-safe", length(.guards)))
}

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
cat("  PATCH 2c-fix2 COMPLETE\n")
cat(strrep("=", 78), "\n\n")
cat("  All six guards in .bd_warn_diagnostics() now use isTRUE(), so a\n")
cat("  diagnostics list missing any field is tolerated rather than fatal.\n")
cat("  Two new tests cover the degeneracy warnings and the legacy shape.\n\n")
cat("  Expected: 0 errors, 0 warnings, 1 note (the clock note).\n")
cat("  This completes Patch 2. Patch 3 is the last one.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
