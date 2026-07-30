## =============================================================================
##  BetaDanish  --  PATCH 2c-fix4 : same code, self-test now reads the disk
## =============================================================================
##
##  fix3 wrote the right file and then tested the wrong function.
##
##  Step 02 wrote R/fit_models.R. Step 04 called
##  BetaDanish:::.bd_warn_diagnostics(), which resolves against the INSTALLED
##  namespace, not the source just written. The package had not been reloaded,
##  so the self-test exercised the previous version and aborted on its faults.
##
##  The proof is in fix3's own output. The five shapes that failed --
##  "converged flag only", "NULL fields", "character field", "empty list",
##  "not a list" -- are exactly the five with no length-one numeric
##  ac_correlation, which is precisely what the fix2 guard could not handle.
##  The four that passed are exactly the four it could. Nine out of nine match
##  the old function, not the new one.
##
##  THIS PATCH
##    Identical file contents. The only change is devtools::load_all() before
##    the self-test, so it exercises the code on disk. It also aborts if
##    load_all() fails, rather than silently testing a stale namespace again.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch2c_fix4.R")
##  IDEMPOTENT   Yes -- re-writing the same content is harmless.
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
cat("  BetaDanish  --  Patch 2c-fix4\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!any(grepl(".bd_make_accept", readLines("R/fit_models.R", warn = FALSE), fixed = TRUE)))
  .die("R/fit_models.R is not guarded -- run Patch 2c-fix first.")
.ok("guarded fit_models.R detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch2cfix4"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

.step("Writing R/fit_models.R with normalised diagnostics")

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

test_that("an incomplete diagnostics list degrades quietly", {
  ## Objects saved before a later diagnostic existed carry no starts_rejected
  ## or loglik_spread. Every field is normalised before it is tested, so an
  ## absent one is simply not reported.
  legacy <- list(converged = TRUE, convergence_code = 1L,
                 vcov_singular = FALSE, near_b_ridge = FALSE,
                 b_distance_se = 5, ac_correlation = 0.2)
  expect_silent(BetaDanish:::.bd_warn_diagnostics(legacy))

  ## Only the convergence flag survives.
  expect_silent(BetaDanish:::.bd_warn_diagnostics(list(converged = TRUE)))

  ## Fields present but of the wrong shape or type.
  expect_silent(BetaDanish:::.bd_warn_diagnostics(
    list(converged = TRUE, ac_correlation = NULL,
         starts_rejected = character(0), loglik_spread = "n/a")))
  expect_silent(BetaDanish:::.bd_warn_diagnostics(
    list(converged = TRUE, ac_correlation = c(0.1, 0.2))))

  ## An empty list has no convergence flag, so it warns -- but it must warn,
  ## not error.
  expect_warning(BetaDanish:::.bd_warn_diagnostics(list()), "optimiser")

  ## A non-list is ignored outright.
  expect_silent(BetaDanish:::.bd_warn_diagnostics(NULL))
}))---")


.step("Loading the package from source, so the self-test reads the disk")
.loaded <- tryCatch({ devtools::load_all(".", quiet = TRUE); TRUE },
                    error = function(e) conditionMessage(e))
if (!isTRUE(.loaded))
  .die("load_all() failed:\n  ", .loaded,
       "\nWithout it the self-test would exercise the installed namespace ",
       "rather than the file just written.\n\nBackups: ", BACKUP_DIR)
.ok("source loaded")

.step("Exercising the warner against every awkward shape")
.warner <- get(".bd_warn_diagnostics", envir = asNamespace("BetaDanish"))
.shapes <- list(
  "complete"          = list(converged = TRUE, convergence_code = 1L,
                             vcov_singular = FALSE, near_b_ridge = FALSE,
                             b_distance_se = 5, ac_correlation = 0.2,
                             starts_ok = 12L, starts_rejected = 0L,
                             loglik_spread = 0.1),
  "legacy (no guard fields)" = list(converged = TRUE, convergence_code = 1L,
                             vcov_singular = FALSE, near_b_ridge = FALSE,
                             b_distance_se = 5, ac_correlation = 0.2),
  "converged flag only" = list(converged = TRUE),
  "NULL fields"        = list(converged = TRUE, ac_correlation = NULL,
                              starts_rejected = NULL),
  "zero-length fields" = list(converged = TRUE, ac_correlation = numeric(0),
                              loglik_spread = character(0)),
  "character field"    = list(converged = TRUE, loglik_spread = "n/a"),
  "length-two field"   = list(converged = TRUE, ac_correlation = c(0.1, 0.9)),
  "empty list"         = list(),
  "not a list"         = NULL
)
.fail <- character(0)
for (nm in names(.shapes)) {
  r <- tryCatch({ suppressWarnings(.warner(.shapes[[nm]])); "ok" },
                error = function(e) conditionMessage(e))
  if (identical(r, "ok")) .ok(nm) else { .warn(paste0(nm, ": ", r)); .fail <- c(.fail, nm) }
}
if (length(.fail))
  .die("The warner errors on: ", paste(.fail, collapse = ", "),
       "\n\nBackups: ", BACKUP_DIR)
.ok("all nine shapes handled without error")

.step("Confirming the degeneracy warnings still fire when they should")
## A function that never errors because it never does anything would sail
## through the shape test above, so check the warnings are still raised.
.catch_warning <- function(expr) {
  msg <- NULL
  tryCatch(
    withCallingHandlers(expr,
                        warning = function(w) {
                          msg <<- conditionMessage(w)
                          invokeRestart("muffleWarning")
                        }),
    error = function(e) NULL)
  msg
}

.probe <- list(
  "degenerate starts" = list(
    d = list(converged = TRUE, starts_rejected = 4L), want = "degenerate ridge"),
  "local optima"      = list(
    d = list(converged = TRUE, loglik_spread = 7.5), want = "several local maxima"),
  "b = 1 ridge"       = list(
    d = list(converged = TRUE, near_b_ridge = TRUE, b_distance_se = 0.6),
    want = "non-identifiability ridge")
)
for (nm in names(.probe)) {
  got <- .catch_warning(.warner(.probe[[nm]]$d))
  if (!is.null(got) && grepl(.probe[[nm]]$want, got, fixed = TRUE)) {
    .ok(paste(nm, "warning fires"))
  } else {
    .warn(paste0(nm, ": expected a warning matching '",
                 .probe[[nm]]$want, "'; got ",
                 if (is.null(got)) "none" else shQuote(got)))
  }
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
cat("  PATCH 2c-fix4 COMPLETE\n")
cat(strrep("=", 78), "\n\n")
cat("  Expected: 0 errors, 0 warnings, 1 note (the clock note).\n")
cat("  That closes Patch 2. Patch 3 is the last one.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
