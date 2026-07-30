## =============================================================================
##  BetaDanish  --  PHASE 2, PATCH 1 of 3 : CORRECTNESS AND DOCUMENTATION
## =============================================================================
##
##  Implements approved recommendations 1-18, 43, 44 and the dev-version part
##  of 45 from the Phase 1 audit.
##
##  HOW TO RUN
##    1. Open the BetaDanish *git working copy* (not the extracted tarball) in
##       RStudio, so that Patch 3 can commit and tag.
##    2. setwd() to the package root -- the directory containing DESCRIPTION.
##    3. source("BetaDanish_Patch1_correctness.R")
##       (or paste the whole file into the console; source() is safer)
##
##  REQUIREMENTS   R >= 4.0 (raw strings), devtools, roxygen2
##  IDEMPOTENT     Yes. Every write is a full-file replacement or a guarded
##                 targeted substitution, so re-running changes nothing.
##  BACKUP         A timestamped copy of every file it touches is written to
##                 .betadanish_backup/<timestamp>/ before any change.
##  SAFETY         Nothing is deleted. No dataset is touched (that is Patch 2).
##                 No git operation is performed.
## =============================================================================

if (getRversion() < "4.0.0")
  stop("This patch needs R >= 4.0 (it uses raw string literals).")

PATCH_ID <- "patch1-correctness"

## ---------------------------------------------------------------- helpers ----

.step_n <- 0L
.step <- function(msg) {
  .step_n <<- .step_n + 1L
  cat(sprintf("\n[%02d] %s\n", .step_n, msg))
}
.ok   <- function(msg) cat("     OK   ", msg, "\n", sep = "")
.info <- function(msg) cat("     ..   ", msg, "\n", sep = "")
.warn <- function(msg) cat("     WARN ", msg, "\n", sep = "")

.die <- function(...) stop("\n\n*** PATCH ABORTED ***\n", ..., "\n", call. = FALSE)

BACKUP_DIR <- NULL

.backup <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  dest <- file.path(BACKUP_DIR, path)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(path, dest, overwrite = TRUE))
    .die("Could not back up ", path, ". Check directory permissions.")
  invisible(TRUE)
}

## Replace a file wholesale with LF endings. Idempotent by construction.
.put <- function(path, content) {
  .backup(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- strsplit(content, "\n", fixed = TRUE)[[1]]
  while (length(lines) && !nzchar(lines[length(lines)]))
    lines <- lines[-length(lines)]
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeLines(lines, con = con, sep = "\n")
  .ok(paste("wrote", path))
  invisible(TRUE)
}

## Guarded targeted substitution. Reports, and never silently no-ops.
.sub_in <- function(path, from, to, fixed = TRUE, required = TRUE, label = NULL) {
  label <- if (is.null(label)) path else label
  if (!file.exists(path)) {
    if (required) .warn(paste(path, "not found; skipped")) else .info(paste(path, "absent; skipped"))
    return(invisible(FALSE))
  }
  txt <- readLines(path, warn = FALSE)
  hits <- sum(grepl(from, txt, fixed = fixed))
  if (hits == 0L) {
    already <- any(grepl(to, txt, fixed = fixed))
    if (already) .info(paste(label, "-- already applied"))
    else if (required) .warn(paste(label, "-- pattern not found; verify manually"))
    return(invisible(already))
  }
  .backup(path)
  out <- gsub(from, to, txt, fixed = fixed)
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(out, con = con, sep = "\n")
  .ok(sprintf("%s -- %d occurrence(s) replaced", label, hits))
  invisible(TRUE)
}

.drop_lines <- function(path, pattern, label) {
  if (!file.exists(path)) return(invisible(FALSE))
  txt <- readLines(path, warn = FALSE)
  keep <- !grepl(pattern, txt)
  if (all(keep)) { .info(paste(label, "-- already applied")); return(invisible(FALSE)) }
  .backup(path)
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(txt[keep], con = con, sep = "\n")
  .ok(sprintf("%s -- %d line(s) removed", label, sum(!keep)))
  invisible(TRUE)
}

## =============================================================================
##  STEP 0  --  PRE-FLIGHT
## =============================================================================

cat(strrep("=", 78), "\n")
cat("  BetaDanish  --  Phase 2, Patch 1 of 3 : correctness and documentation\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight checks")

if (!file.exists("DESCRIPTION"))
  .die("No DESCRIPTION in ", getwd(), ".\n",
       "setwd() to the BetaDanish package root and re-run.")

.desc <- read.dcf("DESCRIPTION")
if (!"Package" %in% colnames(.desc) || .desc[1, "Package"] != "BetaDanish")
  .die("DESCRIPTION says Package: ",
       if ("Package" %in% colnames(.desc)) .desc[1, "Package"] else "<missing>",
       ". This patch is for BetaDanish only.")

.ok(paste("package root:", getwd()))
.ok(paste("current version:", .desc[1, "Version"]))

for (pkg in c("devtools", "roxygen2")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    .die("Package '", pkg, "' is required. install.packages('", pkg, "')")
}
.ok("devtools and roxygen2 available")

if (!dir.exists("R")) .die("No R/ directory found.")

if (!dir.exists(".git"))
  .warn("No .git directory -- you appear to be in an extracted tarball. Patch 3 needs the git working copy.")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch1"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup directory:", BACKUP_DIR))

## =============================================================================
##  STEP 1  --  R/dist_functions.R
##  Recommendations 2, 3, 4, 5, 6, 16
## =============================================================================

.step("Rewriting R/dist_functions.R (recs 2-6, 16: numerics + reference)")

.put("R/dist_functions.R", r"---(#' The Beta-Danish Distribution
#'
#' Density, distribution function, quantile function, hazard function, survival
#' function and random generation for the four-parameter Beta-Danish
#' distribution.
#'
#' @param x,q Vector of quantiles (time points).
#' @param p Vector of probabilities.
#' @param n Number of observations to generate.
#' @param a Shape parameter (beta generator). Set `a = 1` for the
#'   three-parameter Exponentiated Danish (ED) submodel.
#' @param b Shape parameter (beta generator). Governs the upper tail: the
#'   survival function is regularly varying with index `-b`.
#' @param c Shape parameter (baseline shape).
#' @param k Scale parameter (baseline scale).
#' @param log,log.p Logical; if `TRUE`, densities/probabilities are returned on
#'   the log scale.
#' @param lower.tail Logical; if `TRUE` (default) probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @details
#' With baseline \eqn{G(t) = \{kt/(1+kt)\}^{c}}, the Beta-Danish CDF is the
#' regularised incomplete beta function \eqn{F(t) = I_{G(t)}(a, b)} and the
#' density is
#' \deqn{f(t) = \frac{c\,k^{ca}\,t^{ca-1}}{B(a,b)\,(1+kt)^{ca+1}}
#'              \bigl\{1 - G(t)\bigr\}^{b-1}.}
#' The family accommodates decreasing, increasing, unimodal and bathtub-shaped
#' hazard rates.
#'
#' @section Numerical notes:
#' All three of the density, distribution and quantile functions are evaluated
#' so as to retain full relative precision in the upper tail, which is where
#' this family is distinctive and where naive implementations fail:
#'
#' * \eqn{\log\{1 - G(t)\}} is formed as
#'   `log(-expm1(c * -log1p(1/(k*t))))`. Writing \eqn{\log G} as
#'   \eqn{-c\,\log(1 + 1/kt)} avoids the cancellation in
#'   \eqn{\log(kt) - \log(1+kt)}, which loses all significant digits once
#'   \eqn{kt \gtrsim 10^{15}}.
#' * The survival function uses the beta mirror identity
#'   \eqn{1 - I_{y}(a,b) = I_{1-y}(b,a)}, so no probability near one is ever
#'   subtracted from one.
#' * `qbetadanish` obtains \eqn{1 - u} directly via
#'   \eqn{1 - q\beta(p; a, b) = q\beta(p; b, a)} with the tail flag reversed,
#'   so `1 - p` is never formed.
#'
#' Parameters recycle element-wise against `x`, `q` or `p` following the usual
#' convention for R distribution functions, so a per-observation scale (as
#' produced by an AFT link) may be supplied directly.
#'
#' As \eqn{t \to \infty}, \eqn{S(t) \propto t^{-b}}; consequently
#' \eqn{E(X^{r})} is finite if and only if \eqn{b > r}.
#'
#' @return
#' `dbetadanish` gives the density, `pbetadanish` the distribution function,
#' `qbetadanish` the quantile function, `sbetadanish` the survival function,
#' `hbetadanish` the hazard function, and `rbetadanish` generates random
#' deviates. Length is the maximum of the lengths of the numeric arguments.
#'
#' @references
#' Ahmad, B., & Danish, M. Y. (2025). Development and characterization of a
#' flexible three-parameter lifetime distribution: theoretical properties and
#' real-world applications. *Journal of Applied Mathematics, Statistics and
#' Informatics*, 21(1). \doi{10.2478/jamsi-2025-0010}
#'
#' @examples
#' dbetadanish(x = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#' pbetadanish(q = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#' qbetadanish(p = 0.5, a = 1.5, b = 2, c = 3, k = 0.5)
#' hbetadanish(x = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#' rbetadanish(n = 10, a = 1.5, b = 2, c = 3, k = 0.5)
#'
#' # Per-observation scale, as used by the AFT link
#' dbetadanish(c(1, 2, 3), a = 1, b = 2, c = 1.5, k = c(0.5, 1, 2))
#'
#' # The survival tail is regularly varying with index -b
#' round(log(sbetadanish(c(1e10, 1e12), 1.5, 3, 2, 1)), 4)
#'
#' @name BetaDanish
NULL

## Recycle every numeric argument to the common length. Returns NULL if the
## common length is zero.
.bd_conform <- function(x, a, b, c, k) {
  n <- max(length(x), length(a), length(b), length(c), length(k))
  if (n == 0L) return(NULL)
  list(n = n,
       x = rep_len(as.numeric(x), n),
       a = rep_len(as.numeric(a), n),
       b = rep_len(as.numeric(b), n),
       c = rep_len(as.numeric(c), n),
       k = rep_len(as.numeric(k), n))
}

## Element-wise parameter validity.
.bd_bad <- function(p) {
  with(p, is.na(a) | is.na(b) | is.na(c) | is.na(k) |
          a <= 0 | b <= 0 | c <= 0 | k <= 0)
}

## log G(t) = -c * log1p(1/(k t)).  Full relative precision for every k*t > 0.
.bd_logG <- function(kt, c) -c * log1p(1 / kt)

## log{1 - G(t)}, stable for G arbitrarily close to either 0 or 1.
.bd_log1mG <- function(kt, c) log(-expm1(.bd_logG(kt, c)))

#' @rdname BetaDanish
#' @export
dbetadanish <- function(x, a, b, c, k, log = FALSE) {
  p <- .bd_conform(x, a, b, c, k)
  if (is.null(p)) return(numeric(0))
  out <- rep(-Inf, p$n)

  bad <- .bd_bad(p)
  if (any(bad)) {
    out[bad] <- NaN
    warning("Parameters a, b, c and k must be strictly positive; NaN returned ",
            "for ", sum(bad), " element(s).", call. = FALSE)
  }
  out[!bad & is.na(p$x)] <- NA_real_

  ok <- !bad & is.finite(p$x) & p$x > 0
  if (any(ok)) {
    xx <- p$x[ok]; aa <- p$a[ok]; bb <- p$b[ok]
    cc <- p$c[ok]; kk <- p$k[ok]
    kt <- kk * xx
    out[ok] <- log(cc) + cc * aa * log(kk) +
               (cc * aa - 1) * log(xx) -
               (cc * aa + 1) * log1p(kt) +
               (bb - 1) * .bd_log1mG(kt, cc) -
               lbeta(aa, bb)
  }
  if (log) out else exp(out)
}

#' @rdname BetaDanish
#' @export
pbetadanish <- function(q, a, b, c, k, lower.tail = TRUE, log.p = FALSE) {
  p <- .bd_conform(q, a, b, c, k)
  if (is.null(p)) return(numeric(0))

  out <- rep(if (lower.tail) 0 else 1, p$n)             # value at q <= 0
  out[is.infinite(p$x) & p$x > 0] <- if (lower.tail) 1 else 0

  bad <- .bd_bad(p)
  out[bad] <- NaN
  out[!bad & is.na(p$x)] <- NA_real_

  ok <- !bad & is.finite(p$x) & p$x > 0
  if (any(ok)) {
    xx <- p$x[ok]; aa <- p$a[ok]; bb <- p$b[ok]
    cc <- p$c[ok]; kk <- p$k[ok]
    kt <- kk * xx
    out[ok] <- if (lower.tail) {
      stats::pbeta(exp(.bd_logG(kt, cc)), aa, bb,
                   lower.tail = TRUE, log.p = log.p)
    } else {
      ## Beta mirror: S(t) = 1 - I_G(a,b) = I_{1-G}(b,a). No cancellation.
      stats::pbeta(-expm1(.bd_logG(kt, cc)), bb, aa,
                   lower.tail = TRUE, log.p = log.p)
    }
  }
  if (log.p) {
    seeded <- !ok & !bad & !is.na(p$x)
    out[seeded] <- log(out[seeded])
  }
  out
}

#' @rdname BetaDanish
#' @export
qbetadanish <- function(p, a, b, c, k, lower.tail = TRUE, log.p = FALSE) {
  cf <- .bd_conform(p, a, b, c, k)
  if (is.null(cf)) return(numeric(0))
  pp <- cf$x
  out <- rep(NA_real_, cf$n)

  bad <- .bd_bad(cf)
  out[bad] <- NaN

  ok <- !bad & !is.na(pp)
  if (any(ok)) {
    aa <- cf$a[ok]; bb <- cf$b[ok]; cc <- cf$c[ok]; kk <- cf$k[ok]
    ## w = 1 - u obtained directly; qbeta signals NaN for p outside its range.
    w  <- stats::qbeta(pp[ok], bb, aa, lower.tail = !lower.tail, log.p = log.p)
    lv <- log1p(-w) / cc                    # log v, v = u^(1/c)
    out[ok] <- exp(lv - log(-expm1(lv))) / kk
  }
  out
}

#' @rdname BetaDanish
#' @export
rbetadanish <- function(n, a, b, c, k) {
  if (length(n) > 1L) n <- length(n)
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0)
    stop("'n' must be a single non-negative number.", call. = FALSE)
  qbetadanish(stats::runif(n), a, b, c, k)
}

#' @rdname BetaDanish
#' @export
sbetadanish <- function(x, a, b, c, k, log = FALSE) {
  pbetadanish(x, a, b, c, k, lower.tail = FALSE, log.p = log)
}

#' @rdname BetaDanish
#' @export
hbetadanish <- function(x, a, b, c, k, log = FALSE) {
  log_f <- dbetadanish(x, a, b, c, k, log = TRUE)
  log_s <- pbetadanish(x, a, b, c, k, lower.tail = FALSE, log.p = TRUE)
  log_h <- log_f - log_s
  ## An exhausted survival with positive density is a divergent hazard, not a
  ## large finite one. 0/0 stays NaN rather than being silently invented.
  log_h[is.finite(log_f) & is.infinite(log_s) & log_s < 0] <- Inf
  if (log) log_h else exp(log_h)
}
)---")

## =============================================================================
##  STEP 2  --  R/utils-validation.R
##  Recommendation 6
## =============================================================================

.step("Rewriting R/utils-validation.R (rec 6: validation helper)")

.put("R/utils-validation.R", r"---(#' Extract Survival Data from Formula
#'
#' @param formula A survival formula (e.g. `Surv(time, status) ~ 1`).
#' @param data A data frame.
#'
#' @return A list with `time`, `status`, the design matrix `X`, and the model
#'   frame.
#' @noRd
extract_surv_data <- function(formula, data) {
  if (missing(data)) data <- environment(formula)

  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  Y  <- stats::model.extract(mf, "response")

  if (!survival::is.Surv(Y))
    stop("The response must be a 'Surv' object. Example: Surv(time, status) ~ 1",
         call. = FALSE)

  time   <- Y[, 1]
  status <- Y[, 2]

  if (any(!is.finite(time)) || any(time <= 0))
    stop("All survival times must be finite and strictly positive.", call. = FALSE)
  if (!all(status %in% c(0, 1)))
    stop("The event indicator must be coded 0 (censored) or 1 (event).",
         call. = FALSE)
  if (sum(status == 1) < 2L)
    stop("At least two events are needed to fit the model.", call. = FALSE)

  X <- stats::model.matrix(stats::terms(formula), mf)

  list(time = time, status = status, X = X, data_frame = mf)
}

#' Validate Beta-Danish Parameters
#'
#' Scalar validity check retained for the structural-property functions, which
#' take scalar parameters. The vectorised distribution functions validate
#' element-wise internally; see `.bd_bad()`.
#'
#' @param a,b,c,k Numeric parameters.
#' @return `TRUE` if every value is a finite strictly positive number.
#' @noRd
check_positive_params <- function(a, b, c, k) {
  v <- c(a, b, c, k)
  if (!is.numeric(v)) return(FALSE)
  if (anyNA(v)) return(FALSE)
  if (any(!is.finite(v))) return(FALSE)
  if (any(v <= 0)) return(FALSE)
  TRUE
}
)---")

## =============================================================================
##  STEP 3  --  R/fit_models.R
##  Recommendations 7, 10, 12
## =============================================================================

.step("Rewriting R/fit_models.R (recs 7, 10, 12: IC fields, diagnostics, docs)")

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

  start_list <- vector("list", n_starts)
  for (i in seq_len(n_starts)) {
    core <- c(log_b = log(stats::runif(1, 0.5, 5)),
              log_c = log(stats::runif(1, 0.5, 5)),
              log_k = log(k_base * stats::runif(1, 0.5, 2)))
    start_list[[i]] <- if (submodel) core else
      c(log_a = log(stats::runif(1, 0.5, 5)), core)
  }

  fit <- optim_multistart(ll_fun, start_list, method = method)
  if (is.null(fit))
    stop("Optimisation failed to converge from any starting point. Try a ",
         "larger n_starts, or method = \"NM\".", call. = FALSE)

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
    ac_correlation   = NA_real_
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

  invisible(NULL)
}
)---")

## =============================================================================
##  STEP 4  --  R/plotting_extras.R
##  Recommendation 1  (the critical fix)
## =============================================================================

.step("Rewriting R/plotting_extras.R (rec 1: restore plot.bd_aft / plot.bd_cure)")

.put("R/plotting_extras.R", r"---(#' Cox-Snell Residual Plot for AFT and Cure Fits
#'
#' Diagnostic Cox-Snell residual plot for a fitted AFT or cure model. Under a
#' correctly specified model the residuals behave like a unit exponential
#' sample, so the Kaplan-Meier estimate of their cumulative hazard should
#' follow the 45-degree line.
#'
#' @param x A fitted `"bd_aft"` or `"bd_cure"` object.
#' @param ... Further graphical parameters passed to `plot`.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 300
#' xcov  <- stats::rnorm(n)
#' k_i   <- exp(-0.5 - 0.3 * xcov)
#' t_sim <- rbetadanish(n, a = 1, b = 2, c = 1.5, k = k_i)
#' dat   <- data.frame(time = t_sim,
#'                     status = stats::rbinom(n, 1, 0.85),
#'                     x = xcov)
#' fit <- fit_bd_aft(survival::Surv(time, status) ~ x, data = dat, n_starts = 5)
#' plot(fit)
#' }
#'
#' @export
plot.bd_aft <- function(x, ...) {
  .bd_aft_or_cure_coxsnell(x, ...)
}

#' @rdname plot.bd_aft
#' @export
plot.bd_cure <- function(x, ...) {
  .bd_aft_or_cure_coxsnell(x, ...)
}

#' Natural-Scale Shape Parameters from an AFT or Cure Fit
#'
#' `fit_bd_aft()` and `fit_bd_cure()` optimise over `log_b` and `log_c`, so the
#' shape parameters must be looked up under those names and exponentiated.
#' Reading `coefficients["b"]` returns NA, which is what previously disabled
#' the Cox-Snell plots.
#'
#' @param fit A `"bd_aft"` or `"bd_cure"` object.
#' @return A list with numeric `b`, `c` and `delta`.
#' @noRd
.bd_shape_natural <- function(fit) {
  cf <- fit$coefficients
  pick <- function(nm) {
    if (paste0("log_", nm) %in% names(cf)) return(exp(as.numeric(cf[[paste0("log_", nm)]])))
    if (nm %in% names(cf))                 return(as.numeric(cf[[nm]]))
    NA_real_
  }
  delta_names <- grep("^delta_", names(cf), value = TRUE)
  list(b     = pick("b"),
       c     = pick("c"),
       delta = as.numeric(cf[delta_names]))
}

.bd_aft_or_cure_coxsnell <- function(x, ...) {
  d      <- x$data
  time   <- d$time
  status <- d$status
  X      <- d$X

  pars  <- .bd_shape_natural(x)
  b     <- pars$b
  cpar  <- pars$c
  delta <- pars$delta

  if (!is.finite(b) || !is.finite(cpar) || length(delta) == 0L ||
      any(!is.finite(delta))) {
    warning("Cox-Snell plot: the fit contains non-finite coefficients, so ",
            "residuals cannot be computed.", call. = FALSE)
    return(invisible(x))
  }
  if (length(delta) != ncol(X))
    stop("The coefficient layout does not match the design matrix.",
         call. = FALSE)

  k_i <- exp(as.numeric(X %*% delta))

  ## Cox-Snell residual is the fitted cumulative hazard, -log S(t).
  r <- -pbetadanish(time, a = 1, b = b, c = cpar, k = k_i,
                    lower.tail = FALSE, log.p = TRUE)

  keep <- is.finite(r) & r > 0
  if (!any(keep)) {
    warning("Cox-Snell plot: every residual is non-finite; nothing to plot.",
            call. = FALSE)
    return(invisible(x))
  }
  if (sum(keep) < length(r))
    warning(sprintf("Cox-Snell plot: dropped %d non-finite residual(s).",
                    sum(!keep)), call. = FALSE)

  r      <- r[keep]
  status <- status[keep]

  km_r <- survival::survfit(survival::Surv(r, status) ~ 1)
  H_r  <- -log(pmax(km_r$surv, 1e-12))

  graphics::plot(km_r$time, H_r, type = "s",
                 xlab = "Cox-Snell residual",
                 ylab = "Estimated cumulative hazard",
                 main = "Cox-Snell residuals", lwd = 2, ...)
  graphics::abline(0, 1, col = "red", lwd = 2, lty = 2)
  graphics::legend("topleft",
                   legend = c("KM cumulative hazard of residuals", "y = x"),
                   col = c("black", "red"), lwd = 2, lty = c(1, 2), bty = "n")
  invisible(x)
}
)---")

## =============================================================================
##  STEP 5  --  R/advanced_methods.R
##  Recommendation 8
## =============================================================================

.step("Rewriting R/advanced_methods.R (rec 8: honest parameter scales)")

.put("R/advanced_methods.R", r"---(## S3 methods for the AFT, cure and competing-risks fits.
##
## Shape parameters are estimated on the log scale. They are reported here on
## the natural scale with delta-method standard errors and an exponentiated
## log-scale confidence interval, which respects positivity. No Wald test is
## reported for them: the hypothesis "b = 0" is outside the parameter space, so
## a z statistic against zero would be meaningless. Regression coefficients are
## already on their natural scale and are reported with the usual Wald test.

#' Coefficient Table for a Shape/Regression Split
#' @noRd
.bd_split_coef <- function(est, se) {
  shape_idx <- grep("^log_", names(est))
  reg_idx   <- setdiff(seq_along(est), shape_idx)

  shape <- NULL
  if (length(shape_idx)) {
    l <- est[shape_idx]; s <- se[shape_idx]
    shape <- cbind(Estimate     = exp(l),
                   `Std. Error` = exp(l) * s,
                   `Lower 95%`  = exp(l - stats::qnorm(0.975) * s),
                   `Upper 95%`  = exp(l + stats::qnorm(0.975) * s))
    rownames(shape) <- sub("^log_", "", names(l))
  }

  reg <- NULL
  if (length(reg_idx)) {
    e <- est[reg_idx]; s <- se[reg_idx]
    z <- e / s
    reg <- cbind(Estimate     = e,
                 `Std. Error` = s,
                 `z value`    = z,
                 `Pr(>|z|)`   = 2 * stats::pnorm(abs(z), lower.tail = FALSE))
    rownames(reg) <- names(e)
  }
  list(shape = shape, regression = reg)
}

.bd_print_split <- function(x) {
  if (!is.null(x$shape)) {
    cat("Shape parameters (natural scale, delta-method SE):\n")
    print(round(x$shape, 4))
    cat("\n")
  }
  if (!is.null(x$regression)) {
    cat("Regression coefficients (log-scale link):\n")
    stats::printCoefmat(x$regression, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  }
  invisible(NULL)
}

.bd_se <- function(object) {
  se <- sqrt(pmax(diag(object$vcov), 0))
  names(se) <- names(object$coefficients)
  se
}

#' @export
print.bd_aft <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish AFT Model (Exponentiated Danish kernel, a = 1)\n")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n\n")
  cat("Coefficients (optimisation scale):\n")
  print(round(x$coefficients, 4)); cat("\n")
  invisible(x)
}

#' @export
summary.bd_aft <- function(object, ...) {
  res <- list(call = object$call,
              tables = .bd_split_coef(object$coefficients, .bd_se(object)),
              logLik = object$logLik)
  class(res) <- "summary.bd_aft"
  res
}

#' @export
print.summary.bd_aft <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish AFT Model (Exponentiated Danish kernel, a = 1)\n\n")
  .bd_print_split(x$tables)
  cat("---\nLog-Likelihood:", round(x$logLik, 4), "\n")
  invisible(x)
}

#' @export
print.bd_cure <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Cure Model (", x$type, ")\n", sep = "")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n\n")
  cat("Coefficients (optimisation scale):\n")
  print(round(x$coefficients, 4)); cat("\n")
  invisible(x)
}

#' @export
summary.bd_cure <- function(object, ...) {
  res <- list(call = object$call, type = object$type,
              tables = .bd_split_coef(object$coefficients, .bd_se(object)),
              logLik = object$logLik)
  class(res) <- "summary.bd_cure"
  res
}

#' @export
print.summary.bd_cure <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Cure Model (", x$type, ")\n\n", sep = "")
  .bd_print_split(x$tables)
  cat("---\nLog-Likelihood:", round(x$logLik, 4), "\n")
  invisible(x)
}

#' @export
print.bd_competing <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Competing Risks Model\n")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n\n")
  cat("Cause-specific estimates:\n")
  print(round(x$coefficients, 4)); cat("\n")
  invisible(x)
}

#' @export
summary.bd_competing <- function(object, ...) {
  res <- list(call = object$call, coefficients = object$coefficients,
              se = object$se, logLik = object$logLik,
              causes = object$causes)
  class(res) <- "summary.bd_competing"
  res
}

#' @export
print.summary.bd_competing <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Competing Risks Model\n\n")
  cat("Cause-specific estimates (natural scale):\n"); print(round(x$coefficients, 4))
  cat("\nStandard errors:\n"); print(round(x$se, 4))
  cat("---\nLog-Likelihood:", round(x$logLik, 4), "\n")
  cat("\nNote: the cause-specific marginals rest on an assumption of\n")
  cat("independent latent failure times, which is an identifying\n")
  cat("assumption and not testable from the observed data. See\n")
  cat("?fit_bd_competing.\n")
  invisible(x)
}
)---")

## =============================================================================
##  STEP 6  --  R/report_betadanish.R
##  Recommendation 7
## =============================================================================

.step("Rewriting R/report_betadanish.R (rec 7: AIC and BIC actually reported)")

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
    if (length(flags))
      cat("\nDiagnostic flags: ", paste(flags, collapse = "; "),
          "\n  See the Identifiability section of ?fit_betadanish.\n", sep = "")
  }
  invisible(x)
}
)---")

## =============================================================================
##  STEP 7  --  R/data_helpers.R
##  Recommendation 9
## =============================================================================

.step("Rewriting R/data_helpers.R (rec 9: name-based column selection)")

.put("R/data_helpers.R", r"---(#' Read and Prepare Survival Data
#'
#' Reads survival data from a CSV or Excel file and returns a clean data frame
#' ready for [fit_betadanish()]. Columns are selected by name, and covariates
#' keep their original type, so factors are not silently coerced to numbers.
#'
#' @param file Path to a `.csv`, `.xls` or `.xlsx` file.
#' @param time_col Name of the time column.
#' @param status_col Name of the event indicator column (1 = event,
#'   0 = censored). If `NULL` (default), all observations are treated as
#'   uncensored.
#' @param covar_cols Character vector of covariate columns to retain, or `NULL`.
#'
#' @return A data frame with columns `time`, `status` and any retained
#'   covariates. A `"bd_data_report"` attribute records rows read, rows
#'   dropped and the censoring proportion.
#'
#' @details
#' A covariate that happens to be named `time` or `status` is renamed with a
#' `_cov` suffix and a warning, rather than colliding with the response.
#'
#' Excel input requires the `readxl` package.
#'
#' @export
#'
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' write.csv(data.frame(survival_time = c(5, 8, 12, 16),
#'                      status = c(1, 1, 0, 1)),
#'           tmp, row.names = FALSE)
#' dat <- read_survival_data(tmp, time_col = "survival_time",
#'                           status_col = "status")
#' attr(dat, "bd_data_report")
#' unlink(tmp)
read_survival_data <- function(file, time_col, status_col = NULL,
                               covar_cols = NULL) {

  if (!is.character(file) || length(file) != 1L)
    stop("'file' must be a single file path.", call. = FALSE)
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)

  ext <- tolower(tools::file_ext(file))
  dat <- if (ext == "csv") {
    utils::read.csv(file, stringsAsFactors = FALSE, check.names = TRUE)
  } else if (ext %in% c("xls", "xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Reading Excel files requires the 'readxl' package.", call. = FALSE)
    as.data.frame(readxl::read_excel(file))
  } else {
    stop("Unsupported file extension '", ext,
         "'. Supply a .csv, .xls or .xlsx file.", call. = FALSE)
  }

  if (!nrow(dat)) stop("The file contains no data rows.", call. = FALSE)

  if (!time_col %in% names(dat))
    stop("Time column '", time_col, "' not found. Available columns: ",
         paste(names(dat), collapse = ", "), call. = FALSE)

  ## ---- build the response by name, never by position -----------------------
  clean <- data.frame(time = as.numeric(dat[[time_col]]))

  if (is.null(status_col)) {
    message("No status column supplied; treating all observations as ",
            "uncensored (status = 1).")
    clean$status <- 1
  } else {
    if (!status_col %in% names(dat))
      stop("Status column '", status_col, "' not found.", call. = FALSE)
    clean$status <- as.numeric(dat[[status_col]])
  }

  if (!is.null(covar_cols)) {
    missing_cov <- setdiff(covar_cols, names(dat))
    if (length(missing_cov))
      stop("Covariate column(s) not found: ",
           paste(missing_cov, collapse = ", "), call. = FALSE)

    extra <- dat[, covar_cols, drop = FALSE]
    clash <- intersect(names(extra), c("time", "status"))
    if (length(clash)) {
      warning("Covariate(s) named ", paste(clash, collapse = ", "),
              " collide with the response; renamed with a '_cov' suffix.",
              call. = FALSE)
      names(extra)[names(extra) %in% clash] <-
        paste0(names(extra)[names(extra) %in% clash], "_cov")
    }
    clean <- cbind(clean, extra)
  }

  n_read <- nrow(clean)
  clean  <- clean[stats::complete.cases(clean), , drop = FALSE]
  n_kept <- nrow(clean)

  if (n_kept < n_read)
    warning("Dropped ", n_read - n_kept, " row(s) with missing values.",
            call. = FALSE)
  if (!n_kept) stop("No complete rows remain after removing missing values.",
                    call. = FALSE)

  if (any(clean$time <= 0))
    warning("Some times are <= 0; survival models require strictly positive ",
            "times.", call. = FALSE)
  if (!all(clean$status %in% c(0, 1)))
    stop("The status column contains values other than 0 and 1. Recode so ",
         "that 1 = event and 0 = censored.", call. = FALSE)

  row.names(clean) <- NULL
  attr(clean, "bd_data_report") <- list(
    file            = normalizePath(file, mustWork = FALSE),
    rows_read       = n_read,
    rows_kept       = n_kept,
    rows_dropped    = n_read - n_kept,
    n_events        = sum(clean$status == 1),
    censoring_prop  = mean(clean$status == 0),
    covariates      = setdiff(names(clean), c("time", "status"))
  )
  clean
}
)---")

## =============================================================================
##  STEP 8  --  targeted documentation fixes
##  Recommendations 13, 14, 15, 17, 18
## =============================================================================

.step("Applying targeted documentation fixes (recs 13, 14, 15, 17, 18)")

## rec 14 -- melanoma column count
.sub_in("R/data.R",
        "#' @format A data frame with 205 rows and 6 columns:",
        "#' @format A data frame with 205 rows and 7 columns:",
        label = "rec 14: melanoma column count 6 -> 7")

## rec 15 -- carbon fibre unit
.sub_in("R/data.R",  "Gba", "GPa", label = "rec 15: data.R  Gba -> GPa")
.sub_in("README.md", "Gba", "GPa", required = FALSE,
        label = "rec 15: README  Gba -> GPa")
.sub_in("inst/WORDLIST", "Gba", "GPa", required = FALSE,
        label = "rec 15: WORDLIST  Gba -> GPa")

## rec 17 -- CED -> ED
.sub_in("R/aft_models.R",
        "Complementary Exponentiated Danish (CED) baseline (Beta-Danish with a=1)",
        "Exponentiated Danish (ED) kernel, that is the Beta-Danish distribution with a = 1",
        label = "rec 17: aft_models.R  CED -> ED")
.drop_lines("inst/WORDLIST", "^CED\\s*$", "rec 17: WORDLIST drop CED")

## rec 18 -- WORDLIST hygiene. NORI/Comorbid/comorbidities stay until Patch 2
##           removes the brain-cancer documentation that uses them.
.drop_lines("inst/WORDLIST", "^MRL\\s*$",
            "rec 18: WORDLIST drop MRL (no such function yet)")

## rec 13 -- Tsiatis caveat on the competing-risks documentation
if (file.exists("R/competing_risks.R")) {
  .cr <- readLines("R/competing_risks.R", warn = FALSE)
  if (!any(grepl("Tsiatis", .cr, fixed = TRUE))) {
    .anchor <- "#' then optimizes the joint likelihood."
    .idx <- grep(.anchor, .cr, fixed = TRUE)
    if (length(.idx) == 1L) {
      .backup("R/competing_risks.R")
      .note <- c(
        "#'",
        "#' @section Identifiability of the cause-specific marginals:",
        "#' Mutual independence of the latent failure times is an *identifying*",
        "#' assumption, not an empirically testable property. Tsiatis (1975)",
        "#' showed that without it infinitely many joint distributions for",
        "#' \\eqn{(T_1, \\ldots, T_m)} generate the same observable law of",
        "#' \\eqn{(T, \\delta)}, so the marginals cannot be recovered from the",
        "#' observed data alone without further structure such as a copula or a",
        "#' shared frailty.",
        "#'",
        "#' The direction of the resulting bias is known qualitatively. Under",
        "#' positive latent dependence the working independence model overstates",
        "#' each cause-specific survival at moderate times, so the fitted",
        "#' cumulative incidence functions are biased downward; negative",
        "#' dependence reverses both. The overall survival",
        "#' \\eqn{S(t) = \\prod_j S_j(t)} is more robust than the individual",
        "#' marginals, because cause-attribution error partly cancels across",
        "#' causes. Where substantive conclusions rest on the absolute",
        "#' cause-specific CIFs rather than on model selection, supplement this",
        "#' fit with a copula-based sensitivity analysis.",
        "#'",
        "#' @references",
        "#' Tsiatis, A. (1975). A nonidentifiability aspect of the problem of",
        "#' competing risks. *Proceedings of the National Academy of Sciences*,",
        "#' 72(1), 20-22. \\doi{10.1073/pnas.72.1.20}")
      .cr <- append(.cr, .note, after = .idx)
      .con <- file("R/competing_risks.R", open = "wb")
      writeLines(.cr, con = .con, sep = "\n"); close(.con)
      .ok("rec 13: Tsiatis (1975) identifiability caveat added")
    } else {
      .warn("rec 13: anchor comment not found in competing_risks.R; add the caveat by hand")
    }
  } else {
    .info("rec 13: already applied")
  }
}

## rec 18 -- normalise every text file in the package to LF
.step("Normalising line endings to LF (rec 18)")
.crlf_targets <- c(
  list.files("R", pattern = "[.]R$", full.names = TRUE),
  list.files("tests", pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE),
  Filter(file.exists, c("DESCRIPTION", "NEWS.md", "README.md",
                        "inst/WORDLIST", "inst/CITATION")))
.n_fixed <- 0L
for (f in .crlf_targets) {
  raw <- readBin(f, "raw", file.size(f))
  if (!any(raw == as.raw(13L))) next
  .backup(f)
  txt <- readLines(f, warn = FALSE)
  con <- file(f, open = "wb"); writeLines(txt, con = con, sep = "\n"); close(con)
  .n_fixed <- .n_fixed + 1L
}
.ok(sprintf("%d file(s) converted from CRLF to LF", .n_fixed))

## =============================================================================
##  STEP 9  --  DESCRIPTION
##  Recommendation 45 (development version only; 0.3.0 is set in Patch 3)
## =============================================================================

.step("Updating DESCRIPTION (rec 45: development version)")

.d <- readLines("DESCRIPTION", warn = FALSE)
.backup("DESCRIPTION")

.d <- sub("^Version:.*$", "Version: 0.2.0.9000", .d)

## Strip fields that R CMD build generates. They belong in the tarball, not the
## source DESCRIPTION, and a stale Packaged: date is a check NOTE waiting to
## happen.
.generated <- "^(Packaged|Author|Maintainer|NeedsCompilation|Built):"
if (any(grepl(.generated, .d))) {
  .d <- .d[!grepl(.generated, .d)]
  ## Drop continuation lines of a removed multi-line Author: field.
  .d <- .d[!grepl("^\\s+[A-Z][a-z]+ [A-Z][a-z]+ \\[", .d)]
  .info("removed build-generated fields (Packaged/Author/Maintainer/NeedsCompilation)")
}

.con <- file("DESCRIPTION", open = "wb")
writeLines(.d, con = .con, sep = "\n"); close(.con)
.ok("Version: 0.2.0.9000")

## =============================================================================
##  STEP 10  --  NEWS.md
##  Recommendation 11
## =============================================================================

.step("Rewriting NEWS.md (rec 11: remove unsupported 0.2.0 claims)")

.put("NEWS.md", r"---(# BetaDanish 0.2.0.9000 (development version)

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
)---")

## =============================================================================
##  STEP 11  --  TESTS
##  Recommendations 43, 44
## =============================================================================

.step("Writing regression and smoke tests (recs 43, 44)")

.put("tests/testthat/test-dist_functions.R", r"---(test_that("PDF integrates to one", {
  pars <- list(a = 1.5, b = 2.5, c = 2, k = 1)
  f <- function(t) dbetadanish(t, pars$a, pars$b, pars$c, pars$k)
  expect_equal(stats::integrate(f, 0, Inf, rel.tol = 1e-9)$value, 1,
               tolerance = 1e-6)
})

test_that("CDF matches numerical integration of the PDF", {
  a <- 1.2; b <- 3; c <- 2; k <- 0.8
  for (q in c(0.25, 1, 4)) {
    num <- stats::integrate(function(t) dbetadanish(t, a, b, c, k),
                            0, q, rel.tol = 1e-10)$value
    expect_equal(pbetadanish(q, a, b, c, k), num, tolerance = 1e-7)
  }
})

test_that("survival and CDF are complementary in the well-conditioned range", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t <- c(0.1, 0.5, 1, 2, 5)
  expect_equal(pbetadanish(t, a, b, c, k) +
                 sbetadanish(t, a, b, c, k),
               rep(1, length(t)), tolerance = 1e-12)
})

test_that("quantile function inverts the CDF, including far into the tail", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  p <- c(1e-10, 1e-4, 0.25, 0.5, 0.75, 1 - 1e-4, 1 - 1e-9, 1 - 1e-13)
  q <- qbetadanish(p, a, b, c, k)
  expect_true(all(is.finite(q)))
  expect_true(all(diff(q) > 0))
  expect_equal(pbetadanish(q, a, b, c, k), p, tolerance = 1e-8)
})

test_that("quantile boundaries are exact", {
  expect_equal(qbetadanish(0, 1.5, 3, 2, 1), 0)
  expect_identical(qbetadanish(1, 1.5, 3, 2, 1), Inf)
})

test_that("hazard equals density over survival", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t <- c(0.2, 1, 3)
  expect_equal(hbetadanish(t, a, b, c, k),
               dbetadanish(t, a, b, c, k) / sbetadanish(t, a, b, c, k),
               tolerance = 1e-12)
})

test_that("parameters recycle element-wise", {
  k <- c(0.5, 1, 2)
  got <- dbetadanish(c(1, 1, 1), a = 1, b = 2, c = 1.5, k = k)
  want <- vapply(k, function(ki) dbetadanish(1, 1, 2, 1.5, ki), numeric(1))
  expect_equal(got, want)
  expect_length(dbetadanish(1, a = 1, b = 2, c = 1.5, k = k), 3L)
})

test_that("invalid parameters give NaN and a warning", {
  expect_warning(res <- dbetadanish(1:3, a = -1, b = 2, c = 1, k = 1))
  expect_true(all(is.nan(res)))
  expect_true(all(is.nan(pbetadanish(1:3, a = 0, b = 2, c = 1, k = 1))))
  expect_true(all(is.nan(qbetadanish(0.5, a = 1, b = -2, c = 1, k = 1))))
})

test_that("boundary and missing values behave", {
  expect_equal(dbetadanish(c(-1, 0), 1.5, 3, 2, 1), c(0, 0))
  expect_equal(pbetadanish(c(-1, 0), 1.5, 3, 2, 1), c(0, 0))
  expect_equal(pbetadanish(Inf, 1.5, 3, 2, 1), 1)
  expect_equal(sbetadanish(Inf, 1.5, 3, 2, 1), 0)
  expect_true(is.na(dbetadanish(NA_real_, 1.5, 3, 2, 1)))
  expect_length(dbetadanish(numeric(0), 1.5, 3, 2, 1), 0L)
})

test_that("rbetadanish returns the requested length and respects the seed", {
  set.seed(11); x1 <- rbetadanish(50, 1.5, 3, 2, 1)
  set.seed(11); x2 <- rbetadanish(50, 1.5, 3, 2, 1)
  expect_length(x1, 50L)
  expect_identical(x1, x2)
  expect_true(all(x1 > 0))
})
)---")

.put("tests/testthat/test-numerics-tail.R", r"---(## Upper-tail regression tests.
##
## These pin the defects fixed in 0.2.0.9000. The Beta-Danish survival function
## is regularly varying with index -b, so on a log-log scale it must approach a
## straight line of slope -b. An implementation that floors log(1 - G) or that
## forms the survival by subtracting a near-one probability from one cannot
## satisfy this, which is exactly how the earlier defects escaped notice.

test_that("the survival tail is regularly varying with index -b", {
  a <- 1.5; c <- 2; k <- 1
  for (b in c(1.5, 3, 5)) {
    t1 <- 1e8; t2 <- 1e12
    ls1 <- sbetadanish(t1, a, b, c, k, log = TRUE)
    ls2 <- sbetadanish(t2, a, b, c, k, log = TRUE)
    expect_true(is.finite(ls1) && is.finite(ls2))
    slope <- (ls2 - ls1) / (log(t2) - log(t1))
    expect_equal(slope, -b, tolerance = 1e-4)
  }
})

test_that("log-survival stays finite far beyond double-precision saturation", {
  ls <- sbetadanish(c(1e14, 1e16, 1e18), a = 1.5, b = 3, c = 2, k = 1,
                    log = TRUE)
  expect_true(all(is.finite(ls)))
  expect_true(all(diff(ls) < 0))
})

test_that("the log-density tail has index -(b+1)", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t1 <- 1e8; t2 <- 1e12
  ld1 <- dbetadanish(t1, a, b, c, k, log = TRUE)
  ld2 <- dbetadanish(t2, a, b, c, k, log = TRUE)
  expect_true(is.finite(ld1) && is.finite(ld2))
  slope <- (ld2 - ld1) / (log(t2) - log(t1))
  expect_equal(slope, -(b + 1), tolerance = 1e-4)
})

test_that("density and survival match their analytic asymptotes", {
  ## As t -> Inf, with 1 - G(t) ~ c/(kt):
  ##   log f(t) -> b log(c/k) - (b+1) log t - log B(a,b)
  ##   log S(t) -> b log(c/k) -  b    log t - log b - log B(a,b)
  ## Both were verified against a 60-digit computation. At t = 1e18 with
  ## (a,b,c,k) = (1.5,3,2,1) the density asymptote is -161.8253135; the
  ## previous implementation, which floored log(1 - G) at about log(1e-16),
  ## returned -154.0013 instead, so this test separates the two decisively.
  a <- 1.5; b <- 3; c <- 2; k <- 1; t <- 1e18

  expect_equal(dbetadanish(t, a, b, c, k, log = TRUE),
               b * log(c / k) - (b + 1) * log(t) - lbeta(a, b),
               tolerance = 1e-6)

  expect_equal(sbetadanish(t, a, b, c, k, log = TRUE),
               b * log(c / k) - b * log(t) - log(b) - lbeta(a, b),
               tolerance = 1e-6)
})

test_that("the hazard is asymptotically 1/t and does not saturate", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t <- c(1e8, 1e10, 1e12)
  h <- hbetadanish(t, a, b, c, k)
  expect_true(all(is.finite(h)))
  expect_equal(unname(h * t), rep(b, length(t)), tolerance = 1e-3)
})

test_that("quantile round-trip survives the extreme upper tail", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  p <- 1 - 10^-(4:13)
  q <- qbetadanish(p, a, b, c, k)
  expect_true(all(is.finite(q)))
  expect_true(all(diff(q) > 0))
  expect_equal(pbetadanish(q, a, b, c, k), p, tolerance = 1e-9)
})

test_that("upper-tail probabilities agree with the lower-tail complement", {
  ## Only where the complement is itself well conditioned; the point of the
  ## mirror identity is that the survival stays right after this range.
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t <- c(0.5, 1, 5, 20)
  expect_equal(sbetadanish(t, a, b, c, k),
               1 - pbetadanish(t, a, b, c, k), tolerance = 1e-11)
})

test_that("log.p is consistent with the probability scale", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t <- c(0.5, 2, 10)
  expect_equal(pbetadanish(t, a, b, c, k, log.p = TRUE),
               log(pbetadanish(t, a, b, c, k)), tolerance = 1e-12)
  expect_equal(sbetadanish(t, a, b, c, k, log = TRUE),
               log(sbetadanish(t, a, b, c, k)), tolerance = 1e-12)
})
)---")

.put("tests/testthat/test-advanced-models.R", r"---(## Smoke tests for the AFT and cure paths, which had no coverage at all. The
## Cox-Snell plot methods are the specific reason: they were broken from the
## start and nothing exercised them.

make_aft_data <- function(n = 200, seed = 7) {
  set.seed(seed)
  xcov  <- stats::rnorm(n)
  k_i   <- exp(-0.5 - 0.3 * xcov)
  t_sim <- rbetadanish(n, a = 1, b = 2, c = 1.5, k = k_i)
  cens  <- stats::rexp(n, rate = 0.05)
  data.frame(time   = pmin(t_sim, cens),
             status = as.integer(t_sim <= cens),
             x      = xcov,
             group  = rep(c(0, 1), length.out = n))
}

test_that("fit_bd_aft returns a usable object", {
  skip_on_cran()
  dat <- make_aft_data()
  fit <- fit_bd_aft(survival::Surv(time, status) ~ x, data = dat, n_starts = 3)

  expect_s3_class(fit, "bd_aft")
  expect_true(all(c("log_b", "log_c") %in% names(fit$coefficients)))
  expect_true(is.finite(fit$logLik))
  expect_output(print(fit), "AFT")
})

test_that("summary.bd_aft reports shapes on the natural scale", {
  skip_on_cran()
  dat <- make_aft_data()
  fit <- fit_bd_aft(survival::Surv(time, status) ~ x, data = dat, n_starts = 3)
  s   <- summary(fit)

  expect_s3_class(s, "summary.bd_aft")
  expect_true(all(c("b", "c") %in% rownames(s$tables$shape)))
  expect_true(all(s$tables$shape[, "Estimate"] > 0))
  expect_true(all(s$tables$shape[, "Lower 95%"] > 0))
  expect_output(print(s), "natural scale")
})

test_that("plot.bd_aft draws instead of warning about NA coefficients", {
  skip_on_cran()
  dat <- make_aft_data()
  fit <- fit_bd_aft(survival::Surv(time, status) ~ x, data = dat, n_starts = 3)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  ## The defect was a warning about non-finite coefficients followed by an
  ## early return. Target that specifically rather than demanding total
  ## silence, which survfit could break for unrelated reasons.
  msg <- tryCatch({ plot(fit); character(0) },
                  warning = function(w) conditionMessage(w))
  expect_false(any(grepl("non-finite coefficients", msg)))
  expect_false(any(grepl("cannot be computed", msg)))
})

test_that("the natural-scale shape helper reads log_b and log_c", {
  fake <- list(coefficients = c(log_b = log(2), log_c = log(1.5),
                                `delta_(Intercept)` = -0.5))
  got <- BetaDanish:::.bd_shape_natural(fake)
  expect_equal(got$b, 2)
  expect_equal(got$c, 1.5)
  expect_equal(got$delta, -0.5)
})

test_that("fit_bd_cure runs for both formulations", {
  skip_on_cran()
  dat <- make_aft_data(n = 250, seed = 21)
  for (ty in c("mixture", "promotion")) {
    fit <- fit_bd_cure(survival::Surv(time, status) ~ 1,
                       formula_cure = ~ group, data = dat,
                       type = ty, n_starts = 3)
    expect_s3_class(fit, "bd_cure")
    expect_true(is.finite(fit$logLik))
    expect_identical(fit$type, ty)
  }
})

test_that("fit_betadanish stores information criteria and diagnostics", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 3)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3))

  expect_true(is.finite(fit$AIC))
  expect_true(is.finite(fit$BIC))
  expect_equal(fit$nobs, 150L)
  expect_equal(fit$npar, 3L)
  expect_true(is.list(fit$diagnostics))
  expect_equal(unname(stats::AIC(fit)), fit$AIC, tolerance = 1e-8)

  rep <- report_betadanish(fit)
  expect_output(print(rep), "AIC")
})

test_that("read_survival_data selects by name and reports on the data", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(
    data.frame(t = c(1, 2, 3, NA), ev = c(1, 0, 1, 1), grp = c("a","b","a","b")),
    tmp, row.names = FALSE)

  dat <- suppressWarnings(
    read_survival_data(tmp, time_col = "t", status_col = "ev",
                       covar_cols = "grp"))

  expect_named(dat, c("time", "status", "grp"))
  expect_equal(nrow(dat), 3L)
  expect_type(dat$grp, "character")
  expect_equal(attr(dat, "bd_data_report")$rows_dropped, 1L)
})

test_that("a covariate named time is renamed rather than colliding", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(dur = c(1, 2), ev = c(1, 1), time = c(9, 9)),
                   tmp, row.names = FALSE)

  expect_warning(
    dat <- read_survival_data(tmp, time_col = "dur", status_col = "ev",
                              covar_cols = "time"),
    "collide")
  expect_named(dat, c("time", "status", "time_cov"))
  expect_equal(dat$time, c(1, 2))
})
)---")

## =============================================================================
##  STEP 12  --  VERIFY
## =============================================================================

.step("Parsing every modified R file")

.parse_targets <- c(list.files("R", pattern = "[.]R$", full.names = TRUE),
                    list.files("tests", pattern = "[.]R$", recursive = TRUE,
                               full.names = TRUE))
.bad <- character(0)
for (f in .parse_targets) {
  e <- tryCatch({ parse(f); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(e)) .bad <- c(.bad, paste0("  ", f, ": ", e))
}
if (length(.bad))
  .die("These files do not parse:\n", paste(.bad, collapse = "\n"),
       "\n\nRestore from ", BACKUP_DIR, " and report the error.")
.ok(sprintf("%d file(s) parse cleanly", length(.parse_targets)))

.step("devtools::document()")
.res <- tryCatch({ devtools::document(roclets = c("rd", "collate", "namespace")); TRUE },
                 error = function(e) conditionMessage(e))
if (!isTRUE(.res))
  .die("document() failed:\n  ", .res,
       "\n\nBackups are in ", BACKUP_DIR)
.ok("documentation and NAMESPACE regenerated")

.step("Confirming sbetadanish is exported")
.ns <- readLines("NAMESPACE", warn = FALSE)
if (!any(grepl("export(sbetadanish)", .ns, fixed = TRUE)))
  .warn("sbetadanish is not in NAMESPACE; check the roxygen block in R/dist_functions.R")
else .ok("export(sbetadanish) present")

.step("devtools::check() -- this takes a few minutes")
.chk <- tryCatch(
  devtools::check(document = FALSE, args = c("--as-cran"), error_on = "never"),
  error = function(e) { .warn(paste("check() errored:", conditionMessage(e))); NULL })

cat("\n"); cat(strrep("=", 78), "\n")
if (!is.null(.chk)) {
  cat("  CHECK RESULT\n")
  cat(strrep("=", 78), "\n")
  cat("  errors:  ", length(.chk$errors), "\n", sep = "")
  cat("  warnings:", length(.chk$warnings), "\n", sep = "")
  cat("  notes:   ", length(.chk$notes), "\n", sep = "")
  for (nm in c("errors", "warnings", "notes")) {
    v <- .chk[[nm]]
    if (length(v)) { cat("\n-- ", toupper(nm), " --\n", sep = ""); cat(v, sep = "\n\n") }
  }
} else {
  cat("  check() did not complete; run devtools::check() manually.\n")
}

## =============================================================================
##  SUMMARY
## =============================================================================

cat("\n"); cat(strrep("=", 78), "\n")
cat("  PATCH 1 OF 3 COMPLETE  --  correctness and documentation\n")
cat(strrep("=", 78), "\n\n")
cat("  Recommendations implemented\n")
cat("    1   plot.bd_aft / plot.bd_cure restored (log_b, log_c + exp)\n")
cat("    2   dbetadanish far-tail log(1 - G) via -log1p(1/kt) and expm1\n")
cat("    3   qbetadanish rewritten on the beta mirror identity\n")
cat("    4   hbetadanish -700 floor removed\n")
cat("    5   pbetadanish survival via I_{1-G}(b, a)\n")
cat("    6   element-wise parameter recycling; validation rewritten\n")
cat("    7   npar/nobs/AIC/BIC on the fit; report_betadanish fixed\n")
cat("    8   AFT and cure summaries split shape vs regression, natural scale\n")
cat("    9   read_survival_data selects by name; data report attribute\n")
cat("    10  convergence, singular-information and identifiability warnings\n")
cat("    11  NEWS.md rewritten; false 0.2.0 claims corrected in the open\n")
cat("    12  b = 1 ridge and (a, c) confounding documented\n")
cat("    13  Tsiatis (1975) caveat on the competing-risks page\n")
cat("    14  melanoma column count 6 -> 7\n")
cat("    15  Gba -> GPa in data.R, README and WORDLIST\n")
cat("    16  reference year and title aligned with CITATION\n")
cat("    17  CED -> ED throughout\n")
cat("    18  WORDLIST tidied; all text files normalised to LF\n")
cat("    43  smoke tests for AFT, cure, report and CSV paths\n")
cat("    44  upper-tail regression tests anchored on S(t) ~ t^-b\n")
cat("    45  Version: 0.2.0.9000 (0.3.0 is set in Patch 3)\n\n")
cat("  Deliberately NOT done here\n")
cat("    - No dataset was touched. brain_cancer removal is Patch 2.\n")
cat("    - NORI / Comorbid / comorbidities stay in WORDLIST until Patch 2\n")
cat("      removes the documentation that uses them.\n")
cat("    - No git operation. No CRAN submission.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n", sep = "")
cat("  To undo:  file.copy(list.files('", BACKUP_DIR,
    "', recursive = TRUE, full.names = TRUE), '.', overwrite = TRUE)\n\n", sep = "")
cat("  Next: confirm 0 errors and 0 warnings above, skim the NEWS.md diff,\n")
cat("  then ask for Patch 2 (dataset removal, guinea_pig, CSV pipeline).\n\n")
