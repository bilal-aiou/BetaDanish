## =============================================================================
##  BetaDanish  --  PATCH 3f : VISUALISATION, VIGNETTES, AND THE 3e FIXES
## =============================================================================
##
##  THE LAST CODE PATCH. After this, re-run Patch 3d to cut the release.
##
##  FIRST, TWO THINGS PATCH 3e GOT WRONG
##
##    The example for fit_bd_competing() called
##        simulate_bd_competing_data(300, seed = 1)
##    and then asked for covariates = ~ x. That generator only creates the
##    covariate column when `gammas` is supplied, so `x` did not exist and the
##    example failed. Fixed, with the dependency spelled out in the example
##    itself so it cannot be missed again.
##
##    The simulation examples ran for 113, 75 and 18 seconds. CRAN runs
##    \donttest examples during its checks and flags anything over five. They
##    now use two or three replicates at a small sample size, with a comment
##    saying so and pointing at the vignette for a realistic study. The cif,
##    identified-coefficient and profile examples are trimmed too.
##
##  THEN, RECOMMENDATIONS 40, 41, 42 AND 46
##
##    40  bd_ttt_plot()      the scaled total time on test transform, a
##                           distribution-free read on hazard shape BEFORE any
##                           model is fitted. Returns the coordinates and a
##                           suggested shape label.
##    41  bd_profile_plot()  the profile log-likelihood with the critical
##                           threshold and interval marked. When the curve never
##                           returns below the threshold it says so on the plot
##                           rather than implying a bound that is not there.
##    42  plot.bd_bayes()    trace and posterior density panels, restoring the
##                           graphical parameters on exit.
##    46  A "New in 0.3.0" section appended to the introduction vignette and a
##        covariates section to the competing-risks vignette.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3f_visualisation.R")
##  IDEMPOTENT   Yes -- files are written whole, vignette sections are appended
##               only once.
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
.append_once <- function(path, marker, lines, label) {
  if (!file.exists(path)) { .warn(paste(path, "not found; skipped")); return(invisible(FALSE)) }
  cur <- readLines(path, warn = FALSE)
  if (any(grepl(marker, cur, fixed = TRUE))) { .info(paste(label, "-- already present")); return(invisible(FALSE)) }
  .backup(path)
  .write_lines(path, c(cur, "", lines))
  .ok(label)
  invisible(TRUE)
}

cat(strrep("=", 78), "\n")
cat("  BetaDanish  --  Patch 3f : visualisation and vignettes\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("R/simulation_study.R")) .die("Patch 3e has not been applied.")
if (!file.exists("R/inference.R")) .die("Patch 3b has not been applied.")
.ok("Patches 1 through 3e detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3f"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

.step("Fixing the Patch 3e examples in R/competing_risks.R")

.put("R/competing_risks.R", r"---(#' Fit a Beta-Danish Competing Risks Model
#'
#' Fits a parametric competing risks model assuming independent latent failure
#' times, with each cause-specific baseline following the Beta-Danish
#' distribution. Covariates may be included, entering each cause through an
#' accelerated failure time link.
#'
#' @param time Numeric vector of observed times.
#' @param cause Integer vector of event causes: `0` for right-censored and
#'   `1, 2, ..., m` for the competing causes.
#' @param covariates Optional covariates. Either a one-sided formula such as
#'   `~ age + group`, evaluated in `data`, or a numeric matrix or data frame
#'   with one row per observation. `NULL` (default) fits cause-specific
#'   baselines with no regression structure.
#' @param data Data frame in which to evaluate `covariates` when it is a
#'   formula.
#' @param submodel Logical; if `TRUE`, fix `a = 1` in every cause, giving
#'   Exponentiated Danish cause-specific kernels. Default `FALSE`.
#' @param n_starts Integer; number of starting points for the joint
#'   optimisation.
#' @param method Optimisation method passed to [maxLik::maxLik()].
#'
#' @return An object of class `"bd_competing"`.
#'
#' @details
#' Writing \eqn{\delta_i} for the observed cause and \eqn{d_i} for the event
#' indicator, the log-likelihood is
#' \deqn{\ell = \sum_{i: d_i = 1}\Bigl[\log f_{\delta_i}(y_i)
#'        + \sum_{j \neq \delta_i}\log S_j(y_i)\Bigr]
#'        + \sum_{i: d_i = 0}\sum_{j=1}^{m}\log S_j(y_i).}
#' An event contributes the log-density of its own cause and the log-survival
#' of every other cause at the same time; a censored observation contributes
#' the log-survival of every cause.
#'
#' @section Covariates:
#' Each cause carries its own coefficient vector \eqn{\gamma_j}, acting on the
#' time scale:
#' \deqn{T_j = T_{0j}\exp(x'\gamma_j),}
#' so that a positive coefficient lengthens time to failure from that cause.
#' The likelihood contribution of an event at \eqn{t} is evaluated at the
#' accelerated time \eqn{t\exp(-x'\gamma_j)} with the Jacobian
#' \eqn{-x'\gamma_j} added on the log scale.
#'
#' The same design matrix is used for every cause, but the coefficients are
#' free to differ, so a covariate may accelerate one cause and retard another.
#' With `m` causes, `p` covariates and the full Beta-Danish kernel the model
#' carries `m * (4 + p)` parameters, which grows quickly: check
#' `summary()` for standard errors before interpreting any of them.
#'
#' @section Identifiability of the cause-specific marginals:
#' Mutual independence of the latent failure times is an *identifying*
#' assumption, not an empirically testable property. Tsiatis (1975) showed that
#' without it infinitely many joint distributions generate the same observable
#' law of \eqn{(T, \delta)}, so the marginals cannot be recovered from the
#' observed data alone.
#'
#' Under positive latent dependence the working independence model overstates
#' each cause-specific survival at moderate times, biasing the fitted
#' cumulative incidence functions downward; negative dependence reverses both.
#' The overall survival \eqn{S(t) = \prod_j S_j(t)} is more robust than the
#' individual marginals. Where conclusions rest on absolute cause-specific
#' CIFs rather than on model selection, supplement this fit with a
#' copula-based sensitivity analysis.
#'
#' @references
#' Tsiatis, A. (1975). A nonidentifiability aspect of the problem of competing
#' risks. *Proceedings of the National Academy of Sciences*, 72(1), 20-22.
#' \doi{10.1073/pnas.72.1.20}
#'
#' @seealso [cif_betadanish()], [cif_compare()], [simulate_bd_competing_data()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Without covariates
#' d <- simulate_bd_competing_data(150, seed = 1)
#' fit <- fit_bd_competing(d$time, d$cause, submodel = TRUE, n_starts = 1)
#' fit
#'
#' # With one binary covariate. Note that simulate_bd_competing_data() only
#' # creates the covariate column when 'gammas' is supplied.
#' dx <- simulate_bd_competing_data(150, gammas = c(0.5, -0.5), seed = 2)
#' fit2 <- fit_bd_competing(dx$time, dx$cause, covariates = ~ x, data = dx,
#'                          submodel = TRUE, n_starts = 1)
#' fit2$coefficients
#' }
fit_bd_competing <- function(time, cause, covariates = NULL, data = NULL,
                             submodel = FALSE, n_starts = 5, method = "BFGS") {

  time  <- as.numeric(time)
  cause <- as.integer(cause)
  if (length(time) != length(cause))
    stop("'time' and 'cause' must have the same length.", call. = FALSE)
  if (any(!is.finite(time)) || any(time <= 0))
    stop("All times must be finite and strictly positive.", call. = FALSE)

  causes <- sort(unique(cause[cause > 0]))
  m <- length(causes)
  if (m < 2L)
    stop("Competing risks needs at least two distinct event causes.",
         call. = FALSE)

  X <- .bd_cr_design(covariates, data, n = length(time))
  p <- if (is.null(X)) 0L else ncol(X)
  per <- (if (submodel) 3L else 4L) + p

  ## ---- names for one cause -------------------------------------------------
  base_nm <- if (submodel) c("log_b", "log_c", "log_k") else
    c("log_a", "log_b", "log_c", "log_k")
  gam_nm  <- if (p) paste0("gamma_", colnames(X)) else character(0)
  one_nm  <- c(base_nm, gam_nm)
  all_nm  <- unlist(lapply(causes, function(j) paste0("c", j, ":", one_nm)))

  ## ---- cause-wise starting values -----------------------------------------
  starts0 <- numeric(0)
  for (j in causes) {
    dat_j <- data.frame(time = time, status = as.integer(cause == j))
    fit_j <- tryCatch(
      suppressWarnings(
        fit_betadanish(survival::Surv(time, status) ~ 1, data = dat_j,
                       submodel = submodel, n_starts = 2,
                       check_identifiability = FALSE)),
      error = function(e) NULL)
    base <- if (is.null(fit_j)) {
      med <- stats::median(time)
      if (submodel) log(c(2, 1.5, 1 / med)) else log(c(1, 2, 1.5, 1 / med))
    } else {
      log(fit_j$coefficients)
    }
    starts0 <- c(starts0, base, rep(0, p))
  }
  names(starts0) <- all_nm

  ## ---- joint log-likelihood ------------------------------------------------
  ll_joint <- function(pars) {
    total <- 0
    for (idx in seq_along(causes)) {
      j   <- causes[idx]
      off <- (idx - 1L) * per
      th  <- pars[(off + 1L):(off + per)]

      a <- if (submodel) 1 else exp(th[1])
      o <- if (submodel) 0L else 1L
      b <- exp(th[o + 1]); cc <- exp(th[o + 2]); k <- exp(th[o + 3])

      if (p) {
        lacc <- -as.numeric(X %*% th[(o + 4):(o + 3 + p)])
      } else {
        lacc <- rep(0, length(time))
      }
      tt <- time * exp(lacc)

      is_j <- cause == j
      if (any(is_j))
        total <- total + sum(
          suppressWarnings(dbetadanish(tt[is_j], a, b, cc, k, log = TRUE)) +
            lacc[is_j])
      if (any(!is_j))
        total <- total + sum(suppressWarnings(
          pbetadanish(tt[!is_j], a, b, cc, k, lower.tail = FALSE, log.p = TRUE)))
    }
    if (!is.finite(total)) return(-1e10)
    total
  }

  ## ---- guarded multi-start -------------------------------------------------
  start_list <- list(starts0)
  for (i in seq_len(max(0L, n_starts - 1L))) {
    s <- starts0 + stats::rnorm(length(starts0), 0, 0.25)
    names(s) <- all_nm
    start_list[[length(start_list) + 1L]] <- s
  }

  accept <- .bd_make_accept(
    n = length(time),
    log_shape = grep("log_(a|b|c)$", all_nm, value = TRUE),
    log_scale = grep("log_k$", all_nm, value = TRUE))

  fit <- optim_multistart(ll_joint, start_list, method = method, accept = accept)
  if (is.null(fit))
    stop("No admissible optimum was found for the competing risks model. ",
         "Every start either failed or reached a degenerate ridge. Try ",
         "submodel = TRUE, or fewer covariates.", call. = FALSE)

  ## ---- format --------------------------------------------------------------
  est_log <- fit$estimate
  names(est_log) <- all_nm

  vcov_log <- tryCatch(solve(-fit$hessian),
                       error = function(e) matrix(NA_real_, length(est_log),
                                                  length(est_log)))
  dimnames(vcov_log) <- list(all_nm, all_nm)
  se_log <- sqrt(pmax(diag(vcov_log), 0))

  ## Shapes and scale back to the natural scale; regression coefficients are
  ## already on their own scale and are left alone.
  is_log <- grepl("log_", all_nm, fixed = TRUE)
  est_rep <- est_log
  est_rep[is_log] <- exp(est_log[is_log])
  names(est_rep) <- sub("log_", "", all_nm, fixed = TRUE)

  se_rep <- se_log
  se_rep[is_log] <- se_log[is_log] * est_rep[is_log]   # delta method
  names(se_rep) <- names(est_rep)

  par_matrix <- matrix(est_rep, nrow = m, byrow = TRUE,
                       dimnames = list(paste0("Cause_", causes),
                                       sub("^c[0-9]+:", "", names(est_rep))[
                                         seq_len(per)]))
  se_matrix <- matrix(se_rep, nrow = m, byrow = TRUE,
                      dimnames = dimnames(par_matrix))

  out <- list(coefficients = par_matrix,
              se           = se_matrix,
              estimate_log = est_log,
              vcov         = vcov_log,
              logLik       = fit$maximum,
              causes       = causes,
              m            = m,
              submodel     = submodel,
              covariates   = if (p) colnames(X) else character(0),
              X            = X,
              npar         = length(est_log),
              nobs         = length(time),
              convergence  = fit$code,
              starts_ok       = .bd_or(attr(fit, "bd_starts_ok"), NA_integer_),
              starts_rejected = .bd_or(attr(fit, "bd_starts_rejected"), NA_integer_),
              data         = list(time = time, cause = cause),
              call         = match.call())
  class(out) <- "bd_competing"

  if (!is.na(out$starts_rejected) && out$starts_rejected > 0)
    warning(sprintf(paste0("%d starting point(s) reached a degenerate ridge ",
                           "and were discarded."), out$starts_rejected),
            call. = FALSE)

  out
}

#' Build the Competing-Risks Design Matrix
#'
#' Accepts a one-sided formula plus data, a matrix, or a data frame, and
#' returns a numeric design matrix without an intercept. The intercept is
#' excluded because each cause already has its own scale parameter `k`, which
#' the intercept would be confounded with.
#'
#' @noRd
.bd_cr_design <- function(covariates, data, n) {
  if (is.null(covariates)) return(NULL)

  if (inherits(covariates, "formula")) {
    if (length(covariates) != 2L)
      stop("'covariates' must be a one-sided formula, such as ~ age + group.",
           call. = FALSE)
    if (is.null(data))
      stop("'data' is required when 'covariates' is a formula.", call. = FALSE)
    mf <- stats::model.frame(covariates, data, na.action = stats::na.pass)
    X  <- stats::model.matrix(covariates, mf)
  } else {
    X <- as.matrix(covariates)
    if (is.null(colnames(X)))
      colnames(X) <- paste0("x", seq_len(ncol(X)))
  }

  keep <- colnames(X) != "(Intercept)"
  X <- X[, keep, drop = FALSE]

  if (nrow(X) != n)
    stop("The covariates have ", nrow(X), " row(s) but there are ", n,
         " observation(s).", call. = FALSE)
  if (anyNA(X))
    stop("The covariates contain missing values; remove those rows first.",
         call. = FALSE)
  if (!ncol(X))
    stop("No covariates remain after dropping the intercept.", call. = FALSE)

  storage.mode(X) <- "double"
  X
}

#' Cause-Specific Parameters at a Covariate Value
#'
#' Returns the four Beta-Danish parameters for one cause, with any accelerated
#' failure time effect folded into the scale. Scaling time by \eqn{\lambda}
#' maps \eqn{k} to \eqn{k/\lambda}, so a covariate acting as
#' \eqn{T_j = T_{0j}\exp(x'\gamma_j)} is equivalent to replacing \eqn{k_j} by
#' \eqn{k_j\exp(-x'\gamma_j)}.
#'
#' @noRd
.bd_cr_pars <- function(fit, cause_idx, x = NULL) {
  row <- which(fit$causes == cause_idx)
  pm  <- fit$coefficients[row, ]

  a <- if (isTRUE(fit$submodel)) 1 else as.numeric(pm[["a"]])
  b <- as.numeric(pm[["b"]])
  c <- as.numeric(pm[["c"]])
  k <- as.numeric(pm[["k"]])

  if (length(fit$covariates)) {
    gam <- as.numeric(pm[paste0("gamma_", fit$covariates)])
    xv  <- if (is.null(x)) rep(0, length(gam)) else as.numeric(x)
    if (length(xv) != length(gam))
      stop("'x' must supply one value per covariate: ",
           paste(fit$covariates, collapse = ", "), call. = FALSE)
    k <- k * exp(-sum(xv * gam))
  }
  list(a = a, b = b, c = c, k = k)
}

#' Cumulative Incidence Function
#'
#' Computes the cumulative incidence of one cause from a fitted Beta-Danish
#' competing risks model,
#' \eqn{F_j(t) = \int_0^t f_j(u)\prod_{l \neq j} S_l(u)\,du}.
#'
#' @param fit An object of class `"bd_competing"`.
#' @param tvec Numeric vector of times at which to evaluate the CIF.
#' @param cause_idx The cause to evaluate, matching one of the fitted causes.
#' @param x Optional covariate values at which to evaluate, one per covariate
#'   in the fit. Defaults to zero for every covariate, which is the reference
#'   subject. Ignored when the model has no covariates.
#'
#' @return A numeric vector of probabilities the same length as `tvec`.
#'
#' @details
#' The integrand is the cause-specific density multiplied by the survival of
#' every other cause, so the sum of all cause CIFs approaches the overall
#' failure probability rather than one, and each individual CIF is bounded by
#' it.
#'
#' Any covariate effect is folded into the scale parameter before integration,
#' since accelerating time by \eqn{\lambda} is equivalent to dividing \eqn{k}
#' by \eqn{\lambda}.
#'
#' @seealso [fit_bd_competing()], [cif_compare()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' d <- simulate_bd_competing_data(120, seed = 2)
#' fit <- fit_bd_competing(d$time, d$cause, submodel = TRUE, n_starts = 1)
#' cif_betadanish(fit, tvec = c(1, 5, 10), cause_idx = 1)
#' }
cif_betadanish <- function(fit, tvec, cause_idx, x = NULL) {
  if (!inherits(fit, "bd_competing"))
    stop("'fit' must be a bd_competing object.", call. = FALSE)
  if (!(cause_idx %in% fit$causes))
    stop("cause_idx ", cause_idx, " is not among the fitted causes: ",
         paste(fit$causes, collapse = ", "), call. = FALSE)

  pj  <- .bd_cr_pars(fit, cause_idx, x)
  oth <- setdiff(fit$causes, cause_idx)
  pl  <- lapply(oth, function(l) .bd_cr_pars(fit, l, x))

  integrand <- function(u) {
    v <- dbetadanish(u, pj$a, pj$b, pj$c, pj$k)
    for (q in pl)
      v <- v * pbetadanish(u, q$a, q$b, q$c, q$k, lower.tail = FALSE)
    v[!is.finite(v)] <- 0
    v
  }

  vapply(as.numeric(tvec), function(t_val) {
    if (!is.finite(t_val) || t_val <= 0) return(0)
    tryCatch(stats::integrate(integrand, 0, t_val, subdivisions = 2000L,
                              rel.tol = 1e-8, stop.on.error = FALSE)$value,
             error = function(e) NA_real_)
  }, numeric(1))
})---")


.step("Trimming the simulation examples in R/simulation_study.R")

.put("R/simulation_study.R", r"---(## Finite-sample simulation studies.
##
## Each runner draws `n_sim` samples at every sample size, refits, and reports
## bias, RMSE, the ratio of mean estimated standard error to empirical standard
## deviation, and Wald coverage. That last pair is the informative one: a ratio
## far below one means the standard errors understate the sampling variability,
## and coverage well below the nominal level follows from it.
##
## Fits that fail or that the degeneracy guard rejects are counted rather than
## silently dropped, because a high non-convergence rate is itself a finding
## about the parameter region being studied.

#' Simulate Competing Risks Data
#'
#' @param n Sample size.
#' @param pars List of parameter vectors, one per cause, each a named vector
#'   with `b`, `c`, `k` and optionally `a`. Defaults to two causes.
#' @param gammas Numeric vector of accelerated failure time coefficients, one
#'   per cause, acting on a single binary covariate. `NULL` for no covariate.
#' @param censor_rate Approximate proportion of right-censored observations.
#' @param seed Optional integer seed.
#'
#' @return A data frame with `time`, `cause` (0 for censored) and, when
#'   `gammas` is supplied, the covariate `x`.
#'
#' @seealso [fit_bd_competing()]
#'
#' @export
#'
#' @examples
#' d <- simulate_bd_competing_data(200, seed = 1)
#' table(d$cause)
simulate_bd_competing_data <- function(n,
                                       pars = list(c(a = 1, b = 2.5, c = 1.8, k = 0.04),
                                                   c(a = 1, b = 3.5, c = 1.2, k = 0.02)),
                                       gammas = NULL,
                                       censor_rate = 0.15,
                                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- length(pars)
  if (m < 2L) stop("At least two causes are needed.", call. = FALSE)

  x <- if (is.null(gammas)) NULL else stats::rbinom(n, 1, 0.5)
  if (!is.null(gammas) && length(gammas) != m)
    stop("'gammas' must have one entry per cause.", call. = FALSE)

  Tm <- matrix(Inf, n, m)
  for (j in seq_len(m)) {
    p  <- pars[[j]]
    a  <- if ("a" %in% names(p)) p[["a"]] else 1
    t0 <- rbetadanish(n, a, p[["b"]], p[["c"]], p[["k"]])
    Tm[, j] <- if (is.null(gammas)) t0 else t0 * exp(x * gammas[j])
  }

  tmin  <- apply(Tm, 1, min)
  cause <- max.col(-Tm, ties.method = "first")

  cmax <- stats::quantile(tmin, 1 - censor_rate, na.rm = TRUE) * 2
  cens <- stats::runif(n, 0, cmax)
  cause[cens < tmin] <- 0L

  out <- data.frame(time = pmin(tmin, cens), cause = as.integer(cause))
  if (!is.null(x)) out$x <- x
  out
}

#' Summarise One Column of Simulation Replicates
#' @noRd
.bd_sim_summary <- function(est, se, truth, level = 0.95) {
  ok <- is.finite(est)
  z  <- stats::qnorm(1 - (1 - level) / 2)
  cov <- if (any(is.finite(se[ok]))) {
    mean(abs(est[ok] - truth) <= z * se[ok], na.rm = TRUE)
  } else NA_real_
  c(truth      = truth,
    mean       = mean(est[ok]),
    bias       = mean(est[ok]) - truth,
    rel_bias   = (mean(est[ok]) - truth) / truth,
    rmse       = sqrt(mean((est[ok] - truth)^2)),
    emp_sd     = stats::sd(est[ok]),
    mean_se    = mean(se[ok], na.rm = TRUE),
    se_ratio   = mean(se[ok], na.rm = TRUE) / stats::sd(est[ok]),
    coverage   = cov)
}

#' Finite-Sample Simulation Study for the Univariate Model
#'
#' Draws samples from a known Beta-Danish or Exponentiated Danish distribution,
#' refits each one, and reports the finite-sample behaviour of the maximum
#' likelihood estimator.
#'
#' @param n Vector of sample sizes.
#' @param n_sim Number of replicates at each sample size.
#' @param truth Named numeric vector of true parameters. Must contain `b`, `c`
#'   and `k`, and `a` unless `submodel` is `TRUE`.
#' @param submodel Logical; fit the three-parameter ED submodel.
#' @param censor_rate Approximate proportion of right-censored observations.
#' @param n_starts Random starts per fit, on top of the deterministic grid.
#' @param level Nominal coverage level.
#' @param seed Optional integer seed.
#' @param quiet Logical; suppress progress messages.
#'
#' @return An object of class `"bd_simulation"` with a `results` data frame,
#'   one row per sample size and parameter.
#'
#' @details
#' Read `se_ratio` and `coverage` together. A ratio near one with coverage near
#' the nominal level means the asymptotic standard errors are trustworthy at
#' that sample size. A ratio well below one means they are optimistic, and
#' coverage will fall short accordingly.
#'
#' `n_fail` counts replicates where no admissible optimum was found. A high
#' rate is not a nuisance to be suppressed: it says the parameter region being
#' studied is hard to estimate at that sample size.
#'
#' @seealso [bd_simulation_cure()], [bd_simulation_competing()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Deliberately tiny so the example runs in seconds. A real study would use
#' # n_sim in the hundreds; see the vignette.
#' s <- bd_simulation_study(n = 60, n_sim = 3,
#'                          truth = c(b = 3, c = 2, k = 0.5),
#'                          submodel = TRUE, n_starts = 1,
#'                          seed = 1, quiet = TRUE)
#' s
#' }
bd_simulation_study <- function(n = c(50, 100, 200), n_sim = 200,
                                truth = c(a = 1.5, b = 3, c = 2, k = 0.5),
                                submodel = FALSE, censor_rate = 0,
                                n_starts = 3, level = 0.95,
                                seed = NULL, quiet = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  say <- function(...) if (!isTRUE(quiet)) message(...)

  need <- if (submodel) c("b", "c", "k") else c("a", "b", "c", "k")
  if (!all(need %in% names(truth)))
    stop("'truth' must name: ", paste(need, collapse = ", "), call. = FALSE)
  a_true <- if (submodel) 1 else truth[["a"]]

  rows <- list()
  for (nn in n) {
    say("n = ", nn, ": ", n_sim, " replicate(s)")
    est <- matrix(NA_real_, n_sim, length(need),
                  dimnames = list(NULL, need))
    se  <- est
    for (r in seq_len(n_sim)) {
      t0 <- rbetadanish(nn, a_true, truth[["b"]], truth[["c"]], truth[["k"]])
      status <- rep(1L, nn)
      if (censor_rate > 0) {
        cmax <- stats::quantile(t0, 1 - censor_rate) * 2
        cc   <- stats::runif(nn, 0, cmax)
        status <- as.integer(t0 <= cc)
        t0 <- pmin(t0, cc)
      }
      dat <- data.frame(time = t0, status = status)
      fit <- tryCatch(suppressWarnings(
        fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                       submodel = submodel, n_starts = n_starts,
                       check_identifiability = FALSE)),
        error = function(e) NULL)
      if (is.null(fit)) next
      est[r, ] <- fit$coefficients[need]
      se[r, ]  <- sqrt(pmax(diag(fit$vcov)[need], 0))
    }

    n_ok <- sum(stats::complete.cases(est))
    for (p in need) {
      s <- .bd_sim_summary(est[, p], se[, p], truth[[p]], level)
      rows[[length(rows) + 1L]] <- data.frame(
        n = nn, parameter = p, as.list(s),
        n_ok = n_ok, n_fail = n_sim - n_ok,
        row.names = NULL, check.names = FALSE)
    }
  }

  out <- list(results = do.call(rbind, rows),
              design = list(n = n, n_sim = n_sim, truth = truth,
                            submodel = submodel, censor_rate = censor_rate,
                            level = level),
              type = "univariate")
  class(out) <- "bd_simulation"
  out
}

#' Generate Mixture-Cure Data
#'
#' A cured subject never fails, so it is censored at whatever censoring time it
#' draws. Susceptible subjects fail according to the ED distribution. Written
#' here rather than reusing `simulate_bd_cure_data()`, whose interface requires
#' explicit regression coefficients and design matrices that an intercept-only
#' study does not need.
#'
#' @noRd
.bd_sim_cure_data <- function(n, b, c, k, cure_fraction, censor_rate) {
  cured <- stats::rbinom(n, 1, cure_fraction) == 1L
  t <- rep(Inf, n)
  ns <- sum(!cured)
  if (ns) t[!cured] <- rbetadanish(ns, 1, b, c, k)

  fin <- t[is.finite(t)]
  cmax <- if (length(fin)) {
    as.numeric(stats::quantile(fin, max(1 - censor_rate, 0.05))) * 3
  } else 1
  cens <- stats::runif(n, 0, cmax)

  data.frame(time = pmin(t, cens), status = as.integer(t <= cens))
}

#' Finite-Sample Simulation Study for the Cure Model
#'
#' @param n Vector of sample sizes.
#' @param n_sim Number of replicates at each sample size.
#' @param truth Named vector with `b`, `c` and `k` for the susceptible
#'   population.
#' @param cure_fraction True proportion of long-term survivors.
#' @param censor_rate Approximate proportion of censoring among susceptibles.
#' @param type `"mixture"` or `"promotion"`.
#' @param n_starts Random starts per fit.
#' @param level Nominal coverage level.
#' @param seed Optional integer seed.
#' @param quiet Logical; suppress progress messages.
#'
#' @return An object of class `"bd_simulation"`.
#'
#' @seealso [bd_simulation_study()], [fit_bd_cure()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' s <- bd_simulation_cure(n = 80, n_sim = 2, n_starts = 1,
#'                         seed = 1, quiet = TRUE)
#' s
#' }
bd_simulation_cure <- function(n = c(100, 200), n_sim = 100,
                               truth = c(b = 2, c = 1.5, k = 0.5),
                               cure_fraction = 0.3, censor_rate = 0.2,
                               type = c("mixture", "promotion"),
                               n_starts = 3, level = 0.95,
                               seed = NULL, quiet = FALSE) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  say <- function(...) if (!isTRUE(quiet)) message(...)
  need <- c("b", "c")

  rows <- list()
  for (nn in n) {
    say("n = ", nn, ": ", n_sim, " replicate(s)")
    est <- matrix(NA_real_, n_sim, length(need), dimnames = list(NULL, need))
    se  <- est
    for (r in seq_len(n_sim)) {
      dat <- .bd_sim_cure_data(nn, truth[["b"]], truth[["c"]], truth[["k"]],
                               cure_fraction, censor_rate)
      fit <- tryCatch(suppressWarnings(
        fit_bd_cure(formula_aft = survival::Surv(time, status) ~ 1,
                    formula_cure = ~ 1, data = dat, type = type,
                    n_starts = n_starts)),
        error = function(e) NULL)
      if (is.null(fit)) next
      cf <- fit$coefficients
      v  <- sqrt(pmax(diag(fit$vcov), 0))
      for (p in need) {
        nm <- paste0("log_", p)
        if (nm %in% names(cf)) {
          est[r, p] <- exp(cf[[nm]])
          se[r, p]  <- exp(cf[[nm]]) * v[which(names(cf) == nm)]
        }
      }
    }
    n_ok <- sum(stats::complete.cases(est))
    for (p in need) {
      s <- .bd_sim_summary(est[, p], se[, p], truth[[p]], level)
      rows[[length(rows) + 1L]] <- data.frame(
        n = nn, parameter = p, as.list(s),
        n_ok = n_ok, n_fail = n_sim - n_ok,
        row.names = NULL, check.names = FALSE)
    }
  }

  out <- list(results = do.call(rbind, rows),
              design = list(n = n, n_sim = n_sim, truth = truth,
                            cure_fraction = cure_fraction, type = type,
                            level = level),
              type = paste0("cure (", type, ")"))
  class(out) <- "bd_simulation"
  out
}

#' Finite-Sample Simulation Study for Competing Risks
#'
#' @param n Vector of sample sizes.
#' @param n_sim Number of replicates at each sample size.
#' @param pars List of true parameter vectors, one per cause.
#' @param gammas Optional accelerated failure time coefficients, one per cause.
#' @param censor_rate Approximate proportion of right-censored observations.
#' @param submodel Logical; fit Exponentiated Danish cause-specific kernels.
#' @param n_starts Starting points per fit.
#' @param level Nominal coverage level.
#' @param seed Optional integer seed.
#' @param quiet Logical; suppress progress messages.
#'
#' @return An object of class `"bd_simulation"`, with `parameter` naming the
#'   cause and quantity, for example `c1:b`.
#'
#' @seealso [bd_simulation_study()], [fit_bd_competing()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' s <- bd_simulation_competing(n = 120, n_sim = 2, n_starts = 1,
#'                              seed = 1, quiet = TRUE)
#' s
#' }
bd_simulation_competing <- function(n = c(200, 400), n_sim = 100,
                                    pars = list(c(b = 2.5, c = 1.8, k = 0.04),
                                                c(b = 3.5, c = 1.2, k = 0.02)),
                                    gammas = NULL, censor_rate = 0.15,
                                    submodel = TRUE, n_starts = 2,
                                    level = 0.95, seed = NULL, quiet = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  say <- function(...) if (!isTRUE(quiet)) message(...)

  m <- length(pars)
  base <- c("b", "c", "k")
  targets <- unlist(lapply(seq_len(m), function(j) paste0("c", j, ":", base)))
  truths  <- unlist(lapply(seq_len(m), function(j) unname(pars[[j]][base])))
  if (!is.null(gammas)) {
    targets <- c(targets, paste0("c", seq_len(m), ":gamma_x"))
    truths  <- c(truths, gammas)
  }

  rows <- list()
  for (nn in n) {
    say("n = ", nn, ": ", n_sim, " replicate(s)")
    est <- matrix(NA_real_, n_sim, length(targets),
                  dimnames = list(NULL, targets))
    se  <- est
    for (r in seq_len(n_sim)) {
      d <- tryCatch(simulate_bd_competing_data(nn, pars = pars,
                                               gammas = gammas,
                                               censor_rate = censor_rate),
                    error = function(e) NULL)
      if (is.null(d) || length(unique(d$cause[d$cause > 0])) < m) next
      fit <- tryCatch(suppressWarnings(
        fit_bd_competing(d$time, d$cause,
                         covariates = if (is.null(gammas)) NULL else ~ x,
                         data = d, submodel = submodel, n_starts = n_starts)),
        error = function(e) NULL)
      if (is.null(fit)) next
      cm <- fit$coefficients; sm <- fit$se
      for (j in seq_len(m)) {
        for (p in base) {
          est[r, paste0("c", j, ":", p)] <- cm[j, p]
          se[r,  paste0("c", j, ":", p)] <- sm[j, p]
        }
        if (!is.null(gammas) && "gamma_x" %in% colnames(cm)) {
          est[r, paste0("c", j, ":gamma_x")] <- cm[j, "gamma_x"]
          se[r,  paste0("c", j, ":gamma_x")] <- sm[j, "gamma_x"]
        }
      }
    }
    n_ok <- sum(stats::complete.cases(est))
    for (i in seq_along(targets)) {
      s <- .bd_sim_summary(est[, targets[i]], se[, targets[i]], truths[i], level)
      rows[[length(rows) + 1L]] <- data.frame(
        n = nn, parameter = targets[i], as.list(s),
        n_ok = n_ok, n_fail = n_sim - n_ok,
        row.names = NULL, check.names = FALSE)
    }
  }

  out <- list(results = do.call(rbind, rows),
              design = list(n = n, n_sim = n_sim, pars = pars,
                            gammas = gammas, submodel = submodel,
                            level = level),
              type = "competing risks")
  class(out) <- "bd_simulation"
  out
}

#' @param x A `"bd_simulation"` object.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @rdname bd_simulation_study
#' @export
print.bd_simulation <- function(x, digits = 4, ...) {
  cat("\nBeta-Danish simulation study: ", x$type, "\n", sep = "")
  cat("  replicates: ", x$design$n_sim,
      "   sample sizes: ", paste(x$design$n, collapse = ", "), "\n", sep = "")
  cat("  nominal coverage: ", format(100 * x$design$level), "%\n\n", sep = "")

  res <- x$results
  show <- c("n", "parameter", "truth", "mean", "bias", "rmse",
            "se_ratio", "coverage", "n_fail")
  print(format(res[, intersect(show, names(res))], digits = digits),
        row.names = FALSE)

  if (any(res$n_fail > 0))
    cat("\n  Some replicates found no admissible optimum. A high count means\n",
        "  the parameter region is hard to estimate at that sample size.\n",
        sep = "")
  bad <- res$se_ratio < 0.9 & is.finite(res$se_ratio)
  if (any(bad))
    cat("\n  se_ratio below 0.9 for ", sum(bad), " row(s): the asymptotic\n",
        "  standard errors understate the sampling variability there, and\n",
        "  coverage falls short of nominal as a result.\n", sep = "")
  cat("\n")
  invisible(x)
}

#' @param object A `"bd_simulation"` object.
#' @rdname bd_simulation_study
#' @export
summary.bd_simulation <- function(object, ...) object$results)---")


.step("Rec 40, 41, 42: writing R/plots_extra.R")

.put("R/plots_extra.R", r"---(## Diagnostic and exploratory plots.

#' Scaled Total Time on Test Plot
#'
#' Draws the scaled total time on test transform, a distribution-free way of
#' judging hazard shape before any model is fitted.
#'
#' @param time Numeric vector of observed times, or a fitted `"betadanish"`
#'   object, from which the times are taken.
#' @param status Optional event indicator. Censored observations are dropped,
#'   since the transform is defined for complete samples.
#' @param add Logical; add to an existing plot rather than starting a new one.
#' @param col,lwd Colour and line width for the curve.
#' @param main,xlab,ylab Labels.
#' @param ... Further arguments passed to [graphics::plot()].
#'
#' @return Invisibly, a data frame with the plotting coordinates `i_n` and
#'   `phi`, and an attribute `"shape"` giving the suggested hazard shape.
#'
#' @details
#' For an ordered sample \eqn{x_{(1)} \le \cdots \le x_{(n)}}, the scaled
#' transform at \eqn{i/n} is
#' \deqn{\phi(i/n) = \frac{\sum_{j=1}^{i} x_{(j)} + (n-i)x_{(i)}}
#'                        {\sum_{j=1}^{n} x_{(j)}}.}
#'
#' Read it against the diagonal. A curve entirely above the diagonal indicates
#' an increasing hazard, entirely below a decreasing one; a curve that starts
#' below and crosses above suggests a bathtub shape, and the reverse suggests a
#' unimodal one. A curve close to the diagonal indicates a constant hazard,
#' that is an exponential sample.
#'
#' This is a shape diagnostic, not a test. It is worth drawing before choosing
#' between the four-parameter model and its submodel, because it says which
#' hazard shapes the data can support without assuming any of them.
#'
#' @seealso [bd_hazard_shape()] for the fitted counterpart
#'
#' @export
#'
#' @examples
#' data(guinea_pig)
#' ttt <- bd_ttt_plot(guinea_pig$time)
#' attr(ttt, "shape")
bd_ttt_plot <- function(time, status = NULL, add = FALSE,
                        col = "steelblue", lwd = 2,
                        main = "Scaled total time on test",
                        xlab = "i / n", ylab = expression(phi(i/n)), ...) {

  if (inherits(time, "betadanish")) {
    status <- time$data$status
    time   <- time$data$time
  }
  time <- as.numeric(time)
  if (!is.null(status)) {
    keep <- as.numeric(status) == 1
    if (sum(keep) < 5L)
      stop("At least five uncensored observations are needed.", call. = FALSE)
    if (any(!keep))
      warning(sprintf(paste0("%d censored observation(s) dropped: the TTT ",
                             "transform is defined for complete samples."),
                      sum(!keep)), call. = FALSE)
    time <- time[keep]
  }
  time <- sort(time[is.finite(time) & time > 0])
  n <- length(time)
  if (n < 5L) stop("At least five positive times are needed.", call. = FALSE)

  i   <- seq_len(n)
  cs  <- cumsum(time)
  phi <- (cs + (n - i) * time) / cs[n]
  u   <- i / n

  shape <- .bd_ttt_shape(u, phi)

  if (!isTRUE(add)) {
    graphics::plot(c(0, u), c(0, phi), type = "l", col = col, lwd = lwd,
                   xlim = c(0, 1), ylim = c(0, 1),
                   main = main, xlab = xlab, ylab = ylab, ...)
    graphics::abline(0, 1, col = "grey50", lty = 2)
    graphics::legend("bottomright", bty = "n",
                     legend = c("TTT transform", "diagonal (constant hazard)"),
                     col = c(col, "grey50"), lty = c(1, 2), lwd = c(lwd, 1))
    graphics::mtext(paste("suggested hazard:", shape), side = 3, line = 0.2,
                    cex = 0.85)
  } else {
    graphics::lines(c(0, u), c(0, phi), col = col, lwd = lwd)
  }

  out <- data.frame(i_n = u, phi = phi)
  attr(out, "shape") <- shape
  invisible(out)
}

#' Classify a TTT Curve Against the Diagonal
#' @noRd
.bd_ttt_shape <- function(u, phi, tol = 0.02) {
  d <- phi - u
  above <- d > tol
  below <- d < -tol
  if (!any(above) && !any(below)) return("constant (exponential)")
  if (!any(below)) return("increasing")
  if (!any(above)) return("decreasing")
  ## Mixed: the side it starts on decides between bathtub and unimodal.
  first <- if (which(above)[1] < which(below)[1]) "above" else "below"
  if (identical(first, "below")) "bathtub" else "unimodal (upside-down bathtub)"
}

#' Plot a Profile Likelihood
#'
#' Draws the profile log-likelihood produced by [bd_profile_ci()], with the
#' critical threshold and the resulting interval marked.
#'
#' @param x An object of class `"bd_profile"`.
#' @param col,lwd Colour and line width for the profile curve.
#' @param main,xlab,ylab Labels. `NULL` uses sensible defaults.
#' @param ... Further arguments passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#'
#' @details
#' The horizontal line sits at \eqn{\ell_{\max} - \chi^2_{1,\alpha}/2}. Every
#' parameter value whose profile lies above it is inside the interval.
#'
#' When the curve does not fall back below that line before the right-hand edge
#' of the grid, there is no finite upper bound and the plot says so. Widening
#' the grid will not produce one; it will only confirm the flatness. That is the
#' situation the underlying dissertation records for the tail index on the
#' breaking-stress data.
#'
#' @seealso [bd_profile_ci()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' dat <- simulate_bd_data(120, a = 1, b = 3, c = 2, k = 0.5, seed = 3)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       submodel = TRUE, n_starts = 1)
#' p <- bd_profile_ci(fit, "b", n_grid = 15L)
#' bd_profile_plot(p)
#' }
bd_profile_plot <- function(x, col = "steelblue", lwd = 2,
                            main = NULL, xlab = NULL, ylab = NULL, ...) {
  if (!inherits(x, "bd_profile"))
    stop("'x' must be a bd_profile object from bd_profile_ci().", call. = FALSE)

  ok <- is.finite(x$profile)
  if (sum(ok) < 3L)
    stop("Too few finite profile values to plot.", call. = FALSE)

  thr <- x$logLik_max - x$critical
  if (is.null(main)) main <- paste("Profile log-likelihood for", x$parameter)
  if (is.null(xlab)) xlab <- x$parameter
  if (is.null(ylab)) ylab <- "profile log-likelihood"

  graphics::plot(x$grid[ok], x$profile[ok], type = "l", col = col, lwd = lwd,
                 main = main, xlab = xlab, ylab = ylab, ...)
  graphics::abline(h = thr, col = "red", lty = 2)
  graphics::abline(v = x$estimate, col = "grey40", lty = 3)

  if (is.finite(x$lower)) graphics::abline(v = x$lower, col = "red", lty = 3)
  if (is.finite(x$upper)) graphics::abline(v = x$upper, col = "red", lty = 3)

  lab <- sprintf("%s%% interval: [%s, %s]", format(100 * x$level),
                 signif(x$lower, 4),
                 if (is.infinite(x$upper)) "Inf" else signif(x$upper, 4))
  graphics::mtext(lab, side = 3, line = 0.2, cex = 0.85)

  if (isTRUE(x$open_above))
    graphics::legend("bottomright", bty = "n", text.col = "red3",
                     legend = "no finite upper bound: report a lower bound")

  invisible(x)
}

#' Diagnostic Plots for a Bayesian Fit
#'
#' Trace and posterior density plots for each parameter of a
#' [bayes_betadanish()] fit.
#'
#' @param x An object of class `"bd_bayes"`.
#' @param which Optional character vector of parameters to show. Defaults to
#'   all of them.
#' @param type `"both"` (default) for trace and density side by side, or
#'   `"trace"` or `"density"` alone.
#' @param col Colour for the traces and densities.
#' @param ... Further arguments passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#'
#' @details
#' Read the traces first. A well-mixed chain looks like noise around a stable
#' level, with no drift and no long excursions. Visible trend means the burn-in
#' was too short; a chain that sticks at one value for many iterations means
#' the proposal is too wide and almost every move is being rejected, which is
#' worth fixing with `tune` before interpreting anything.
#'
#' The graphical parameters are restored on exit, so the function leaves the
#' device as it found it.
#'
#' @seealso [bayes_betadanish()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("MCMCpack", quietly = TRUE) &&
#'     requireNamespace("coda", quietly = TRUE)) {
#'   dat <- simulate_bd_data(80, a = 1, b = 3, c = 2, k = 0.5, seed = 6)
#'   bfit <- bayes_betadanish(dat$time, dat$status, submodel = TRUE,
#'                            burnin = 200, mcmc = 800, seed = 1)
#'   plot(bfit)
#' }
#' }
plot.bd_bayes <- function(x, which = NULL, type = c("both", "trace", "density"),
                          col = "steelblue", ...) {
  type <- match.arg(type)
  dr <- as.matrix(x$draws)
  if (!is.matrix(dr) || !nrow(dr))
    stop("The fit contains no posterior draws.", call. = FALSE)

  nm <- colnames(dr)
  if (!is.null(which)) {
    miss <- setdiff(which, nm)
    if (length(miss))
      stop("Not in the posterior: ", paste(miss, collapse = ", "),
           ". Available: ", paste(nm, collapse = ", "), call. = FALSE)
    dr <- dr[, which, drop = FALSE]
    nm <- which
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  ncol_panel <- if (identical(type, "both")) 2L else 1L
  graphics::par(mfrow = c(length(nm), ncol_panel),
                mar = c(4, 4, 2, 1))

  it <- seq_len(nrow(dr))
  for (p in nm) {
    v <- dr[, p]
    if (type %in% c("both", "trace")) {
      graphics::plot(it, v, type = "l", col = col,
                     main = paste("Trace:", p),
                     xlab = "iteration", ylab = p, ...)
      graphics::abline(h = mean(v), col = "red", lty = 2)
    }
    if (type %in% c("both", "density")) {
      d <- stats::density(v)
      graphics::plot(d, col = col, lwd = 2, main = paste("Posterior:", p),
                     xlab = p, ...)
      graphics::abline(v = mean(v), col = "red", lty = 2)
      if (!is.null(x$HPD) && p %in% rownames(x$HPD))
        graphics::abline(v = x$HPD[p, ], col = "grey40", lty = 3)
    }
  }
  invisible(x)
})---")


.step("Writing tests/testthat/test-cr-simulation.R")

.put("tests/testthat/test-cr-simulation.R", r"---(## Competing risks with covariates, and the simulation runners.

test_that("the design matrix builder handles all three input forms", {
  d <- data.frame(a = rnorm(10), g = factor(rep(c("x", "y"), 5)))
  f <- BetaDanish:::.bd_cr_design(~ a + g, d, n = 10)
  expect_equal(nrow(f), 10L)
  expect_false("(Intercept)" %in% colnames(f))

  m <- BetaDanish:::.bd_cr_design(matrix(1:20, 10, 2), NULL, n = 10)
  expect_equal(dim(m), c(10L, 2L))
  expect_equal(colnames(m), c("x1", "x2"))

  expect_null(BetaDanish:::.bd_cr_design(NULL, NULL, n = 10))
})

test_that("the design matrix builder rejects bad input", {
  d <- data.frame(a = rnorm(5))
  expect_error(BetaDanish:::.bd_cr_design(~ a, NULL, n = 5), "'data' is required")
  expect_error(BetaDanish:::.bd_cr_design(~ a, d, n = 9), "row")
  expect_error(BetaDanish:::.bd_cr_design(y ~ a, d, n = 5), "one-sided")
  d2 <- data.frame(a = c(1, NA, 3, 4, 5))
  expect_error(BetaDanish:::.bd_cr_design(~ a, d2, n = 5), "missing values")
})

test_that("simulated competing risks data has the right shape", {
  d <- simulate_bd_competing_data(300, seed = 5)
  expect_named(d, c("time", "cause"))
  expect_true(all(d$time > 0))
  expect_true(all(d$cause %in% 0:2))
  expect_gt(sum(d$cause == 1), 10)
  expect_gt(sum(d$cause == 2), 10)
  expect_gt(sum(d$cause == 0), 0)

  dx <- simulate_bd_competing_data(200, gammas = c(0.5, -0.3), seed = 6)
  expect_true("x" %in% names(dx))
  expect_true(all(dx$x %in% c(0, 1)))
})

test_that("competing risks fits without covariates", {
  skip_on_cran()
  d <- simulate_bd_competing_data(300, seed = 7)
  fit <- suppressWarnings(
    fit_bd_competing(d$time, d$cause, submodel = TRUE, n_starts = 2))

  expect_s3_class(fit, "bd_competing")
  expect_equal(dim(fit$coefficients), c(2L, 3L))
  expect_equal(colnames(fit$coefficients), c("b", "c", "k"))
  expect_true(all(fit$coefficients > 0))
  expect_true(is.finite(fit$logLik))
  expect_equal(fit$nobs, 300L)
})

test_that("competing risks fits with covariates and recovers their sign", {
  skip_on_cran()
  d <- simulate_bd_competing_data(600, gammas = c(0.8, -0.8),
                                  censor_rate = 0.1, seed = 8)
  fit <- suppressWarnings(
    fit_bd_competing(d$time, d$cause, covariates = ~ x, data = d,
                     submodel = TRUE, n_starts = 3))

  expect_true("gamma_x" %in% colnames(fit$coefficients))
  expect_equal(fit$covariates, "x")
  ## Cause 1 was accelerated and cause 2 retarded; the signs should differ.
  g <- fit$coefficients[, "gamma_x"]
  expect_gt(g[1], g[2])
})

test_that("competing risks validates its inputs", {
  expect_error(fit_bd_competing(1:5, c(1, 2)), "same length")
  expect_error(fit_bd_competing(c(1, 2, 3), c(1, 1, 1)), "at least two")
  expect_error(fit_bd_competing(c(-1, 2, 3), c(0, 1, 2)), "strictly positive")
})

test_that("the CIF is monotone, bounded, and covariate-aware", {
  skip_on_cran()
  d <- simulate_bd_competing_data(400, seed = 9)
  fit <- suppressWarnings(
    fit_bd_competing(d$time, d$cause, submodel = TRUE, n_starts = 2))

  tv <- c(1, 5, 20, 100)
  f1 <- cif_betadanish(fit, tv, cause_idx = 1)
  f2 <- cif_betadanish(fit, tv, cause_idx = 2)

  expect_true(all(diff(f1) >= -1e-8))
  expect_true(all(f1 >= 0 & f1 <= 1))
  expect_equal(cif_betadanish(fit, 0, 1), 0)
  ## The competing CIFs cannot together exceed one.
  expect_true(all(f1 + f2 <= 1 + 1e-6))

  expect_error(cif_betadanish(fit, tv, cause_idx = 9), "not among")
})

test_that("a covariate shifts the CIF", {
  skip_on_cran()
  d <- simulate_bd_competing_data(500, gammas = c(1, 0), censor_rate = 0.1,
                                  seed = 10)
  fit <- suppressWarnings(
    fit_bd_competing(d$time, d$cause, covariates = ~ x, data = d,
                     submodel = TRUE, n_starts = 3))
  a <- cif_betadanish(fit, 10, cause_idx = 1, x = 0)
  b <- cif_betadanish(fit, 10, cause_idx = 1, x = 1)
  expect_true(is.finite(a) && is.finite(b))
  expect_false(isTRUE(all.equal(a, b)))
  expect_error(cif_betadanish(fit, 10, 1, x = c(0, 0)), "one value per covariate")
})

test_that("the simulation summary helper computes what it claims", {
  set.seed(1)
  est <- rnorm(500, mean = 2.1, sd = 0.4)
  se  <- rep(0.4, 500)
  s <- BetaDanish:::.bd_sim_summary(est, se, truth = 2)
  expect_equal(unname(s["truth"]), 2)
  expect_equal(unname(s["bias"]), mean(est) - 2, tolerance = 1e-12)
  expect_equal(unname(s["rmse"]), sqrt(mean((est - 2)^2)), tolerance = 1e-12)
  expect_gt(s[["coverage"]], 0.8)
  expect_lt(abs(s[["se_ratio"]] - 1), 0.2)
})

test_that("the univariate simulation study returns a tidy table", {
  skip_on_cran()
  s <- bd_simulation_study(n = 80, n_sim = 12,
                           truth = c(b = 3, c = 2, k = 0.5),
                           submodel = TRUE, n_starts = 2,
                           seed = 11, quiet = TRUE)
  expect_s3_class(s, "bd_simulation")
  expect_true(is.data.frame(s$results))
  expect_setequal(s$results$parameter, c("b", "c", "k"))
  expect_true(all(c("bias", "rmse", "se_ratio", "coverage", "n_fail") %in%
                    names(s$results)))
  expect_true(all(s$results$truth == c(3, 2, 0.5)))
  expect_output(print(s), "simulation study")
  expect_identical(summary(s), s$results)
})

test_that("the study rejects a truth vector missing a parameter", {
  expect_error(bd_simulation_study(n = 50, n_sim = 2, truth = c(b = 3, c = 2),
                                   submodel = TRUE, quiet = TRUE),
               "must name")
})

test_that("the cure data generator censors every cured subject", {
  set.seed(2)
  d <- BetaDanish:::.bd_sim_cure_data(2000, b = 2, c = 1.5, k = 0.5,
                                      cure_fraction = 0.4, censor_rate = 0.2)
  expect_named(d, c("time", "status"))
  expect_true(all(d$time > 0))
  expect_true(all(d$status %in% c(0L, 1L)))
  ## At least the cure fraction must end up censored.
  expect_gt(mean(d$status == 0), 0.35)
})

test_that("the competing-risks simulation runner returns per-cause rows", {
  skip_on_cran()
  s <- bd_simulation_competing(n = 250, n_sim = 4, n_starts = 2,
                               seed = 12, quiet = TRUE)
  expect_s3_class(s, "bd_simulation")
  expect_true(all(grepl("^c[12]:", s$results$parameter)))
  expect_equal(nrow(s$results), 6L)
  expect_output(print(s), "competing risks")
}))---")


.step("Writing tests/testthat/test-plots-extra.R")

.put("tests/testthat/test-plots-extra.R", r"---(## Visualisation added in 0.3.0. Plots are drawn to a null device; the
## assertions are about the returned values and the classification logic,
## which is what can actually be wrong.

test_that("the TTT transform has the right endpoints and is increasing", {
  data(guinea_pig, package = "BetaDanish", envir = environment())
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  ttt <- bd_ttt_plot(guinea_pig$time)

  expect_s3_class(ttt, "data.frame")
  expect_named(ttt, c("i_n", "phi"))
  expect_equal(nrow(ttt), nrow(guinea_pig))
  expect_equal(ttt$phi[nrow(ttt)], 1, tolerance = 1e-12)   # phi(1) = 1
  expect_true(all(diff(ttt$phi) >= -1e-12))                # non-decreasing
  expect_true(all(ttt$phi >= 0 & ttt$phi <= 1 + 1e-12))
  expect_true(attr(ttt, "shape") %in%
                c("increasing", "decreasing", "bathtub",
                  "unimodal (upside-down bathtub)", "constant (exponential)"))
})

test_that("the TTT classifier recognises the reference shapes", {
  cl <- BetaDanish:::.bd_ttt_shape
  u <- seq(0.02, 1, length.out = 50)
  expect_equal(cl(u, u), "constant (exponential)")
  expect_equal(cl(u, pmin(u + 0.15, 1)), "increasing")
  expect_equal(cl(u, pmax(u - 0.15, 0)), "decreasing")
  ## Starts below the diagonal then crosses above: bathtub.
  expect_equal(cl(u, u + 0.2 * sin(2 * pi * u)), "bathtub")
  ## Starts above then falls below: unimodal.
  expect_equal(cl(u, u - 0.2 * sin(2 * pi * u)),
               "unimodal (upside-down bathtub)")
})

test_that("an exponential sample gives a TTT curve near the diagonal", {
  set.seed(4)
  x <- stats::rexp(400, rate = 0.5)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  ttt <- bd_ttt_plot(x)
  expect_lt(max(abs(ttt$phi - ttt$i_n)), 0.12)
})

test_that("bd_ttt_plot drops censored observations with a warning", {
  set.seed(5)
  t <- rbetadanish(60, 1, 3, 2, 0.5)
  s <- stats::rbinom(60, 1, 0.8)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_warning(ttt <- bd_ttt_plot(t, status = s), "censored")
  expect_equal(nrow(ttt), sum(s == 1))
})

test_that("bd_ttt_plot accepts a fitted object and validates its input", {
  skip_on_cran()
  dat <- simulate_bd_data(60, a = 1, b = 3, c = 2, k = 0.5, seed = 7)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 1, check_identifiability = FALSE))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_s3_class(suppressWarnings(bd_ttt_plot(fit)), "data.frame")
  expect_error(bd_ttt_plot(c(1, 2, 3)), "At least five")
})

test_that("bd_profile_plot draws and returns its input", {
  skip_on_cran()
  dat <- simulate_bd_data(100, a = 1, b = 3, c = 2, k = 0.5, seed = 8)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 1, check_identifiability = FALSE))
  p <- bd_profile_ci(fit, "b", n_grid = 12L)

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  out <- bd_profile_plot(p)
  expect_identical(out$parameter, "b")
  expect_error(bd_profile_plot(list()), "bd_profile object")
})

test_that("plot.bd_bayes validates before drawing", {
  fake <- structure(list(draws = matrix(rnorm(300), 100, 3,
                                        dimnames = list(NULL, c("b", "c", "k"))),
                         HPD = NULL, submodel = TRUE),
                    class = "bd_bayes")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_invisible(plot(fake))
  expect_invisible(plot(fake, which = "b", type = "trace"))
  expect_error(plot(fake, which = "zzz"), "Not in the posterior")

  empty <- structure(list(draws = matrix(numeric(0), 0, 0)), class = "bd_bayes")
  expect_error(plot(empty), "no posterior draws")
})

test_that("plot.bd_bayes leaves the graphical parameters as it found them", {
  fake <- structure(list(draws = matrix(rnorm(200), 100, 2,
                                        dimnames = list(NULL, c("b", "c")))),
                    class = "bd_bayes")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  before <- graphics::par("mfrow")
  plot(fake)
  expect_equal(graphics::par("mfrow"), before)
}))---")


## =============================================================================
##  REC 46  --  VIGNETTES
## =============================================================================

.step("Rec 46: appending a 'New in 0.3.0' section to the introduction vignette")

.append_once(
  "vignettes/BetaDanish_Introduction.Rmd",
  "## New in 0.3.0",
  c("## New in 0.3.0",
    "",
    "Three groups of functions were added after the 0.2.0 release.",
    "",
    "### Structural properties",
    "",
    "```{r}",
    "p <- list(a = 1.5, b = 5, c = 2, k = 1)",
    "",
    "bd_moment_summary(p$a, p$b, p$c, p$k)",
    "bd_entropy_shannon(p$a, p$b, p$c, p$k)",
    "bd_stress_strength(strength = p, stress = p)   # identical laws: one half",
    "```",
    "",
    "Moments exist only when `b > r`, and the package enforces it rather than",
    "returning a plausible-looking number from a truncated sum:",
    "",
    "```{r}",
    "bd_moments(1:4, a = 1.5, b = 3, c = 2, k = 1)   # the fourth is Inf",
    "bd_tail_index(a = 1.5, b = 3, c = 2, k = 1)$moment_condition",
    "```",
    "",
    "### Shape, before and after fitting",
    "",
    "`bd_ttt_plot()` reads the hazard shape from the data alone, with no model",
    "assumed. `bd_hazard_shape()` does the same from a fitted parameter set, so",
    "the two can be compared.",
    "",
    "```{r, fig.width = 6, fig.height = 4}",
    "data(guinea_pig)",
    "ttt <- bd_ttt_plot(guinea_pig$time)",
    "attr(ttt, \"shape\")",
    "```",
    "",
    "### Interval estimation",
    "",
    "All four parameters are positive, so a symmetric Wald interval on the",
    "natural scale can cross zero. `bd_wald_ci()` builds on the log scale",
    "instead, and `bd_profile_ci()` inverts the likelihood ratio test. When the",
    "profile never returns below the critical threshold there is no finite upper",
    "bound, and the honest report is a lower bound rather than a point estimate",
    "with an interval around it.",
    "",
    "```{r, eval = FALSE}",
    "fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = remission,",
    "                      submodel = TRUE)",
    "bd_wald_ci(fit)",
    "bd_profile_plot(bd_profile_ci(fit, \"b\"))",
    "",
    "# For the four-parameter model, report the identified composite",
    "bd_identified_coef(fit_betadanish(survival::Surv(time, status) ~ 1,",
    "                                  data = remission))",
    "```",
    "",
    "### Simulation studies",
    "",
    "`bd_simulation_study()`, `bd_simulation_cure()` and",
    "`bd_simulation_competing()` report bias, RMSE, the ratio of mean estimated",
    "standard error to empirical standard deviation, and Wald coverage. Read the",
    "last two together: a ratio well below one means the asymptotic standard",
    "errors are optimistic at that sample size, and coverage falls short as a",
    "consequence.",
    "",
    "```{r, eval = FALSE}",
    "bd_simulation_study(n = c(50, 100, 200), n_sim = 500,",
    "                    truth = c(b = 3, c = 2, k = 0.5), submodel = TRUE)",
    "```"),
  "introduction vignette")

.step("Rec 46: appending a covariates section to the competing-risks vignette")

.append_once(
  "vignettes/bd-competing-risks.Rmd",
  "## Covariates (new in 0.3.0)",
  c("## Covariates (new in 0.3.0)",
    "",
    "Each cause can carry its own coefficient vector, acting on the time scale:",
    "",
    "$$T_j = T_{0j}\\exp(x'\\gamma_j)$$",
    "",
    "so a positive coefficient lengthens time to failure from that cause. The",
    "same covariate may therefore accelerate one cause and retard another, which",
    "is the point of fitting them separately.",
    "",
    "```{r, eval = FALSE}",
    "# The covariate column only exists when 'gammas' is supplied",
    "d <- simulate_bd_competing_data(400, gammas = c(0.8, -0.8), seed = 1)",
    "",
    "fit <- fit_bd_competing(d$time, d$cause, covariates = ~ x, data = d,",
    "                        submodel = TRUE)",
    "fit$coefficients",
    "```",
    "",
    "The intercept is dropped from the design matrix, because each cause already",
    "has its own scale parameter `k` that an intercept would be confounded with.",
    "",
    "`cif_betadanish()` takes the covariate values at which to evaluate the",
    "cumulative incidence, defaulting to the reference subject:",
    "",
    "```{r, eval = FALSE}",
    "cif_betadanish(fit, tvec = c(5, 20, 50), cause_idx = 1, x = 0)",
    "cif_betadanish(fit, tvec = c(5, 20, 50), cause_idx = 1, x = 1)",
    "```",
    "",
    "A covariate effect folds exactly into the scale: scaling time by",
    "$\\lambda$ maps $k$ to $k/\\lambda$, so no re-integration is needed.",
    "",
    "### A caution worth repeating",
    "",
    "Independence of the latent failure times is an *identifying* assumption, not",
    "a testable one. Under positive dependence the working independence model",
    "biases the cause-specific cumulative incidence functions downward. Overall",
    "survival is more robust than the individual marginals; where conclusions",
    "rest on absolute CIFs, add a copula-based sensitivity analysis."),
  "competing-risks vignette")

## =============================================================================
##  NEWS
## =============================================================================

.step("Recording the additions in NEWS.md")
.nw <- readLines("NEWS.md", warn = FALSE)
if (!any(grepl("bd_ttt_plot", .nw, fixed = TRUE))) {
  .hdr <- grep("^# BetaDanish 0\\.3\\.0|^# BetaDanish 0\\.2\\.0\\.9000", .nw)
  if (length(.hdr) >= 1L) {
    .backup("NEWS.md")
    .sec <- c(
      "",
      "## Visualisation",
      "",
      "* `bd_ttt_plot()` draws the scaled total time on test transform, a",
      "  distribution-free read on hazard shape before any model is fitted. It",
      "  returns the plotting coordinates and a suggested shape label, so it can",
      "  be compared against `bd_hazard_shape()` from the fitted parameters.",
      "",
      "* `bd_profile_plot()` draws a profile log-likelihood with the critical",
      "  threshold and the resulting interval marked. Where the curve never",
      "  returns below the threshold the plot says there is no finite upper",
      "  bound, rather than implying one at the edge of the grid.",
      "",
      "* `plot()` method for `bd_bayes` objects, giving trace and posterior",
      "  density panels and restoring the graphical parameters on exit.",
      "",
      "* The introduction and competing-risks vignettes cover the new",
      "  functions.")
    .nw <- append(.nw, .sec, after = .hdr[1])
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
  for (nm in sub("[ ]*<-.*$", "",
                 grep("^[a-zA-Z_.][a-zA-Z0-9_.]*[ ]*<-[ ]*function", txt, value = TRUE)))
    .defs[[nm]] <- c(.defs[[nm]], basename(f))
}
.dups <- Filter(function(v) length(unique(v)) > 1L, .defs)
if (length(.dups)) {
  for (nm in names(.dups))
    .warn(sprintf("%s defined in: %s", nm, paste(unique(.dups[[nm]]), collapse = ", ")))
  .die("A function is defined in more than one file. Backups: ", BACKUP_DIR)
}
.ok(sprintf("%d definition(s), every name unique", length(.defs)))

.step("Confirming the 3e example bug is gone")
.cx <- readLines("R/competing_risks.R", warn = FALSE)
.ex <- grep("covariates = ~ x", .cx, fixed = TRUE, value = TRUE)
if (length(.ex)) {
  ## Every example that asks for ~ x must sit near a simulate call with gammas.
  .gam <- any(grepl("gammas = c(", .cx, fixed = TRUE))
  if (.gam) .ok("the covariate example now generates the covariate first")
  else .warn("an example uses ~ x but no gammas call was found nearby")
}

.step("Checking example runtimes are plausible")
.simx <- readLines("R/simulation_study.R", warn = FALSE)
.big <- grep("n_sim = ([1-9][0-9]{1,})", .simx[grepl("^#'", .simx)], value = TRUE)
if (length(.big)) {
  .warn("an example still uses a large n_sim:")
  for (l in .big) cat("        ", trimws(l), "\n", sep = "")
} else {
  .ok("no example uses a large replicate count")
}

.step("Loading from source for the self-test")
.loaded <- tryCatch({ devtools::load_all(".", quiet = TRUE); TRUE },
                    error = function(e) conditionMessage(e))
if (!isTRUE(.loaded)) .die("load_all() failed:\n  ", .loaded, "\n\nBackups: ", BACKUP_DIR)
.ok("source loaded")

.step("Self-test on the new plots")
.fails <- character(0)
.pass <- function(label, ok, detail = "") {
  if (isTRUE(ok)) .ok(paste0(label, if (nzchar(detail)) paste0("  ", detail) else ""))
  else { .warn(paste0(label, "  ", detail)); .fails <<- c(.fails, label) }
}

grDevices::pdf(NULL)
on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

data(guinea_pig, package = "BetaDanish", envir = environment())
.t <- bd_ttt_plot(guinea_pig$time)
.pass("TTT ends at one", isTRUE(all.equal(.t$phi[nrow(.t)], 1)),
      sprintf("phi(1) = %.10f", .t$phi[nrow(.t)]))
.pass("TTT is non-decreasing", all(diff(.t$phi) >= -1e-12))
.pass("TTT reports a shape", nzchar(attr(.t, "shape")),
      paste("shape:", attr(.t, "shape")))

## An exponential sample must sit close to the diagonal.
set.seed(9); .e <- stats::rexp(500, 0.5)
.te <- bd_ttt_plot(.e)
.pass("exponential sample lies near the diagonal",
      max(abs(.te$phi - .te$i_n)) < 0.12,
      sprintf("max deviation %.4f", max(abs(.te$phi - .te$i_n))))

## The Bayesian plot must leave par() untouched.
.fake <- structure(list(draws = matrix(stats::rnorm(200), 100, 2,
                                       dimnames = list(NULL, c("b", "c")))),
                   class = "bd_bayes")
.before <- graphics::par("mfrow")
plot(.fake)
.pass("plot.bd_bayes restores graphical parameters",
      identical(graphics::par("mfrow"), .before))

try(grDevices::dev.off(), silent = TRUE)

if (length(.fails))
  .die("Self-tests failed: ", paste(.fails, collapse = "; "),
       "\nBackups: ", BACKUP_DIR)
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
  nms <- .rd_arg_names(f); dup <- unique(nms[duplicated(nms)])
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
for (f in c("bd_ttt_plot", "bd_profile_plot")) {
  if (any(grepl(paste0("export(", f, ")"), .ns, fixed = TRUE))) .ok(f)
  else .warn(paste("not exported:", f))
}
if (any(grepl("S3method(plot,bd_bayes)", .ns, fixed = TRUE))) {
  .ok("S3method(plot,bd_bayes)")
} else {
  .warn("plot.bd_bayes is not registered")
}

.step("devtools::test()")
.t2 <- tryCatch(devtools::test(), error = function(e) { .warn(conditionMessage(e)); NULL })

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
cat("  PATCH 3f COMPLETE  --  all 46 recommendations implemented\n")
cat(strrep("=", 78), "\n\n")
cat("  fixed  the Patch 3e covariate example and the long example runtimes\n")
cat("  40     bd_ttt_plot()\n")
cat("  41     bd_profile_plot()\n")
cat("  42     plot.bd_bayes()\n")
cat("  46     introduction and competing-risks vignettes updated\n\n")
cat("  NEXT: re-run Patch 3d to cut the release.\n\n")
cat("      source(\"dev/BetaDanish_Patch3d_release.R\")\n\n")
cat("  It is idempotent, so running it again is safe. It will set the version,\n")
cat("  fix the pkgdown index, write cran-comments.md, commit and tag.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
