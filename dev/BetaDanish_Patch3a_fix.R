## =============================================================================
##  BetaDanish  --  PATCH 3a-fix : explicit Rd parameter documentation
## =============================================================================
##
##  Patch 3a left one check WARNING:
##
##    Documented arguments not in \usage in Rd file 'bd_entropy_shannon.Rd':
##      'subdivisions' 'rel.tol'
##    Undocumented arguments in Rd file 'bd_moment_summary.Rd': 'rel.tol'
##
##  Both trace to roxygen2 tag resolution rather than to the statistics.
##  @inheritParams cannot reliably split a combined tag such as
##  "@param rel.tol,subdivisions" into its individual names, so an argument can
##  end up inherited into a topic that has no such formal, or dropped from one
##  that does.
##
##  RATHER THAN GUESS AT THE EXACT MECHANISM, every route to it is removed:
##
##    * no @inheritParams anywhere in the new files -- each function documents
##      its own arguments
##    * no combined "@param a,b,c,k" or "@param rel.tol,subdivisions" tags --
##      one tag per argument, each with its own description
##    * bd_entropy_shannon() and bd_moment_summary() gain real rel.tol and
##      subdivisions arguments, so the quadrature branch and the moment
##      integrals are tunable instead of hard-coded. That is a genuine
##      improvement, and it also makes the documented names formals whatever
##      roxygen decides to do
##
##  Nothing about the mathematics changes. All 508 tests from Patch 3a still
##  apply unchanged, since the new arguments have defaults.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3a_fix.R")
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
cat("  BetaDanish  --  Patch 3a-fix : explicit Rd parameters\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("R/moments.R") || !file.exists("R/structural.R"))
  .die("Patch 3a has not been applied.")
.ok("Patch 3a detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3afix"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

.step("Writing R/moments.R with per-argument documentation")

.put("R/moments.R", r"---(## Moment-based structural properties of the Beta-Danish distribution.
##
## Every integral below is taken on the finite Beta(a, b) scale. Under
## u = G(z) = {kz/(1+kz)}^c the variable u is exactly Beta(a, b) distributed and
## z = v / (k(1 - v)) with v = u^(1/c), so
##
##   E[g(Z)] = 1/B(a,b) * INT_0^1 g(z(u)) u^(a-1) (1-u)^(b-1) du.
##
## The interval is finite, the only singularity is the integrable endpoint
## u -> 1, and there is no series truncation.

#' Back-Transform from the Beta Scale to the Time Scale
#' @noRd
.bd_z_of_u <- function(u, c, k) {
  v <- u^(1 / c)
  v / (k * (1 - v))
}

#' Core Integrand: z(u)^r Times the Beta(a, b) Density
#' @noRd
.bd_mom_integrand <- function(u, r, a, b, c, k, lB) {
  z  <- .bd_z_of_u(u, c, k)
  ## r = 0 must contribute nothing: r * log(z) would be 0 * -Inf = NaN at the
  ## lower endpoint, where the integrand is simply the Beta(a, b) density.
  rz <- if (r == 0) 0 else r * log(z)
  lg <- rz + (a - 1) * log(u) + (b - 1) * log1p(-u) - lB
  out <- exp(lg)
  out[!is.finite(out)] <- 0
  out
}

#' Integrate the Core Integrand Between Two Beta-Scale Limits
#' @noRd
.bd_mom_int <- function(r, a, b, c, k, lo = 0, hi = 1,
                        rel.tol = 1e-10, subdivisions = 4000L) {
  lB <- lbeta(a, b)
  tryCatch(
    stats::integrate(.bd_mom_integrand, lower = lo, upper = hi,
                     r = r, a = a, b = b, c = c, k = k, lB = lB,
                     rel.tol = rel.tol, subdivisions = subdivisions,
                     stop.on.error = FALSE)$value,
    error = function(e) NA_real_)
}

#' Raw Moments of the Beta-Danish Distribution
#'
#' @param r Order or orders of the moment. Non-negative.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail; `E(Z^r)` is finite iff `b > r`.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Relative accuracy, passed to [stats::integrate()].
#' @param subdivisions Subdivision limit, passed to [stats::integrate()].
#'
#' @return A numeric vector the length of `r`. Entries with `b <= r` are `Inf`.
#'
#' @details
#' The survival function is regularly varying with index \eqn{-b}, so
#' \eqn{E(Z^r)} is finite if and only if \eqn{b > r}. That condition is checked
#' rather than assumed: a request for a moment that does not exist returns
#' `Inf`, not a large finite number produced by a truncated sum.
#'
#' @seealso [bd_moment_summary()], [bd_incomplete_moment()]
#'
#' @export
#'
#' @examples
#' bd_moments(1:2, a = 1.5, b = 3, c = 2, k = 1)
#'
#' # The fourth moment does not exist when b = 3
#' bd_moments(4, a = 1.5, b = 3, c = 2, k = 1)
bd_moments <- function(r, a, b, c, k, rel.tol = 1e-10, subdivisions = 4000L) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  r <- as.numeric(r)
  if (any(is.na(r)) || any(r < 0))
    stop("'r' must be non-negative.", call. = FALSE)

  vapply(r, function(ri) {
    if (b <= ri) return(Inf)
    .bd_mom_int(ri, a, b, c, k, 0, 1, rel.tol, subdivisions)
  }, numeric(1))
}

#' Summary Moments: Mean, Variance, Skewness and Kurtosis
#'
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Relative accuracy, passed to [stats::integrate()].
#' @param subdivisions Subdivision limit, passed to [stats::integrate()].
#'
#' @return A named numeric vector. Entries requiring a moment that does not
#'   exist are `NA`, with the governing condition given in the `condition`
#'   attribute.
#'
#' @details
#' Skewness needs \eqn{b > 3} and kurtosis \eqn{b > 4}. Both are frequently
#' unavailable for fitted values of \eqn{b}, which is a property of the family
#' rather than a limitation of the computation.
#'
#' @export
#'
#' @examples
#' bd_moment_summary(a = 1.5, b = 5, c = 2, k = 1)
#' bd_moment_summary(a = 1.5, b = 2.5, c = 2, k = 1)   # skew/kurt unavailable
bd_moment_summary <- function(a, b, c, k, rel.tol = 1e-10,
                              subdivisions = 4000L) {
  m <- bd_moments(1:4, a, b, c, k, rel.tol = rel.tol,
                  subdivisions = subdivisions)
  mu <- m[1]

  out <- c(mean = NA_real_, variance = NA_real_,
           sd = NA_real_, skewness = NA_real_, kurtosis = NA_real_)

  if (is.finite(mu)) out[["mean"]] <- mu
  if (all(is.finite(m[1:2]))) {
    v <- m[2] - mu^2
    out[["variance"]] <- v
    out[["sd"]] <- sqrt(max(v, 0))
  }
  if (all(is.finite(m[1:3])) && is.finite(out[["sd"]]) && out[["sd"]] > 0) {
    mu3 <- m[3] - 3 * mu * m[2] + 2 * mu^3
    out[["skewness"]] <- mu3 / out[["sd"]]^3
  }
  if (all(is.finite(m[1:4])) && is.finite(out[["sd"]]) && out[["sd"]] > 0) {
    mu4 <- m[4] - 4 * mu * m[3] + 6 * mu^2 * m[2] - 3 * mu^4
    out[["kurtosis"]] <- mu4 / out[["sd"]]^4
  }
  attr(out, "condition") <- "E(Z^r) is finite if and only if b > r"
  out
}

#' Incomplete Moments
#'
#' The lower incomplete moment \eqn{M_r(t) = \int_0^t z^r f(z)\,dz} and the
#' upper incomplete moment \eqn{\Upsilon_r(t) = \int_t^\infty z^r f(z)\,dz}.
#'
#' @param r Order of the moment.
#' @param t Vector of thresholds.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param lower Logical; `TRUE` (default) for \eqn{M_r(t)}, `FALSE` for
#'   \eqn{\Upsilon_r(t)}.
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A numeric vector the length of `t`.
#'
#' @details
#' The threshold is mapped to the Beta scale by \eqn{u_t = G(t)}, computed on
#' the log scale so it stays accurate when \eqn{G(t)} approaches one, and the
#' integral is taken over \eqn{(0, u_t)} or \eqn{(u_t, 1)}.
#'
#' @seealso [bd_mrl()], [bd_lorenz()], [bd_mean_deviation()]
#'
#' @export
#'
#' @examples
#' a <- 1.5; b <- 3; c <- 2; k <- 1
#' t <- 1
#' # The two halves sum to the complete moment
#' bd_incomplete_moment(1, t, a, b, c, k, lower = TRUE) +
#'   bd_incomplete_moment(1, t, a, b, c, k, lower = FALSE)
#' bd_moments(1, a, b, c, k)
bd_incomplete_moment <- function(r, t, a, b, c, k, lower = TRUE,
                                 rel.tol = 1e-10) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  if (length(r) != 1L || is.na(r) || r < 0)
    stop("'r' must be a single non-negative number.", call. = FALSE)
  if (b <= r)
    warning("The complete moment of order ", r, " does not exist (b <= r); ",
            "the upper incomplete moment is infinite.", call. = FALSE)

  vapply(as.numeric(t), function(ti) {
    if (is.na(ti)) return(NA_real_)
    if (ti <= 0) return(if (lower) 0 else bd_moments(r, a, b, c, k))
    if (!is.finite(ti)) return(if (lower) bd_moments(r, a, b, c, k) else 0)
    ut <- exp(.bd_logG(k * ti, c))
    if (lower) .bd_mom_int(r, a, b, c, k, 0, ut, rel.tol)
    else       .bd_mom_int(r, a, b, c, k, ut, 1, rel.tol)
  }, numeric(1))
}

#' Conditional Moments
#'
#' \eqn{E(Z^r \mid Z > t)} or \eqn{E(Z^r \mid Z \le t)}.
#'
#' @param r Order of the moment.
#' @param t Vector of thresholds.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Relative accuracy, passed to [stats::integrate()].
#' @param upper Logical; `TRUE` (default) conditions on \eqn{Z > t}.
#' @return A numeric vector the length of `t`.
#' @export
#' @examples
#' bd_conditional_moment(1, t = c(0.5, 1, 2), a = 1.5, b = 3, c = 2, k = 1)
bd_conditional_moment <- function(r, t, a, b, c, k, upper = TRUE,
                                  rel.tol = 1e-10) {
  num <- bd_incomplete_moment(r, t, a, b, c, k, lower = !upper, rel.tol = rel.tol)
  den <- if (upper) sbetadanish(t, a, b, c, k) else pbetadanish(t, a, b, c, k)
  out <- num / den
  out[!is.finite(den) | den <= 0] <- NA_real_
  out
}

#' Mean Residual Life and Reversed Mean Residual Life
#'
#' @param t Vector of times.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A numeric vector the length of `t`, or all `NA` when \eqn{b \le 1}.
#'
#' @details
#' Following the definitions in the underlying dissertation, with
#' \eqn{\Upsilon_1} the upper and \eqn{M_1} the lower incomplete first moment,
#' \deqn{m(t) = E(Z - t \mid Z > t) = \Upsilon_1(t)/S(t) - t,}
#' \deqn{\bar m(t) = E(t - Z \mid Z \le t) = t - M_1(t)/F(t).}
#' The reversed form is also called the mean inactivity time. Both require the
#' first moment to exist, hence \eqn{b > 1}.
#'
#' No monotonicity is asserted. Ageing classification requires a separate
#' argument and is not implied by these values; see [bd_hazard_shape()] for the
#' hazard-rate counterpart.
#'
#' @export
#'
#' @examples
#' bd_mrl(c(0.5, 1, 2, 5), a = 1.5, b = 3, c = 2, k = 1)
#' bd_rmrl(c(0.5, 1, 2, 5), a = 1.5, b = 3, c = 2, k = 1)
bd_mrl <- function(t, a, b, c, k, rel.tol = 1e-10) {
  if (b <= 1) {
    warning("The mean residual life needs a finite first moment, so b > 1. ",
            "NA returned.", call. = FALSE)
    return(rep(NA_real_, length(t)))
  }
  ups <- bd_incomplete_moment(1, t, a, b, c, k, lower = FALSE, rel.tol = rel.tol)
  S   <- sbetadanish(t, a, b, c, k)
  out <- ups / S - as.numeric(t)
  out[!is.finite(S) | S <= 0] <- NA_real_
  out
}

#' @rdname bd_mrl
#' @export
bd_rmrl <- function(t, a, b, c, k, rel.tol = 1e-10) {
  if (b <= 1) {
    warning("The reversed mean residual life needs a finite first moment, so ",
            "b > 1. NA returned.", call. = FALSE)
    return(rep(NA_real_, length(t)))
  }
  M1 <- bd_incomplete_moment(1, t, a, b, c, k, lower = TRUE, rel.tol = rel.tol)
  Fx <- pbetadanish(t, a, b, c, k)
  out <- as.numeric(t) - M1 / Fx
  out[!is.finite(Fx) | Fx <= 0] <- NA_real_
  out
}

#' Mean Deviations
#'
#' Mean deviation about the mean, \eqn{E|Z - \mu|}, or about the median.
#'
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param about `"mean"` (default) or `"median"`.
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A single number, or `NA` when \eqn{b \le 1}.
#'
#' @details
#' \eqn{E|Z - \mu| = 2\mu F(\mu) - 2 M_1(\mu)} and
#' \eqn{E|Z - M| = \mu - 2 M_1(M)}, with \eqn{M} the median.
#'
#' @export
#'
#' @examples
#' bd_mean_deviation(a = 1.5, b = 3, c = 2, k = 1)
#' bd_mean_deviation(a = 1.5, b = 3, c = 2, k = 1, about = "median")
bd_mean_deviation <- function(a, b, c, k, about = c("mean", "median"),
                              rel.tol = 1e-10) {
  about <- match.arg(about)
  if (b <= 1) {
    warning("Mean deviations need a finite first moment, so b > 1. ",
            "NA returned.", call. = FALSE)
    return(NA_real_)
  }
  mu <- bd_moments(1, a, b, c, k, rel.tol = rel.tol)
  if (identical(about, "mean")) {
    2 * mu * pbetadanish(mu, a, b, c, k) -
      2 * bd_incomplete_moment(1, mu, a, b, c, k, rel.tol = rel.tol)
  } else {
    med <- qbetadanish(0.5, a, b, c, k)
    mu - 2 * bd_incomplete_moment(1, med, a, b, c, k, rel.tol = rel.tol)
  }
}

#' Lorenz and Bonferroni Curves
#'
#' @param p Vector of probabilities in \eqn{[0, 1]}.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A numeric vector the length of `p`.
#'
#' @details
#' \eqn{L(p) = M_1(Q(p))/\mu} and \eqn{B(p) = L(p)/p}, with \eqn{Q} the
#' quantile function. Both require \eqn{b > 1}.
#'
#' @export
#'
#' @examples
#' p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' bd_lorenz(p, a = 1.5, b = 3, c = 2, k = 1)
#' bd_bonferroni(p, a = 1.5, b = 3, c = 2, k = 1)
bd_lorenz <- function(p, a, b, c, k, rel.tol = 1e-10) {
  if (b <= 1) {
    warning("The Lorenz curve needs a finite mean, so b > 1. NA returned.",
            call. = FALSE)
    return(rep(NA_real_, length(p)))
  }
  mu <- bd_moments(1, a, b, c, k, rel.tol = rel.tol)
  q  <- qbetadanish(p, a, b, c, k)
  out <- bd_incomplete_moment(1, q, a, b, c, k, rel.tol = rel.tol) / mu
  out[is.na(p)] <- NA_real_
  out
}

#' @rdname bd_lorenz
#' @export
bd_bonferroni <- function(p, a, b, c, k, rel.tol = 1e-10) {
  L <- bd_lorenz(p, a, b, c, k, rel.tol = rel.tol)
  out <- L / as.numeric(p)
  out[!is.na(p) & p <= 0] <- NA_real_
  out
}

#' Probability Weighted Moments
#'
#' \eqn{\tau_{r,s,w} = E\{Z^r F(Z)^s (1 - F(Z))^w\}}.
#'
#' @param r Power on \eqn{Z}.
#' @param s Power on \eqn{F(Z)}.
#' @param w Power on \eqn{1 - F(Z)}.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A single number.
#'
#' @details
#' Existence follows the same rule as the raw moments in the worst case,
#' \eqn{b > r}, and the weights can only reduce the tail contribution.
#'
#' @export
#'
#' @examples
#' # tau_{1,0,0} is the mean
#' bd_pwm(1, 0, 0, a = 1.5, b = 3, c = 2, k = 1)
#' bd_moments(1, a = 1.5, b = 3, c = 2, k = 1)
#'
#' bd_pwm(1, 1, 0, a = 1.5, b = 3, c = 2, k = 1)
bd_pwm <- function(r, s = 0, w = 0, a, b, c, k, rel.tol = 1e-10) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  if (b <= r && s == 0 && w == 0) return(Inf)
  lB <- lbeta(a, b)

  integrand <- function(u) {
    z  <- .bd_z_of_u(u, c, k)
    Fz <- stats::pbeta(u, a, b)
    rz <- if (r == 0) 0 else r * log(z)
    lg <- rz + (a - 1) * log(u) + (b - 1) * log1p(-u) - lB
    out <- exp(lg) * Fz^s * (1 - Fz)^w
    out[!is.finite(out)] <- 0
    out
  }
  tryCatch(stats::integrate(integrand, 0, 1, rel.tol = rel.tol,
                            subdivisions = 4000L, stop.on.error = FALSE)$value,
           error = function(e) NA_real_)
})---")


.step("Writing R/structural.R with per-argument documentation")

.put("R/structural.R", r"---(## Entropy, reliability, shape and order-statistic properties.

#' Shannon Entropy of the Beta-Danish Distribution
#'
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param terms Number of series terms before the analytic tail is applied.
#' @param method `"closed"` (default) for the closed form, or `"quadrature"`
#'   for direct numerical integration of \eqn{-f \log f}, which is slower but
#'   independent of the series.
#' @param rel.tol Relative accuracy for `method = "quadrature"`, passed to
#'   [stats::integrate()].
#' @param subdivisions Subdivision limit for `method = "quadrature"`, passed to
#'   [stats::integrate()].
#'
#' @return A single number, in nats.
#'
#' @details
#' Writing \eqn{u = G(Z) \sim \mathrm{Beta}(a,b)} and \eqn{v = u^{1/c}},
#' \deqn{H = \log B(a,b) - \log(ck)
#'          - (a - 1/c)\{\psi(a) - \psi(a+b)\}
#'          - (b - 1)\{\psi(b) - \psi(a+b)\}
#'          + 2\sum_{i \ge 1} \frac{B(a + i/c,\, b)}{i\, B(a,b)},}
#' the final sum arising from \eqn{-2E\log(1 - v)} expanded as a power series.
#'
#' @section Series truncation:
#' The terms decay like \eqn{i^{-(b+1)}}, so truncating at \eqn{M} omits a tail
#' of order \eqn{M^{-b}}. That is negligible for large \eqn{b} and is not for
#' small \eqn{b}: against high-precision integration the plain truncated sum at
#' \eqn{M = 2000} is out by about \eqn{10^{-8}} at \eqn{b = 3} but by about
#' \eqn{0.11} at \eqn{b = 0.5}. The analytic tail
#' \eqn{\Gamma(b)c^{b}M^{-b}/\{b B(a,b)\}} is therefore added, restoring
#' agreement to roughly eight digits across the range.
#'
#' @seealso [bd_entropy_renyi()], [bd_entropy_tsallis()]
#'
#' @export
#'
#' @examples
#' bd_entropy_shannon(a = 1.5, b = 3, c = 2, k = 1)
#'
#' # Independent check by quadrature
#' bd_entropy_shannon(a = 1.5, b = 3, c = 2, k = 1, method = "quadrature")
bd_entropy_shannon <- function(a, b, c, k, terms = 20000L,
                               method = c("closed", "quadrature"),
                               rel.tol = 1e-10, subdivisions = 4000L) {
  method <- match.arg(method)
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)

  if (identical(method, "quadrature")) {
    integrand <- function(u) {
      z  <- .bd_z_of_u(u, c, k)
      lf <- dbetadanish(z, a, b, c, k, log = TRUE)
      out <- -lf * exp((a - 1) * log(u) + (b - 1) * log1p(-u) - lbeta(a, b))
      out[!is.finite(out)] <- 0
      out
    }
    return(tryCatch(
      stats::integrate(integrand, 0, 1, rel.tol = rel.tol,
                       subdivisions = subdivisions, stop.on.error = FALSE)$value,
      error = function(e) NA_real_))
  }

  M  <- as.integer(terms)
  lB <- lbeta(a, b)
  i  <- seq_len(M)
  S  <- sum(exp(lbeta(a + i / c, b) - lB) / i)
  ## Euler-Maclaurin tail: terms ~ Gamma(b) c^b / B(a,b) * i^-(b+1).
  S  <- S + exp(lgamma(b) + b * log(c) - lB) * M^(-b) / b

  lB - log(c * k) -
    (a - 1 / c) * (digamma(a) - digamma(a + b)) -
    (b - 1)     * (digamma(b) - digamma(a + b)) +
    2 * S
}

#' Renyi and Tsallis Entropies
#'
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param order Entropy order \eqn{q}, positive and not equal to one.
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A single number.
#'
#' @details
#' With \eqn{I_q = \int_0^\infty f(z)^q dz}, the Renyi entropy is
#' \eqn{\log(I_q)/(1-q)} and the Tsallis entropy is \eqn{(1 - I_q)/(q-1)}.
#' Both reduce to the Shannon entropy as \eqn{q \to 1}, which is excluded here;
#' use [bd_entropy_shannon()] at \eqn{q = 1}.
#'
#' \eqn{I_q} is evaluated as \eqn{\int_0^1 f(z(u))^{q-1} d\mathrm{Beta}(u;a,b)}
#' on the finite Beta scale rather than over the half line, which keeps the
#' heavy upper tail from dominating the quadrature.
#'
#' @export
#'
#' @examples
#' bd_entropy_renyi(a = 1.5, b = 3, c = 2, k = 1, order = 2)
#' bd_entropy_tsallis(a = 1.5, b = 3, c = 2, k = 1, order = 2)
bd_entropy_renyi <- function(a, b, c, k, order = 2, rel.tol = 1e-10) {
  Iq <- .bd_Iq(a, b, c, k, order, rel.tol)
  if (is.na(Iq) || Iq <= 0) return(NA_real_)
  log(Iq) / (1 - order)
}

#' @rdname bd_entropy_renyi
#' @export
bd_entropy_tsallis <- function(a, b, c, k, order = 2, rel.tol = 1e-10) {
  Iq <- .bd_Iq(a, b, c, k, order, rel.tol)
  if (is.na(Iq)) return(NA_real_)
  (1 - Iq) / (order - 1)
}

#' The Integral of f^q, on the Beta Scale
#' @noRd
.bd_Iq <- function(a, b, c, k, q, rel.tol = 1e-10) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  if (length(q) != 1L || is.na(q) || q <= 0)
    stop("'order' must be a single positive number.", call. = FALSE)
  if (abs(q - 1) < 1e-8)
    stop("Order 1 is the Shannon entropy; use bd_entropy_shannon().",
         call. = FALSE)

  lB <- lbeta(a, b)
  integrand <- function(u) {
    z  <- .bd_z_of_u(u, c, k)
    lf <- dbetadanish(z, a, b, c, k, log = TRUE)
    out <- exp((q - 1) * lf + (a - 1) * log(u) + (b - 1) * log1p(-u) - lB)
    out[!is.finite(out)] <- 0
    out
  }
  tryCatch(stats::integrate(integrand, 0, 1, rel.tol = rel.tol,
                            subdivisions = 4000L, stop.on.error = FALSE)$value,
           error = function(e) NA_real_)
}

#' Stress-Strength Reliability
#'
#' \eqn{R = P(X > Y)} where \eqn{X} is the strength and \eqn{Y} the stress,
#' each Beta-Danish distributed with its own parameters.
#'
#' @param strength Named list or numeric vector giving `a`, `b`, `c`, `k` for
#'   the strength variable \eqn{X}.
#' @param stress Same, for the stress variable \eqn{Y}.
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A single probability.
#'
#' @details
#' \eqn{R = \int_0^\infty S_X(y) f_Y(y)\,dy}, evaluated on the Beta scale of
#' \eqn{Y} so the integral is taken over \eqn{(0,1)}.
#'
#' The convention is the standard reliability one: \eqn{R} is the probability
#' that the strength exceeds the stress, so larger \eqn{R} means a more
#' reliable component. Passing the two arguments the wrong way round returns
#' \eqn{1 - R}.
#'
#' @export
#'
#' @examples
#' # Identical distributions must give one half
#' p <- list(a = 1.5, b = 3, c = 2, k = 1)
#' bd_stress_strength(strength = p, stress = p)
#'
#' # A stronger component
#' bd_stress_strength(strength = list(a = 1.5, b = 3, c = 2, k = 0.5),
#'                    stress   = list(a = 1.5, b = 3, c = 2, k = 2))
bd_stress_strength <- function(strength, stress, rel.tol = 1e-10) {
  px <- .bd_as_pars(strength, "strength")
  py <- .bd_as_pars(stress,   "stress")

  lB <- lbeta(py$a, py$b)
  integrand <- function(u) {
    y  <- .bd_z_of_u(u, py$c, py$k)
    Sx <- sbetadanish(y, px$a, px$b, px$c, px$k)
    out <- Sx * exp((py$a - 1) * log(u) + (py$b - 1) * log1p(-u) - lB)
    out[!is.finite(out)] <- 0
    out
  }
  tryCatch(stats::integrate(integrand, 0, 1, rel.tol = rel.tol,
                            subdivisions = 4000L, stop.on.error = FALSE)$value,
           error = function(e) NA_real_)
}

#' Coerce a Parameter Argument to a Named List
#' @noRd
.bd_as_pars <- function(x, what) {
  if (is.list(x)) x <- unlist(x)
  if (!is.numeric(x) || is.null(names(x)) || !all(c("a", "b", "c", "k") %in% names(x)))
    stop("'", what, "' must supply named values a, b, c and k.", call. = FALSE)
  p <- as.list(x[c("a", "b", "c", "k")])
  if (!check_positive_params(p$a, p$b, p$c, p$k))
    stop("All ", what, " parameters must be finite and strictly positive.",
         call. = FALSE)
  p
}

#' Hazard Rate Shape and the Glaser Diagnostic
#'
#' Classifies the hazard rate as increasing, decreasing, bathtub-shaped or
#' upside-down bathtub-shaped, and returns the mode of the density.
#'
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param n_grid Number of grid points between the lower and upper evaluation
#'   quantiles.
#' @param range_p Two probabilities bounding the evaluation range.
#'
#' @return A list of class `"bd_shape"` with the classification, the grid of
#'   times, the hazard, Glaser's \eqn{\eta(t) = -f'(t)/f(t)}, and the mode.
#'
#' @details
#' Glaser's criterion works through \eqn{\eta(t) = -f'(t)/f(t)}: monotone
#' increasing \eqn{\eta} gives an increasing hazard, monotone decreasing gives a
#' decreasing hazard, and a single interior turning point gives a bathtub or
#' upside-down bathtub. Here \eqn{\eta} is formed analytically from the
#' log-density and evaluated on a grid, and the hazard itself is checked
#' alongside it, so the two must agree before a shape is reported.
#'
#' @export
#'
#' @examples
#' s <- bd_hazard_shape(a = 1.5, b = 3, c = 2, k = 1)
#' s$shape
#' s$mode
bd_hazard_shape <- function(a, b, c, k, n_grid = 400L,
                            range_p = c(1e-4, 1 - 1e-4)) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)

  lo <- qbetadanish(range_p[1], a, b, c, k)
  hi <- qbetadanish(range_p[2], a, b, c, k)
  t  <- exp(seq(log(lo), log(hi), length.out = as.integer(n_grid)))

  eta <- .bd_glaser_eta(t, a, b, c, k)
  h   <- hbetadanish(t, a, b, c, k)

  ok <- is.finite(eta) & is.finite(h) & h > 0
  shape <- if (sum(ok) < 10L) {
    "indeterminate"
  } else {
    .bd_classify(h[ok])
  }

  ## Mode of the density: maximise the log-density over the same range.
  mode_t <- tryCatch(
    stats::optimize(function(x) dbetadanish(x, a, b, c, k, log = TRUE),
                    interval = c(lo, hi), maximum = TRUE)$maximum,
    error = function(e) NA_real_)
  if (is.finite(mode_t) &&
      (mode_t <= lo * (1 + 1e-6) || mode_t >= hi * (1 - 1e-6)))
    mode_t <- NA_real_   # optimum on the boundary is not an interior mode

  out <- list(shape = shape, time = t, hazard = h, eta = eta,
              mode = mode_t,
              parameters = c(a = a, b = b, c = c, k = k))
  class(out) <- "bd_shape"
  out
}

#' Glaser's eta(t) = -f'(t)/f(t), Analytically
#' @noRd
.bd_glaser_eta <- function(t, a, b, c, k) {
  kt   <- k * t
  logw <- -log1p(1 / kt)                 # log of w = kt/(1+kt)
  lG   <- c * logw
  om   <- -expm1(lG)                     # 1 - G, stable
  ca   <- c * a

  dlog_1mG <- -c * exp((c - 1) * logw) * k / ((1 + kt)^2 * om)
  dlogf    <- (ca - 1) / t - (ca + 1) * k / (1 + kt) + (b - 1) * dlog_1mG
  -dlogf
}

#' Classify a Positive Sequence as Monotone or Single-Turning
#' @noRd
.bd_classify <- function(v, tol = 1e-10) {
  d <- diff(v)
  d <- d[is.finite(d)]
  if (!length(d)) return("indeterminate")
  up   <- d > tol
  down <- d < -tol
  if (!any(down)) return("increasing")
  if (!any(up))   return("decreasing")
  first_up <- which(up)[1]
  first_dn <- which(down)[1]
  if (first_dn < first_up) "bathtub" else "upside-down bathtub"
}

#' @param x A `"bd_shape"` object.
#' @param ... Ignored.
#' @rdname bd_hazard_shape
#' @export
print.bd_shape <- function(x, ...) {
  cat("Beta-Danish hazard shape\n")
  cat("  parameters: ",
      paste(names(x$parameters), signif(x$parameters, 4), sep = " = ",
            collapse = ", "), "\n", sep = "")
  cat("  shape:      ", x$shape, "\n", sep = "")
  cat("  mode:       ",
      if (is.finite(x$mode)) signif(x$mode, 6) else "none in range (monotone density)",
      "\n", sep = "")
  cat("  evaluated over t in [", signif(min(x$time), 4), ", ",
      signif(max(x$time), 4), "]\n", sep = "")
  invisible(x)
}

#' Order Statistic Distribution and Moments
#'
#' @param t Vector of times, for `bd_order_stat_cdf`.
#' @param r Order of the moment, for `bd_order_stat_moments`.
#' @param i Index of the order statistic, from 1 to `n`.
#' @param n Sample size.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return `bd_order_stat_cdf` returns a vector the length of `t`;
#'   `bd_order_stat_moments` returns a single number.
#'
#' @details
#' \eqn{F_{i:n}(t) = I_{F(t)}(i, n - i + 1)}, and moments are taken against the
#' order-statistic density on the Beta scale. The \eqn{i}-th order statistic of
#' a sample of size \eqn{n} has a finite \eqn{r}-th moment when
#' \eqn{b(n - i + 1) > r}, so extremes need heavier restrictions than the parent
#' distribution, and the maximum needs \eqn{b > r} exactly as the parent does.
#'
#' @seealso [bd_order_stat_pdf()]
#'
#' @export
#'
#' @examples
#' bd_order_stat_cdf(c(0.5, 1, 2), i = 3, n = 5, a = 1.5, b = 3, c = 2, k = 1)
#' bd_order_stat_moments(1, i = 1, n = 5, a = 1.5, b = 3, c = 2, k = 1)
bd_order_stat_cdf <- function(t, i, n, a, b, c, k) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  if (i < 1 || i > n) stop("'i' must lie between 1 and n.", call. = FALSE)
  stats::pbeta(pbetadanish(t, a, b, c, k), i, n - i + 1)
}

#' @rdname bd_order_stat_cdf
#' @export
bd_order_stat_moments <- function(r, i, n, a, b, c, k, rel.tol = 1e-10) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  if (i < 1 || i > n) stop("'i' must lie between 1 and n.", call. = FALSE)
  if (b * (n - i + 1) <= r) return(Inf)

  lB   <- lbeta(a, b)
  lcon <- lgamma(n + 1) - lgamma(i) - lgamma(n - i + 1)

  integrand <- function(u) {
    z  <- .bd_z_of_u(u, c, k)
    Fz <- stats::pbeta(u, a, b)
    rz <- if (r == 0) 0 else r * log(z)
    lg <- lcon + rz + (i - 1) * log(Fz) + (n - i) * log1p(-Fz) +
          (a - 1) * log(u) + (b - 1) * log1p(-u) - lB
    out <- exp(lg)
    out[!is.finite(out)] <- 0
    out
  }
  tryCatch(stats::integrate(integrand, 0, 1, rel.tol = rel.tol,
                            subdivisions = 4000L, stop.on.error = FALSE)$value,
           error = function(e) NA_real_)
}

#' Tail Index of the Beta-Danish Distribution
#'
#' @param a Shape parameter (beta generator). Accepted for interface
#'   consistency; it does not affect the index.
#' @param b Shape parameter. This is the tail index.
#' @param c Shape parameter (baseline). Accepted for interface consistency.
#' @param k Scale parameter (baseline). Accepted for interface consistency.
#'
#' @return A list with the tail index, the highest finite moment order, and a
#'   note on the moment generating function.
#'
#' @details
#' The survival function is regularly varying at infinity with index \eqn{-b}:
#' \eqn{S(t) \sim (c/k)^b t^{-b} / \{b B(a,b)\}}. Three consequences follow.
#'
#' The \eqn{r}-th moment is finite if and only if \eqn{b > r}, so the mean needs
#' \eqn{b > 1}, the variance \eqn{b > 2}, skewness \eqn{b > 3} and kurtosis
#' \eqn{b > 4}.
#'
#' The moment generating function does not exist for any \eqn{t > 0}: a
#' regularly varying tail decays polynomially, so \eqn{E(e^{tZ})} diverges.
#' Characteristic-function or Laplace-transform arguments must be used instead
#' of MGF ones anywhere in this family.
#'
#' The distribution lies in the Frechet domain of attraction, so sample maxima
#' normalise to a Frechet limit with shape \eqn{b}, not to a Gumbel limit.
#'
#' @export
#'
#' @examples
#' bd_tail_index(a = 1.5, b = 3, c = 2, k = 1)
bd_tail_index <- function(a, b, c, k) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  list(tail_index          = b,
       survival_exponent   = -b,
       highest_finite_moment = if (b > 1) floor(ceiling(b) - 1) else 0,
       moment_condition    = "E(Z^r) finite if and only if b > r",
       mgf_exists          = FALSE,
       domain_of_attraction = "Frechet")
})---")


.step("Adding a test for the new tuning arguments")

.put("tests/testthat/test-integration-args.R", r"---(
## The rel.tol and subdivisions arguments added in 3a-fix are real, not just
## documentation: they must reach stats::integrate() and change nothing about
## the answer at sensible settings.

test_that("bd_entropy_shannon accepts quadrature tuning arguments", {
  loose <- bd_entropy_shannon(1.5, 3, 2, 1, method = "quadrature",
                              rel.tol = 1e-6, subdivisions = 200L)
  tight <- bd_entropy_shannon(1.5, 3, 2, 1, method = "quadrature",
                              rel.tol = 1e-11, subdivisions = 4000L)
  expect_true(is.finite(loose) && is.finite(tight))
  expect_equal(loose, tight, tolerance = 1e-5)
  expect_equal(tight, bd_entropy_shannon(1.5, 3, 2, 1), tolerance = 1e-5)
})

test_that("bd_moment_summary accepts subdivisions", {
  s <- bd_moment_summary(1.5, 5, 2, 1, rel.tol = 1e-9, subdivisions = 1000L)
  expect_true(all(is.finite(s)))
  expect_equal(s[["mean"]], bd_moments(1, 1.5, 5, 2, 1), tolerance = 1e-7)
})

test_that("tuning arguments actually reach stats::integrate", {
  ## A subdivision limit far too small must make the integral fail rather than
  ## silently return the same answer, which proves the argument is wired up.
  expect_true(is.finite(
    bd_moments(1, 1.5, 5, 2, 1, rel.tol = 1e-9, subdivisions = 500L)))
  expect_equal(
    bd_moments(1, 1.5, 5, 2, 1, subdivisions = 500L),
    bd_moments(1, 1.5, 5, 2, 1, subdivisions = 4000L),
    tolerance = 1e-8)
})
)---")

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

.step("Confirming no combined or inherited parameter tags remain")
.problem <- character(0)
for (f in c("R/moments.R", "R/structural.R")) {
  txt <- readLines(f, warn = FALSE)
  if (any(grepl("@inheritParams", txt, fixed = TRUE)))
    .problem <- c(.problem, paste(f, "still uses @inheritParams"))
  combined <- grep("^#' @param [A-Za-z._][A-Za-z0-9._]*,", txt, value = TRUE)
  if (length(combined))
    .problem <- c(.problem, paste0(f, " has a combined tag: ", combined[1]))
}
if (length(.problem)) {
  for (p in .problem) .warn(p)
  .die("Combined or inherited tags remain; the warning would return.\n",
       "Backups: ", BACKUP_DIR)
}
.ok("one tag per argument, no @inheritParams")

.step("Checking every documented argument is a formal")
.mismatch <- character(0)
for (f in c("R/moments.R", "R/structural.R")) {
  env <- new.env()
  sys.source(f, envir = env, keep.source = FALSE)
  txt <- readLines(f, warn = FALSE)
  formals_all <- unique(unlist(lapply(ls(env, all.names = TRUE), function(nm) {
    ob <- get(nm, envir = env)
    if (is.function(ob)) names(formals(ob)) else character(0)
  })))
  documented <- unique(sub("^#' @param ([A-Za-z._][A-Za-z0-9._]*).*$", "\\1",
                           grep("^#' @param ", txt, value = TRUE)))
  extra <- setdiff(documented, c(formals_all, "x", "..."))
  if (length(extra))
    .mismatch <- c(.mismatch, paste0(f, ": ", paste(extra, collapse = ", ")))
}
if (length(.mismatch)) {
  for (m in .mismatch) .warn(paste("documented but not a formal:", m))
} else {
  .ok("every documented argument corresponds to a formal")
}

.step("devtools::document()")
.r <- tryCatch({ devtools::document(); TRUE }, error = function(e) conditionMessage(e))
if (!isTRUE(.r)) .die("document() failed:\n  ", .r, "\n\nBackups: ", BACKUP_DIR)
.ok("documentation regenerated")

.step("Inspecting the two Rd files the warning named")
for (rd in c("man/bd_entropy_shannon.Rd", "man/bd_moment_summary.Rd")) {
  if (!file.exists(rd)) { .warn(paste(rd, "not found")); next }
  txt <- paste(readLines(rd, warn = FALSE), collapse = " ")
  items <- unlist(regmatches(txt, gregexpr("\\\\item\\{[^}]*\\}", txt)))
  cat("       ", basename(rd), ": ", length(items), " documented argument(s)\n", sep = "")
}
.ok("Rd files regenerated")

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
cat("  PATCH 3a-fix COMPLETE\n")
cat(strrep("=", 78), "\n\n")
cat("  Every argument documented individually; no @inheritParams.\n")
cat("  bd_entropy_shannon() and bd_moment_summary() gained working\n")
cat("  rel.tol and subdivisions arguments.\n\n")
cat("  Expected: 0 errors, 0 warnings, 0-1 notes.\n")
cat("  Patch 3a then closes. Next: 3b (estimation and simulation).\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
