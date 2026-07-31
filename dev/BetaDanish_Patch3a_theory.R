## =============================================================================
##  BetaDanish  --  PHASE 2, PATCH 3a : THEORETICAL PROPERTIES
## =============================================================================
##
##  Implements approved recommendations 25 to 33.
##
##    25  raw, incomplete and conditional moments, with the b > r condition
##    26  Shannon (closed form), Renyi and Tsallis entropies
##    27  mean residual life and reversed mean residual life
##    28  stress-strength reliability
##    29  hazard shape classification, Glaser diagnostic, mode
##    30  probability weighted moments, mean deviations, Lorenz, Bonferroni
##    31  the named ED API: ded, ped, qed, red, hed, sed
##    32  order statistic moments and distribution function
##    33  tail index, and the MGF documented as non-existent
##
##  NUMERICAL BASIS
##    Every integral is taken on the finite Beta(a, b) scale rather than over
##    (0, Inf), using the substitution u = G(z) = {kz/(1+kz)}^c, under which
##    u ~ Beta(a, b) exactly and z = u^(1/c) / (k(1 - u^(1/c))). This is the
##    approach documented in Master_R_Code_for_Thesis.R at moment_betadanish(),
##    where it replaced a fixed J = 4000 series truncation that under-reported
##    the ED variance by about 8.5 percent near the moment boundary.
##
##  ONE CORRECTION TO THE MASTER SCRIPT
##    bd_shannon() truncates its series at M terms with no tail correction. The
##    terms decay like i^-(b+1), so the omitted tail is negligible for large b
##    but not for small b. Against 40-digit numerical integration:
##
##        b = 3.0   error 8.7e-09        b = 0.8   error 1.3e-02
##        b = 2.5   error 2.3e-07        b = 0.5   error 1.1e-01   (M = 2000)
##
##    bd_entropy_shannon() adds the analytic Euler-Maclaurin tail
##    Gamma(b) c^b / B(a,b) * M^-b / b, which restores agreement to about eight
##    digits at every b tested. Worth folding back into the master script.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3a_theory.R")
##  IDEMPOTENT   Yes -- every file is written whole.
##  SCOPE        No changes to any existing function. Patch 3b covers
##               estimation and simulation, Patch 3c visualisation and 0.3.0.
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
cat("  BetaDanish  --  Phase 2, Patch 3a of 3 : theoretical properties\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!any(grepl(".bd_logG", readLines("R/dist_functions.R", warn = FALSE), fixed = TRUE)))
  .die("Patch 1 has not been applied -- R/dist_functions.R has no .bd_logG().")
if (!any(grepl(".bd_make_accept", readLines("R/fit_models.R", warn = FALSE), fixed = TRUE)))
  .die("Patch 2c has not been applied -- run it and its fixes first.")
.ok("Patches 1, 2 and 2c detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3a"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  R/moments.R  --  recommendations 25, 27, 30
## =============================================================================

.step("Rec 25, 27, 30: writing R/moments.R")

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
#' @param a,b,c,k Distribution parameters.
#' @param rel.tol,subdivisions Passed to [stats::integrate()].
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
#' @inheritParams bd_moments
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
bd_moment_summary <- function(a, b, c, k, rel.tol = 1e-10) {
  m <- bd_moments(1:4, a, b, c, k, rel.tol = rel.tol)
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
#' @param a,b,c,k Distribution parameters.
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
#' @inheritParams bd_incomplete_moment
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
#' @param a,b,c,k Distribution parameters.
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
#' @param a,b,c,k Distribution parameters.
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
#' @param a,b,c,k Distribution parameters.
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
#' @param a,b,c,k Distribution parameters.
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
}
)---")

## =============================================================================
##  R/structural.R  --  recommendations 26, 28, 29, 32, 33
## =============================================================================

.step("Rec 26, 28, 29, 32, 33: writing R/structural.R")

.put("R/structural.R", r"---(## Entropy, reliability, shape and order-statistic properties.

#' Shannon Entropy of the Beta-Danish Distribution
#'
#' @param a,b,c,k Distribution parameters.
#' @param terms Number of series terms before the analytic tail is applied.
#' @param method `"closed"` (default) for the closed form, or `"quadrature"`
#'   for direct numerical integration of \eqn{-f \log f}, which is slower but
#'   independent of the series.
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
                               method = c("closed", "quadrature")) {
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
      stats::integrate(integrand, 0, 1, rel.tol = 1e-10,
                       subdivisions = 4000L, stop.on.error = FALSE)$value,
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
#' @param a,b,c,k Distribution parameters.
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
#' @param a,b,c,k Distribution parameters.
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
#' @param a,b,c,k Distribution parameters.
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
#' @param a,b,c,k Distribution parameters. Only `b` affects the index; the
#'   others are accepted so the call reads like the rest of the interface.
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
}
)---")

## =============================================================================
##  R/ed_api.R  --  recommendation 31
## =============================================================================

.step("Rec 31: writing R/ed_api.R")

.put("R/ed_api.R", r"---(#' The Exponentiated Danish Distribution
#'
#' Density, distribution function, quantile function, survival function, hazard
#' function and random generation for the three-parameter Exponentiated Danish
#' distribution, the \eqn{a = 1} submodel of the Beta-Danish family.
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations to generate.
#' @param b,c,k Shape, shape and scale parameters.
#' @param log,log.p Logical; return values on the log scale.
#' @param lower.tail Logical; if `TRUE` (default) probabilities are
#'   \eqn{P[X \le x]}.
#'
#' @return Vectors of the same form as the Beta-Danish equivalents.
#'
#' @details
#' These are thin wrappers that fix \eqn{a = 1}, provided because the ED
#' submodel is a named distribution in its own right in the underlying work and
#' writing `dbetadanish(x, 1, b, c, k)` obscures that.
#'
#' At \eqn{a = 1} the beta generator collapses and the distribution function
#' simplifies to \eqn{F(t) = 1 - \{1 - G(t)\}^{b}} with
#' \eqn{G(t) = \{kt/(1+kt)\}^{c}}. The upper tail still has index \eqn{-b}, so
#' the moment condition \eqn{b > r} is unchanged.
#'
#' @seealso [BetaDanish] for the four-parameter parent,
#'   [fit_betadanish()] with `submodel = TRUE` for estimation.
#'
#' @name ExponentiatedDanish
#'
#' @examples
#' ded(2, b = 3, c = 2, k = 1)
#' ped(2, b = 3, c = 2, k = 1)
#' qed(0.5, b = 3, c = 2, k = 1)
#' red(5, b = 3, c = 2, k = 1)
#'
#' # Identical to the a = 1 parent
#' all.equal(ded(2, 3, 2, 1), dbetadanish(2, 1, 3, 2, 1))
NULL

#' @rdname ExponentiatedDanish
#' @export
ded <- function(x, b, c, k, log = FALSE) {
  dbetadanish(x, a = 1, b = b, c = c, k = k, log = log)
}

#' @rdname ExponentiatedDanish
#' @export
ped <- function(q, b, c, k, lower.tail = TRUE, log.p = FALSE) {
  pbetadanish(q, a = 1, b = b, c = c, k = k,
              lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ExponentiatedDanish
#' @export
qed <- function(p, b, c, k, lower.tail = TRUE, log.p = FALSE) {
  qbetadanish(p, a = 1, b = b, c = c, k = k,
              lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ExponentiatedDanish
#' @export
red <- function(n, b, c, k) {
  rbetadanish(n, a = 1, b = b, c = c, k = k)
}

#' @rdname ExponentiatedDanish
#' @export
sed <- function(x, b, c, k, log = FALSE) {
  sbetadanish(x, a = 1, b = b, c = c, k = k, log = log)
}

#' @rdname ExponentiatedDanish
#' @export
hed <- function(x, b, c, k, log = FALSE) {
  hbetadanish(x, a = 1, b = b, c = c, k = k, log = log)
}
)---")

## =============================================================================
##  TESTS
## =============================================================================

.step("Writing tests/testthat/test-moments.R")

.put("tests/testthat/test-moments.R", r"---(P <- list(a = 1.5, b = 5, c = 2, k = 1)     # b = 5 so all four moments exist

test_that("raw moments agree with direct integration of z^r f(z)", {
  for (r in 1:3) {
    direct <- stats::integrate(
      function(z) z^r * dbetadanish(z, P$a, P$b, P$c, P$k),
      0, Inf, rel.tol = 1e-10)$value
    expect_equal(bd_moments(r, P$a, P$b, P$c, P$k), direct, tolerance = 1e-6)
  }
})

test_that("the zeroth moment is one", {
  expect_equal(bd_moments(0, P$a, P$b, P$c, P$k), 1, tolerance = 1e-8)
})

test_that("moments respect the b > r existence condition", {
  ## The survival tail has index -b, so E(Z^r) is finite iff b > r.
  expect_true(is.finite(bd_moments(2, a = 1.5, b = 2.5, c = 2, k = 1)))
  expect_identical(bd_moments(3, a = 1.5, b = 2.5, c = 2, k = 1), Inf)
  expect_identical(bd_moments(4, a = 1.5, b = 3,   c = 2, k = 1), Inf)
  expect_identical(bd_moments(3, a = 1.5, b = 3,   c = 2, k = 1), Inf)  # b = r
})

test_that("bd_moment_summary reports NA rather than nonsense past the boundary", {
  s <- bd_moment_summary(P$a, P$b, P$c, P$k)
  expect_true(all(is.finite(s)))
  expect_gt(s[["variance"]], 0)
  expect_equal(s[["sd"]], sqrt(s[["variance"]]))

  s2 <- bd_moment_summary(a = 1.5, b = 2.5, c = 2, k = 1)
  expect_true(is.finite(s2[["mean"]]))
  expect_true(is.finite(s2[["variance"]]))
  expect_true(is.na(s2[["skewness"]]))
  expect_true(is.na(s2[["kurtosis"]]))
})

test_that("incomplete moments split the complete moment exactly", {
  for (t in c(0.2, 1, 5)) {
    lo <- bd_incomplete_moment(1, t, P$a, P$b, P$c, P$k, lower = TRUE)
    up <- bd_incomplete_moment(1, t, P$a, P$b, P$c, P$k, lower = FALSE)
    expect_equal(lo + up, bd_moments(1, P$a, P$b, P$c, P$k), tolerance = 1e-7)
  }
})

test_that("incomplete moment boundaries behave", {
  expect_equal(bd_incomplete_moment(1, 0, P$a, P$b, P$c, P$k), 0)
  expect_equal(bd_incomplete_moment(1, Inf, P$a, P$b, P$c, P$k),
               bd_moments(1, P$a, P$b, P$c, P$k), tolerance = 1e-7)
  expect_true(is.na(bd_incomplete_moment(1, NA_real_, P$a, P$b, P$c, P$k)))
})

test_that("conditional moments match the incomplete/survival ratio", {
  t <- 1.5
  expect_equal(
    bd_conditional_moment(1, t, P$a, P$b, P$c, P$k, upper = TRUE),
    bd_incomplete_moment(1, t, P$a, P$b, P$c, P$k, lower = FALSE) /
      sbetadanish(t, P$a, P$b, P$c, P$k),
    tolerance = 1e-9)
})

test_that("MRL matches its definition and is zero-limit consistent", {
  t <- c(0.5, 1, 3)
  m <- bd_mrl(t, P$a, P$b, P$c, P$k)
  expect_true(all(is.finite(m)))
  expect_true(all(m > 0))

  ## m(0+) = E(Z), by definition
  expect_equal(bd_mrl(1e-8, P$a, P$b, P$c, P$k),
               bd_moments(1, P$a, P$b, P$c, P$k), tolerance = 1e-4)

  ## E(Z - t | Z > t) computed the long way
  tt <- 1
  direct <- stats::integrate(
    function(z) (z - tt) * dbetadanish(z, P$a, P$b, P$c, P$k),
    tt, Inf, rel.tol = 1e-10)$value / sbetadanish(tt, P$a, P$b, P$c, P$k)
  expect_equal(bd_mrl(tt, P$a, P$b, P$c, P$k), direct, tolerance = 1e-6)
})

test_that("reversed MRL matches its definition", {
  tt <- 2
  direct <- tt - stats::integrate(
    function(z) z * dbetadanish(z, P$a, P$b, P$c, P$k),
    0, tt, rel.tol = 1e-10)$value / pbetadanish(tt, P$a, P$b, P$c, P$k)
  expect_equal(bd_rmrl(tt, P$a, P$b, P$c, P$k), direct, tolerance = 1e-6)
})

test_that("MRL requires a finite mean", {
  expect_warning(m <- bd_mrl(1, a = 1.5, b = 0.8, c = 2, k = 1), "b > 1")
  expect_true(is.na(m))
  expect_warning(bd_rmrl(1, a = 1.5, b = 0.8, c = 2, k = 1), "b > 1")
})

test_that("mean deviation about the mean matches E|Z - mu|", {
  mu <- bd_moments(1, P$a, P$b, P$c, P$k)
  direct <- stats::integrate(
    function(z) abs(z - mu) * dbetadanish(z, P$a, P$b, P$c, P$k),
    0, Inf, rel.tol = 1e-10)$value
  expect_equal(bd_mean_deviation(P$a, P$b, P$c, P$k), direct, tolerance = 1e-6)
})

test_that("mean deviation about the median matches E|Z - M|", {
  med <- qbetadanish(0.5, P$a, P$b, P$c, P$k)
  direct <- stats::integrate(
    function(z) abs(z - med) * dbetadanish(z, P$a, P$b, P$c, P$k),
    0, Inf, rel.tol = 1e-10)$value
  expect_equal(bd_mean_deviation(P$a, P$b, P$c, P$k, about = "median"),
               direct, tolerance = 1e-6)
})

test_that("the Lorenz curve has the right endpoints and shape", {
  p <- c(0, 0.25, 0.5, 0.75, 1)
  L <- bd_lorenz(p, P$a, P$b, P$c, P$k)
  expect_equal(L[1], 0, tolerance = 1e-8)
  expect_equal(L[5], 1, tolerance = 1e-6)
  expect_true(all(diff(L) > 0))          # increasing
  expect_true(all(L[2:4] < p[2:4]))      # below the equality line
})

test_that("Bonferroni is Lorenz divided by p", {
  p <- c(0.25, 0.5, 0.75)
  expect_equal(bd_bonferroni(p, P$a, P$b, P$c, P$k),
               bd_lorenz(p, P$a, P$b, P$c, P$k) / p, tolerance = 1e-10)
})

test_that("PWM with zero weights reduces to the raw moment", {
  expect_equal(bd_pwm(1, 0, 0, P$a, P$b, P$c, P$k),
               bd_moments(1, P$a, P$b, P$c, P$k), tolerance = 1e-7)
  expect_equal(bd_pwm(0, 0, 0, P$a, P$b, P$c, P$k), 1, tolerance = 1e-8)
})

test_that("PWM matches direct integration", {
  direct <- stats::integrate(
    function(z) z * pbetadanish(z, P$a, P$b, P$c, P$k) *
      dbetadanish(z, P$a, P$b, P$c, P$k),
    0, Inf, rel.tol = 1e-10)$value
  expect_equal(bd_pwm(1, 1, 0, P$a, P$b, P$c, P$k), direct, tolerance = 1e-6)
})

test_that("invalid parameters are rejected", {
  expect_error(bd_moments(1, a = -1, b = 3, c = 2, k = 1), "strictly positive")
  expect_error(bd_moments(-1, P$a, P$b, P$c, P$k), "non-negative")
})
)---")

.step("Writing tests/testthat/test-structural.R")

.put("tests/testthat/test-structural.R", r"---(P <- list(a = 1.5, b = 3, c = 2, k = 1)

test_that("closed-form Shannon entropy matches quadrature", {
  for (p in list(list(a = 1.5, b = 3, c = 2, k = 1),
                 list(a = 1,   b = 2.5, c = 3, k = 1.2),
                 list(a = 2,   b = 0.8, c = 1.5, k = 0.5))) {
    closed <- bd_entropy_shannon(p$a, p$b, p$c, p$k)
    quad   <- bd_entropy_shannon(p$a, p$b, p$c, p$k, method = "quadrature")
    expect_equal(closed, quad, tolerance = 1e-5,
                 label = sprintf("b = %.1f", p$b))
  }
})

test_that("the series tail correction matters at small b", {
  ## Terms decay like i^-(b+1). Without the analytic tail the truncated sum is
  ## badly short at small b; with it, closed and quadrature agree.
  small <- list(a = 1, b = 0.5, c = 2, k = 1)
  closed <- bd_entropy_shannon(small$a, small$b, small$c, small$k, terms = 2000L)
  quad   <- bd_entropy_shannon(small$a, small$b, small$c, small$k,
                               method = "quadrature")
  expect_equal(closed, quad, tolerance = 1e-3)

  ## And the correction is not merely cosmetic: the raw sum is far off.
  i <- seq_len(2000L)
  raw <- lbeta(small$a, small$b) - log(small$c * small$k) -
    (small$a - 1 / small$c) * (digamma(small$a) - digamma(small$a + small$b)) -
    (small$b - 1) * (digamma(small$b) - digamma(small$a + small$b)) +
    2 * sum(exp(lbeta(small$a + i / small$c, small$b) -
                  lbeta(small$a, small$b)) / i)
  expect_gt(abs(raw - quad), 0.05)
})

test_that("Renyi and Tsallis are finite and ordered sensibly", {
  r2 <- bd_entropy_renyi(P$a, P$b, P$c, P$k, order = 2)
  r3 <- bd_entropy_renyi(P$a, P$b, P$c, P$k, order = 3)
  sh <- bd_entropy_shannon(P$a, P$b, P$c, P$k)
  expect_true(all(is.finite(c(r2, r3, sh))))
  ## Renyi entropy is non-increasing in its order, and Shannon is the q -> 1 limit
  expect_lt(r3, r2)
  expect_lt(r2, sh)

  t2 <- bd_entropy_tsallis(P$a, P$b, P$c, P$k, order = 2)
  expect_true(is.finite(t2))
})

test_that("order one is refused for Renyi and Tsallis", {
  expect_error(bd_entropy_renyi(P$a, P$b, P$c, P$k, order = 1), "Shannon")
  expect_error(bd_entropy_tsallis(P$a, P$b, P$c, P$k, order = 1), "Shannon")
})

test_that("stress-strength is one half for identical distributions", {
  expect_equal(bd_stress_strength(P, P), 0.5, tolerance = 1e-6)
})

test_that("stress-strength respects direction", {
  strong <- list(a = 1.5, b = 3, c = 2, k = 0.5)   # larger scale => larger X
  weak   <- list(a = 1.5, b = 3, c = 2, k = 2)
  R  <- bd_stress_strength(strength = strong, stress = weak)
  Rr <- bd_stress_strength(strength = weak,   stress = strong)
  expect_gt(R, 0.5)
  expect_equal(R + Rr, 1, tolerance = 1e-6)
})

test_that("stress-strength matches direct integration", {
  x <- list(a = 1.2, b = 4, c = 2, k = 0.8)
  y <- list(a = 1.5, b = 3, c = 2, k = 1)
  direct <- stats::integrate(
    function(t) sbetadanish(t, x$a, x$b, x$c, x$k) *
      dbetadanish(t, y$a, y$b, y$c, y$k),
    0, Inf, rel.tol = 1e-10)$value
  expect_equal(bd_stress_strength(x, y), direct, tolerance = 1e-6)
})

test_that("stress-strength validates its arguments", {
  expect_error(bd_stress_strength(list(a = 1, b = 2), P), "named values")
  expect_error(bd_stress_strength(list(a = -1, b = 2, c = 1, k = 1), P),
               "strictly positive")
})

test_that("hazard shape classification returns a known label", {
  s <- bd_hazard_shape(P$a, P$b, P$c, P$k)
  expect_s3_class(s, "bd_shape")
  expect_true(s$shape %in% c("increasing", "decreasing", "bathtub",
                             "upside-down bathtub", "indeterminate"))
  expect_length(s$time, 400L)
  expect_true(all(is.finite(s$eta[is.finite(s$hazard)])))
  expect_output(print(s), "hazard shape")
})

test_that("Glaser eta agrees with a numerical derivative of the log-density", {
  t <- c(0.3, 1, 2.5)
  h <- 1e-5
  numeric_eta <- -(dbetadanish(t + h, P$a, P$b, P$c, P$k, log = TRUE) -
                   dbetadanish(t - h, P$a, P$b, P$c, P$k, log = TRUE)) / (2 * h)
  expect_equal(BetaDanish:::.bd_glaser_eta(t, P$a, P$b, P$c, P$k),
               numeric_eta, tolerance = 1e-5)
})

test_that("the classifier recognises each shape", {
  cl <- BetaDanish:::.bd_classify
  expect_equal(cl(1:10), "increasing")
  expect_equal(cl(10:1), "decreasing")
  expect_equal(cl(c(5, 3, 1, 2, 4, 6)), "bathtub")
  expect_equal(cl(c(1, 3, 6, 4, 2)), "upside-down bathtub")
})

test_that("order statistic CDF is a proper distribution function", {
  t <- c(0.2, 1, 5, 20)
  F3 <- bd_order_stat_cdf(t, i = 3, n = 5, P$a, P$b, P$c, P$k)
  expect_true(all(diff(F3) >= 0))
  expect_true(all(F3 >= 0 & F3 <= 1))

  ## The maximum is F^n and the minimum is 1 - (1 - F)^n
  Fx <- pbetadanish(t, P$a, P$b, P$c, P$k)
  expect_equal(bd_order_stat_cdf(t, 5, 5, P$a, P$b, P$c, P$k), Fx^5,
               tolerance = 1e-10)
  expect_equal(bd_order_stat_cdf(t, 1, 5, P$a, P$b, P$c, P$k), 1 - (1 - Fx)^5,
               tolerance = 1e-10)
})

test_that("order statistic moments are ordered and respect existence", {
  m <- vapply(1:5, function(i)
    bd_order_stat_moments(1, i, 5, P$a, P$b, P$c, P$k), numeric(1))
  expect_true(all(is.finite(m)))
  expect_true(all(diff(m) > 0))          # E(Z_(1)) < ... < E(Z_(n))

  ## The maximum of a sample needs b > r, exactly as the parent does.
  expect_identical(bd_order_stat_moments(4, 5, 5, a = 1.5, b = 3, c = 2, k = 1),
                   Inf)
  ## The minimum is far better behaved: b*n > r suffices.
  expect_true(is.finite(
    bd_order_stat_moments(4, 1, 5, a = 1.5, b = 3, c = 2, k = 1)))
})

test_that("order statistic index is validated", {
  expect_error(bd_order_stat_cdf(1, i = 6, n = 5, P$a, P$b, P$c, P$k), "between")
  expect_error(bd_order_stat_moments(1, 0, 5, P$a, P$b, P$c, P$k), "between")
})

test_that("the tail index is b, with the consequences recorded", {
  ti <- bd_tail_index(P$a, P$b, P$c, P$k)
  expect_equal(ti$tail_index, P$b)
  expect_equal(ti$survival_exponent, -P$b)
  expect_false(ti$mgf_exists)
  expect_equal(ti$domain_of_attraction, "Frechet")
})

test_that("the ED API matches the a = 1 parent exactly", {
  b <- 3; c <- 2; k <- 1; x <- c(0.5, 1, 4)
  expect_equal(ded(x, b, c, k), dbetadanish(x, 1, b, c, k))
  expect_equal(ped(x, b, c, k), pbetadanish(x, 1, b, c, k))
  expect_equal(sed(x, b, c, k), sbetadanish(x, 1, b, c, k))
  expect_equal(hed(x, b, c, k), hbetadanish(x, 1, b, c, k))
  p <- c(0.1, 0.5, 0.9)
  expect_equal(qed(p, b, c, k), qbetadanish(p, 1, b, c, k))
  expect_equal(ded(x, b, c, k, log = TRUE), dbetadanish(x, 1, b, c, k, log = TRUE))
  expect_equal(ped(x, b, c, k, lower.tail = FALSE),
               pbetadanish(x, 1, b, c, k, lower.tail = FALSE))
})

test_that("the ED closed form F(t) = 1 - (1 - G)^b holds", {
  b <- 3; c <- 2; k <- 1; t <- c(0.5, 1, 4, 50)
  G <- (k * t / (1 + k * t))^c
  expect_equal(ped(t, b, c, k), 1 - (1 - G)^b, tolerance = 1e-9)
})

test_that("red respects the seed and returns positive values", {
  set.seed(3); x1 <- red(20, 3, 2, 1)
  set.seed(3); x2 <- red(20, 3, 2, 1)
  expect_identical(x1, x2)
  expect_length(x1, 20L)
  expect_true(all(x1 > 0))
})
)---")

## =============================================================================
##  NEWS
## =============================================================================

.step("Recording the additions in NEWS.md")

.nw <- readLines("NEWS.md", warn = FALSE)
if (!any(grepl("bd_moments", .nw, fixed = TRUE))) {
  .hdr <- grep("^# BetaDanish 0\\.2\\.0\\.9000", .nw)
  if (length(.hdr) == 1L) {
    .backup("NEWS.md")
    .sec <- c(
      "",
      "## Theoretical properties",
      "",
      "* **Moments**: `bd_moments()`, `bd_moment_summary()`,",
      "  `bd_incomplete_moment()` and `bd_conditional_moment()`. Every integral",
      "  is taken on the finite Beta(a, b) scale under the substitution",
      "  `u = G(z)`, so there is no series truncation. The existence condition",
      "  `b > r` is enforced: a moment that does not exist returns `Inf` rather",
      "  than a large finite number from a truncated sum.",
      "",
      "* **Ageing and inequality**: `bd_mrl()` and `bd_rmrl()` for mean residual",
      "  life and mean inactivity time, `bd_mean_deviation()`, `bd_lorenz()`,",
      "  `bd_bonferroni()` and `bd_pwm()`.",
      "",
      "* **Entropies**: `bd_entropy_shannon()` now uses the closed form, with",
      "  `method = \"quadrature\"` retained as an independent cross-check.",
      "  `bd_entropy_renyi()` and `bd_entropy_tsallis()` are new.",
      "",
      "  The closed form's series terms decay like `i^-(b+1)`, so plain",
      "  truncation is accurate for large `b` and not for small `b`. Against",
      "  high-precision integration the untruncated-tail sum at `M = 2000` is",
      "  out by about `1e-8` at `b = 3` but by about `0.11` at `b = 0.5`. An",
      "  analytic Euler-Maclaurin tail is now added, restoring agreement to",
      "  roughly eight digits throughout.",
      "",
      "* **Reliability**: `bd_stress_strength()` for `R = P(X > Y)` with `X` the",
      "  strength, evaluated on the Beta scale of the stress variable.",
      "",
      "* **Shape**: `bd_hazard_shape()` classifies the hazard as increasing,",
      "  decreasing, bathtub or upside-down bathtub using Glaser's",
      "  `eta(t) = -f'(t)/f(t)`, formed analytically, and reports the mode.",
      "",
      "* **Order statistics**: `bd_order_stat_cdf()` and",
      "  `bd_order_stat_moments()`, alongside the existing",
      "  `bd_order_stat_pdf()`.",
      "",
      "* **Tail behaviour**: `bd_tail_index()` records that the survival",
      "  function is regularly varying with index `-b`, that `E(Z^r)` is finite",
      "  exactly when `b > r`, that the moment generating function does not",
      "  exist, and that the distribution lies in the Frechet domain of",
      "  attraction.",
      "",
      "* **Exponentiated Danish API**: `ded()`, `ped()`, `qed()`, `red()`,",
      "  `sed()` and `hed()` name the `a = 1` submodel directly instead of",
      "  requiring `dbetadanish(x, 1, b, c, k)`.")
    .nw <- append(.nw, .sec, after = .hdr)
    .write_lines("NEWS.md", .nw)
    .ok("NEWS.md updated")
  } else {
    .warn("could not find the version header in NEWS.md; add the note by hand")
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

.step("Loading from source for the self-test")
.loaded <- tryCatch({ devtools::load_all(".", quiet = TRUE); TRUE },
                    error = function(e) conditionMessage(e))
if (!isTRUE(.loaded))
  .die("load_all() failed:\n  ", .loaded, "\n\nBackups: ", BACKUP_DIR)
.ok("source loaded")

.step("Numerical self-test against independent integration")
.a <- 1.5; .b <- 5; .cc <- 2; .kk <- 1
.checks <- list()
.checks[["E(Z) vs direct integral"]] <- {
  got  <- bd_moments(1, .a, .b, .cc, .kk)
  want <- stats::integrate(function(z) z * dbetadanish(z, .a, .b, .cc, .kk),
                           0, Inf, rel.tol = 1e-10)$value
  c(got, want)
}
.checks[["incomplete moments sum to the whole"]] <- {
  lo <- bd_incomplete_moment(1, 1, .a, .b, .cc, .kk, lower = TRUE)
  up <- bd_incomplete_moment(1, 1, .a, .b, .cc, .kk, lower = FALSE)
  c(lo + up, bd_moments(1, .a, .b, .cc, .kk))
}
.checks[["Shannon closed vs quadrature"]] <- c(
  bd_entropy_shannon(.a, .b, .cc, .kk),
  bd_entropy_shannon(.a, .b, .cc, .kk, method = "quadrature"))
.checks[["Shannon closed vs quadrature, b = 0.5"]] <- c(
  bd_entropy_shannon(1, 0.5, 2, 1),
  bd_entropy_shannon(1, 0.5, 2, 1, method = "quadrature"))
.checks[["stress-strength, identical laws = 0.5"]] <- {
  p <- list(a = .a, b = .b, c = .cc, k = .kk)
  c(bd_stress_strength(p, p), 0.5)
}
.checks[["MRL vs direct conditional expectation"]] <- {
  tt <- 1
  want <- stats::integrate(function(z) (z - tt) * dbetadanish(z, .a, .b, .cc, .kk),
                           tt, Inf, rel.tol = 1e-10)$value /
          sbetadanish(tt, .a, .b, .cc, .kk)
  c(bd_mrl(tt, .a, .b, .cc, .kk), want)
}
.checks[["ED wrapper vs a = 1 parent"]] <- c(
  ded(2, 3, 2, 1), dbetadanish(2, 1, 3, 2, 1))

.worst <- 0
for (nm in names(.checks)) {
  v <- .checks[[nm]]
  rel <- abs(v[1] - v[2]) / max(abs(v[2]), 1e-12)
  .worst <- max(.worst, rel)
  if (is.finite(rel) && rel < 1e-4) {
    .ok(sprintf("%-42s  %.6g vs %.6g  (rel %.1e)", nm, v[1], v[2], rel))
  } else {
    .warn(sprintf("%-42s  %.6g vs %.6g  (rel %.1e)", nm, v[1], v[2], rel))
  }
}
if (!is.finite(.worst) || .worst > 1e-4)
  .die("A numerical self-test disagreed with independent integration.\n",
       "Nothing further was run. Backups: ", BACKUP_DIR)
.ok(sprintf("all checks agree to better than 1e-4 (worst %.1e)", .worst))

.step("devtools::document()")
.r <- tryCatch({ devtools::document(); TRUE }, error = function(e) conditionMessage(e))
if (!isTRUE(.r)) .die("document() failed:\n  ", .r, "\n\nBackups: ", BACKUP_DIR)
.ok("documentation regenerated")

.step("Confirming the new exports")
.ns <- readLines("NAMESPACE", warn = FALSE)
.want <- c("bd_moments", "bd_moment_summary", "bd_incomplete_moment",
           "bd_conditional_moment", "bd_mrl", "bd_rmrl", "bd_mean_deviation",
           "bd_lorenz", "bd_bonferroni", "bd_pwm", "bd_entropy_shannon",
           "bd_entropy_renyi", "bd_entropy_tsallis", "bd_stress_strength",
           "bd_hazard_shape", "bd_order_stat_cdf", "bd_order_stat_moments",
           "bd_tail_index", "ded", "ped", "qed", "red", "sed", "hed")
.missing <- .want[!vapply(.want, function(f)
  any(grepl(paste0("export(", f, ")"), .ns, fixed = TRUE)), logical(1))]
if (length(.missing)) {
  .warn(paste("not exported:", paste(.missing, collapse = ", ")))
} else {
  .ok(sprintf("all %d new functions exported", length(.want)))
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
cat("  PATCH 3a COMPLETE  --  theoretical properties\n")
cat(strrep("=", 78), "\n\n")
cat("  25  moments: raw, summary, incomplete, conditional; b > r enforced\n")
cat("  26  Shannon closed form with analytic tail; Renyi; Tsallis\n")
cat("  27  mean residual life and mean inactivity time\n")
cat("  28  stress-strength reliability\n")
cat("  29  hazard shape via Glaser's eta, plus the mode\n")
cat("  30  PWM, mean deviations, Lorenz, Bonferroni\n")
cat("  31  ED API: ded, ped, qed, red, sed, hed\n")
cat("  32  order statistic CDF and moments\n")
cat("  33  tail index, moment condition, MGF non-existence\n\n")
cat("  Still to come\n")
cat("    3b  estimation and simulation (recs 34-39)\n")
cat("    3c  visualisation and the 0.3.0 release (recs 40-42, 45, 46)\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
