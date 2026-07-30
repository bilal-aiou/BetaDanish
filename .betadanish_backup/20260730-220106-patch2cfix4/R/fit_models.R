#' Fit the Beta-Danish Distribution to Survival Data
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
}
