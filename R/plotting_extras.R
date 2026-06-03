#' Cox-Snell Residual Plot for AFT and Cure Fits
#'
#' Diagnostic Cox-Snell residual plot for a fitted AFT or cure model.
#'
#' @param x A fitted \code{"bd_aft"} or \code{"bd_cure"} object.
#' @param ... Further graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 200; x <- stats::rnorm(n)
#' k <- exp(-0.5 - 0.3 * x)
#' t_sim <- rbetadanish(n, a = 1, b = 2, c = 1.5, k = k)
#' dat <- data.frame(time = t_sim, status = 1, x = x)
#' fit <- fit_bd_aft(survival::Surv(time, status) ~ x, data = dat,
#'                   n_starts = 2)
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

.bd_aft_or_cure_coxsnell <- function(x, ...) {
  d <- x$data
  time   <- d$time
  status <- d$status
  X      <- d$X
  b <- as.numeric(x$coefficients["b"])
  c <- as.numeric(x$coefficients["c"])
  delta_names <- grep("^delta_", names(x$coefficients), value = TRUE)
  delta <- as.numeric(x$coefficients[delta_names])
  if (length(delta) != ncol(X))
    stop("AFT/cure coefficient layout does not match the design matrix.")
  k_i <- exp(as.numeric(X %*% delta))
  r <- -pbetadanish(time, a = 1, b = b, c = c, k = k_i,
                    lower.tail = FALSE, log.p = TRUE)
  r[!is.finite(r)] <- NA
  km_r <- survival::survfit(survival::Surv(r, status) ~ 1)
  H_r  <- -log(pmax(km_r$surv, 1e-12))
  graphics::plot(km_r$time, H_r, type = "s",
                 xlab = "Cox-Snell residual",
                 ylab = "Estimated cumulative hazard",
                 main = "Cox-Snell residuals", lwd = 2, ...)
  graphics::abline(0, 1, col = "red", lwd = 2, lty = 2)
  graphics::legend("topleft",
                   legend = c("KM cum-hazard of residuals", "y = x"),
                   col = c("black", "red"), lwd = 2,
                   lty = c(1, 2), bty = "n")
  invisible(x)
}
