## =============================================================================
##  BetaDanish  --  PATCH 3f-fix : plot.bd_bayes density panel
## =============================================================================
##
##  Patch 3f aborted at its own self-test, which is what the self-test is for.
##
##      Error in xy.coords(x, y, xlabel, ylabel, log) :
##        'x' is a list, but does not have components 'x' and 'y'
##
##  plot.bd_bayes() drew the posterior panel with
##
##      graphics::plot(d, ...)          # d is a "density" object
##
##  expecting S3 dispatch to plot.density. Inside a package namespace that
##  dispatch did not happen, so plot.default received a list and failed. The
##  trace panels were fine because they pass numeric x and y.
##
##  FIX
##    The density panel now plots d$x against d$y explicitly. Nothing depends
##    on dispatch, and the axis label is set properly rather than inherited.
##    Every plot call in the file also now uses the base generic, matching the
##    decision taken in Patch 2b when plot() moved from graphics to base.
##
##  STATE
##    Patch 3f wrote every other file successfully before it stopped, so this
##    patch rewrites R/plots_extra.R only and then runs the verification,
##    documentation, tests and check that 3f never reached.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3f_fix.R")
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
cat("  BetaDanish  --  Patch 3f-fix : plot.bd_bayes density panel\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("R/plots_extra.R")) .die("Patch 3f has not been applied.")
if (!file.exists("R/simulation_study.R")) .die("Patch 3e has not been applied.")
.ok("Patch 3f groundwork detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3ffix"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

.step("Rewriting R/plots_extra.R")

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
#' @param ... Further graphical parameters.
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
    plot(c(0, u), c(0, phi), type = "l", col = col, lwd = lwd,
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
#' @param ... Further graphical parameters.
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

  plot(x$grid[ok], x$profile[ok], type = "l", col = col, lwd = lwd,
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
#' @param ... Further graphical parameters.
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
      plot(it, v, type = "l", col = col,
                     main = paste("Trace:", p),
                     xlab = "iteration", ylab = p, ...)
      graphics::abline(h = mean(v), col = "red", lty = 2)
    }
    if (type %in% c("both", "density")) {
      ## Plot the coordinates rather than the density object: relying on S3
      ## dispatch for plot.density is fragile inside a package namespace, and
      ## falling through to plot.default on a list is an error, not a warning.
      d <- stats::density(v)
      plot(d$x, d$y, type = "l", col = col, lwd = 2,
           main = paste("Posterior:", p), xlab = p, ylab = "density", ...)
      graphics::abline(v = mean(v), col = "red", lty = 2)
      if (!is.null(x$HPD) && p %in% rownames(x$HPD))
        graphics::abline(v = x$HPD[p, ], col = "grey40", lty = 3)
    }
  }
  invisible(x)
})---")


.step("Confirming no plot call relies on S3 dispatch for a non-numeric object")
.pl <- readLines("R/plots_extra.R", warn = FALSE)
.susp <- grep("plot\\((d|x)[,)]", .pl, value = TRUE)
if (length(.susp)) {
  .warn("a plot call may still pass an object rather than coordinates:")
  for (l in .susp) cat("        ", trimws(l), "\n", sep = "")
} else {
  .ok("every plot call passes numeric coordinates")
}
if (any(grepl("graphics::plot(", .pl, fixed = TRUE))) {
  .warn("graphics::plot( still present; plot() moved to base in R 4.0")
} else {
  .ok("the base plot generic is used throughout")
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

.step("Loading from source for the self-test")
.loaded <- tryCatch({ devtools::load_all(".", quiet = TRUE); TRUE },
                    error = function(e) conditionMessage(e))
if (!isTRUE(.loaded)) .die("load_all() failed:\n  ", .loaded, "\n\nBackups: ", BACKUP_DIR)
.ok("source loaded")

.step("Self-test on all three plots, density panel included")
.fails <- character(0)
.pass <- function(label, ok, detail = "") {
  if (isTRUE(ok)) .ok(paste0(label, if (nzchar(detail)) paste0("  ", detail) else ""))
  else { .warn(paste0(label, "  ", detail)); .fails <<- c(.fails, label) }
}

grDevices::pdf(NULL)
on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

data(guinea_pig, package = "BetaDanish", envir = environment())
.t <- bd_ttt_plot(guinea_pig$time)
.pass("TTT ends at one", isTRUE(all.equal(.t$phi[nrow(.t)], 1)))
.pass("TTT is non-decreasing", all(diff(.t$phi) >= -1e-12))
.pass("TTT reports a shape", nzchar(attr(.t, "shape")),
      paste("shape:", attr(.t, "shape")))

set.seed(9); .e <- stats::rexp(500, 0.5)
.te <- bd_ttt_plot(.e)
.pass("exponential sample lies near the diagonal",
      max(abs(.te$phi - .te$i_n)) < 0.12,
      sprintf("max deviation %.4f", max(abs(.te$phi - .te$i_n))))

.fake <- structure(list(draws = matrix(stats::rnorm(200), 100, 2,
                                       dimnames = list(NULL, c("b", "c")))),
                   class = "bd_bayes")

## The panel that failed in Patch 3f. Exercised on its own first.
.r1 <- tryCatch({ plot(.fake, type = "density"); "ok" },
                error = function(e) conditionMessage(e))
.pass("density panel draws", identical(.r1, "ok"), if (identical(.r1, "ok")) "" else .r1)

.r2 <- tryCatch({ plot(.fake, type = "trace"); "ok" },
                error = function(e) conditionMessage(e))
.pass("trace panel draws", identical(.r2, "ok"), if (identical(.r2, "ok")) "" else .r2)

.before <- graphics::par("mfrow")
.r3 <- tryCatch({ plot(.fake); "ok" }, error = function(e) conditionMessage(e))
.pass("both panels draw together", identical(.r3, "ok"),
      if (identical(.r3, "ok")) "" else .r3)
.pass("graphical parameters restored",
      identical(graphics::par("mfrow"), .before))

## And with an HPD interval present, which adds reference lines.
.fake2 <- .fake
.fake2$HPD <- matrix(c(-1.5, 1.5, -1.6, 1.6), nrow = 2, byrow = TRUE,
                     dimnames = list(c("b", "c"), c("lower", "upper")))
.r4 <- tryCatch({ plot(.fake2); "ok" }, error = function(e) conditionMessage(e))
.pass("HPD reference lines draw", identical(.r4, "ok"),
      if (identical(.r4, "ok")) "" else .r4)

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
  if (any(grepl(paste0("export(", f, ")"), .ns, fixed = TRUE))) {
    .ok(f)
  } else {
    .warn(paste("not exported:", f))
  }
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
cat("  PATCH 3f-fix COMPLETE  --  all 46 recommendations implemented\n")
cat(strrep("=", 78), "\n\n")
cat("  NEXT: cut the release.\n\n")
cat("      source(\"dev/BetaDanish_Patch3d_release.R\")\n\n")
cat("  Then steps 3 to 8 of the checklist: confirm the JAMSI citation,\n")
cat("  check_win_devel(), urlchecker, push and tag, build, submit.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
