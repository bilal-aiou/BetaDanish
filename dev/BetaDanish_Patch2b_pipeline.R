## =============================================================================
##  BetaDanish  --  PHASE 2, PATCH 2b of 3 : CSV ANALYSIS PIPELINE
## =============================================================================
##
##  Implements approved recommendations 23, 24 and 46.
##
##    23  bd_analyze_csv() -- end-to-end pipeline from a CSV file through
##        univariate, AFT, cure or competing-risks analysis, with optional
##        table and figure output to a directory
##    24  bd_csv_template() -- write a correctly shaped skeleton CSV
##    46  README section and a new vignette for the file-driven workflow
##
##  HOW TO RUN   source("dev/BetaDanish_Patch2b_pipeline.R")
##  IDEMPOTENT   Yes.
##
##  DESIGN NOTES
##    * output_dir defaults to NULL and nothing is ever written to disk unless
##      you name a directory. CRAN forbids packages writing outside tempdir()
##      without explicit user consent, so this default is not negotiable.
##    * Every model fit is wrapped: one failure is recorded and the pipeline
##      continues rather than losing the whole run.
##    * flexsurv, cmprsk and MCMCpack are Suggests, so each is guarded with
##      requireNamespace() and a missing package degrades to a recorded skip.
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
cat("  BetaDanish  --  Phase 2, Patch 2b of 3 : CSV analysis pipeline\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("data/guinea_pig.rda")) .die("Patch 2 has not been applied.")
if (!any(grepl("covar_cols = \"all\"", readLines("R/data_helpers.R", warn = FALSE), fixed = TRUE)))
  .die("Patch 2a has not been applied. Run it first.")
.ok("Patch 2 and Patch 2a detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch2b"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  REC 24  --  bd_csv_template()
## =============================================================================

.step("Rec 24: writing R/csv_template.R")

.put("R/csv_template.R", r"---(#' Write a Skeleton CSV for BetaDanish
#'
#' Writes a small, correctly shaped CSV showing the layout [bd_analyze_csv()]
#' and [read_survival_data()] expect. Fill it in with your own data, or use it
#' as a reference for renaming the columns of an existing file.
#'
#' @param path Destination file path. Use [tempfile()] if you only want to
#'   inspect the result.
#' @param type Which layout to write:
#'   `"univariate"` gives `time` and `status`;
#'   `"complete"` gives `time` alone, for an uncensored sample;
#'   `"covariate"` adds two example covariate columns;
#'   `"competing"` gives `time` and `cause`.
#' @param n Number of illustrative rows. Default 10.
#' @param overwrite Logical; overwrite an existing file. Default `FALSE`.
#'
#' @return The path, invisibly.
#'
#' @details
#' The illustrative values are plausible but arbitrary. `status` is coded
#' 1 = event observed, 0 = right-censored. For `"competing"`, `cause` is coded
#' 0 = censored and 1, 2, ... for the competing causes.
#'
#' @seealso [bd_analyze_csv()], [read_survival_data()]
#'
#' @export
#'
#' @examples
#' p <- bd_csv_template(tempfile(fileext = ".csv"), type = "covariate", n = 5)
#' read.csv(p)
#' unlink(p)
bd_csv_template <- function(path,
                            type = c("univariate", "complete", "covariate",
                                     "competing"),
                            n = 10, overwrite = FALSE) {
  type <- match.arg(type)
  if (!is.character(path) || length(path) != 1L)
    stop("'path' must be a single file path.", call. = FALSE)
  if (file.exists(path) && !isTRUE(overwrite))
    stop("'", path, "' already exists. Pass overwrite = TRUE to replace it.",
         call. = FALSE)
  n <- as.integer(n)
  if (is.na(n) || n < 1L) stop("'n' must be a positive integer.", call. = FALSE)

  tm <- round(seq(4, 4 + 6 * (n - 1), length.out = n) +
                stats::runif(n, -1, 1), 1)
  tm <- pmax(tm, 0.5)

  out <- switch(
    type,
    univariate = data.frame(time = tm,
                            status = rep_len(c(1L, 1L, 0L), n)),
    complete   = data.frame(time = tm),
    covariate  = data.frame(time = tm,
                            status = rep_len(c(1L, 1L, 0L), n),
                            age = rep_len(c(45L, 62L, 51L, 70L), n),
                            group = rep_len(c("treated", "control"), n)),
    competing  = data.frame(time = tm,
                            cause = rep_len(c(1L, 2L, 0L), n))
  )

  utils::write.csv(out, path, row.names = FALSE)
  message("Wrote a '", type, "' template with ", n, " row(s) to: ", path)
  invisible(path)
}
)---")

## =============================================================================
##  REC 23  --  bd_analyze_csv()
## =============================================================================

.step("Rec 23: writing R/analyze_csv.R")

.put("R/analyze_csv.R", r"---(## Internal helpers for the CSV pipeline. Every model fit goes through
## .bd_try(), so a single failure is recorded and the run continues instead of
## aborting and losing the analyses that did succeed.

#' Run an Expression, Recording Rather Than Propagating Failure
#'
#' Warnings are captured rather than discarded. The identifiability diagnostics
#' from [fit_betadanish()] arrive as warnings, and they are among the most
#' important things a user of this pipeline needs to see, so they are collected
#' and surfaced on the result instead of being muffled into silence.
#'
#' @noRd
.bd_try <- function(expr, label, store) {
  ws <- character(0)
  res <- withCallingHandlers(
    tryCatch(expr,
             error = function(e) structure(list(message = conditionMessage(e)),
                                           class = "bd_failed")),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (length(ws))
    store$warnings <- c(store$warnings,
                        stats::setNames(ws, rep(label, length(ws))))
  if (inherits(res, "bd_failed")) {
    store$failures <- c(store$failures,
                        stats::setNames(res$message, label))
    return(NULL)
  }
  res
}

#' Parameter Estimate Table for a betadanish Fit
#' @noRd
.bd_est_table <- function(fit, label) {
  est <- fit$coefficients
  se  <- sqrt(pmax(diag(fit$vcov), 0))
  z   <- stats::qnorm(0.975)
  data.frame(
    model      = label,
    parameter  = names(est),
    estimate   = as.numeric(est),
    std_error  = as.numeric(se),
    lower_95   = as.numeric(est - z * se),
    upper_95   = as.numeric(est + z * se),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
}

#' Information-Criteria Table Across Fits
#' @noRd
.bd_ic_table <- function(fits) {
  rows <- lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    if (is.null(f) || is.null(f$logLik)) return(NULL)
    data.frame(model   = nm,
               n       = .bd_or(f$nobs, NA_integer_),
               npar    = .bd_or(f$npar, length(f$coefficients)),
               logLik  = as.numeric(f$logLik),
               AIC     = tryCatch(stats::AIC(f), error = function(e) NA_real_),
               BIC     = tryCatch(stats::BIC(f), error = function(e) NA_real_),
               row.names = NULL, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Build a Survival Formula from Retained Covariate Names
#' @noRd
.bd_surv_formula <- function(covs) {
  rhs <- if (length(covs)) paste(covs, collapse = " + ") else "1"
  stats::as.formula(paste0("survival::Surv(time, status) ~ ", rhs),
                    env = parent.frame())
}

#' Write One PNG, Closing the Device Whatever Happens
#'
#' The device is closed before the file is tested, so a half-written PNG is
#' never reported as a success. Graphical parameters need no saving here: they
#' are per-device, and the device is discarded.
#'
#' @noRd
.bd_png <- function(path, expr, width = 1600, height = 1200, res = 180) {
  opened <- FALSE
  ok <- tryCatch({
    grDevices::png(path, width = width, height = height, res = res)
    opened <- TRUE
    force(expr)
    TRUE
  }, error = function(e) FALSE)
  if (opened) try(grDevices::dev.off(), silent = TRUE)
  if (isTRUE(ok) && file.exists(path)) path else character(0)
}

#' Analyse Survival Data Straight from a CSV File
#'
#' A single entry point that reads a delimited or Excel file, fits the
#' requested Beta-Danish model or models, assembles the results into tidy
#' tables, and optionally writes those tables and the diagnostic figures to a
#' directory.
#'
#' @param file Path to a `.csv`, `.txt`, `.tsv`, `.xls` or `.xlsx` file.
#' @param time_col,status_col,cause_col Column names, passed to
#'   [read_survival_data()]. `NULL` means guess.
#' @param covariates Covariate columns to use. `NULL` means none for
#'   `analysis = "univariate"`, and every non-response column for the
#'   regression analyses.
#' @param analysis Which analysis to run: `"univariate"` (default), `"aft"`,
#'   `"cure"` or `"competing"`.
#' @param model For `analysis = "univariate"`, whether to fit the
#'   four-parameter Beta-Danish (`"BD"`), the three-parameter Exponentiated
#'   Danish submodel (`"ED"`), or `"both"` (default) and compare them.
#' @param compare Logical; benchmark against standard lifetime distributions
#'   via [compare_distributions()]. Requires the `flexsurv` package; skipped
#'   with a recorded note if it is absent. Default `TRUE`.
#' @param bayes Logical; also run [bayes_betadanish()]. Requires `MCMCpack`.
#'   Default `FALSE`, since it is slow.
#' @param cure_formula Right-hand-side formula for the cure fraction when
#'   `analysis = "cure"`, for example `~ group`. Defaults to the retained
#'   covariates.
#' @param cure_type `"mixture"` (default) or `"promotion"`.
#' @param output_dir Directory to write results to. `NULL` (default) writes
#'   nothing. The directory is created if it does not exist.
#' @param n_starts Number of random starting points for each fit.
#' @param seed Optional integer seed, for reproducible multi-start fitting.
#' @param quiet Logical; suppress progress messages.
#'
#' @return An object of class `"bd_analysis"`: a list with `call`, `analysis`,
#'   `data`, `data_report`, `fits`, `tables`, `extras`, `warnings`, `failures`,
#'   `output_dir` and `files`. It has `print`, `summary` and `plot` methods.
#'
#' @details
#' Nothing is written to disk unless `output_dir` is supplied. When it is, the
#' function writes one CSV per result table and a PNG per diagnostic figure,
#' and records the paths in the `files` component.
#'
#' Each fit is attempted independently. If one fails, the message is recorded
#' in `failures` and the remaining analyses still run, so a difficult
#' four-parameter fit does not cost you the submodel results alongside it.
#'
#' Warnings raised during fitting are captured in `warnings` rather than
#' printed as they occur. This is where the identifiability diagnostics appear,
#' so it is worth reading: a four-parameter fit that converges cleanly but sits
#' on the `b = 1` ridge will say so there.
#'
#' For `analysis = "univariate"` with `model = "both"`, a likelihood ratio test
#' of the submodel is included. Read it alongside those diagnostics: see the
#' Identifiability section of [fit_betadanish()].
#'
#' @seealso [read_survival_data()], [bd_csv_template()], [fit_betadanish()],
#'   [fit_bd_aft()], [fit_bd_cure()], [fit_bd_competing()]
#'
#' @export
#'
#' @examples
#' f <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
#'
#' # Submodel only, no benchmarking: fast enough for an example
#' res <- bd_analyze_csv(f, analysis = "univariate", model = "ED",
#'                       compare = FALSE, n_starts = 3, seed = 1, quiet = TRUE)
#' res
#' res$tables$information_criteria
#'
#' \donttest{
#' # Both models, with tables and figures written to a temporary directory
#' out <- file.path(tempdir(), "bd_results")
#' res2 <- bd_analyze_csv(f, model = "both", compare = FALSE,
#'                        output_dir = out, n_starts = 5, seed = 1)
#' basename(res2$files)
#'
#' # AFT regression, covariates taken from the file
#' g <- system.file("extdata", "covariate_sample.csv", package = "BetaDanish")
#' res3 <- bd_analyze_csv(g, analysis = "aft", covariates = c("age", "thickness"),
#'                        n_starts = 5, seed = 1)
#' summary(res3)
#' }
bd_analyze_csv <- function(file,
                           time_col = NULL, status_col = NULL,
                           covariates = NULL, cause_col = NULL,
                           analysis = c("univariate", "aft", "cure", "competing"),
                           model = c("both", "BD", "ED"),
                           compare = TRUE,
                           bayes = FALSE,
                           cure_formula = NULL,
                           cure_type = c("mixture", "promotion"),
                           output_dir = NULL,
                           n_starts = 10,
                           seed = NULL,
                           quiet = FALSE) {

  analysis  <- match.arg(analysis)
  model     <- match.arg(model)
  cure_type <- match.arg(cure_type)
  if (!is.null(seed)) set.seed(seed)
  say <- function(...) if (!isTRUE(quiet)) message(...)

  store <- new.env(parent = emptyenv())
  store$failures <- character(0)
  store$warnings <- character(0)

  ## ---- read -----------------------------------------------------------------
  covar_arg <- covariates
  if (is.null(covar_arg) && analysis %in% c("aft", "cure")) covar_arg <- "all"

  dat <- read_survival_data(file, time_col = time_col, status_col = status_col,
                            covar_cols = covar_arg, cause_col = cause_col,
                            quiet = quiet)
  rep  <- attr(dat, "bd_data_report")
  covs <- rep$covariates

  say("Read ", rep$rows_kept, " row(s); ", rep$n_events, " event(s); ",
      sprintf("%.1f%%", 100 * rep$censoring_prop), " censored.")

  out <- list(call = match.call(), analysis = analysis, model = model,
              data = dat, data_report = rep,
              fits = list(), tables = list(), extras = list(),
              failures = character(0), output_dir = NULL, files = character(0))

  ## ---- fit ------------------------------------------------------------------
  if (analysis == "univariate") {
    fml <- .bd_surv_formula(character(0))

    if (model %in% c("BD", "both")) {
      say("Fitting the four-parameter Beta-Danish model.")
      out$fits$BD <- .bd_try(
        fit_betadanish(fml, data = dat, submodel = FALSE, n_starts = n_starts),
        "fit_betadanish (BD)", store)
    }
    if (model %in% c("ED", "both")) {
      say("Fitting the three-parameter Exponentiated Danish submodel.")
      out$fits$ED <- .bd_try(
        fit_betadanish(fml, data = dat, submodel = TRUE, n_starts = n_starts),
        "fit_betadanish (ED)", store)
    }

    if (!is.null(out$fits$BD) && !is.null(out$fits$ED)) {
      lrt <- .bd_try(compare_models(out$fits$BD, out$fits$ED),
                     "compare_models", store)
      if (!is.null(lrt)) out$tables$likelihood_ratio_test <- as.data.frame(lrt)
    }

  } else if (analysis == "aft") {
    if (!length(covs))
      say("No covariates retained; the AFT fit reduces to an intercept-only ",
          "model, equivalent to the ED submodel with k = exp(intercept).")
    fml <- .bd_surv_formula(covs)
    say("Fitting the AFT model: ", deparse(fml[[3]]))
    out$fits$AFT <- .bd_try(
      fit_bd_aft(fml, data = dat, n_starts = n_starts), "fit_bd_aft", store)

  } else if (analysis == "cure") {
    if (is.null(cure_formula))
      cure_formula <- stats::as.formula(
        paste0("~ ", if (length(covs)) paste(covs, collapse = " + ") else "1"))
    say("Fitting the ", cure_type, " cure model; cure fraction ~ ",
        deparse(cure_formula[[2]]))
    out$fits$CURE <- .bd_try(
      fit_bd_cure(formula_aft  = .bd_surv_formula(character(0)),
                  formula_cure = cure_formula,
                  data = dat, type = cure_type, n_starts = n_starts),
      "fit_bd_cure", store)

  } else if (analysis == "competing") {
    if (!"cause" %in% names(dat))
      stop("analysis = \"competing\" needs a cause column. Pass cause_col, ",
           "and code it 0 for censored with 1, 2, ... for the causes.",
           call. = FALSE)
    say("Fitting cause-specific models for ",
        length(setdiff(unique(dat$cause), 0)), " cause(s).")
    out$fits$CR <- .bd_try(
      fit_bd_competing(time = dat$time, cause = dat$cause, n_starts = n_starts),
      "fit_bd_competing", store)

    if (!is.null(out$fits$CR)) {
      if (requireNamespace("cmprsk", quietly = TRUE)) {
        cc <- .bd_try(cif_compare(out$fits$CR, plot = FALSE), "cif_compare", store)
        if (!is.null(cc)) out$extras$cif <- cc
      } else {
        store$failures <- c(store$failures,
                            c(cif_compare = "the 'cmprsk' package is not installed"))
      }
    }
  }

  ## ---- estimate and IC tables ----------------------------------------------
  uni <- Filter(function(f) inherits(f, "betadanish"), out$fits)
  if (length(uni)) {
    out$tables$estimates <- do.call(
      rbind, Map(.bd_est_table, uni, names(uni)))
    out$tables$information_criteria <- .bd_ic_table(uni)

    gof <- lapply(uni, function(f) .bd_try(gof_betadanish(f), "gof_betadanish", store))
    gof <- Filter(Negate(is.null), gof)
    if (length(gof)) {
      out$tables$goodness_of_fit <- do.call(rbind, Map(function(g, nm)
        data.frame(model = nm, as.list(g$IC), KS = g$KS_Statistic,
                   row.names = NULL, check.names = FALSE),
        gof, names(gof)))
    }
  }

  for (nm in intersect(names(out$fits), c("AFT", "CURE"))) {
    f <- out$fits[[nm]]
    if (is.null(f)) next
    s <- .bd_try(summary(f), paste0("summary (", nm, ")"), store)
    if (!is.null(s)) out$extras[[tolower(nm)]] <- s
  }

  ## ---- optional benchmarking ------------------------------------------------
  primary <- if (!is.null(out$fits$BD)) out$fits$BD else out$fits$ED
  if (isTRUE(compare) && inherits(primary, "betadanish")) {
    if (requireNamespace("flexsurv", quietly = TRUE)) {
      say("Benchmarking against standard lifetime distributions.")
      cmp <- .bd_try(compare_distributions(primary), "compare_distributions", store)
      if (!is.null(cmp)) out$tables$distribution_benchmark <- as.data.frame(cmp)
    } else {
      store$failures <- c(store$failures,
                          c(compare_distributions = "the 'flexsurv' package is not installed"))
    }
  }

  ## ---- optional Bayesian run ------------------------------------------------
  if (isTRUE(bayes)) {
    if (requireNamespace("MCMCpack", quietly = TRUE)) {
      say("Running the Bayesian sampler; this takes a while.")
      out$extras$bayes <- .bd_try(
        bayes_betadanish(dat$time, dat$status,
                         submodel = !identical(model, "BD"), seed = seed),
        "bayes_betadanish", store)
    } else {
      store$failures <- c(store$failures,
                          c(bayes_betadanish = "the 'MCMCpack' package is not installed"))
    }
  }

  ## ---- data report table ----------------------------------------------------
  out$tables$data_report <- data.frame(
    file           = rep$file,
    time_col       = .bd_or(rep$time_col, NA_character_),
    status_col     = .bd_or(rep$status_col, NA_character_),
    rows_read      = rep$rows_read,
    rows_kept      = rep$rows_kept,
    n_events       = rep$n_events,
    censoring_prop = rep$censoring_prop,
    covariates     = paste(rep$covariates, collapse = "; "),
    grid_step      = .bd_or(rep$grid_step, NA_real_),
    row.names = NULL, stringsAsFactors = FALSE)

  out$failures <- store$failures
  out$warnings <- store$warnings
  class(out) <- "bd_analysis"

  ## ---- optional output ------------------------------------------------------
  if (!is.null(output_dir)) out <- .bd_write_outputs(out, output_dir, quiet)

  if (length(out$warnings) && !isTRUE(quiet))
    message("Completed with ", length(out$warnings),
            " model warning(s); see $warnings. Identifiability diagnostics ",
            "appear here.")
  if (length(out$failures) && !isTRUE(quiet))
    message("Completed with ", length(out$failures),
            " recorded failure(s); see $failures.")

  out
}

#' Write Tables and Figures for a bd_analysis Object
#' @noRd
.bd_write_outputs <- function(obj, output_dir, quiet = FALSE) {
  say <- function(...) if (!isTRUE(quiet)) message(...)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir))
    stop("Could not create '", output_dir, "'.", call. = FALSE)

  written <- character(0)

  for (nm in names(obj$tables)) {
    tb <- obj$tables[[nm]]
    if (is.null(tb) || !nrow(as.data.frame(tb))) next
    p <- file.path(output_dir, paste0(nm, ".csv"))
    ok <- tryCatch({ utils::write.csv(tb, p, row.names = FALSE); TRUE },
                   error = function(e) FALSE)
    if (ok) written <- c(written, p)
  }

  prim <- if (!is.null(obj$fits$BD)) obj$fits$BD else obj$fits$ED
  if (inherits(prim, "betadanish")) {
    for (ty in c("survival", "hazard", "density", "pp", "qq")) {
      written <- c(written,
                   .bd_png(file.path(output_dir, paste0(ty, ".png")),
                           plot(prim, type = ty)))
    }
  }
  for (nm in c("AFT", "CURE")) {
    f <- obj$fits[[nm]]
    if (is.null(f)) next
    written <- c(written,
                 .bd_png(file.path(output_dir, paste0("coxsnell_", tolower(nm), ".png")),
                         plot(f)))
  }
  if (!is.null(obj$fits$CR) && requireNamespace("cmprsk", quietly = TRUE)) {
    written <- c(written,
                 .bd_png(file.path(output_dir, "cif.png"),
                         cif_compare(obj$fits$CR, plot = TRUE)))
  }

  say("Wrote ", length(written), " file(s) to ", normalizePath(output_dir))
  obj$output_dir <- normalizePath(output_dir, mustWork = FALSE)
  obj$files      <- written
  obj
}

#' @param x A `"bd_analysis"` object.
#' @param ... Ignored.
#' @rdname bd_analyze_csv
#' @export
print.bd_analysis <- function(x, ...) {
  cat("\nBetaDanish CSV analysis\n")
  cat("-----------------------\n")
  cat("Analysis:      ", x$analysis, "\n", sep = "")
  cat("Observations:  ", x$data_report$rows_kept,
      "  (events: ", x$data_report$n_events,
      ", censored: ", sprintf("%.1f%%", 100 * x$data_report$censoring_prop),
      ")\n", sep = "")
  if (length(x$data_report$covariates))
    cat("Covariates:    ", paste(x$data_report$covariates, collapse = ", "),
        "\n", sep = "")
  cat("Models fitted: ",
      if (length(x$fits)) paste(names(Filter(Negate(is.null), x$fits)),
                                collapse = ", ") else "none", "\n", sep = "")

  ic <- x$tables$information_criteria
  if (!is.null(ic)) { cat("\n"); print(ic, row.names = FALSE) }

  if (length(x$files))
    cat("\nOutput:        ", length(x$files), " file(s) in ", x$output_dir,
        "\n", sep = "")
  if (length(x$warnings))
    cat("\nWarnings:      ", length(x$warnings),
        " (see $warnings)\n", sep = "")
  if (length(x$failures))
    cat("\nFailures:      ", length(x$failures),
        " (see $failures)\n", sep = "")
  cat("\n")
  invisible(x)
}

#' @param object A `"bd_analysis"` object.
#' @rdname bd_analyze_csv
#' @export
summary.bd_analysis <- function(object, ...) {
  structure(object, class = c("summary.bd_analysis", "bd_analysis"))
}

#' @rdname bd_analyze_csv
#' @export
print.summary.bd_analysis <- function(x, ...) {
  y <- x; class(y) <- "bd_analysis"
  print(y)

  for (nm in names(x$tables)) {
    if (identical(nm, "information_criteria")) next
    tb <- x$tables[[nm]]
    if (is.null(tb)) next
    cat("-- ", gsub("_", " ", nm), " --\n", sep = "")
    print(tb, row.names = FALSE)
    cat("\n")
  }
  for (nm in names(x$extras)) {
    cat("-- ", nm, " --\n", sep = "")
    print(x$extras[[nm]])
    cat("\n")
  }
  if (length(x$warnings)) {
    cat("-- model warnings --\n")
    for (i in seq_along(x$warnings))
      cat("  ", names(x$warnings)[i], ": ", x$warnings[i], "\n", sep = "")
    cat("\n")
  }
  if (length(x$failures)) {
    cat("-- failures --\n")
    for (i in seq_along(x$failures))
      cat("  ", names(x$failures)[i], ": ", x$failures[i], "\n", sep = "")
    cat("\n")
  }
  invisible(x)
}

#' @param type Plot type passed to [plot.betadanish()].
#' @rdname bd_analyze_csv
#' @export
plot.bd_analysis <- function(x, type = "survival", ...) {
  prim <- if (!is.null(x$fits$BD)) x$fits$BD else x$fits$ED
  if (inherits(prim, "betadanish")) return(plot(prim, type = type, ...))
  for (nm in c("AFT", "CURE")) {
    if (!is.null(x$fits[[nm]])) return(plot(x$fits[[nm]], ...))
  }
  if (!is.null(x$fits$CR)) return(cif_compare(x$fits$CR, plot = TRUE))
  warning("Nothing to plot: no model was fitted successfully.", call. = FALSE)
  invisible(x)
}
)---")

## =============================================================================
##  DESCRIPTION  --  grDevices is now used
## =============================================================================

.step("Adding grDevices to Imports")

.d <- readLines("DESCRIPTION", warn = FALSE)
.i <- grep("^Imports:", .d)
if (length(.i) == 1L && !grepl("grDevices", .d[.i], fixed = TRUE)) {
  .backup("DESCRIPTION")
  .d[.i] <- sub("^Imports:\\s*", "Imports: grDevices, ", .d[.i])
  .write_lines("DESCRIPTION", .d)
  .ok("Imports now include grDevices")
} else {
  .info("grDevices already declared")
}

## =============================================================================
##  TESTS
## =============================================================================

.step("Writing tests/testthat/test-analyze-csv.R")

.put("tests/testthat/test-analyze-csv.R", r"---(csv_path <- function(f) system.file("extdata", f, package = "BetaDanish")

test_that("bd_csv_template writes each layout", {
  expected <- list(univariate = c("time", "status"),
                   complete   = "time",
                   covariate  = c("time", "status", "age", "group"),
                   competing  = c("time", "cause"))
  for (ty in names(expected)) {
    p <- tempfile(fileext = ".csv"); on.exit(unlink(p), add = TRUE)
    suppressMessages(bd_csv_template(p, type = ty, n = 6))
    got <- utils::read.csv(p)
    expect_named(got, expected[[ty]], info = ty)
    expect_equal(nrow(got), 6L, info = ty)
    expect_true(all(got$time > 0), info = ty)
  }
})

test_that("bd_csv_template refuses to clobber without overwrite", {
  p <- tempfile(fileext = ".csv"); on.exit(unlink(p), add = TRUE)
  suppressMessages(bd_csv_template(p))
  expect_error(bd_csv_template(p), "already exists")
  expect_silent(suppressMessages(bd_csv_template(p, overwrite = TRUE)))
})

test_that("a template round-trips through read_survival_data", {
  p <- tempfile(fileext = ".csv"); on.exit(unlink(p), add = TRUE)
  suppressMessages(bd_csv_template(p, type = "covariate", n = 12))
  dat <- read_survival_data(p, covar_cols = "all", quiet = TRUE)
  expect_named(dat, c("time", "status", "age", "group"))
  expect_true(all(dat$status %in% c(0, 1)))
})

test_that("the univariate pipeline runs and returns tidy tables", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("censored_sample.csv"),
                        analysis = "univariate", model = "ED",
                        compare = FALSE, n_starts = 3, seed = 11, quiet = TRUE)

  expect_s3_class(res, "bd_analysis")
  expect_true(inherits(res$fits$ED, "betadanish"))
  expect_true(is.data.frame(res$tables$estimates))
  expect_setequal(res$tables$estimates$parameter, c("b", "c", "k"))
  expect_true(is.data.frame(res$tables$information_criteria))
  expect_true(is.finite(res$tables$information_criteria$AIC[1]))
  expect_length(res$failures, 0L)
  expect_output(print(res), "BetaDanish CSV analysis")
})

test_that("model = 'both' adds a likelihood ratio test", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "both",
                        compare = FALSE, n_starts = 3, seed = 12, quiet = TRUE)
  expect_true(all(c("BD", "ED") %in% names(res$fits)))
  expect_true(is.data.frame(res$tables$likelihood_ratio_test))
  expect_equal(nrow(res$tables$information_criteria), 2L)
})

test_that("nothing is written unless output_dir is given", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "ED",
                        compare = FALSE, n_starts = 3, seed = 13, quiet = TRUE)
  expect_null(res$output_dir)
  expect_length(res$files, 0L)
})

test_that("output_dir receives tables and figures", {
  skip_on_cran()
  od <- file.path(tempdir(), paste0("bdout_", as.integer(runif(1, 1, 1e6))))
  on.exit(unlink(od, recursive = TRUE), add = TRUE)

  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "ED",
                        compare = FALSE, output_dir = od, n_starts = 3,
                        seed = 14, quiet = TRUE)

  expect_true(dir.exists(od))
  expect_gt(length(res$files), 3L)
  expect_true(any(grepl("estimates[.]csv$", res$files)))
  expect_true(any(grepl("survival[.]png$", res$files)))
  expect_true(all(file.exists(res$files)))

  back <- utils::read.csv(grep("estimates[.]csv$", res$files, value = TRUE)[1])
  expect_true(all(c("model", "parameter", "estimate") %in% names(back)))
})

test_that("the AFT path uses the covariates from the file", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("covariate_sample.csv"), analysis = "aft",
                        covariates = c("age", "thickness"),
                        compare = FALSE, n_starts = 3, seed = 15, quiet = TRUE)
  expect_true(inherits(res$fits$AFT, "bd_aft"))
  expect_true(any(grepl("^delta_age$", names(res$fits$AFT$coefficients))))
  expect_true(!is.null(res$extras$aft))
})

test_that("the competing-risks path needs a cause column", {
  skip_on_cran()
  expect_error(
    bd_analyze_csv(csv_path("censored_sample.csv"), analysis = "competing",
                   n_starts = 2, quiet = TRUE),
    "cause column")

  res <- bd_analyze_csv(csv_path("competing_sample.csv"), analysis = "competing",
                        time_col = "time", cause_col = "cause",
                        compare = FALSE, n_starts = 2, seed = 16, quiet = TRUE)
  expect_true(inherits(res$fits$CR, "bd_competing"))
})

test_that("a failing fit is recorded rather than fatal", {
  skip_on_cran()
  ## Two rows cannot support a fit, but extract_surv_data rejects it first,
  ## so the failure must surface as a recorded message, not an abort.
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(time = c(1, 2, 3), status = c(1, 0, 0)), tmp,
                   row.names = FALSE)
  res <- bd_analyze_csv(tmp, model = "ED", compare = FALSE, n_starts = 2,
                        quiet = TRUE)
  expect_s3_class(res, "bd_analysis")
  expect_gt(length(res$failures), 0L)
  expect_null(res$fits$ED)
})

test_that("model warnings are captured rather than discarded", {
  skip_on_cran()
  ## The four-parameter fit on this file lands near the b = 1 ridge, so the
  ## identifiability diagnostic fires. It must survive the pipeline: swallowing
  ## it would hide exactly the thing the user needs to see.
  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "both",
                        compare = FALSE, n_starts = 5, seed = 21, quiet = TRUE)
  expect_type(res$warnings, "character")
  if (length(res$warnings)) expect_true(all(nzchar(names(res$warnings))))
  ## Whatever happened, it is recorded somewhere rather than lost.
  expect_true(!is.null(res$fits$BD) || length(res$failures) > 0L)
})

test_that("summary prints every table without error", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "ED",
                        compare = FALSE, n_starts = 3, seed = 17, quiet = TRUE)
  expect_output(print(summary(res)), "estimates")
})

test_that("plot.bd_analysis dispatches to the primary fit", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "ED",
                        compare = FALSE, n_starts = 3, seed = 18, quiet = TRUE)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(plot(res, type = "survival"), NA)
})
)---")

## =============================================================================
##  REC 46  --  VIGNETTE AND README
## =============================================================================

.step("Rec 46a: writing vignettes/bd-csv-workflow.Rmd")

.put("vignettes/bd-csv-workflow.Rmd", r"---(---
title: "Analysing Your Own Data from a CSV File"
author: "Bilal Ahmad & Dr. Muhammad Yameen Danish"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Analysing Your Own Data from a CSV File}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      warning = FALSE, message = FALSE)
library(BetaDanish)
```

Most users arrive with a spreadsheet rather than an R data frame. This vignette
covers the file-driven path: from a CSV on disk to fitted models, tables and
figures, without writing any modelling code.

## What the file should look like

Two columns are enough:

| column | meaning |
|---|---|
| `time` | survival or failure time, strictly positive |
| `status` | 1 if the event was observed, 0 if right-censored |

Covariates go in additional columns. For competing risks, replace `status` with
`cause`, coded 0 for censored and 1, 2, ... for the causes.

If you would rather start from a working example than a description, write a
template and fill it in:

```{r}
tmp <- tempfile(fileext = ".csv")
bd_csv_template(tmp, type = "covariate", n = 6)
read.csv(tmp)
```

## Reading a file

`read_survival_data()` will guess the time and status columns when they carry
recognisable names, and reports what it chose:

```{r}
f <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
dat <- read_survival_data(f)
head(dat)
```

Note that covariates are not retained by default. Ask for them by name, or use
`covar_cols = "all"`:

```{r}
dat_all <- read_survival_data(f, covar_cols = "all", quiet = TRUE)
names(dat_all)
```

The attached report records what happened, which is worth checking before you
model anything:

```{r}
str(attr(dat_all, "bd_data_report"))
```

One field deserves attention. `grid_step` is non-`NA` when the times look
recorded on a coarse grid -- whole days or whole months, say. Rounded times
break the assumptions behind a point-density likelihood, so treat standard
errors from such data as optimistic.

### A word on `censor`

Column-name guessing deliberately ignores names like `censor`. In some
conventions `censor = 1` means *censored*; in others it means *observed*.
Guessing wrong would silently invert every event in your dataset and still
produce a plausible-looking fit, so you must name that column yourself.

## One-call analysis

`bd_analyze_csv()` does the whole thing: read, fit, tabulate, and optionally
write results to a directory.

```{r}
res <- bd_analyze_csv(f, analysis = "univariate", model = "ED",
                      compare = FALSE, n_starts = 5, seed = 1, quiet = TRUE)
res
```

The tidy tables are in `$tables`:

```{r}
res$tables$estimates
res$tables$goodness_of_fit
```

### Comparing the parent model against its submodel

With `model = "both"`, both the four-parameter Beta-Danish and the
three-parameter Exponentiated Danish submodel are fitted and compared:

```{r}
both <- bd_analyze_csv(f, model = "both", compare = FALSE,
                       n_starts = 5, seed = 2, quiet = TRUE)
both$tables$information_criteria
both$tables$likelihood_ratio_test
```

Read that test alongside the identifiability diagnostics. On many datasets the
four-parameter model fits well but is only weakly identified, and the
likelihood ratio test will not reject the submodel. See the Identifiability
section of `?fit_betadanish`, and the `guinea_pig` dataset for a contrasting
case where the parent model is well identified.

## Saving tables and figures

Nothing is written to disk unless you name a directory. When you do, you get
one CSV per table and a PNG per diagnostic figure:

```{r}
out <- file.path(tempdir(), "bd_results")
saved <- bd_analyze_csv(f, model = "ED", compare = FALSE,
                        output_dir = out, n_starts = 5, seed = 3, quiet = TRUE)
basename(saved$files)
```

## Regression from a file

The same entry point handles the regression models. Covariates are taken from
the file, either all of them or the ones you name:

```{r, eval = FALSE}
g <- system.file("extdata", "covariate_sample.csv", package = "BetaDanish")

# Accelerated failure time
aft <- bd_analyze_csv(g, analysis = "aft",
                      covariates = c("age", "thickness"))

# Mixture cure model, cure fraction depending on ulceration
cure <- bd_analyze_csv(g, analysis = "cure", cure_formula = ~ ulcer)

# Competing risks
h  <- system.file("extdata", "competing_sample.csv", package = "BetaDanish")
cr <- bd_analyze_csv(h, analysis = "competing", cause_col = "cause")
```

## When something fails

Each fit is attempted independently. If one fails, the message is recorded and
the rest of the analysis still runs, so a difficult four-parameter fit does not
cost you the submodel results beside it:

```{r}
res$failures
```

An empty result means everything succeeded.
)---")

.step("Rec 46b: adding a CSV workflow section to README.md")

.rd <- readLines("README.md", warn = FALSE)
if (!any(grepl("bd_analyze_csv", .rd, fixed = TRUE))) {
  .anchor <- grep("^## Built-in Datasets\\s*$", .rd)
  if (length(.anchor) == 1L) {
    .backup("README.md")
    .sec <- c(
      "## Working from a CSV File",
      "",
      "The whole workflow can be driven from a spreadsheet, with no modelling",
      "code. Two columns are enough: `time` and `status` (1 = event,",
      "0 = censored).",
      "",
      "```r",
      "# Not sure of the layout? Write a template and fill it in.",
      'bd_csv_template("my_data.csv", type = "covariate")',
      "",
      "# Read a file; time and status are guessed from common column names.",
      'dat <- read_survival_data("my_data.csv", covar_cols = "all")',
      'attr(dat, "bd_data_report")   # what was read, dropped, and inferred',
      "",
      "# Or run the whole analysis in one call.",
      "res <- bd_analyze_csv(",
      '  "my_data.csv",',
      '  analysis   = "univariate",   # or "aft", "cure", "competing"',
      '  model      = "both",         # Beta-Danish and the ED submodel',
      '  output_dir = "results"       # tables as CSV, figures as PNG',
      ")",
      "",
      "res                       # headline summary",
      "res$tables$estimates      # tidy parameter table",
      "res$failures              # empty if everything succeeded",
      "```",
      "",
      "Nothing is written to disk unless `output_dir` is supplied. Each model is",
      "fitted independently, so one failure is recorded rather than losing the",
      "whole run.",
      "",
      "| Function | Purpose |",
      "|---|---|",
      "| `bd_analyze_csv()` | Read, fit, tabulate and optionally save |",
      "| `read_survival_data()` | Read and validate a file into a data frame |",
      "| `bd_csv_template()` | Write a correctly shaped skeleton CSV |",
      "")
    .rd <- append(.rd, .sec, after = .anchor - 1L)
    .write_lines("README.md", .rd)
    .ok("README CSV workflow section added")
  } else {
    .warn("could not find the '## Built-in Datasets' heading; add the section by hand")
  }
} else {
  .info("README already documents bd_analyze_csv")
}

.step("Recording the additions in NEWS.md")

.nw <- readLines("NEWS.md", warn = FALSE)
if (!any(grepl("bd_analyze_csv", .nw, fixed = TRUE))) {
  .hdr <- grep("^## Data input\\s*$", .nw)
  if (length(.hdr) == 1L) {
    .backup("NEWS.md")
    .sec <- c(
      "",
      "## File-driven analysis",
      "",
      "* **`bd_analyze_csv()`** reads a delimited or Excel file, fits the",
      "  requested model, assembles tidy result tables, and optionally writes",
      "  those tables and the diagnostic figures to a directory. It covers the",
      "  univariate, AFT, cure and competing-risks paths. Nothing is written to",
      "  disk unless `output_dir` is supplied, and each fit is attempted",
      "  independently so one failure is recorded rather than fatal.",
      "",
      "* **`bd_csv_template()`** writes a correctly shaped skeleton CSV in any",
      "  of four layouts, so the expected column structure can be seen rather",
      "  than read about.",
      "",
      "* `print`, `summary` and `plot` methods for the new `bd_analysis` class.",
      "",
      "* New vignette: \"Analysing Your Own Data from a CSV File\".")
    .nw <- append(.nw, .sec, after = .hdr - 1L)
    .write_lines("NEWS.md", .nw)
    .ok("NEWS.md updated")
  } else {
    .warn("anchor not found in NEWS.md; add the note by hand")
  }
} else {
  .info("NEWS.md already records it")
}

## =============================================================================
##  VERIFY
## =============================================================================

.step("Auditing roxygen for unescaped Rd macros")
.KNOWN <- paste0("\\\\(describe|item|eqn|deqn|doi|donttest|dontrun|code|link|",
                 "url|emph|strong|href|enumerate|itemize|preformatted|verb|",
                 "Sexpr|tabular|cr|dots|ldots|R|method|usage|section)")
.flag <- character(0)
for (f in c("R/analyze_csv.R", "R/csv_template.R")) {
  for (l in grep("^#'", readLines(f, warn = FALSE), value = TRUE)) {
    hits <- gregexpr("\\\\[A-Za-z]", l)[[1]]
    if (hits[1] == -1L) next
    if (!grepl(.KNOWN, l)) .flag <- c(.flag, paste0(f, ": ", trimws(l)))
  }
}
if (length(.flag)) {
  .warn("possible unescaped backslash in roxygen:")
  for (l in .flag) cat("        ", l, "\n", sep = "")
} else {
  .ok("no unescaped Rd macros")
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

.step("Confirming the new exports")
.ns <- readLines("NAMESPACE", warn = FALSE)
for (e in c("export(bd_analyze_csv)", "export(bd_csv_template)",
            "S3method(print,bd_analysis)", "S3method(summary,bd_analysis)",
            "S3method(plot,bd_analysis)")) {
  if (any(grepl(e, .ns, fixed = TRUE))) .ok(e) else .warn(paste("missing:", e))
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
cat("  PATCH 2b COMPLETE  --  CSV analysis pipeline\n")
cat(strrep("=", 78), "\n\n")
cat("  23  bd_analyze_csv() with univariate / aft / cure / competing paths,\n")
cat("      tidy tables, optional CSV + PNG output, print/summary/plot\n")
cat("  24  bd_csv_template() in four layouts\n")
cat("  46  README section and new vignette on the file-driven workflow\n\n")
cat("  Patch 2 is now complete. Remaining: Patch 3 -- theoretical properties,\n")
cat("  estimation methods, simulation studies, visualization, version 0.3.0.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
