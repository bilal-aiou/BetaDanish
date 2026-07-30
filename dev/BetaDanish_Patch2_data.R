## =============================================================================
##  BetaDanish  --  PHASE 2, PATCH 2 of 3 : DATA LAYER
## =============================================================================
##
##  Implements approved recommendations 19, 20, 21, 22.
##
##    19  Remove the brain_cancer dataset completely (Part C of the audit)
##    20  Add guinea_pig (Bjerkedal 1960, n = 72) as the identified-interior
##        exemplar that brain_cancer's removal would otherwise leave missing
##    21  Create inst/extdata/ with four example CSVs
##    22  Complete the read_survival_data() hardening begun in Patch 1
##
##  HOW TO RUN   From the package root, as before:
##                 source("dev/BetaDanish_Patch2_data.R")
##
##  IDEMPOTENT   Yes.
##  BACKUP       Everything it touches, including the deleted files, goes to
##               .betadanish_backup/<timestamp>-patch2/ first. The brain_cancer
##               data file is recoverable from there if anything goes wrong.
##  SCOPE        No new analysis functions. bd_analyze_csv() and
##               bd_csv_template() are Patch 2b; the theory, estimation,
##               simulation and visualization work is Patch 3.
## =============================================================================

if (getRversion() < "4.0.0") stop("This patch needs R >= 4.0.")

.step_n <- 0L
.step <- function(m) { .step_n <<- .step_n + 1L; cat(sprintf("\n[%02d] %s\n", .step_n, m)) }
.ok   <- function(m) cat("     OK   ", m, "\n", sep = "")
.info <- function(m) cat("     ..   ", m, "\n", sep = "")
.warn <- function(m) cat("     WARN ", m, "\n", sep = "")
.die  <- function(...) stop("\n\n*** PATCH ABORTED ***\n", ..., "\n", call. = FALSE)

BACKUP_DIR <- NULL
.backup <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  dest <- file.path(BACKUP_DIR, path)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(path, dest, overwrite = TRUE)) .die("Could not back up ", path)
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
cat("  BetaDanish  --  Phase 2, Patch 2 of 3 : data layer\n")
cat(strrep("=", 78), "\n")

## ------------------------------------------------------------- pre-flight ----

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
.d <- read.dcf("DESCRIPTION")
if (.d[1, "Package"] != "BetaDanish") .die("This is not the BetaDanish package.")

if (!any(grepl(".bd_conform", readLines("R/dist_functions.R", warn = FALSE), fixed = TRUE)))
  .die("Patch 1 has not been applied. Run Patch 1 and Patch 1a first.")
if (any(grepl("n <- max(length(x), length(a)",
              readLines("R/dist_functions.R", warn = FALSE), fixed = TRUE)))
  .die("Patch 1a has not been applied. Run it before this patch.")
.ok("Patch 1 and Patch 1a detected")
.ok(paste("version:", .d[1, "Version"]))

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch2"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  REC 19  --  REMOVE brain_cancer
## =============================================================================

.step("Rec 19a: removing data/brain_cancer.rda")

if (file.exists("data/brain_cancer.rda")) {
  .backup("data/brain_cancer.rda")
  if (!file.remove("data/brain_cancer.rda")) .die("Could not delete data/brain_cancer.rda")
  .ok("deleted (a copy is in the backup directory)")
} else {
  .info("already removed")
}

.step("Rec 19b: removing the roxygen block from R/data.R")

.dat <- readLines("R/data.R", warn = FALSE)
.start <- grep("^#' Brain Cancer Survival Data\\s*$", .dat)
.end   <- grep('^"brain_cancer"\\s*$', .dat)

if (length(.start) == 1L && length(.end) == 1L && .end > .start) {
  .backup("R/data.R")
  .write_lines("R/data.R", .dat[-(.start:.end)])
  .ok(sprintf("removed lines %d-%d (%d lines)", .start, .end, .end - .start + 1L))
} else if (length(.start) == 0L && length(.end) == 0L) {
  .info("already removed")
} else {
  .die("Could not bracket the brain_cancer block in R/data.R.\n",
       "  '#' Brain Cancer Survival Data' matches: ", length(.start), "\n",
       '  \'"brain_cancer"\' matches: ', length(.end), "\n",
       "Remove that roxygen block by hand, then re-run.")
}

.step("Rec 19c: removing man/brain_cancer.Rd")

if (file.exists("man/brain_cancer.Rd")) {
  .backup("man/brain_cancer.Rd")
  file.remove("man/brain_cancer.Rd")
  .ok("deleted")
} else {
  .info("already absent")
}

.step("Rec 19d: replacing the README AFT example and dataset table row")

.rd <- readLines("README.md", warn = FALSE)
.backup("README.md")

## The AFT example block, lines beginning data("brain_cancer") through the fence.
.i <- grep('^data\\("brain_cancer"\\)\\s*$', .rd)
if (length(.i) == 1L) {
  ## Find the closing fence after it.
  .fence <- .i + which(grepl("^```\\s*$", .rd[(.i + 1L):min(length(.rd), .i + 20L)]))[1]
  if (is.na(.fence)) .die("Could not find the closing fence of the README AFT example.")
  .replacement <- c(
    'data("melanoma")',
    'melanoma$event <- ifelse(melanoma$status == 1, 1, 0)',
    '',
    'fit_aft <- fit_bd_aft(',
    '  survival::Surv(time, event) ~ age + thickness,',
    '  data = melanoma',
    ')',
    'summary(fit_aft)',
    'plot(fit_aft)   # Cox-Snell residual diagnostic',
    '```')
  .rd <- append(.rd[-(.i:.fence)], .replacement, after = .i - 1L)
  .ok("AFT example now uses melanoma")
} else if (any(grepl('data("melanoma")', .rd, fixed = TRUE))) {
  .info("AFT example already replaced")
} else {
  .warn("README AFT example not found in the expected form; check it by hand")
}

## Dataset table: drop the brain_cancer row, add guinea_pig.
.row <- grep("^\\| `brain_cancer` \\|", .rd)
if (length(.row)) {
  .rd <- .rd[-.row]
  .ok("brain_cancer row removed from the dataset table")
}
if (!any(grepl("`guinea_pig`", .rd, fixed = TRUE))) {
  .mel <- grep("^\\| `melanoma` \\|", .rd)
  if (length(.mel) == 1L) {
    .rd <- append(
      .rd,
      "| `guinea_pig` | 72 | Guinea pig survival, virulent tubercle bacilli (days) |",
      after = .mel)
    .ok("guinea_pig row added to the dataset table")
  } else {
    .warn("could not locate the melanoma table row; add the guinea_pig row by hand")
  }
} else {
  .info("guinea_pig row already present")
}

.write_lines("README.md", .rd)

.step("Rec 19e: removing orphaned brain-cancer vocabulary from inst/WORDLIST")

if (file.exists("inst/WORDLIST")) {
  .wl <- readLines("inst/WORDLIST", warn = FALSE)
  .drop <- c("NORI", "Comorbid", "comorbidities")
  .keep <- !(trimws(.wl) %in% .drop)
  if (!all(.keep)) {
    .backup("inst/WORDLIST")
    .wl <- .wl[.keep]
    ## Add the vocabulary the new dataset needs.
    .add <- setdiff(c("Bjerkedal", "bacilli", "tubercle"), trimws(.wl))
    if (length(.add)) .wl <- sort(c(.wl[nzchar(.wl)], .add))
    .write_lines("inst/WORDLIST", .wl)
    .ok(sprintf("removed %d term(s), added %d", sum(!.keep), length(.add)))
  } else {
    .info("already tidied")
  }
}

## =============================================================================
##  REC 20  --  ADD guinea_pig
## =============================================================================

.step("Rec 20a: writing data-raw/guinea_pig.R and generating data/guinea_pig.rda")

dir.create("data-raw", showWarnings = FALSE)

.put("data-raw/guinea_pig.R", r"---(## Guinea pig survival data (Bjerkedal 1960).
##
## Survival times in days of 72 guinea pigs injected with virulent tubercle
## bacilli, from the principal regimen of the study. This is a complete sample:
## every animal was observed to death, so status is 1 throughout.
##
## Included because it is the one dataset in the thesis on which the
## four-parameter Beta-Danish model attains a genuine finite interior optimum,
## with b-hat = 3.64 (SE 1.20), placing b about 2.2 standard errors clear of the
## b = 1 identifiability ridge. Every other application sits on the flat
## (a, c) direction with a-hat < 1.
##
## Re-run with:  source("data-raw/guinea_pig.R")

gp <- c( 12,  15,  22,  24,  24,  32,  32,  33,  34,  38,
         38,  43,  44,  48,  52,  53,  54,  54,  55,  56,
         57,  58,  58,  59,  60,  60,  60,  60,  61,  62,
         63,  65,  65,  67,  68,  70,  70,  72,  73,  75,
         76,  76,  81,  83,  84,  85,  87,  91,  95,  96,
         98,  99, 109, 110, 121, 127, 129, 131, 143, 146,
        146, 175, 175, 211, 233, 258, 258, 263, 297, 341,
        341, 376)

stopifnot(length(gp) == 72L, !anyNA(gp), all(gp > 0),
          identical(gp, sort(gp)))

guinea_pig <- data.frame(time = as.numeric(gp), status = 1L)

save(guinea_pig, file = "data/guinea_pig.rda", version = 2, compress = "xz")
)---")

.res <- tryCatch({ source("data-raw/guinea_pig.R", local = new.env()); TRUE },
                 error = function(e) conditionMessage(e))
if (!isTRUE(.res)) .die("Generating guinea_pig.rda failed:\n  ", .res)
if (!file.exists("data/guinea_pig.rda")) .die("data/guinea_pig.rda was not created.")

.chkenv <- new.env(); load("data/guinea_pig.rda", envir = .chkenv)
.gp <- .chkenv$guinea_pig
if (!is.data.frame(.gp) || nrow(.gp) != 72L || !all(c("time", "status") %in% names(.gp)))
  .die("data/guinea_pig.rda does not have the expected shape.")
.ok(sprintf("guinea_pig: %d rows, mean time %.2f days, all uncensored",
            nrow(.gp), mean(.gp$time)))

.step("Rec 20b: documenting guinea_pig in R/data.R")

.dat <- readLines("R/data.R", warn = FALSE)
if (!any(grepl('^"guinea_pig"\\s*$', .dat))) {
  .backup("R/data.R")
  .doc <- c(
    "#' Guinea Pig Survival Times",
    "#'",
    "#' Survival times, in days, of 72 guinea pigs injected with virulent",
    "#' tubercle bacilli, taken from the principal regimen of the study of",
    "#' Bjerkedal (1960). A complete sample: every animal was observed to death,",
    "#' so `status` is 1 throughout.",
    "#'",
    "#' @format A data frame with 72 rows and 2 columns:",
    "#' \\describe{",
    "#'   \\item{time}{Survival time in days}",
    "#'   \\item{status}{Event indicator, 1 for all observations (no censoring)}",
    "#' }",
    "#'",
    "#' @details",
    "#' These data are the reference case for a well-identified four-parameter",
    "#' fit. The scaled total-time-on-test transform is unimodal and the upper",
    "#' tail is determinate rather than heavy, and the maximum-likelihood fit of",
    "#' the full Beta-Danish model attains a finite interior optimum: every",
    "#' estimate lies strictly inside the parameter space, with",
    "#' \\eqn{\\hat b = 3.64} (Wald standard error 1.20), so that",
    "#' \\eqn{(\\hat b - 1)/\\mathrm{SE} = 2.20} places \\eqn{b} more than two",
    "#' standard errors clear of the \\eqn{b = 1} identifiability ridge.",
    "#'",
    "#' By contrast, on `remission` and `carbon_fibres` the four-parameter fit",
    "#' sits on the flat \\eqn{(a, c)} direction of the likelihood with",
    "#' \\eqn{\\hat a < 1} and is only weakly identified. Use this dataset when",
    "#' you want to see the parent model behaving well; see the Identifiability",
    "#' section of [fit_betadanish()] for what to watch for elsewhere.",
    "#'",
    "#' @source Bjerkedal, T. (1960). Acquisition of resistance in guinea pigs",
    "#'   infected with different doses of virulent tubercle bacilli.",
    "#'   *American Journal of Epidemiology*, 72(1), 130-148.",
    "#'   \\doi{10.1093/oxfordjournals.aje.a120129}",
    "#'",
    "#' @examples",
    "#' data(guinea_pig)",
    "#' summary(guinea_pig$time)",
    "#' \\donttest{",
    "#' fit_full <- fit_betadanish(survival::Surv(time, status) ~ 1,",
    "#'                            data = guinea_pig)",
    "#' fit_sub  <- fit_betadanish(survival::Surv(time, status) ~ 1,",
    "#'                            data = guinea_pig, submodel = TRUE)",
    "#' compare_models(fit_full, fit_sub)",
    "#' }",
    '"guinea_pig"')
  .write_lines("R/data.R", c(.dat, .doc))
  .ok("guinea_pig documented")
} else {
  .info("guinea_pig already documented")
}

## =============================================================================
##  REC 21  --  EXAMPLE CSVs IN inst/extdata
## =============================================================================

.step("Rec 21: generating inst/extdata example CSVs")

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

.load1 <- function(nm) {
  e <- new.env(); load(file.path("data", paste0(nm, ".rda")), envir = e); e[[nm]]
}

## 1. Complete uncensored sample, deliberately using a non-standard column name
##    so the examples exercise explicit column selection.
.gp <- .load1("guinea_pig")
utils::write.csv(data.frame(survival_days = .gp$time),
                 "inst/extdata/complete_sample.csv", row.names = FALSE)
.ok("complete_sample.csv  (72 rows, no status column)")

## 2. Right-censored sample with a grouping covariate.
.tr <- .load1("transplant")
utils::write.csv(.tr, "inst/extdata/censored_sample.csv", row.names = FALSE)
.ok(sprintf("censored_sample.csv  (%d rows, %.0f%% censored)",
            nrow(.tr), 100 * mean(.tr$status == 0)))

## 3. Covariate sample for the AFT and cure paths.
.mel <- .load1("melanoma")
.mel_csv <- data.frame(time      = .mel$time,
                       event     = ifelse(.mel$status == 1, 1L, 0L),
                       age       = .mel$age,
                       thickness = .mel$thickness,
                       ulcer     = .mel$ulcer)
utils::write.csv(.mel_csv, "inst/extdata/covariate_sample.csv", row.names = FALSE)
.ok(sprintf("covariate_sample.csv  (%d rows, 3 covariates)", nrow(.mel_csv)))

## 4. Two-cause competing-risks sample, simulated with a fixed seed so the file
##    is reproducible from this script.
set.seed(20260730)
.n  <- 200
.t1 <- qbetadanish(stats::runif(.n), a = 1, b = 2.5, c = 1.8, k = 0.04)
.t2 <- qbetadanish(stats::runif(.n), a = 1, b = 3.5, c = 1.2, k = 0.02)
.cn <- stats::rexp(.n, rate = 1 / 150)
.tm <- pmin(.t1, .t2, .cn)
.cs <- ifelse(.tm == .cn, 0L, ifelse(.t1 < .t2, 1L, 2L))
utils::write.csv(data.frame(time = round(.tm, 3), cause = .cs),
                 "inst/extdata/competing_sample.csv", row.names = FALSE)
.ok(sprintf("competing_sample.csv  (%d rows; cause 0/1/2 = %d/%d/%d)",
            .n, sum(.cs == 0), sum(.cs == 1), sum(.cs == 2)))

## =============================================================================
##  REC 22  --  COMPLETE THE read_survival_data() HARDENING
## =============================================================================

.step("Rec 22: rewriting R/data_helpers.R with delimiter, encoding and cause support")

.put("R/data_helpers.R", r"---(## Candidate column names used when time_col or status_col is left NULL.
## Deliberately excludes "censor" and "cens": those names are used for both
## codings in practice, and guessing wrong silently inverts every event.
.BD_TIME_NAMES   <- c("time", "times", "t", "survtime", "surv_time",
                      "survival_time", "survival_days", "duration", "lifetime",
                      "futime", "days", "months", "years")
.BD_STATUS_NAMES <- c("status", "event", "died", "death", "delta", "survstatus",
                      "failure", "observed")
.BD_CAUSE_NAMES  <- c("cause", "causes", "failtype", "fail_type", "event_type",
                      "risk", "failcause")

#' Match a Column Name Case-Insensitively Against a Candidate List
#' @noRd
.bd_guess_col <- function(nms, candidates) {
  lower <- tolower(nms)
  for (cand in candidates) {
    hit <- which(lower == cand)
    if (length(hit)) return(nms[hit[1]])
  }
  NULL
}

#' Detect Whether Times Look Recorded on a Coarse Grid
#'
#' Returns the inferred recording increment, or `NA` if the times do not look
#' grid-recorded. Whole-month or whole-day recording matters because the
#' point-density likelihood is not appropriate for coarsely rounded times; the
#' grouped likelihood added in a later release is.
#'
#' @noRd
.bd_grid_step <- function(t) {
  t <- t[is.finite(t) & t > 0]
  if (length(t) < 5L) return(NA_real_)
  u <- sort(unique(t))
  if (length(u) < 3L) return(NA_real_)
  d <- diff(u); d <- d[d > 0]
  if (!length(d)) return(NA_real_)
  g <- min(d)
  if (!is.finite(g) || g <= 0) return(NA_real_)
  r <- t / g
  if (max(abs(r - round(r))) < 1e-8) g else NA_real_
}

#' Read and Prepare Survival Data
#'
#' Reads survival data from a delimited text file or an Excel workbook and
#' returns a clean data frame ready for [fit_betadanish()] and the regression
#' fitters. Columns are selected by name, covariates keep their original type,
#' and the result carries a report describing what was read.
#'
#' @param file Path to a `.csv`, `.txt`, `.tsv`, `.xls` or `.xlsx` file.
#' @param time_col Name of the time column. If `NULL` (default), the column is
#'   guessed from a list of common names and the choice is reported.
#' @param status_col Name of the event indicator (1 = event, 0 = censored). If
#'   `NULL`, guessed the same way; if no candidate is found, all observations
#'   are treated as uncensored.
#' @param covar_cols Character vector of covariate columns to retain, or `NULL`.
#' @param cause_col Name of a competing-risks cause column, or `NULL`. When
#'   supplied, code 0 means censored and positive integers index the causes.
#' @param sep Field separator. Defaults to `","`; use `"\t"` for
#'   tab-delimited files.
#' @param dec Decimal mark. Use `","` for European-format numerics.
#' @param encoding File encoding passed to [utils::read.table()], for example
#'   `"UTF-8"` or `"latin1"`.
#' @param na.strings Strings to treat as missing.
#' @param drop_na Logical; drop incomplete rows. Default `TRUE`.
#' @param quiet Logical; suppress the informational messages.
#'
#' @return A data frame with columns `time`, `status`, an optional `cause`, and
#'   any retained covariates. The `"bd_data_report"` attribute is a list
#'   recording the file, rows read and kept, number of events, censoring
#'   proportion, covariate names, and the inferred recording grid.
#'
#' @details
#' Column guessing is deliberately conservative. Names meaning "censoring
#' indicator" are excluded from the candidate list, because `censor = 1` means
#' *censored* in some conventions and *observed* in others, and guessing wrong
#' would silently invert every event. Name the column explicitly if in doubt.
#'
#' A covariate that happens to be called `time`, `status` or `cause` is renamed
#' with a `_cov` suffix and a warning rather than colliding with the response.
#'
#' Excel input requires the `readxl` package.
#'
#' @seealso [fit_betadanish()], [fit_bd_aft()], [fit_bd_competing()]
#'
#' @export
#'
#' @examples
#' # A complete uncensored sample with a non-standard column name
#' f <- system.file("extdata", "complete_sample.csv", package = "BetaDanish")
#' dat <- read_survival_data(f, time_col = "survival_days")
#' attr(dat, "bd_data_report")$rows_kept
#'
#' # Column names guessed automatically
#' f2 <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
#' dat2 <- read_survival_data(f2, quiet = TRUE)
#' mean(dat2$status == 0)
#'
#' # Competing risks
#' f3 <- system.file("extdata", "competing_sample.csv", package = "BetaDanish")
#' dat3 <- read_survival_data(f3, time_col = "time", cause_col = "cause",
#'                            quiet = TRUE)
#' table(dat3$cause)
read_survival_data <- function(file, time_col = NULL, status_col = NULL,
                               covar_cols = NULL, cause_col = NULL,
                               sep = ",", dec = ".", encoding = "unknown",
                               na.strings = c("NA", "", ".", "#N/A"),
                               drop_na = TRUE, quiet = FALSE) {

  if (!is.character(file) || length(file) != 1L)
    stop("'file' must be a single file path.", call. = FALSE)
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)

  say <- function(...) if (!isTRUE(quiet)) message(...)

  ## ---- read -----------------------------------------------------------------
  ext <- tolower(tools::file_ext(file))
  dat <- if (ext %in% c("csv", "txt", "tsv", "dat")) {
    if (ext == "tsv" && identical(sep, ",")) sep <- "\t"
    utils::read.table(file, header = TRUE, sep = sep, dec = dec,
                      na.strings = na.strings, encoding = encoding,
                      stringsAsFactors = FALSE, check.names = TRUE,
                      comment.char = "")
  } else if (ext %in% c("xls", "xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Reading Excel files requires the 'readxl' package.", call. = FALSE)
    as.data.frame(readxl::read_excel(file))
  } else {
    stop("Unsupported file extension '", ext, "'. Supply a .csv, .txt, .tsv, ",
         ".xls or .xlsx file.", call. = FALSE)
  }

  if (!nrow(dat)) stop("The file contains no data rows.", call. = FALSE)
  if (ncol(dat) == 1L && identical(sep, ","))
    warning("Only one column was read. If the file is not comma separated, ",
            "pass sep = \"\\t\" or sep = \";\".", call. = FALSE)

  nms <- names(dat)

  ## ---- resolve the time column ---------------------------------------------
  if (is.null(time_col)) {
    time_col <- .bd_guess_col(nms, .BD_TIME_NAMES)
    if (is.null(time_col))
      stop("Could not identify a time column. Available columns: ",
           paste(nms, collapse = ", "),
           ". Pass time_col explicitly.", call. = FALSE)
    say("Using '", time_col, "' as the time column.")
  }
  if (!time_col %in% nms)
    stop("Time column '", time_col, "' not found. Available columns: ",
         paste(nms, collapse = ", "), call. = FALSE)

  clean <- data.frame(time = suppressWarnings(as.numeric(dat[[time_col]])))
  if (all(is.na(clean$time)))
    stop("Column '", time_col, "' contains no numeric values. Check 'dec' if ",
         "the file uses a comma as the decimal mark.", call. = FALSE)

  ## ---- resolve the status column -------------------------------------------
  guessed_status <- FALSE
  if (is.null(status_col)) {
    status_col <- .bd_guess_col(nms, .BD_STATUS_NAMES)
    guessed_status <- !is.null(status_col)
  }
  if (is.null(status_col)) {
    say("No status column found; treating all observations as uncensored.")
    clean$status <- 1
  } else {
    if (!status_col %in% nms)
      stop("Status column '", status_col, "' not found.", call. = FALSE)
    if (guessed_status)
      say("Using '", status_col, "' as the event indicator ",
          "(1 = event, 0 = censored).")
    clean$status <- suppressWarnings(as.numeric(dat[[status_col]]))
  }

  ## ---- optional cause column ------------------------------------------------
  if (is.null(cause_col)) {
    cc <- .bd_guess_col(nms, .BD_CAUSE_NAMES)
    if (!is.null(cc) && !identical(cc, status_col)) {
      cause_col <- cc
      say("Using '", cause_col, "' as the competing-risks cause column.")
    }
  }
  if (!is.null(cause_col)) {
    if (!cause_col %in% nms)
      stop("Cause column '", cause_col, "' not found.", call. = FALSE)
    clean$cause <- suppressWarnings(as.integer(dat[[cause_col]]))
    if (anyNA(clean$cause) && !all(is.na(dat[[cause_col]])))
      warning("Some cause codes could not be read as integers.", call. = FALSE)
  }

  ## ---- covariates -----------------------------------------------------------
  if (!is.null(covar_cols)) {
    missing_cov <- setdiff(covar_cols, nms)
    if (length(missing_cov))
      stop("Covariate column(s) not found: ", paste(missing_cov, collapse = ", "),
           call. = FALSE)
    extra   <- dat[, covar_cols, drop = FALSE]
    reserved <- c("time", "status", "cause")
    clash   <- intersect(names(extra), reserved)
    if (length(clash)) {
      warning("Covariate(s) named ", paste(clash, collapse = ", "),
              " collide with the response; renamed with a '_cov' suffix.",
              call. = FALSE)
      names(extra)[names(extra) %in% clash] <-
        paste0(names(extra)[names(extra) %in% clash], "_cov")
    }
    clean <- cbind(clean, extra)
  }

  ## ---- validate -------------------------------------------------------------
  n_read <- nrow(clean)
  if (isTRUE(drop_na)) {
    clean <- clean[stats::complete.cases(clean), , drop = FALSE]
    if (nrow(clean) < n_read)
      warning("Dropped ", n_read - nrow(clean), " row(s) with missing values.",
              call. = FALSE)
  }
  n_kept <- nrow(clean)
  if (!n_kept)
    stop("No complete rows remain after removing missing values.", call. = FALSE)

  if (any(clean$time <= 0, na.rm = TRUE))
    warning("Some times are <= 0; survival models require strictly positive ",
            "times.", call. = FALSE)
  if (!all(clean$status %in% c(0, 1)))
    stop("The status column contains values other than 0 and 1. Recode so ",
         "that 1 = event and 0 = censored. If '", status_col, "' is a ",
         "censoring indicator, invert it first.", call. = FALSE)

  grid <- .bd_grid_step(clean$time)
  if (!is.na(grid))
    say("Times appear recorded on a grid of ", format(grid),
        ". For coarsely rounded data the point-density likelihood ",
        "understates uncertainty.")

  row.names(clean) <- NULL
  attr(clean, "bd_data_report") <- list(
    file           = normalizePath(file, mustWork = FALSE),
    time_col       = time_col,
    status_col     = status_col,
    cause_col      = cause_col,
    rows_read      = n_read,
    rows_kept      = n_kept,
    rows_dropped   = n_read - n_kept,
    n_events       = sum(clean$status == 1),
    censoring_prop = mean(clean$status == 0),
    covariates     = setdiff(names(clean), c("time", "status", "cause")),
    grid_step      = grid
  )
  clean
}
)---")

.step("Adding tests for the data layer")

.put("tests/testthat/test-data-layer.R", r"---(test_that("brain_cancer is gone", {
  expect_false("brain_cancer" %in% utils::data(package = "BetaDanish")$results[, "Item"])
  expect_error(get("brain_cancer", envir = asNamespace("BetaDanish")))
})

test_that("guinea_pig has the documented shape", {
  data(guinea_pig, package = "BetaDanish", envir = environment())
  expect_s3_class(guinea_pig, "data.frame")
  expect_equal(nrow(guinea_pig), 72L)
  expect_named(guinea_pig, c("time", "status"))
  expect_true(all(guinea_pig$status == 1))
  expect_true(all(guinea_pig$time > 0))
  expect_equal(min(guinea_pig$time), 12)
  expect_equal(max(guinea_pig$time), 376)
})

test_that("guinea_pig gives a well-identified four-parameter fit", {
  skip_on_cran()
  data(guinea_pig, package = "BetaDanish", envir = environment())
  set.seed(1)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = guinea_pig,
                   n_starts = 12, check_identifiability = FALSE))

  ## The thesis reports b-hat = 3.64 (SE 1.20), i.e. 2.2 SEs clear of the
  ## b = 1 ridge. Tolerances are loose: the point is that the optimum is
  ## interior and b is separated from one, not that it reproduces exactly.
  expect_true(all(fit$coefficients > 0))
  expect_gt(fit$coefficients[["b"]], 1.5)
  expect_true(is.finite(fit$AIC))
  expect_false(isTRUE(fit$diagnostics$near_b_ridge))
})

test_that("the example CSVs ship and read back", {
  for (f in c("complete_sample.csv", "censored_sample.csv",
              "covariate_sample.csv", "competing_sample.csv")) {
    p <- system.file("extdata", f, package = "BetaDanish")
    expect_true(nzchar(p), info = f)
    expect_gt(nrow(utils::read.csv(p)), 20L)
  }
})

test_that("column names are guessed when not supplied", {
  p   <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
  dat <- read_survival_data(p, quiet = TRUE)
  expect_named(dat, c("time", "status", "group"), ignore.order = TRUE)
  rep <- attr(dat, "bd_data_report")
  expect_equal(rep$time_col, "time")
  expect_equal(rep$status_col, "status")
})

test_that("an unguessable time column is a clear error", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(zzz = 1:5, qqq = 1L), tmp, row.names = FALSE)
  expect_error(read_survival_data(tmp, quiet = TRUE), "time column")
})

test_that("a competing-risks cause column is read", {
  p   <- system.file("extdata", "competing_sample.csv", package = "BetaDanish")
  dat <- read_survival_data(p, time_col = "time", cause_col = "cause",
                            quiet = TRUE)
  expect_true("cause" %in% names(dat))
  expect_true(all(dat$cause %in% 0:2))
})

test_that("semicolon-delimited European-format input works", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  writeLines(c("time;status", "1,5;1", "2,25;0", "3,75;1", "4,5;1", "5,25;1"),
             tmp)
  dat <- read_survival_data(tmp, sep = ";", dec = ",", quiet = TRUE)
  expect_equal(dat$time, c(1.5, 2.25, 3.75, 4.5, 5.25))
  expect_equal(dat$status, c(1, 0, 1, 1, 1))
})

test_that("a censoring indicator is rejected rather than guessed", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(time = c(1, 2, 3, 4, 5), censor = c(0, 1, 0, 0, 1)),
                   tmp, row.names = FALSE)
  ## "censor" is not in the candidate list, so status is not guessed from it.
  dat <- read_survival_data(tmp, quiet = TRUE)
  expect_true(all(dat$status == 1))
})

test_that("grid-recorded times are detected and reported", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(time = c(3, 6, 9, 12, 18, 24, 36),
                              status = c(1, 1, 0, 1, 1, 0, 1)),
                   tmp, row.names = FALSE)
  dat <- read_survival_data(tmp, quiet = TRUE)
  expect_equal(attr(dat, "bd_data_report")$grid_step, 3)
})
)---")

## =============================================================================
##  NEWS
## =============================================================================

.step("Recording the change in NEWS.md")

.nw <- readLines("NEWS.md", warn = FALSE)
if (!any(grepl("brain_cancer", .nw, fixed = TRUE) &
         grepl("removed", .nw, fixed = TRUE))) {
  .hdr <- grep("^# BetaDanish 0\\.2\\.0\\.9000", .nw)
  if (length(.hdr) == 1L) {
    .backup("NEWS.md")
    .sec <- c(
      "",
      "## Datasets",
      "",
      "* **`brain_cancer` has been removed.** The dataset of 500 brain cancer",
      "  patients that shipped in 0.1.0 and 0.2.0 is no longer included, and",
      "  the AFT example in the README now uses `melanoma`. The entry under",
      "  0.1.0 below is left as written, since it records what that release",
      "  actually contained.",
      "",
      "* **`guinea_pig` added** (Bjerkedal 1960, n = 72): survival times in days",
      "  of guinea pigs injected with virulent tubercle bacilli. This is the",
      "  reference case for a well-identified four-parameter fit, with",
      "  \\eqn{\\hat b = 3.64} (SE 1.20) sitting about 2.2 standard errors clear",
      "  of the \\eqn{b = 1} ridge. On `remission` and `carbon_fibres` the",
      "  parent model sits on the flat \\eqn{(a, c)} direction instead, so",
      "  without this dataset no built-in example showed the full model",
      "  behaving well.",
      "",
      "* **`inst/extdata/` added** with four example CSVs -- complete,",
      "  censored, covariate and competing-risks -- so that examples and",
      "  vignettes can demonstrate the file-driven workflow.",
      "",
      "## Data input",
      "",
      "* `read_survival_data()` gains `cause_col`, `sep`, `dec`, `encoding`,",
      "  `na.strings`, `drop_na` and `quiet`, and will guess the time and",
      "  status columns when they are not named. Guessing deliberately excludes",
      "  names such as `censor`, whose polarity differs between conventions.",
      "",
      "* The `bd_data_report` attribute now also records which columns were",
      "  used and whether the times look recorded on a coarse grid.")
    .nw <- append(.nw, .sec, after = .hdr)
    .write_lines("NEWS.md", .nw)
    .ok("NEWS.md updated")
  } else {
    .warn("could not find the 0.2.0.9000 header in NEWS.md; add the note by hand")
  }
} else {
  .info("NEWS.md already records the removal")
}

.step("Build-ignoring data-raw/")

.rbi <- readLines(".Rbuildignore", warn = FALSE)
if (!"^data-raw$" %in% .rbi) {
  .backup(".Rbuildignore")
  .write_lines(".Rbuildignore", c(.rbi[nzchar(.rbi)], "^data-raw$"))
  .ok("added ^data-raw$")
} else {
  .info("already build-ignored")
}

## =============================================================================
##  VERIFY
## =============================================================================

.step("Confirming brain_cancer is gone everywhere")

.residual <- character(0)
for (f in c(list.files("R", pattern = "[.]R$", full.names = TRUE),
            list.files("man", pattern = "[.]Rd$", full.names = TRUE),
            list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE),
            list.files("tests", pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
            Filter(file.exists, c("README.md", "DESCRIPTION", "NAMESPACE",
                                  "inst/WORDLIST", "inst/CITATION")))) {
  if (any(grepl("brain_cancer", readLines(f, warn = FALSE), fixed = TRUE)))
    .residual <- c(.residual, f)
}
## NEWS.md legitimately mentions it: the 0.1.0 history and the removal note.
if (length(.residual)) {
  .warn(paste("brain_cancer still referenced in:", paste(.residual, collapse = ", ")))
} else {
  .ok("no references outside NEWS.md, where the history is kept deliberately")
}
if (file.exists("data/brain_cancer.rda")) .die("data/brain_cancer.rda still present.")
.ok("data/brain_cancer.rda absent")

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
cat("  PATCH 2 OF 3 COMPLETE  --  data layer\n")
cat(strrep("=", 78), "\n\n")
cat("  19  brain_cancer removed: data file, roxygen block, man page,\n")
cat("      README example and table row, WORDLIST vocabulary\n")
cat("  20  guinea_pig added (Bjerkedal 1960, n = 72) with data-raw/ source\n")
cat("  21  inst/extdata/ with four example CSVs\n")
cat("  22  read_survival_data(): cause_col, sep, dec, encoding, na.strings,\n")
cat("      quiet, column guessing, grid-recording detection\n\n")
cat("  Left for Patch 2b\n")
cat("    23  bd_analyze_csv() end-to-end pipeline\n")
cat("    24  bd_csv_template()\n")
cat("    46  README and vignette rewrite around the CSV workflow\n\n")
cat("  Backups: ", BACKUP_DIR, "\n", sep = "")
cat("  brain_cancer.rda is recoverable from there if needed.\n\n")
cat("  Target state: 0 errors, 0 warnings, 1 note (the clock note).\n\n")
