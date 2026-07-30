## =============================================================================
##  BetaDanish  --  PATCH 2a : follow-up fixes after the Patch 2 check
## =============================================================================
##
##  WARNING x2  man/read_survival_data.Rd:36 "unknown macro '\t'".
##              The roxygen @param for `sep` contained a literal backslash-t.
##              Backslash is an escape character in Rd, so it must be written
##              \\t in the roxygen comment. One root cause, reported by both
##              the install step and the Rd check.
##
##  ERROR       test-data-layer.R:46 expected the `group` column to survive
##              read_survival_data(). It did not: columns not named in
##              covar_cols are dropped. The code behaved as designed and the
##              test was wrong -- but the design was also wrong. Silently
##              discarding a user's covariates is poor behaviour for a
##              file-driven tool, so:
##                * covar_cols = "all" now retains every non-response column
##                * columns that are dropped are named in a message
##                * the report records available and dropped columns
##
##  Also fixed  The "Only one column was read" warning fired on any
##              single-column file, including the shipped complete_sample.csv
##              and its own @examples. It now fires only when the header
##              actually contains a tab, semicolon or pipe.
##
##  Not a problem, no action needed:
##    * "unable to verify current time"  -- your clock/network, not the package
##    * "Removed empty directory .../_snaps"  -- routine build tidying
##    * 'Unknown command "TMPDIR=..."' and the quarto warning at the very end
##      -- a devtools/quarto interaction on Windows, after the check finished.
##      It has no bearing on the check result.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch2a_followup.R")
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

cat(strrep("=", 78), "\n")
cat("  BetaDanish  --  Patch 2a : follow-up fixes\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("data/guinea_pig.rda")) .die("Patch 2 has not been applied.")
.ok("Patch 2 detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch2a"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  R/data_helpers.R  --  all three fixes
## =============================================================================

.step("Rewriting R/data_helpers.R (Rd macro, covariate retention, separator warning)")

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
#' point-density likelihood is not appropriate for coarsely rounded times.
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

#' Warn Only When the Header Really Suggests a Different Separator
#'
#' A single-column result is perfectly normal for a file of bare times. It is
#' only suspicious when the header line also contains a tab, semicolon or pipe.
#'
#' @noRd
.bd_check_separator <- function(file, n_col, sep) {
  if (n_col > 1L || !identical(sep, ",")) return(invisible(NULL))
  first <- tryCatch(readLines(file, n = 1L, warn = FALSE),
                    error = function(e) character(0))
  if (!length(first)) return(invisible(NULL))
  if (grepl("[\t;|]", first))
    warning("Only one column was read, but the header contains another ",
            "delimiter. If the file is not comma separated, pass sep = ",
            "'\\t', sep = ';' or sep = '|'.", call. = FALSE)
  invisible(NULL)
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
#' @param covar_cols Covariate columns to retain. Either a character vector of
#'   column names, the string `"all"` to retain every column that is not part
#'   of the response, or `NULL` (default) to retain none. Whatever is dropped
#'   is named in a message unless `quiet = TRUE`.
#' @param cause_col Name of a competing-risks cause column, or `NULL`. When
#'   supplied, code 0 means censored and positive integers index the causes.
#' @param sep Field separator. Defaults to a comma. Use `"\\t"` for
#'   tab-delimited input, `";"` for semicolon-delimited.
#' @param dec Decimal mark. Use `","` for European-format numerics.
#' @param encoding File encoding passed to [utils::read.table()], for example
#'   `"UTF-8"` or `"latin1"`.
#' @param na.strings Strings to treat as missing.
#' @param drop_na Logical; drop incomplete rows. Default `TRUE`.
#' @param quiet Logical; suppress the informational messages.
#'
#' @return A data frame with columns `time`, `status`, an optional `cause`, and
#'   any retained covariates. The `"bd_data_report"` attribute is a list
#'   recording the file, which columns were used, rows read and kept, number of
#'   events, censoring proportion, retained and dropped columns, and the
#'   inferred recording grid.
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
#' dat <- read_survival_data(f, time_col = "survival_days", quiet = TRUE)
#' attr(dat, "bd_data_report")$rows_kept
#'
#' # Column names guessed, and every covariate retained
#' f2 <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
#' dat2 <- read_survival_data(f2, covar_cols = "all", quiet = TRUE)
#' names(dat2)
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
  .bd_check_separator(file, ncol(dat), sep)

  nms <- names(dat)

  ## ---- resolve the time column ---------------------------------------------
  if (is.null(time_col)) {
    time_col <- .bd_guess_col(nms, .BD_TIME_NAMES)
    if (is.null(time_col))
      stop("Could not identify a time column. Available columns: ",
           paste(nms, collapse = ", "), ". Pass time_col explicitly.",
           call. = FALSE)
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
    status_col     <- .bd_guess_col(nms, .BD_STATUS_NAMES)
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
  response_cols <- c(time_col, status_col, cause_col)

  if (identical(covar_cols, "all")) {
    covar_cols <- setdiff(nms, response_cols)
    if (length(covar_cols))
      say("Retaining ", length(covar_cols), " covariate(s): ",
          paste(covar_cols, collapse = ", "), ".")
  }

  if (length(covar_cols)) {
    missing_cov <- setdiff(covar_cols, nms)
    if (length(missing_cov))
      stop("Covariate column(s) not found: ", paste(missing_cov, collapse = ", "),
           call. = FALSE)
    extra    <- dat[, covar_cols, drop = FALSE]
    reserved <- c("time", "status", "cause")
    clash    <- intersect(names(extra), reserved)
    if (length(clash)) {
      warning("Covariate(s) named ", paste(clash, collapse = ", "),
              " collide with the response; renamed with a '_cov' suffix.",
              call. = FALSE)
      names(extra)[names(extra) %in% clash] <-
        paste0(names(extra)[names(extra) %in% clash], "_cov")
    }
    clean <- cbind(clean, extra)
  }

  dropped_cols <- setdiff(nms, c(response_cols, covar_cols))
  if (length(dropped_cols))
    say("Not retained: ", paste(dropped_cols, collapse = ", "),
        ". Pass covar_cols to keep specific columns, or covar_cols = \"all\".")

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
    file              = normalizePath(file, mustWork = FALSE),
    time_col          = time_col,
    status_col        = status_col,
    cause_col         = cause_col,
    available_columns = nms,
    covariates        = setdiff(names(clean), c("time", "status", "cause")),
    dropped_columns   = dropped_cols,
    rows_read         = n_read,
    rows_kept         = n_kept,
    rows_dropped      = n_read - n_kept,
    n_events          = sum(clean$status == 1),
    censoring_prop    = mean(clean$status == 0),
    grid_step         = grid
  )
  clean
}
)---")

## =============================================================================
##  tests/testthat/test-data-layer.R  --  corrected and extended
## =============================================================================

.step("Rewriting tests/testthat/test-data-layer.R")

.put("tests/testthat/test-data-layer.R", r"---(test_that("brain_cancer is gone", {
  items <- utils::data(package = "BetaDanish")$results[, "Item"]
  expect_false("brain_cancer" %in% items)
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
  ## interior and b is separated from one, not exact reproduction.
  expect_true(all(fit$coefficients > 0))
  expect_gt(fit$coefficients[["b"]], 1.5)
  expect_true(is.finite(fit$AIC))
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
  ## Without covar_cols only the response is retained.
  expect_named(dat, c("time", "status"))
  rep <- attr(dat, "bd_data_report")
  expect_equal(rep$time_col, "time")
  expect_equal(rep$status_col, "status")
  expect_true("group" %in% rep$dropped_columns)
})

test_that("covar_cols = 'all' retains every non-response column", {
  p   <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
  dat <- read_survival_data(p, covar_cols = "all", quiet = TRUE)
  expect_true("group" %in% names(dat))
  expect_length(attr(dat, "bd_data_report")$dropped_columns, 0L)

  p2   <- system.file("extdata", "covariate_sample.csv", package = "BetaDanish")
  dat2 <- read_survival_data(p2, covar_cols = "all", quiet = TRUE)
  expect_true(all(c("age", "thickness", "ulcer") %in% names(dat2)))
})

test_that("named covariates are retained and the rest reported as dropped", {
  p   <- system.file("extdata", "covariate_sample.csv", package = "BetaDanish")
  dat <- read_survival_data(p, covar_cols = "age", quiet = TRUE)
  expect_named(dat, c("time", "status", "age"))
  expect_setequal(attr(dat, "bd_data_report")$dropped_columns,
                  c("thickness", "ulcer"))
})

test_that("dropped columns are announced unless quiet", {
  p <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")

  ## Several messages are emitted, so collect them all rather than relying on
  ## expect_message matching whichever comes first.
  msgs <- character(0)
  withCallingHandlers(
    read_survival_data(p),
    message = function(cond) {
      msgs <<- c(msgs, conditionMessage(cond))
      invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("Not retained", msgs)))
  expect_true(any(grepl("group", msgs)))

  expect_silent(read_survival_data(p, quiet = TRUE))
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

test_that("a genuinely single-column file does not trigger a separator warning", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(time = c(1, 2, 3, 4, 5)), tmp, row.names = FALSE)
  expect_warning(read_survival_data(tmp, quiet = TRUE), regexp = NA)
})

test_that("the separator heuristic flags only a mis-delimited header", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)

  ## Single column, but the header betrays a semicolon: flag it.
  writeLines(c("time;status", "1;1"), tmp)
  expect_warning(BetaDanish:::.bd_check_separator(tmp, 1L, ","),
                 "another delimiter")

  ## Single column with a clean header: perfectly normal, stay quiet.
  writeLines(c("time", "1", "2"), tmp)
  expect_warning(BetaDanish:::.bd_check_separator(tmp, 1L, ","), regexp = NA)

  ## Several columns were parsed, so the separator was clearly right.
  writeLines(c("time;status", "1;1"), tmp)
  expect_warning(BetaDanish:::.bd_check_separator(tmp, 3L, ","), regexp = NA)
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

.step("Recording covar_cols = 'all' in NEWS.md")

.nw <- readLines("NEWS.md", warn = FALSE)
if (!any(grepl('covar_cols = `"all"`', .nw, fixed = TRUE)) &&
    !any(grepl('covar_cols = "all"', .nw, fixed = TRUE))) {
  .anchor <- grep("^\\* The `bd_data_report` attribute now also records", .nw)
  if (length(.anchor) == 1L) {
    .backup("NEWS.md")
    .add <- c(
      "",
      "* `read_survival_data(covar_cols = \"all\")` retains every column that is",
      "  not part of the response. Columns that are dropped are now named in a",
      "  message and recorded in the report, rather than disappearing silently.")
    .nw <- append(.nw, .add, after = .anchor)
    con <- file("NEWS.md", open = "wb"); writeLines(.nw, con = con, sep = "\n"); close(con)
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

.step("Checking the roxygen source for unescaped Rd macros")

.src <- readLines("R/data_helpers.R", warn = FALSE)
.rox <- grep("^#'", .src, value = TRUE)
.bad_macro <- grep("(^|[^\\\\])\\\\[a-zA-Z]", .rox, value = TRUE)
.bad_macro <- .bad_macro[!grepl("\\\\\\\\", .bad_macro)]
.known <- "\\\\(describe|item|eqn|deqn|doi|donttest|dontrun|code|link|url|emph|strong|href|enumerate|itemize|preformatted|verb|Sexpr)"
.bad_macro <- .bad_macro[!grepl(.known, .bad_macro)]
if (length(.bad_macro)) {
  .warn("possible unescaped backslash in roxygen:")
  for (l in .bad_macro) cat("        ", l, "\n", sep = "")
} else {
  .ok("no unescaped Rd macros in R/data_helpers.R")
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

.step("Confirming man/read_survival_data.Rd is clean")
if (file.exists("man/read_survival_data.Rd")) {
  .rdl <- readLines("man/read_survival_data.Rd", warn = FALSE)
  .hit <- grep("(^|[^\\\\])\\\\t([^a-zA-Z]|$)", .rdl)
  if (length(.hit)) {
    .warn(sprintf("a bare \\t remains on Rd line(s): %s",
                  paste(.hit, collapse = ", ")))
  } else {
    .ok("no bare backslash-t in the generated Rd")
  }
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
cat("  PATCH 2a COMPLETE\n")
cat(strrep("=", 78), "\n\n")
cat("  Fixed\n")
cat("    WARNING x2  \\t escaped as \\\\t in the sep documentation\n")
cat("    ERROR       test corrected; covar_cols = \"all\" added; dropped\n")
cat("                columns now reported instead of vanishing silently\n")
cat("    (extra)     separator warning no longer fires on single-column files\n\n")
cat("  Target state: 0 errors, 0 warnings, 1 note (the clock note).\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
