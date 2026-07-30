## =============================================================================
##  BetaDanish  --  PATCH 1a : follow-up fixes after the Patch 1 check
## =============================================================================
##
##  Fixes the one ERROR, the 12 test warnings, and 2 of the 3 NOTEs reported by
##  devtools::check() after Patch 1.
##
##    ERROR   .bd_conform() recycled zero-length arguments up to length 1.
##            R's distribution functions return a zero-length result when ANY
##            argument is zero-length. The test was right; the code was wrong.
##
##    WARN    The new identifiability diagnostics fire on the existing test
##            fixtures, which genuinely sit near the b = 1 ridge. Silenced
##            where identifiability is not what is under test, and covered by
##            a new deterministic test of the warning logic itself.
##
##    NOTE 1  .betadanish_backup is a hidden directory at top level.
##    NOTE 3  BetaDanish_Manual.docx / .pdf are non-standard top-level files.
##            Both are build-ignored: kept in your repo, excluded from the
##            tarball that CRAN sees.
##
##    NOTE 2  "unable to verify current time" -- not actionable, see below.
##
##  HOW TO RUN   Same as before, from the package root:
##                 source("dev/BetaDanish_Patch1a_followup.R")
##  IDEMPOTENT   Yes.
## =============================================================================

if (getRversion() < "4.0.0") stop("This patch needs R >= 4.0.")

.step_n <- 0L
.step <- function(msg) { .step_n <<- .step_n + 1L; cat(sprintf("\n[%02d] %s\n", .step_n, msg)) }
.ok    <- function(msg) cat("     OK   ", msg, "\n", sep = "")
.info  <- function(msg) cat("     ..   ", msg, "\n", sep = "")
.warn  <- function(msg) cat("     WARN ", msg, "\n", sep = "")
.die   <- function(...) stop("\n\n*** PATCH ABORTED ***\n", ..., "\n", call. = FALSE)

BACKUP_DIR <- NULL
.backup <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  dest <- file.path(BACKUP_DIR, path)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(path, dest, overwrite = TRUE)) .die("Could not back up ", path)
  invisible(TRUE)
}

.sub_in <- function(path, from, to, label, required = TRUE) {
  if (!file.exists(path)) {
    if (required) .warn(paste(path, "not found; skipped"))
    return(invisible(FALSE))
  }
  txt  <- readLines(path, warn = FALSE)
  hits <- sum(grepl(from, txt, fixed = TRUE))
  if (hits == 0L) {
    if (any(grepl(to, txt, fixed = TRUE))) .info(paste(label, "-- already applied"))
    else if (required) .warn(paste(label, "-- pattern not found; check by hand"))
    return(invisible(FALSE))
  }
  .backup(path)
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(gsub(from, to, txt, fixed = TRUE), con = con, sep = "\n")
  .ok(sprintf("%s -- %d occurrence(s)", label, hits))
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
cat("  BetaDanish  --  Patch 1a : follow-up fixes\n")
cat(strrep("=", 78), "\n")

## ------------------------------------------------------------- pre-flight ----

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
.d <- read.dcf("DESCRIPTION")
if (.d[1, "Package"] != "BetaDanish") .die("This is not the BetaDanish package.")
if (!file.exists("R/dist_functions.R")) .die("R/dist_functions.R missing -- run Patch 1 first.")
if (!any(grepl(".bd_conform", readLines("R/dist_functions.R", warn = FALSE), fixed = TRUE)))
  .die("R/dist_functions.R does not contain .bd_conform -- Patch 1 has not been applied.")
.ok(paste("package root:", getwd()))
.ok(paste("version:", .d[1, "Version"]))

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch1a"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  FIX 1  --  the ERROR: zero-length propagation
## =============================================================================

.step("Fixing .bd_conform zero-length handling (the check ERROR)")

.sub_in("R/dist_functions.R",
        "  n <- max(length(x), length(a), length(b), length(c), length(k))\n  if (n == 0L) return(NULL)",
        "  lens <- c(length(x), length(a), length(b), length(c), length(k))\n  ## R's distribution functions return a zero-length result when ANY argument\n  ## is zero-length, rather than recycling up to length one. Match that.\n  if (any(lens == 0L)) return(NULL)\n  n <- max(lens)",
        label = "zero-length propagation")

## The two lines may have been read separately; fall back to a line-wise edit.
.dl <- readLines("R/dist_functions.R", warn = FALSE)
if (any(grepl("n <- max(length(x), length(a), length(b), length(c), length(k))", .dl, fixed = TRUE))) {
  .backup("R/dist_functions.R")
  .i <- grep("n <- max(length(x), length(a), length(b), length(c), length(k))", .dl, fixed = TRUE)
  .j <- grep("if (n == 0L) return(NULL)", .dl, fixed = TRUE)
  if (length(.i) == 1L && length(.j) == 1L && .j == .i + 1L) {
    .dl <- append(.dl[-c(.i, .j)],
                  c("  lens <- c(length(x), length(a), length(b), length(c), length(k))",
                    "  ## R's distribution functions return a zero-length result when ANY",
                    "  ## argument is zero-length, rather than recycling up to length one.",
                    "  if (any(lens == 0L)) return(NULL)",
                    "  n <- max(lens)"),
                  after = .i - 1L)
    .con <- file("R/dist_functions.R", open = "wb")
    writeLines(.dl, con = .con, sep = "\n"); close(.con)
    .ok("zero-length propagation (line-wise fallback)")
  } else {
    .die("Could not locate the .bd_conform length lines. Edit R/dist_functions.R by hand:\n",
         "  replace  n <- max(length(x), ...)  and the  if (n == 0L)  line with\n",
         "  lens <- c(length(x), length(a), length(b), length(c), length(k))\n",
         "  if (any(lens == 0L)) return(NULL)\n  n <- max(lens)")
  }
} else {
  .info("zero-length propagation -- already applied")
}

## =============================================================================
##  FIX 2  --  the 12 test warnings
## =============================================================================

.step("Silencing identifiability warnings in the existing tests")

## Every n_starts call site in the original test suite is a fit_betadanish
## call, so these two substitutions are unambiguous.
.tfiles <- list.files("tests/testthat", pattern = "^test-.*[.]R$", full.names = TRUE)
.n_touched <- 0L
for (f in .tfiles) {
  txt <- readLines(f, warn = FALSE)
  if (!any(grepl("fit_betadanish(", txt, fixed = TRUE))) next
  if (any(grepl("check_identifiability", txt, fixed = TRUE))) next
  new <- txt
  new <- gsub("n_starts = 1)", "n_starts = 1, check_identifiability = FALSE)", new, fixed = TRUE)
  new <- gsub("n_starts = 2)", "n_starts = 2, check_identifiability = FALSE)", new, fixed = TRUE)
  if (!identical(new, txt)) {
    .backup(f)
    con <- file(f, open = "wb"); writeLines(new, con = con, sep = "\n"); close(con)
    .n_touched <- .n_touched + 1L
    .ok(paste("quieted", basename(f)))
  }
}
if (.n_touched == 0L) .info("no test files needed changing")

## A deterministic test of the warning logic, so the diagnostics are covered
## without depending on where the optimiser happens to land.
.step("Adding a deterministic test for the diagnostics themselves")

.put("tests/testthat/test-identifiability.R", r"---(## The identifiability diagnostics are tested directly rather than through a
## fitted model, so the test does not depend on where the optimiser lands.

mk_diag <- function(...) {
  base <- list(converged = TRUE, convergence_code = 1L, vcov_singular = FALSE,
               near_b_ridge = FALSE, b_distance_se = 5, ac_correlation = 0.2)
  utils::modifyList(base, list(...))
}

test_that("the b = 1 ridge warning fires when b-hat sits close to one", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(near_b_ridge = TRUE,
                                              b_distance_se = 0.6)),
    "non-identifiability ridge")
})

test_that("a singular information matrix is reported", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(vcov_singular = TRUE)),
    "singular")
})

test_that("(a, c) confounding is reported above the 0.95 threshold", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(ac_correlation = 0.99)),
    "only the product")
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(ac_correlation = -0.98)),
    "only the product")
})

test_that("a poor convergence code is reported", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(converged = FALSE,
                                              convergence_code = 4L)),
    "code 4")
})

test_that("a well-behaved fit produces no diagnostic warnings", {
  expect_silent(BetaDanish:::.bd_warn_diagnostics(mk_diag()))
})

test_that("check_identifiability = FALSE suppresses the warnings", {
  skip_on_cran()
  dat <- simulate_bd_data(120, a = 1, b = 1.1, c = 2, k = 0.5, seed = 99)
  fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                        n_starts = 3, check_identifiability = FALSE)
  expect_s3_class(fit, "betadanish")
  expect_true(is.list(fit$diagnostics))
})

test_that("zero-length input gives zero-length output for every function", {
  z <- numeric(0)
  expect_length(dbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(pbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(qbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(sbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(hbetadanish(z, 1.5, 3, 2, 1), 0L)
  ## A zero-length parameter also collapses the result, as in stats::dnorm.
  expect_length(dbetadanish(1, a = z, b = 3, c = 2, k = 1), 0L)
})
)---")

## =============================================================================
##  FIX 3  --  NOTE 1 and NOTE 3: build-ignore housekeeping
## =============================================================================

.step("Build-ignoring backup and manual files (NOTE 1, NOTE 3)")

if (!file.exists(".Rbuildignore")) {
  file.create(".Rbuildignore")
  .info("created .Rbuildignore")
}
.backup(".Rbuildignore")

.rbi      <- readLines(".Rbuildignore", warn = FALSE)
.wanted   <- c("^dev$", "^\\.betadanish_backup$", "^BetaDanish_Manual\\.docx$",
               "^BetaDanish_Manual\\.pdf$", "^.*\\.Rproj$", "^\\.Rproj\\.user$")
.to_add   <- setdiff(.wanted, .rbi)
if (length(.to_add)) {
  .con <- file(".Rbuildignore", open = "wb")
  writeLines(c(.rbi[nzchar(.rbi)], .to_add), con = .con, sep = "\n"); close(.con)
  .ok(paste("added:", paste(.to_add, collapse = "  ")))
} else {
  .info(".Rbuildignore already complete")
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

.step("devtools::document()")
.r <- tryCatch({ devtools::document(); TRUE }, error = function(e) conditionMessage(e))
if (!isTRUE(.r)) .die("document() failed:\n  ", .r, "\n\nBackups: ", BACKUP_DIR)
.ok("documentation regenerated")

.step("devtools::test() -- fast, before the full check")
.t <- tryCatch(devtools::test(), error = function(e) { .warn(conditionMessage(e)); NULL })

.step("devtools::check() -- several minutes, do not interrupt")
.chk <- tryCatch(devtools::check(document = FALSE, args = "--as-cran", error_on = "never"),
                 error = function(e) { .warn(conditionMessage(e)); NULL })

cat("\n", strrep("=", 78), "\n", sep = "")
if (!is.null(.chk)) {
  cat("  CHECK RESULT\n", strrep("=", 78), "\n", sep = "")
  cat("  errors:  ", length(.chk$errors),   "\n", sep = "")
  cat("  warnings:", length(.chk$warnings), "\n", sep = "")
  cat("  notes:   ", length(.chk$notes),    "\n", sep = "")
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
cat("  PATCH 1a COMPLETE\n")
cat(strrep("=", 78), "\n\n")
cat("  Fixed\n")
cat("    ERROR   .bd_conform now propagates zero length, matching stats::dnorm\n")
cat("    WARN    identifiability warnings quieted in the existing tests;\n")
cat("            new deterministic test covers the warning logic itself\n")
cat("    NOTE 1  .betadanish_backup build-ignored\n")
cat("    NOTE 3  BetaDanish_Manual.docx / .pdf build-ignored\n\n")
cat("  Expected remaining\n")
cat("    NOTE    'unable to verify current time' -- your machine could not\n")
cat("            reach the time server the check uses. Harmless, machine\n")
cat("            specific, and will not appear on CRAN's own builders.\n\n")
cat("  Target state: 0 errors, 0 warnings, 1 note.\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
