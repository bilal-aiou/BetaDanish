## =============================================================================
##  BetaDanish  --  PATCH 3c-meta : PUBLIC-FACING DESCRIPTION AND CITATIONS
## =============================================================================
##
##  The DESCRIPTION text on the CRAN page and the pkgdown site is the first
##  thing anyone reads, and it still describes 0.2.0. Since then the package has
##  gained the whole structural-properties layer, the CSV pipeline, penalised
##  and grouped estimation, profile and Wald inference, the identified
##  reparametrisation, the named ED interface, and the degeneracy guard. None of
##  that is mentioned; Bayesian inference was never mentioned even in 0.2.0.
##
##  WHAT THIS PATCH CHANGES
##    DESCRIPTION   Description field rewritten to describe what the package
##                  actually does. Title, authors, URLs and licence untouched --
##                  they are accurate, and changing a published Title invites
##                  avoidable review friction.
##    inst/CITATION Hard-coded "R package version 0.2.0" replaced with
##                  meta$Version, so the citation can never again disagree with
##                  the installed package.
##    README.md     Overview and feature list rewritten to match the code.
##
##  ONE THING FOR YOU TO CHECK, NOT FOR ME TO GUESS
##    The package cites the JAMSI article as volume 21, number 1. Volume 21 is
##    confirmed as 2025. But the neighbouring DOI 10.2478/jamsi-2025-0006 is
##    published in Volume 21, ISSUE 2 (December 2025), pages 5-19, and JAMSI
##    numbers its DOIs sequentially across a volume -- so a higher number
##    landing in issue 1 would be unusual.
##
##    I have NOT changed it. You are the author and will know from your offprint
##    or acceptance letter. The patch prints the current entry and asks. If it
##    should be 21(2), edit inst/CITATION and DESCRIPTION and re-run
##    devtools::document().
##
##    Page numbers are also absent from the CITATION entry. If the article has
##    them, adding `pages = "..."` completes the record.
##
##  A NOTE ON THE 0.2.0 README
##    Its overview already claimed mean residual life, hazard-shape
##    classification and stress-strength reliability. At 0.2.0 none of those
##    existed -- the same class of unsupported claim the Phase 1 audit found in
##    NEWS.md, which I flagged there and missed here. Patch 3a made every one of
##    them true, so the list below is now accurate rather than aspirational.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3c_meta.R")
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
.write_lines <- function(path, lines) {
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(lines, con = con, sep = "\n")
}

cat(strrep("=", 78), "\n")
cat("  BetaDanish  --  Patch 3c-meta : public description and citations\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("R/inference.R")) .die("Patch 3b has not been applied.")
.ok("Patch 3b detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3cmeta"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  DESCRIPTION
## =============================================================================

.step("Rewriting the Description field")

.d <- readLines("DESCRIPTION", warn = FALSE)
.backup("DESCRIPTION")

## Replace the Description block: the field line plus its indented continuations.
.i0 <- grep("^Description:", .d)
if (length(.i0) != 1L) .die("Could not locate a single Description: field.")
.i1 <- .i0
while (.i1 < length(.d) && grepl("^[ \t]", .d[.i1 + 1L])) .i1 <- .i1 + 1L

.new_desc <- c(
"Description: Implements the four-parameter Beta-Danish distribution and its",
"    three-parameter Exponentiated Danish submodel for survival, reliability",
"    and lifetime data analysis, following Ahmad and Danish (2025)",
"    <doi:10.2478/jamsi-2025-0010>. Density, distribution, quantile, survival,",
"    hazard and random generation functions are evaluated so as to retain",
"    accuracy in the heavy upper tail, where the survival function is regularly",
"    varying. Estimation covers maximum likelihood for complete and",
"    right-censored samples, ridge-penalized fitting for weakly identified",
"    regimes, a grouped likelihood for times recorded on a coarse grid, and",
"    Bayesian sampling. Inference provides log-scale Wald and profile",
"    likelihood intervals, together with a reparameterization in terms of the",
"    identified composite of the two shape parameters. Structural properties",
"    include raw, incomplete and conditional moments with their existence",
"    conditions, Shannon, Renyi and Tsallis entropies, mean residual life,",
"    mean deviations, Lorenz and Bonferroni curves, probability weighted",
"    moments, order statistics, stress-strength reliability, hazard shape",
"    classification and the tail index. Regression modules cover accelerated",
"    failure time models, mixture and promotion-time cure models, and competing",
"    risks with Aalen-Johansen comparison and Gray's test. Analyses can be run",
"    directly from a delimited text file or spreadsheet.")

.d <- append(.d[-(.i0:.i1)], .new_desc, after = .i0 - 1L)
.write_lines("DESCRIPTION", .d)
.ok(sprintf("Description rewritten (%d lines)", length(.new_desc)))

.step("Checking spelling conventions against the declared Language")
.lang <- tryCatch(read.dcf("DESCRIPTION")[1, "Language"], error = function(e) NA)
if (!is.na(.lang)) {
  .ok(paste("DESCRIPTION declares Language:", .lang))
  if (identical(.lang, "en-US")) {
    .info("The Description field uses US spelling to match. Note that the")
    .info("roxygen documentation added in Patches 1-3b is predominantly")
    .info("British ('penalised', 'normalised', 'parametrisation'). If you")
    .info("ever enable spelling::spell_check_test(), either switch Language")
    .info("to en-GB or add those forms to inst/WORDLIST.")
  }
}

.step("Checking the Description for non-ASCII characters")
.desc_txt <- paste(.new_desc, collapse = " ")
.bad_chars <- unique(unlist(strsplit(.desc_txt, ""))[
  grepl("[^ -~]", unlist(strsplit(.desc_txt, "")))])
if (length(.bad_chars)) {
  .warn(paste("non-ASCII character(s):", paste(.bad_chars, collapse = " ")))
} else {
  .ok("plain ASCII throughout")
}

.step("Confirming the DOI is still referenced in the required form")
if (any(grepl("<doi:10.2478/jamsi-2025-0010>", .d, fixed = TRUE))) {
  .ok("<doi:10.2478/jamsi-2025-0010> present")
} else {
  .warn("the DOI reference is missing; CRAN expects Authors (year) <doi:...>")
}

## =============================================================================
##  CITATION
## =============================================================================

.step("Making inst/CITATION track the installed version")

if (file.exists("inst/CITATION")) {
  .c <- readLines("inst/CITATION", warn = FALSE)
  if (any(grepl("0.2.0", .c, fixed = TRUE))) {
    .backup("inst/CITATION")
    .c <- gsub('note    = "R package version 0.2.0"',
               'note    = paste("R package version", meta$Version)',
               .c, fixed = TRUE)
    .c <- gsub('"Analysis. R package version 0.2.0.",',
               '"Analysis. R package version ", meta$Version, ".",',
               .c, fixed = TRUE)
    .write_lines("inst/CITATION", .c)
    if (any(grepl("0.2.0", .c, fixed = TRUE))) {
      .warn("a hard-coded 0.2.0 remains in inst/CITATION; check it by hand")
    } else {
      .ok("version now taken from meta$Version")
    }
  } else {
    .info("already version-independent")
  }
}

.step("Reporting the citation entry that needs your confirmation")
cat("\n")
cat("       The article is currently cited as:\n\n")
cat("         Ahmad, B., & Danish, M. Y. (2025). Development and\n")
cat("         characterization of a flexible three-parameter lifetime\n")
cat("         distribution: theoretical properties and real-world\n")
cat("         applications. Journal of Applied Mathematics, Statistics and\n")
cat("         Informatics, 21(1). doi:10.2478/jamsi-2025-0010\n\n")
cat("       Volume 21 = 2025 is confirmed. The ISSUE NUMBER is not.\n")
cat("       DOI 10.2478/jamsi-2025-0006 is in volume 21, issue 2\n")
cat("       (December 2025), pages 5-19, and JAMSI numbers its DOIs\n")
cat("       sequentially across a volume.\n\n")
cat("       Please confirm the issue number and the page range from your\n")
cat("       offprint, then edit inst/CITATION and DESCRIPTION if needed.\n")
cat("       Nothing has been changed automatically.\n\n")

## =============================================================================
##  README
## =============================================================================

.step("Rewriting the README overview and feature list")

.r <- readLines("README.md", warn = FALSE)
.o0 <- grep("^## Overview[ ]*$", .r)
if (length(.o0) != 1L) {
  .warn("could not find a single '## Overview' heading; README left alone")
} else {
  ## Replace everything from Overview up to the next level-2 heading.
  .rest <- grep("^## ", .r)
  .o1 <- .rest[.rest > .o0]
  .o1 <- if (length(.o1)) .o1[1] - 1L else length(.r)

  .new_readme <- c(
"## Overview",
"",
"The **BetaDanish** package implements the four-parameter Beta-Danish",
"distribution and its three-parameter Exponentiated Danish (ED) submodel for",
"survival, reliability and lifetime-data analysis. The distribution was",
"introduced by Ahmad and Danish (2025) and accommodates monotonic, unimodal",
"and bathtub-shaped hazards within a single parametric family.",
"",
"Its upper tail is regularly varying with index `-b`, so the survival function",
"decays polynomially rather than exponentially. That is what makes the family",
"useful for heavy-tailed lifetime data, and it has two consequences the package",
"takes seriously: `E(Z^r)` is finite only when `b > r`, and the moment",
"generating function does not exist at all.",
"",
"### Distribution functions",
"",
"`dbetadanish()`, `pbetadanish()`, `qbetadanish()`, `sbetadanish()`,",
"`hbetadanish()` and `rbetadanish()`, with `ded()`, `ped()`, `qed()`, `sed()`,",
"`hed()` and `red()` naming the ED submodel directly. All are evaluated so as",
"to hold accuracy far into the tail: the quantile function uses the beta mirror",
"identity rather than subtracting a near-one probability from one, and the",
"survival function is never formed by cancellation.",
"",
"### Estimation and inference",
"",
"- Maximum likelihood for complete and right-censored samples",
"- A **grouped likelihood** for times recorded on a coarse grid, where treating",
"  a rounded value as exact would understate every standard error",
"- **Ridge-penalized** fitting for the weakly identified regime",
"- **Bayesian** sampling by random-walk Metropolis",
"- Log-scale Wald and **profile-likelihood** intervals; an unbounded profile is",
"  reported as such rather than truncated at the grid edge",
"- `bd_identified_coef()` reports the fit through the identified composite",
"  `ac`, which is what the data determine when `a` and `c` cannot be separated",
"- Diagnostics that warn when a fit lands on the `b = 1` ridge, when the",
"  information matrix is singular, or when starts were discarded as degenerate",
"",
"### Structural properties",
"",
"Raw, incomplete and conditional moments with their existence conditions;",
"Shannon (closed form), Renyi and Tsallis entropies; mean residual life and",
"mean inactivity time; mean deviations; Lorenz and Bonferroni curves;",
"probability weighted moments; order-statistic densities, distributions and",
"moments; stress-strength reliability; hazard-shape classification via",
"Glaser's criterion; and the tail index.",
"",
"### Regression",
"",
"Accelerated failure time models, mixture and promotion-time cure models, and",
"competing risks with Aalen-Johansen comparison and Gray's test, each with",
"Cox-Snell residual diagnostics.",
"",
"### Working from a file",
"",
"`bd_analyze_csv()` takes a delimited file or spreadsheet through reading,",
"fitting, tabulation and optional figure output in one call. See the",
"\"Analysing Your Own Data from a CSV File\" vignette.",
"")

  .backup("README.md")
  .r <- append(.r[-(.o0:.o1)], .new_readme, after = .o0 - 1L)
  .write_lines("README.md", .r)
  .ok(sprintf("README overview rewritten (%d lines)", length(.new_readme)))
}

## =============================================================================
##  PKGDOWN REFERENCE INDEX
## =============================================================================

.step("Checking the pkgdown reference index against the exports")

if (!file.exists("_pkgdown.yml")) {
  .info("no _pkgdown.yml in the package root; nothing to check")
} else {
  .ns  <- readLines("NAMESPACE", warn = FALSE)
  .exp <- sub("^export\\((.*)\\)$", "\\1",
              grep("^export\\(", .ns, value = TRUE))
  .exp <- gsub('"', "", .exp, fixed = TRUE)

  .yml <- readLines("_pkgdown.yml", warn = FALSE)
  if (!any(grepl("^reference:", .yml))) {
    .info("_pkgdown.yml has no explicit reference: section, so pkgdown will")
    .info("index everything automatically. Nothing to do.")
  } else {
    .listed <- trimws(gsub("^[ \t-]*", "", grep("^[ \t]*-[ ]", .yml, value = TRUE)))
    .listed <- gsub('"', "", .listed, fixed = TRUE)
    .missing <- setdiff(.exp, .listed)
    if (!length(.missing)) {
      .ok("every exported function appears in the reference index")
    } else {
      .warn(sprintf("%d exported function(s) are missing from _pkgdown.yml.",
                    length(.missing)))
      .warn("pkgdown fails the build when a topic is absent from the index.")
      cat("\n       Paste this into the reference: section of _pkgdown.yml:\n\n")
      cat("  - title: Added in 0.3.0\n    contents:\n")
      for (f in sort(.missing)) cat("      - ", f, "\n", sep = "")
      cat("\n")
    }
  }
}

## =============================================================================
##  VERIFY
## =============================================================================

.step("Re-reading DESCRIPTION to confirm it parses")
.chk <- tryCatch(read.dcf("DESCRIPTION"), error = function(e) NULL)
if (is.null(.chk)) .die("DESCRIPTION no longer parses. Backups: ", BACKUP_DIR)
.ok(sprintf("DESCRIPTION parses; %d field(s)", ncol(.chk)))
cat("\n       Title:   ", .chk[1, "Title"], "\n", sep = "")
cat("       Version: ", .chk[1, "Version"], "\n", sep = "")
cat("       Desc:    ", substr(.chk[1, "Description"], 1, 110), "...\n", sep = "")

.step("Sourcing inst/CITATION to confirm it evaluates")
if (file.exists("inst/CITATION")) {
  .env <- new.env()
  assign("meta", list(Version = .chk[1, "Version"]), envir = .env)
  .r2 <- tryCatch({ source("inst/CITATION", local = .env); TRUE },
                  error = function(e) conditionMessage(e))
  if (isTRUE(.r2)) .ok("inst/CITATION evaluates cleanly")
  else .warn(paste("inst/CITATION failed to evaluate:", .r2))
}

.step("devtools::document()")
.r3 <- tryCatch({ devtools::document(); TRUE }, error = function(e) conditionMessage(e))
if (!isTRUE(.r3)) .die("document() failed:\n  ", .r3, "\n\nBackups: ", BACKUP_DIR)
.ok("documentation regenerated")

.step("devtools::check() -- several minutes, do not interrupt")
.chkres <- tryCatch(devtools::check(document = FALSE, args = "--as-cran", error_on = "never"),
                    error = function(e) { .warn(conditionMessage(e)); NULL })

cat("\n", strrep("=", 78), "\n", sep = "")
if (!is.null(.chkres)) {
  cat("  CHECK RESULT\n", strrep("=", 78), "\n", sep = "")
  cat(sprintf("  errors=%d  warnings=%d  notes=%d\n",
              length(.chkres$errors), length(.chkres$warnings), length(.chkres$notes)))
  for (nm in c("errors", "warnings", "notes")) {
    if (length(.chkres[[nm]])) {
      cat("\n---- ", toupper(nm), " ----\n", sep = "")
      cat(.chkres[[nm]], sep = "\n\n")
    }
  }
} else {
  cat("  check() did not complete; run devtools::check() manually.\n")
}

cat("\n", strrep("=", 78), "\n", sep = "")
cat("  PATCH 3c-meta COMPLETE\n")
cat(strrep("=", 78), "\n\n")
cat("  DESCRIPTION   Description field now describes the current package\n")
cat("  inst/CITATION version taken from meta$Version, never hard-coded\n")
cat("  README.md     overview and feature list rewritten to match the code\n\n")
cat("  NEEDS YOUR CONFIRMATION\n")
cat("    The JAMSI issue number, 21(1) or 21(2), and the page range.\n")
cat("    See the note printed above. Nothing was changed automatically.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
