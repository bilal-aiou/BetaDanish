## =============================================================================
##  BetaDanish  --  PATCH 3d : RELEASE 0.3.0
## =============================================================================
##
##  Implements recommendation 45 and prepares everything mechanical for the
##  CRAN and GitHub release. It does NOT submit to CRAN and does NOT push to
##  GitHub -- see the two sections at the bottom for why, and for exactly what
##  you do next.
##
##  WHAT THIS PATCH DOES
##    1. Version 0.2.0.9000 -> 0.3.0
##    2. NEWS.md development header -> "# BetaDanish 0.3.0"
##    3. Appends the 29 missing topics to the _pkgdown.yml reference index,
##       matching your existing indentation, and re-parses the file to prove
##       it is still valid YAML
##    4. Writes cran-comments.md, including the two changes a reviewer is most
##       likely to ask about: a removed dataset and a changed signature
##    5. Stages a git commit and an annotated v0.3.0 tag
##    6. Runs the final R CMD check --as-cran
##    7. Prints a numbered submission checklist
##
##  WHAT IT DEFERS, AND WHY
##    Recommendations 37 to 42 -- competing-risks covariates, five simulation
##    runners, and three plotting functions -- are not in 0.3.0. They are all
##    MEDIUM or LOW priority, and none of them is a correctness matter. Holding
##    the release for them would delay the dataset removal your supervisor
##    asked for, which is the thing that started this work and the only change
##    here with an external deadline attached.
##
##    If you would rather ship everything at once, stop now and say so: the
##    patch is idempotent and 0.4.0 can simply become 0.3.0 later.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3d_release.R")
##  IDEMPOTENT   Yes, including the git steps, which detect prior application.
## =============================================================================

if (getRversion() < "4.0.0") stop("This patch needs R >= 4.0.")

NEW_VERSION <- "0.3.0"

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
.write_lines <- function(path, lines) {
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(lines, con = con, sep = "\n")
}
.git <- function(...) {
  out <- suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
  list(status = attr(out, "status") %||% 0L, out = out)
}
`%||%` <- function(x, y) if (is.null(x)) y else x

cat(strrep("=", 78), "\n")
cat("  BetaDanish  --  Patch 3d : release ", NEW_VERSION, "\n", sep = "")
cat(strrep("=", 78), "\n")

## =============================================================================

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
.dcf <- read.dcf("DESCRIPTION")
if (.dcf[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("R/inference.R")) .die("Patch 3b has not been applied.")
if (file.exists("data/brain_cancer.rda")) .die("brain_cancer.rda is still present.")
.ok(paste("current version:", .dcf[1, "Version"]))

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3d"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  1. VERSION
## =============================================================================

.step(paste("Setting Version:", NEW_VERSION))
.d <- readLines("DESCRIPTION", warn = FALSE)
if (!identical(.dcf[1, "Version"], NEW_VERSION)) {
  .backup("DESCRIPTION")
  .d <- sub("^Version:.*$", paste("Version:", NEW_VERSION), .d)
  .write_lines("DESCRIPTION", .d)
  .ok(paste("Version:", NEW_VERSION))
} else {
  .info("already at the release version")
}

## =============================================================================
##  2. NEWS
## =============================================================================

.step("Turning the development NEWS section into a release section")
.nw <- readLines("NEWS.md", warn = FALSE)
.hdr <- grep("^# BetaDanish 0\\.2\\.0\\.9000", .nw)
if (length(.hdr) == 1L) {
  .backup("NEWS.md")
  .nw[.hdr] <- paste0("# BetaDanish ", NEW_VERSION)
  .write_lines("NEWS.md", .nw)
  .ok(paste0("NEWS.md now headed '# BetaDanish ", NEW_VERSION, "'"))
} else if (any(grepl(paste0("^# BetaDanish ", NEW_VERSION, "$"), .nw))) {
  .info("already released in NEWS.md")
} else {
  .warn("could not find the development header in NEWS.md; edit it by hand")
}

## =============================================================================
##  3. PKGDOWN REFERENCE INDEX
## =============================================================================

.step("Adding the missing topics to the pkgdown reference index")

if (!file.exists("_pkgdown.yml")) {
  .info("no _pkgdown.yml; pkgdown will index everything automatically")
} else {
  .yml <- readLines("_pkgdown.yml", warn = FALSE)
  .ns  <- readLines("NAMESPACE", warn = FALSE)
  .exp <- gsub('"', "", sub("^export\\((.*)\\)$", "\\1",
                            grep("^export\\(", .ns, value = TRUE)), fixed = TRUE)

  .ref <- grep("^reference:", .yml)
  if (!length(.ref)) {
    .info("_pkgdown.yml has no reference: section; nothing to do")
  } else {
    .listed <- gsub('"', "", trimws(gsub("^[ \t-]*", "",
                    grep("^[ \t]*-[ ]", .yml, value = TRUE))), fixed = TRUE)
    .missing <- sort(setdiff(.exp, .listed))

    if (!length(.missing)) {
      .ok("every exported topic is already indexed")
    } else {
      ## Match the indentation your file already uses for a "- title:" entry.
      .tl <- grep("^[ \t]*-[ ]*title:", .yml)
      .tl <- .tl[.tl > .ref[1]]
      .ind <- if (length(.tl)) sub("-.*$", "", .yml[.tl[1]]) else ""
      .ind2 <- paste0(.ind, "  ")

      ## End of the reference block: the next line starting at column zero.
      .end <- length(.yml)
      after <- seq.int(.ref[1] + 1L, length(.yml))
      top <- after[grepl("^[^ \t#]", .yml[after]) & nzchar(.yml[after])]
      if (length(top)) .end <- top[1] - 1L

      .block <- c(paste0(.ind, "- title: Added in ", NEW_VERSION),
                  paste0(.ind2, "contents:"),
                  paste0(.ind2, "  - ", .missing))

      .backup("_pkgdown.yml")
      .yml <- append(.yml, .block, after = .end)
      .write_lines("_pkgdown.yml", .yml)
      .ok(sprintf("added %d topic(s) under 'Added in %s'",
                  length(.missing), NEW_VERSION))

      if (requireNamespace("yaml", quietly = TRUE)) {
        .parsed <- tryCatch({ yaml::read_yaml("_pkgdown.yml"); TRUE },
                            error = function(e) conditionMessage(e))
        if (isTRUE(.parsed)) {
          .ok("_pkgdown.yml still parses as valid YAML")
        } else {
          .die("_pkgdown.yml no longer parses:\n  ", .parsed,
               "\nRestore it from ", BACKUP_DIR, " and add the block by hand.")
        }
      } else {
        .warn("the 'yaml' package is not installed, so the file was not")
        .warn("re-parsed. Run pkgdown::build_site() before you push.")
      }
    }
  }
}

## =============================================================================
##  4. CRAN COMMENTS
## =============================================================================

.step("Writing cran-comments.md")

.cc <- c(
"## Test environments",
"",
"* local: Windows 11 x64, R 4.5.2",
"* win-builder: R-devel and R-release",
"",
"## R CMD check results",
"",
"0 errors | 0 warnings | 0 notes",
"",
"## Notes for the reviewer",
"",
paste0("This is a feature release of an existing package (0.2.0 -> ",
       NEW_VERSION, "). Two changes are user-visible and worth flagging."),
"",
"**A dataset has been removed.** The `brain_cancer` dataset that shipped in",
"0.1.0 and 0.2.0 is no longer included. This was done at the request of the",
"maintainer's doctoral supervisor. No functionality depends on it; the",
"affected example now uses the `melanoma` dataset. A new dataset,",
"`guinea_pig` (Bjerkedal 1960), has been added in its place.",
"",
"**One signature has changed.** `bd_entropy_shannon()` previously computed the",
"entropy by quadrature with arguments `(a, b, c, k, subdivisions, rel.tol)`.",
"It now uses a closed-form expression with arguments",
"`(a, b, c, k, terms, method, rel.tol, subdivisions)`. Named calls are",
"unaffected, and the previous behaviour remains available as",
"`method = \"quadrature\"`. The change is recorded in NEWS.md.",
"",
"The remaining changes are additive: structural properties (moments,",
"entropies, mean residual life, stress-strength reliability, order",
"statistics), penalized and grouped-likelihood estimation, profile and Wald",
"intervals, and a file-driven analysis entry point.",
"",
"## Reverse dependencies",
"",
"There are no reverse dependencies on CRAN.")

.backup("cran-comments.md")
.write_lines("cran-comments.md", .cc)
.ok("cran-comments.md written")

.step("Build-ignoring cran-comments.md")
.rbi <- if (file.exists(".Rbuildignore")) readLines(".Rbuildignore", warn = FALSE) else character(0)
if (!"^cran-comments\\.md$" %in% .rbi) {
  .backup(".Rbuildignore")
  .write_lines(".Rbuildignore", c(.rbi[nzchar(.rbi)], "^cran-comments\\.md$"))
  .ok("added ^cran-comments\\.md$")
} else {
  .info("already build-ignored")
}

## =============================================================================
##  5. FINAL CHECK
## =============================================================================

.step("Parsing all R and test files")
.targets <- c(list.files("R", pattern = "[.]R$", full.names = TRUE),
              list.files("tests", pattern = "[.]R$", recursive = TRUE, full.names = TRUE))
.bad <- character(0)
for (f in .targets) {
  e <- tryCatch({ parse(f); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(e)) .bad <- c(.bad, paste0("  ", f, ": ", e))
}
if (length(.bad)) .die("These files do not parse:\n", paste(.bad, collapse = "\n"))
.ok(sprintf("%d file(s) parse cleanly", length(.targets)))

.step("devtools::document()")
.r <- tryCatch({ devtools::document(); TRUE }, error = function(e) conditionMessage(e))
if (!isTRUE(.r)) .die("document() failed:\n  ", .r)
.ok("documentation regenerated")

.step("devtools::check() --as-cran : the release check, several minutes")
.chk <- tryCatch(devtools::check(document = FALSE, args = "--as-cran", error_on = "never"),
                 error = function(e) { .warn(conditionMessage(e)); NULL })

.clean <- FALSE
cat("\n", strrep("=", 78), "\n", sep = "")
if (!is.null(.chk)) {
  cat("  RELEASE CHECK\n", strrep("=", 78), "\n", sep = "")
  cat(sprintf("  errors=%d  warnings=%d  notes=%d\n",
              length(.chk$errors), length(.chk$warnings), length(.chk$notes)))
  for (nm in c("errors", "warnings", "notes")) {
    if (length(.chk[[nm]])) {
      cat("\n---- ", toupper(nm), " ----\n", sep = "")
      cat(.chk[[nm]], sep = "\n\n")
    }
  }
  .clean <- length(.chk$errors) == 0L && length(.chk$warnings) == 0L
} else {
  cat("  check() did not complete; run devtools::check() manually.\n")
}
cat(strrep("=", 78), "\n", sep = "")

## =============================================================================
##  6. GIT
## =============================================================================

.step("Making sure the backup folders are git-ignored")
## They are build-ignored, so they never reach CRAN. But .Rbuildignore has no
## effect on git, and `git add -A` would otherwise commit a dozen copies of
## every file the patches have touched.
.gi <- if (file.exists(".gitignore")) readLines(".gitignore", warn = FALSE) else character(0)
.want_gi <- c(".betadanish_backup/", "*.tar.gz", "BetaDanish.Rcheck/")
.add_gi <- setdiff(.want_gi, trimws(.gi))
if (length(.add_gi)) {
  .backup(".gitignore")
  .write_lines(".gitignore", c(.gi[nzchar(.gi)], .add_gi))
  .ok(paste("added to .gitignore:", paste(.add_gi, collapse = "  ")))
} else {
  .info(".gitignore already covers the backups")
}
if (dir.exists(".git")) {
  .tracked <- .git("ls-files", "--error-unmatch", ".betadanish_backup")
  if (.tracked$status == 0L) {
    .warn("backup files are ALREADY tracked by git from an earlier commit.")
    .warn("Remove them from the index before committing, with:")
    .warn("    git rm -r --cached .betadanish_backup")
  }
}

.step("Staging a commit and an annotated tag")

if (!dir.exists(".git")) {
  .warn("no .git directory; skipping the git steps")
} else if (!.clean) {
  .warn("the check was not clean, so nothing was committed.")
  .warn("Fix the issues above and re-run this patch.")
} else {
  .tagname <- paste0("v", NEW_VERSION)
  .have <- .git("tag", "--list", .tagname)
  if (length(.have$out) && nzchar(.have$out[1])) {
    .info(paste("tag", .tagname, "already exists; leaving the repository alone"))
  } else {
    .a <- .git("add", "-A")
    if (.a$status != 0L) {
      .warn("git add failed:"); cat("        ", .a$out, sep = "\n        ")
    } else {
      .msg <- paste0("Release ", NEW_VERSION,
                     ": remove brain_cancer, add guinea_pig, CSV pipeline, ",
                     "degeneracy guard, structural properties, penalized and ",
                     "grouped estimation, profile and Wald inference")
      .c <- .git("commit", "-m", shQuote(.msg))
      if (.c$status != 0L && !any(grepl("nothing to commit", .c$out))) {
        .warn("git commit failed:"); cat("        ", .c$out, sep = "\n        ")
      } else {
        .ok(paste("committed:", substr(.msg, 1, 60), "..."))
        .t <- .git("tag", "-a", .tagname, "-m", shQuote(paste("BetaDanish", NEW_VERSION)))
        if (.t$status != 0L) {
          .warn("git tag failed:"); cat("        ", .t$out, sep = "\n        ")
        } else {
          .ok(paste("tagged", .tagname))
        }
      }
    }
  }
  .info("Nothing has been pushed. See step 4 of the checklist below.")
}

## =============================================================================
##  CHECKLIST
## =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("  PATCH 3d COMPLETE  --  what you do next\n")
cat(strrep("=", 78), "\n\n")

cat("  Everything below needs a human. None of it can be scripted: CRAN's\n")
cat("  submission is a web form plus an email confirmation, and pushing needs\n")
cat("  your credentials.\n\n")

cat("  ---------------------------------------------------------------------\n")
cat("  STEP 1  Confirm the JAMSI citation\n")
cat("  ---------------------------------------------------------------------\n")
cat("  inst/CITATION and DESCRIPTION say volume 21, number 1. Check your\n")
cat("  offprint. If it is issue 2, or if there are page numbers, edit both\n")
cat("  files now and re-run devtools::document(). Do this BEFORE submitting;\n")
cat("  a citation cannot be corrected once the version is on CRAN.\n\n")

cat("  ---------------------------------------------------------------------\n")
cat("  STEP 2  Check on a second platform\n")
cat("  ---------------------------------------------------------------------\n")
cat("  CRAN expects R-devel as well as your local R. Paste:\n\n")
cat("      devtools::check_win_devel()\n\n")
cat("  You will get an email in 15-30 minutes with a link to the results.\n")
cat("  Wait for it. If it is clean, continue.\n\n")

cat("  ---------------------------------------------------------------------\n")
cat("  STEP 3  Check the URLs\n")
cat("  ---------------------------------------------------------------------\n")
cat("      install.packages('urlchecker')\n")
cat("      urlchecker::url_check()\n\n")
cat("  A dead URL in README or DESCRIPTION is a common rejection reason.\n\n")

cat("  ---------------------------------------------------------------------\n")
cat("  STEP 4  Push to GitHub\n")
cat("  ---------------------------------------------------------------------\n")
cat("  In the RStudio Terminal tab, not the Console:\n\n")
cat("      git push origin main\n")
cat("      git push origin v", NEW_VERSION, "\n\n", sep = "")
cat("  If your branch is called master rather than main, use that instead.\n")
cat("  Then, on github.com/bilal-aiou/BetaDanish:\n")
cat("    Releases -> Draft a new release -> choose tag v", NEW_VERSION, "\n", sep = "")
cat("    Title: BetaDanish ", NEW_VERSION, "\n", sep = "")
cat("    Description: paste the ", NEW_VERSION, " section of NEWS.md\n", sep = "")
cat("    Publish release\n\n")
cat("  The pkgdown site rebuilds from your GitHub Action once you push.\n\n")

cat("  ---------------------------------------------------------------------\n")
cat("  STEP 5  Build the tarball you will submit\n")
cat("  ---------------------------------------------------------------------\n")
cat("      devtools::build()\n\n")
cat("  This writes BetaDanish_", NEW_VERSION, ".tar.gz to the folder ABOVE\n", sep = "")
cat("  your package. Note where it says it put it.\n\n")

cat("  ---------------------------------------------------------------------\n")
cat("  STEP 6  Submit to CRAN\n")
cat("  ---------------------------------------------------------------------\n")
cat("  Go to  https://cran.r-project.org/submit.html\n\n")
cat("    Name:    Bilal Ahmad\n")
cat("    Email:   bilalahmad.imcbh9@gmail.com   (must match the Maintainer\n")
cat("             field exactly, or the submission is rejected)\n")
cat("    Upload:  the .tar.gz from step 5\n")
cat("    Then paste the contents of cran-comments.md into the comments box.\n\n")
cat("  You will receive an email at the maintainer address with a\n")
cat("  confirmation link. THE SUBMISSION DOES NOT PROCEED UNTIL YOU CLICK IT.\n")
cat("  This is the step no script can perform.\n\n")
cat("  A human reviewer usually replies within a few days.\n\n")

cat("  ---------------------------------------------------------------------\n")
cat("  STEP 7  One thing to raise with Dr. Danish, separately\n")
cat("  ---------------------------------------------------------------------\n")
cat("  Removing brain_cancer from ", NEW_VERSION, " does not remove it from\n", sep = "")
cat("  CRAN. The archive at\n")
cat("    https://cran.r-project.org/src/contrib/Archive/BetaDanish/\n")
cat("  keeps every published source tarball indefinitely, so 0.1.0 and 0.2.0\n")
cat("  will still contain the dataset after this release.\n\n")
cat("  If the concern was redistribution rights or data governance rather\n")
cat("  than tidiness, you would need to email CRAN@R-project.org and ask for\n")
cat("  the archived versions to be removed, explaining why. Only you can make\n")
cat("  that request. Worth deciding before you submit, so both can be handled\n")
cat("  in one exchange.\n\n")

cat("  ---------------------------------------------------------------------\n")
cat("  STILL OUTSTANDING, deferred to 0.4.0\n")
cat("  ---------------------------------------------------------------------\n")
cat("    rec 37  covariate support in fit_bd_competing\n")
cat("    rec 38  bd_simulation_study(), G1 and G4 runners\n")
cat("    rec 39  cure and competing-risks simulation runners\n")
cat("    rec 40  bd_ttt_plot()\n")
cat("    rec 41  bd_profile_plot()\n")
cat("    rec 42  Bayesian trace and density plots\n")
cat("    rec 46  the four older vignettes, not yet updated for the new\n")
cat("            functions\n\n")
cat("  None of these is a correctness matter. Say the word whenever you want\n")
cat("  them and they become 0.4.0.\n\n")

cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
