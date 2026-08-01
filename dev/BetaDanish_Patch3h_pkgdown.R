## =============================================================================
##  BetaDanish  --  PATCH 3h : the pkgdown site build
## =============================================================================
##
##  pkgdown::build_site() reported exactly two faults, and both are mine:
##
##    x In _pkgdown.yml, reference[11].contents[7] (brain_cancer) must be a
##      known topic name or alias.
##    x In _pkgdown.yml, 1 vignette missing from index: "bd-csv-workflow".
##
##  WHY MY EARLIER CHECK MISSED THIS
##    Patch 3d compared the NAMESPACE exports against the index and reported
##    "0 missing". That check ran in one direction only: it found exports
##    absent from the index, never index entries pointing at topics that no
##    longer exist. brain_cancer was removed in Patch 2 -- dataset, roxygen
##    block and man page -- but its line in _pkgdown.yml survived, and my
##    check could not see it. The vignette added in Patch 2b was never indexed
##    at all.
##
##  WHAT THIS PATCH DOES
##    Removes the dangling brain_cancer entry, adds the missing vignette to the
##    articles index, then checks the index in BOTH directions and rebuilds the
##    site locally so you see it succeed before pushing.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3h_pkgdown.R")
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
.write_lines <- function(path, lines) {
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(lines, con = con, sep = "\n")
}

cat(strrep("=", 78), "\n")
cat("  BetaDanish  --  Patch 3h : pkgdown site build\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("_pkgdown.yml")) .die("No _pkgdown.yml in the package root.")
if (!requireNamespace("pkgdown", quietly = TRUE))
  .die("The 'pkgdown' package is required. install.packages('pkgdown')")
.ok("package and _pkgdown.yml found")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3h"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  1. REMOVE THE DANGLING brain_cancer ENTRY
## =============================================================================

.step("Removing the brain_cancer entry from the reference index")

.yml <- readLines("_pkgdown.yml", warn = FALSE)
.hit <- grep("^[ \t]*-[ ]*brain_cancer[ ]*$", .yml)
if (length(.hit)) {
  .backup("_pkgdown.yml")
  cat("       removing line ", paste(.hit, collapse = ", "), ": ",
      trimws(.yml[.hit[1]]), "\n", sep = "")
  .yml <- .yml[-.hit]
  .write_lines("_pkgdown.yml", .yml)
  .ok(sprintf("removed %d line(s)", length(.hit)))
} else {
  .info("no brain_cancer entry found; already removed")
}

## =============================================================================
##  2. ADD THE MISSING VIGNETTE TO THE ARTICLES INDEX
## =============================================================================

.step("Adding bd-csv-workflow to the articles index")

.yml <- readLines("_pkgdown.yml", warn = FALSE)
.vigs <- sub("[.]Rmd$", "", list.files("vignettes", pattern = "[.]Rmd$"))
.art  <- grep("^articles:", .yml)

if (!length(.art)) {
  .info("_pkgdown.yml has no articles: section, so pkgdown indexes them all")
} else {
  ## Everything listed anywhere in the file, so a vignette named under any
  ## articles subsection counts as indexed.
  .listed <- gsub('"', "", trimws(gsub("^[ \t-]*", "",
                  grep("^[ \t]*-[ ]", .yml, value = TRUE))), fixed = TRUE)
  .missing <- setdiff(.vigs, .listed)

  if (!length(.missing)) {
    .ok("every vignette is already indexed")
  } else {
    .backup("_pkgdown.yml")
    ## Match the indentation the file already uses for a contents entry.
    .c <- grep("^[ \t]*-[ ]", .yml)
    .c <- .c[.c > .art[1]]
    .ind <- if (length(.c)) sub("-.*$", "", .yml[.c[1]]) else "    "

    ## End of the articles block: the next line starting at column zero.
    .end <- length(.yml)
    after <- seq.int(.art[1] + 1L, length(.yml))
    top <- after[grepl("^[^ \t#]", .yml[after]) & nzchar(.yml[after])]
    if (length(top)) .end <- top[1] - 1L

    .block <- c(paste0(.ind, "- title: Working from a file"),
                paste0(.ind, "  contents:"),
                paste0(.ind, "    - ", .missing))
    .yml <- append(.yml, .block, after = .end)
    .write_lines("_pkgdown.yml", .yml)
    .ok(sprintf("added: %s", paste(.missing, collapse = ", ")))
  }
}

## =============================================================================
##  3. CHECK THE INDEX IN BOTH DIRECTIONS
## =============================================================================

.step("Checking the reference index in BOTH directions")

.yml <- readLines("_pkgdown.yml", warn = FALSE)
.listed <- gsub('"', "", trimws(gsub("^[ \t-]*", "",
                grep("^[ \t]*-[ ]", .yml, value = TRUE))), fixed = TRUE)
.listed <- .listed[!grepl(":", .listed)]          # drop "- title: ..." lines

## Every topic pkgdown knows about: Rd file names plus their aliases.
.topics <- character(0)
for (f in list.files("man", pattern = "[.]Rd$", full.names = TRUE)) {
  .topics <- c(.topics, sub("[.]Rd$", "", basename(f)))
  a <- grep("^.alias[{]", readLines(f, warn = FALSE), value = TRUE)
  .topics <- c(.topics, trimws(sub("[}].*$", "", sub("^.alias[{]", "", a))))
}
.topics <- unique(c(.topics, .vigs))

## Direction A: listed but not real. This is what broke the build, and what
## my Patch 3d check could not see.
.dangling <- setdiff(.listed, .topics)
if (length(.dangling)) {
  .warn("index entries with no matching topic:")
  for (d in .dangling) cat("        - ", d, "\n", sep = "")
  .die("pkgdown fails the build on a dangling entry. Backups: ", BACKUP_DIR)
}
.ok(sprintf("%d index entr(ies), every one resolves to a real topic",
            length(.listed)))

## Direction B: real but not listed.
.ns  <- readLines("NAMESPACE", warn = FALSE)
.exp <- gsub('"', "", sub("^export\\((.*)\\)$", "\\1",
                          grep("^export\\(", .ns, value = TRUE)), fixed = TRUE)
.unlisted <- sort(setdiff(c(.exp, .vigs), .listed))
if (length(.unlisted)) {
  .warn("exported topics or vignettes not in the index:")
  for (u in .unlisted) cat("        - ", u, "\n", sep = "")
  .warn("pkgdown will also fail on these; add them to _pkgdown.yml")
} else {
  .ok("every export and vignette appears in the index")
}

.step("Confirming _pkgdown.yml is still valid YAML")
if (requireNamespace("yaml", quietly = TRUE)) {
  .p <- tryCatch({ yaml::read_yaml("_pkgdown.yml"); TRUE },
                 error = function(e) conditionMessage(e))
  if (isTRUE(.p)) .ok("parses cleanly")
  else .die("_pkgdown.yml no longer parses:\n  ", .p,
            "\nRestore it from ", BACKUP_DIR, " and edit by hand.")
} else {
  .info("the 'yaml' package is not installed; skipping the parse check")
}

## =============================================================================
##  4. BUILD THE SITE
## =============================================================================

.step("pkgdown::build_site() -- a few minutes, do not interrupt")

.res <- tryCatch({
  pkgdown::build_site(preview = FALSE)
  TRUE
}, error = function(e) conditionMessage(e))

cat("\n", strrep("=", 78), "\n", sep = "")
if (isTRUE(.res)) {
  cat("  SITE BUILT SUCCESSFULLY\n")
  cat(strrep("=", 78), "\n\n")
  .n <- length(list.files("docs", recursive = TRUE))
  cat("  docs/ now contains ", .n, " file(s).\n", sep = "")
  cat("  Open docs/index.html in a browser to preview it.\n\n")
  cat("  TO PUBLISH, in the Terminal tab:\n\n")
  cat("      git add -A\n")
  cat("      git commit -m \"Fix pkgdown reference index\"\n")
  cat("      git push origin main\n\n")
  cat("  The pkgdown Action rebuilds the site from that push.\n\n")
} else {
  cat("  SITE BUILD FAILED\n")
  cat(strrep("=", 78), "\n\n")
  cat("  ", .res, "\n\n", sep = "")
  cat("  The message above names the offending entry. If it is another\n")
  cat("  dangling topic, remove that line from _pkgdown.yml and re-run.\n\n")
}

cat("  Note: 'already exists' warnings about docs/deps are harmless -- they\n")
cat("  mean pkgdown is refreshing a site that was built before.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
