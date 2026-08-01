## =============================================================================
##  BetaDanish  --  PATCH 3h-fix : index by Rd topic, not by NAMESPACE export
## =============================================================================
##
##  Patch 3h fixed the two faults pkgdown originally named, and then its own
##  check said "every export and vignette appears in the index" -- immediately
##  before pkgdown reported two topics missing. My check was wrong, and it was
##  wrong for a reason worth writing down.
##
##  I COMPARED THE INDEX AGAINST THE WRONG UNIVERSE
##
##    My check:  exports parsed out of NAMESPACE
##    pkgdown:   every man/*.Rd topic not marked @keywords internal
##
##  The second set is strictly larger, and the two missing items are exactly
##  the kinds of thing that live in the gap:
##
##    guinea_pig     a DATASET. Documented as "guinea_pig" in data.R and made
##                   available through LazyData, never through export().
##    plot.bd_bayes  an S3 METHOD. Registered as S3method(plot,bd_bayes), not
##                   as export(plot.bd_bayes).
##
##  Neither appears in NAMESPACE as an export, so a setdiff against the export
##  list could not see them however carefully it was written. Patch 3d looked
##  one direction, Patch 3h looked both directions but at the wrong set. This
##  patch uses the set pkgdown actually uses.
##
##  WHAT THIS PATCH DOES
##    Adds the two missing topics under headings that fit the existing index,
##    then re-checks using Rd topics as the universe, in both directions, and
##    rebuilds the site so you see it succeed rather than assume it.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3h_fix.R")
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
cat("  BetaDanish  --  Patch 3h-fix : index by Rd topic\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("_pkgdown.yml")) .die("No _pkgdown.yml in the package root.")
if (!requireNamespace("pkgdown", quietly = TRUE))
  .die("The 'pkgdown' package is required.")
.ok("package and _pkgdown.yml found")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3hfix"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  THE RIGHT UNIVERSE
## =============================================================================

#' Every topic pkgdown will expect to find in the index: one entry per Rd
#' file, skipping those marked internal and the package-level page, which
#' pkgdown handles itself.
.rd_topics <- function() {
  out <- character(0)
  for (f in list.files("man", pattern = "[.]Rd$", full.names = TRUE)) {
    txt <- readLines(f, warn = FALSE)
    if (any(grepl("^.keyword[{]internal[}]", txt))) next
    nm <- sub("[.]Rd$", "", basename(f))
    if (grepl("-package$", nm)) next
    out <- c(out, nm)
  }
  sort(unique(out))
}

#' Aliases let an index entry name any one of a topic's aliases rather than
#' its file name, so resolution must consider all of them.
.rd_aliases <- function() {
  out <- character(0)
  for (f in list.files("man", pattern = "[.]Rd$", full.names = TRUE)) {
    a <- grep("^.alias[{]", readLines(f, warn = FALSE), value = TRUE)
    out <- c(out, trimws(sub("[}].*$", "", sub("^.alias[{]", "", a))))
  }
  sort(unique(out))
}

.index_entries <- function(yml) {
  e <- gsub('"', "", trimws(gsub("^[ \t-]*", "",
            grep("^[ \t]*-[ ]", yml, value = TRUE))), fixed = TRUE)
  e[!grepl(":", e)]
}

## =============================================================================
##  1. ADD WHAT IS MISSING
## =============================================================================

.step("Finding topics pkgdown expects but the index omits")

.yml   <- readLines("_pkgdown.yml", warn = FALSE)
.vigs  <- sub("[.]Rmd$", "", list.files("vignettes", pattern = "[.]Rmd$"))
.listed <- .index_entries(.yml)
.topics <- .rd_topics()
.alias  <- .rd_aliases()

## A topic counts as indexed if the index names its file name or any alias.
.covered <- function(t) {
  if (t %in% .listed) return(TRUE)
  f <- file.path("man", paste0(t, ".Rd"))
  if (!file.exists(f)) return(FALSE)
  a <- grep("^.alias[{]", readLines(f, warn = FALSE), value = TRUE)
  a <- trimws(sub("[}].*$", "", sub("^.alias[{]", "", a)))
  any(a %in% .listed)
}
.missing <- Filter(function(t) !.covered(t), .topics)

if (!length(.missing)) {
  .ok("every Rd topic is already indexed")
} else {
  cat("       missing: ", paste(.missing, collapse = ", "), "\n", sep = "")

  ## Datasets go with the datasets, plot methods with the plotting section,
  ## anything else into a catch-all. Matching the file's own indentation.
  .c <- grep("^[ \t]*-[ ]", .yml)
  .ind <- if (length(.c)) sub("-.*$", "", .yml[.c[1]]) else "    "

  ## Walk to the end of the contents list that the anchor belongs to. The
  ## stop condition must exclude "- title:", which also matches a bare "- "
  ## pattern and would otherwise place the new entry inside the NEXT section.
  .place <- function(topic, after_pattern) {
    hit <- grep(after_pattern, .yml)
    if (!length(hit)) return(FALSE)
    j <- hit[1]
    repeat {
      if (j >= length(.yml)) break
      nxt <- .yml[j + 1L]
      if (!grepl("^[ \t]*-[ ]", nxt)) break
      if (grepl("^[ \t]*-[ ]*title[ ]*:", nxt)) break
      j <- j + 1L
    }
    .yml <<- append(.yml, paste0(.ind, "- ", topic), after = j)
    TRUE
  }

  .leftover <- character(0)
  for (t in .missing) {
    done <- FALSE
    if (grepl("^plot[.]", t)) {
      done <- .place(t, "^[ \t]*-[ ]*plot[.]betadanish[ ]*$")
    } else if (t %in% c("guinea_pig", "remission", "carbon_fibres", "aarset",
                        "leukemia", "melanoma", "transplant")) {
      done <- .place(t, "^[ \t]*-[ ]*(remission|carbon_fibres|transplant)[ ]*$")
    }
    if (!done) .leftover <- c(.leftover, t)
  }

  if (length(.leftover)) {
    ## Anything unplaced goes in its own section at the end of reference:.
    .ref <- grep("^reference:", .yml)
    if (length(.ref)) {
      after <- seq.int(.ref[1] + 1L, length(.yml))
      top <- after[grepl("^[^ \t#]", .yml[after]) & nzchar(.yml[after])]
      end <- if (length(top)) top[1] - 1L else length(.yml)
      blk <- c(paste0(.ind, "- title: Additional topics"),
               paste0(.ind, "  contents:"),
               paste0(.ind, "    - ", .leftover))
      .yml <- append(.yml, blk, after = end)
      .info(paste("placed in a new section:", paste(.leftover, collapse = ", ")))
    }
  }

  .backup("_pkgdown.yml")
  .write_lines("_pkgdown.yml", .yml)
  .ok(sprintf("added %d topic(s) to the index", length(.missing)))
}

## =============================================================================
##  2. RE-CHECK, USING THE RIGHT UNIVERSE, IN BOTH DIRECTIONS
## =============================================================================

.step("Re-checking the index against man/*.Rd, both directions")

.yml    <- readLines("_pkgdown.yml", warn = FALSE)
.listed <- .index_entries(.yml)
.topics <- .rd_topics()
.alias  <- .rd_aliases()

## A: index entries that resolve to nothing.
.dangling <- setdiff(.listed, c(.topics, .alias, .vigs))
if (length(.dangling)) {
  .warn("index entries with no matching topic:")
  for (d in .dangling) cat("        - ", d, "\n", sep = "")
  .die("pkgdown fails on a dangling entry. Backups: ", BACKUP_DIR)
}
.ok(sprintf("%d index entr(ies), every one resolves", length(.listed)))

## B: topics pkgdown expects that the index omits.
.still <- Filter(function(t) !.covered(t), .topics)
if (length(.still)) {
  .warn("Rd topics still missing from the index:")
  for (s in .still) cat("        - ", s, "\n", sep = "")
  .die("pkgdown fails on these too. Backups: ", BACKUP_DIR)
}
.ok(sprintf("%d Rd topic(s), every one indexed", length(.topics)))

## C: vignettes.
.vlisted <- setdiff(.vigs, .listed)
if (length(.vlisted)) {
  .warn(paste("vignettes not indexed:", paste(.vlisted, collapse = ", ")))
} else {
  .ok(sprintf("%d vignette(s), every one indexed", length(.vigs)))
}

.step("Confirming _pkgdown.yml is still valid YAML")
if (requireNamespace("yaml", quietly = TRUE)) {
  .p <- tryCatch({ yaml::read_yaml("_pkgdown.yml"); TRUE },
                 error = function(e) conditionMessage(e))
  if (isTRUE(.p)) {
    .ok("parses cleanly")
  } else {
    .die("_pkgdown.yml no longer parses:\n  ", .p,
         "\nRestore from ", BACKUP_DIR, " and edit by hand.")
  }
} else {
  .info("the 'yaml' package is not installed; skipping the parse check")
}

## =============================================================================
##  3. BUILD
## =============================================================================

.step("pkgdown::build_site() -- a few minutes, do not interrupt")

.res <- tryCatch({ pkgdown::build_site(preview = FALSE); TRUE },
                 error = function(e) conditionMessage(e))

cat("\n", strrep("=", 78), "\n", sep = "")
if (isTRUE(.res)) {
  cat("  SITE BUILT SUCCESSFULLY\n")
  cat(strrep("=", 78), "\n\n")
  cat("  docs/ contains ", length(list.files("docs", recursive = TRUE)),
      " file(s). Open docs/index.html to preview.\n\n", sep = "")
  cat("  TO PUBLISH, in the Terminal tab:\n\n")
  cat("      git add -A\n")
  cat("      git commit -m \"Fix pkgdown reference index\"\n")
  cat("      git push origin main\n\n")
} else {
  cat("  SITE BUILD FAILED\n")
  cat(strrep("=", 78), "\n\n")
  cat("  ", .res, "\n\n", sep = "")
  cat("  If it names further missing topics, send me the message and I will\n")
  cat("  place them; the check above should have caught them, so a residual\n")
  cat("  failure means pkgdown is applying a rule I have not accounted for.\n\n")
}
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
