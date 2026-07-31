## =============================================================================
##  BetaDanish  --  PATCH 3a-fix3 : remove the duplicate bd_entropy_shannon
## =============================================================================
##
##  ROOT CAUSE, finally visible in the fix2 log:
##
##      Writing bd_entropy_shannon.Rd
##      Writing bd_entropy.Rd
##
##  Two Rd files, because there are two functions of that name. R/entropy.R has
##  shipped since 0.2.0 with a quadrature bd_entropy_shannon(a, b, c, k,
##  subdivisions, rel.tol). Patch 3a added a second, closed-form version in
##  R/structural.R and left the first in place. Every warning in this sequence
##  came from the two roxygen blocks contesting one topic:
##
##    fix   -> "documented arguments not in usage: subdivisions, rel.tol"
##             (merged blocks, usage taken from only one of them)
##    fix2  -> "duplicated argument entries: a b c k rel.tol subdivisions"
##             (exactly the six names the two signatures then shared)
##    fix3  -> "duplicated alias" once the new topic was renamed
##
##  R/entropy.R is removed. Its only contents were that function. The
##  closed-form version in R/structural.R has been the one actually in use all
##  along -- collation is alphabetical, so structural.R loaded last and
##  overwrote the definition -- which is why all 515 tests pass either way.
##
##  SIGNATURE CHANGE
##    old: bd_entropy_shannon(a, b, c, k, subdivisions = 2000, rel.tol = 1e-8)
##    new: bd_entropy_shannon(a, b, c, k, terms = 20000L,
##                            method = c("closed", "quadrature"),
##                            rel.tol = 1e-10, subdivisions = 4000L)
##    Named calls are unaffected. A positional call passing subdivisions fifth
##    would now set `terms`. This is recorded in NEWS.
##
##  NEW GUARD
##    The patch scans R/ for any function name defined in more than one file
##    and aborts if it finds one. That check would have caught this at Patch 3a.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3a_fix3.R")
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

cat(strrep("=", 78), "\n")
cat("  BetaDanish  --  Patch 3a-fix3 : remove the duplicate definition\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("R/structural.R")) .die("Patch 3a has not been applied.")
if (!any(grepl("@name bd_entropy", readLines("R/structural.R", warn = FALSE), fixed = TRUE)))
  .die("Patch 3a-fix2 has not been applied -- R/structural.R has no unified entropy topic.")
.ok("Patch 3a and 3a-fix2 detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3afix3"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

.step("Confirming R/entropy.R holds nothing but the superseded function")
if (file.exists("R/entropy.R")) {
  .e <- readLines("R/entropy.R", warn = FALSE)
  .defs_here <- grep("^[a-zA-Z_.][a-zA-Z0-9_.]*[ ]*<-[ ]*function", .e, value = TRUE)
  cat("        R/entropy.R defines: ",
      paste(sub("[ ]*<-.*$", "", .defs_here), collapse = ", "), "\n", sep = "")
  if (length(.defs_here) != 1L ||
      !grepl("^bd_entropy_shannon", .defs_here[1])) {
    .die("R/entropy.R defines something other than bd_entropy_shannon alone.\n",
         "Inspect it by hand rather than letting this patch delete it.")
  }
  .backup("R/entropy.R")
  file.remove("R/entropy.R")
  .ok("R/entropy.R removed (a copy is in the backup directory)")
} else {
  .info("R/entropy.R already absent")
}

.step("Removing the orphaned Rd file")
if (file.exists("man/bd_entropy_shannon.Rd")) {
  .backup("man/bd_entropy_shannon.Rd")
  file.remove("man/bd_entropy_shannon.Rd")
  .ok("man/bd_entropy_shannon.Rd removed")
} else {
  .info("already absent")
}

.step("Scanning R/ for any function defined in more than one file")
.defs <- list()
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  txt <- readLines(f, warn = FALSE)
  hits <- grep("^[a-zA-Z_.][a-zA-Z0-9_.]*[ ]*<-[ ]*function", txt, value = TRUE)
  for (nm in sub("[ ]*<-.*$", "", hits))
    .defs[[nm]] <- c(.defs[[nm]], basename(f))
}
.dups <- Filter(function(v) length(unique(v)) > 1L, .defs)
if (length(.dups)) {
  for (nm in names(.dups))
    .warn(sprintf("%s is defined in: %s", nm, paste(unique(.dups[[nm]]), collapse = ", ")))
  .die("A function is defined in more than one file. Whichever collates last ",
       "silently wins, and roxygen will contest the Rd topic.\nBackups: ", BACKUP_DIR)
}
.ok(sprintf("%d function definition(s) scanned, every name unique", length(.defs)))

.step("Recording the replacement in NEWS.md")
.nw <- readLines("NEWS.md", warn = FALSE)
if (!any(grepl("superseded by the closed form", .nw, fixed = TRUE))) {
  .anchor <- grep("^\\* \\*\\*Entropies\\*\\*", .nw)
  if (length(.anchor) >= 1L) {
    .backup("NEWS.md")
    .add <- c(
      "",
      "  The quadrature `bd_entropy_shannon()` that shipped in 0.2.0 from",
      "  `R/entropy.R` is superseded by the closed form and that file is",
      "  removed. The signature changes from",
      "  `(a, b, c, k, subdivisions, rel.tol)` to",
      "  `(a, b, c, k, terms, method, rel.tol, subdivisions)`. Named calls are",
      "  unaffected; a positional call passing `subdivisions` fifth would now",
      "  set `terms`. The old behaviour is available as",
      "  `method = \"quadrature\"`.")
    .nw <- append(.nw, .add, after = .anchor[1] + 2L)
    con <- file("NEWS.md", open = "wb"); writeLines(.nw, con = con, sep = "\n"); close(con)
    .ok("NEWS.md updated")
  } else {
    .warn("anchor not found in NEWS.md; add the note by hand")
  }
} else {
  .info("NEWS.md already records it")
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

.rd_arg_names <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  i0  <- regexpr("arguments[{]", txt)
  if (i0 < 0) return(character(0))
  rest <- substring(txt, i0)
  i1 <- regexpr("\n[}]\n", rest)
  if (i1 > 0) rest <- substring(rest, 1, i1)
  raw <- unlist(regmatches(rest, gregexpr("item[{][^}]*[}]", rest)))
  if (!length(raw)) return(character(0))
  raw <- sub("^item[{]", "", raw)
  raw <- sub("[}]$", "", raw)
  trimws(unlist(strsplit(raw, ",")))
}
.rd_aliases <- function(path) {
  txt <- readLines(path, warn = FALSE)
  a <- grep("^.alias[{]", txt, value = TRUE)
  trimws(sub("[}].*$", "", sub("^.alias[{]", "", a)))
}

.step("Checking every Rd file for duplicated arguments")
.rds <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
.bad <- character(0)
for (f in .rds) {
  nms <- .rd_arg_names(f)
  dup <- unique(nms[duplicated(nms)])
  if (length(dup)) {
    .warn(sprintf("%s duplicates: %s", basename(f), paste(dup, collapse = ", ")))
    .bad <- c(.bad, basename(f))
  }
}
if (length(.bad)) .die("Duplicated arguments remain. Backups: ", BACKUP_DIR)
.ok(sprintf("%d Rd file(s), no duplicated arguments", length(.rds)))

.step("Checking no alias appears in two Rd files")
.alias_map <- list()
for (f in .rds) for (a in .rd_aliases(f))
  .alias_map[[a]] <- c(.alias_map[[a]], basename(f))
.dupal <- Filter(function(v) length(v) > 1L, .alias_map)
if (length(.dupal)) {
  for (a in names(.dupal))
    .warn(sprintf("alias '%s' in: %s", a, paste(.dupal[[a]], collapse = ", ")))
  .die("A duplicated alias would fail R CMD check. Backups: ", BACKUP_DIR)
}
.ok(sprintf("%d alias(es) across %d file(s), all unique",
            length(.alias_map), length(.rds)))

.step("Checking every documented argument is a formal of its topic")
suppressWarnings(devtools::load_all(".", quiet = TRUE))
.mismatch <- character(0)
for (f in .rds) {
  al <- .rd_aliases(f)
  fns <- Filter(function(x) exists(x, envir = asNamespace("BetaDanish")), al)
  if (!length(fns)) next
  formals_all <- unique(unlist(lapply(fns, function(nm) {
    ob <- get(nm, envir = asNamespace("BetaDanish"))
    if (is.function(ob)) names(formals(ob)) else character(0)
  })))
  if (!length(formals_all)) next
  documented <- .rd_arg_names(f)
  extra <- setdiff(documented, formals_all)
  missing <- setdiff(setdiff(formals_all, "..."), documented)
  if (length(extra) || length(missing))
    .mismatch <- c(.mismatch, sprintf("%s: extra {%s} missing {%s}",
                                      basename(f),
                                      paste(extra, collapse = ", "),
                                      paste(missing, collapse = ", ")))
}
if (length(.mismatch)) {
  for (m in .mismatch) .warn(m)
} else {
  .ok("documented arguments match the formals everywhere")
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
cat("  PATCH 3a-fix3 COMPLETE\n")
cat(strrep("=", 78), "\n\n")
cat("  R/entropy.R removed: it held a second bd_entropy_shannon that had been\n")
cat("  silently overwritten by the closed-form version at load time.\n")
cat("  Three new guards now run before check(): duplicate definitions across\n")
cat("  R/, duplicated Rd arguments, and duplicated Rd aliases.\n\n")
cat("  Expected: 0 errors, 0 warnings, 0-1 notes -- Patch 3a closes.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
