## =============================================================================
##  BetaDanish  --  PATCH 3f-fix3 : vignette chunk style and two test signs
## =============================================================================
##
##  The self-test and the S3 registration are now correct: every plot drew, by
##  direct call and through dispatch, and NAMESPACE registers the method. Two
##  faults remain, both mine, and neither is in the package code.
##
##  1. THE INTRODUCTION VIGNETTE
##
##     BetaDanish_Introduction.Rmd shows its code in plain ```r display blocks,
##     which knitr does not execute. Its library(BetaDanish) call is inside one
##     of those, so the package is never attached while the vignette builds.
##     My appended section used ```{r} chunks, which ARE executed -- and were
##     therefore the first code ever run in that file, with nothing attached:
##
##       Error in bd_moment_summary(): could not find function
##
##     An evaluated setup chunk is inserted at the top of the appended section.
##     The rest of the vignette is untouched: converting its display blocks to
##     live chunks would be a rewrite, and a risky one during a release.
##
##  2. TWO TEST EXPECTATIONS WERE SWAPPED
##
##     Barlow-Campo: the scaled TTT transform is concave for an increasing
##     hazard and convex for a decreasing one. So convex-then-concave -- below
##     the diagonal, then above -- is a BATHTUB, and concave-then-convex is
##     UNIMODAL.
##
##       u + 0.2 sin(2 pi u)  lies above the diagonal on (0, 1/2)  -> unimodal
##       u - 0.2 sin(2 pi u)  lies below the diagonal on (0, 1/2)  -> bathtub
##
##     .bd_ttt_shape() had this right. My two expectations were the wrong way
##     round, which is why it reported "unimodal" where the test said
##     "bathtub" and vice versa. The corrected test records the reasoning so
##     the next reader does not have to re-derive it.
##
##  HOW TO RUN   source("dev/BetaDanish_Patch3f_fix3.R")
##  IDEMPOTENT   Yes -- the vignette chunk is inserted only once.
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
cat("  BetaDanish  --  Patch 3f-fix3\n")
cat(strrep("=", 78), "\n")

.step("Pre-flight")
if (!file.exists("DESCRIPTION")) .die("No DESCRIPTION here. setwd() to the package root.")
if (read.dcf("DESCRIPTION")[1, "Package"] != "BetaDanish") .die("Not the BetaDanish package.")
if (!file.exists("R/plots_extra.R")) .die("Patch 3f has not been applied.")
.ns <- readLines("NAMESPACE", warn = FALSE)
if (!any(grepl("S3method(plot,bd_bayes)", .ns, fixed = TRUE)))
  .die("NAMESPACE does not register plot.bd_bayes -- run Patch 3f-fix2 first.")
.ok("Patch 3f-fix2 detected")

BACKUP_DIR <- file.path(".betadanish_backup",
                        format(Sys.time(), "%Y%m%d-%H%M%S-patch3ffix3"))
dir.create(BACKUP_DIR, recursive = TRUE, showWarnings = FALSE)
.ok(paste("backup:", BACKUP_DIR))

## =============================================================================
##  1. VIGNETTE
## =============================================================================

.step("Attaching the package in the appended vignette section")

.vp <- "vignettes/BetaDanish_Introduction.Rmd"
if (!file.exists(.vp)) {
  .warn("introduction vignette not found; skipped")
} else {
  .v <- readLines(.vp, warn = FALSE)
  .h <- grep("^## New in 0\\.3\\.0[ ]*$", .v)
  if (!length(.h)) {
    .warn("the 'New in 0.3.0' heading was not found; check the vignette by hand")
  } else if (any(grepl("bd_intro_setup", .v, fixed = TRUE))) {
    .info("the setup chunk is already present")
  } else {
    .backup(.vp)
    .chunk <- c("",
                "```{r bd_intro_setup, include = FALSE}",
                "# The blocks above this point are ```r display blocks, which knitr does",
                "# not execute, so the package has not been attached yet.",
                "library(BetaDanish)",
                "```")
    .v <- append(.v, .chunk, after = .h[1])
    .write_lines(.vp, .v)
    .ok("evaluated setup chunk inserted after the '## New in 0.3.0' heading")
  }

  ## Report the chunk styles so the mismatch is visible rather than implicit.
  .v2 <- readLines(.vp, warn = FALSE)
  .live <- sum(grepl("^```[{]r", .v2))
  .disp <- sum(grepl("^```r[ ]*$", .v2))
  .info(sprintf("%d evaluated chunk(s), %d display-only block(s)", .live, .disp))
}

## =============================================================================
##  2. TESTS
## =============================================================================

.step("Correcting the two swapped TTT expectations")

.put("tests/testthat/test-plots-extra.R", r"---(## Visualisation added in 0.3.0. Plots are drawn to a null device; the
## assertions are about the returned values and the classification logic,
## which is what can actually be wrong.

test_that("the TTT transform has the right endpoints and is increasing", {
  data(guinea_pig, package = "BetaDanish", envir = environment())
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  ttt <- bd_ttt_plot(guinea_pig$time)

  expect_s3_class(ttt, "data.frame")
  expect_named(ttt, c("i_n", "phi"))
  expect_equal(nrow(ttt), nrow(guinea_pig))
  expect_equal(ttt$phi[nrow(ttt)], 1, tolerance = 1e-12)   # phi(1) = 1
  expect_true(all(diff(ttt$phi) >= -1e-12))                # non-decreasing
  expect_true(all(ttt$phi >= 0 & ttt$phi <= 1 + 1e-12))
  expect_true(attr(ttt, "shape") %in%
                c("increasing", "decreasing", "bathtub",
                  "unimodal (upside-down bathtub)", "constant (exponential)"))
})

test_that("the TTT classifier recognises the reference shapes", {
  cl <- BetaDanish:::.bd_ttt_shape
  u <- seq(0.02, 1, length.out = 50)
  expect_equal(cl(u, u), "constant (exponential)")
  expect_equal(cl(u, pmin(u + 0.15, 1)), "increasing")
  expect_equal(cl(u, pmax(u - 0.15, 0)), "decreasing")
  ## Barlow-Campo: the TTT is concave for an increasing hazard and convex for
  ## a decreasing one, so convex-then-concave (below the diagonal, then above)
  ## is a bathtub, and concave-then-convex is unimodal.
  ##
  ##   u + 0.2 sin(2 pi u)  is above the diagonal on (0, 1/2)  -> unimodal
  ##   u - 0.2 sin(2 pi u)  is below the diagonal on (0, 1/2)  -> bathtub
  expect_equal(cl(u, u + 0.2 * sin(2 * pi * u)),
               "unimodal (upside-down bathtub)")
  expect_equal(cl(u, u - 0.2 * sin(2 * pi * u)), "bathtub")
})

test_that("an exponential sample gives a TTT curve near the diagonal", {
  set.seed(4)
  x <- stats::rexp(400, rate = 0.5)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  ttt <- bd_ttt_plot(x)
  expect_lt(max(abs(ttt$phi - ttt$i_n)), 0.12)
})

test_that("bd_ttt_plot drops censored observations with a warning", {
  set.seed(5)
  t <- rbetadanish(60, 1, 3, 2, 0.5)
  s <- stats::rbinom(60, 1, 0.8)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_warning(ttt <- bd_ttt_plot(t, status = s), "censored")
  expect_equal(nrow(ttt), sum(s == 1))
})

test_that("bd_ttt_plot accepts a fitted object and validates its input", {
  skip_on_cran()
  dat <- simulate_bd_data(60, a = 1, b = 3, c = 2, k = 0.5, seed = 7)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 1, check_identifiability = FALSE))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_s3_class(suppressWarnings(bd_ttt_plot(fit)), "data.frame")
  expect_error(bd_ttt_plot(c(1, 2, 3)), "At least five")
})

test_that("bd_profile_plot draws and returns its input", {
  skip_on_cran()
  dat <- simulate_bd_data(100, a = 1, b = 3, c = 2, k = 0.5, seed = 8)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 1, check_identifiability = FALSE))
  p <- bd_profile_ci(fit, "b", n_grid = 12L)

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  out <- bd_profile_plot(p)
  expect_identical(out$parameter, "b")
  expect_error(bd_profile_plot(list()), "bd_profile object")
})

test_that("plot.bd_bayes validates before drawing", {
  fake <- structure(list(draws = matrix(rnorm(300), 100, 3,
                                        dimnames = list(NULL, c("b", "c", "k"))),
                         HPD = NULL, submodel = TRUE),
                    class = "bd_bayes")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_invisible(plot(fake))
  expect_invisible(plot(fake, which = "b", type = "trace"))
  expect_error(plot(fake, which = "zzz"), "Not in the posterior")

  empty <- structure(list(draws = matrix(numeric(0), 0, 0)), class = "bd_bayes")
  expect_error(plot(empty), "no posterior draws")
})

test_that("plot.bd_bayes leaves the graphical parameters as it found them", {
  fake <- structure(list(draws = matrix(rnorm(200), 100, 2,
                                        dimnames = list(NULL, c("b", "c")))),
                    class = "bd_bayes")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  before <- graphics::par("mfrow")
  plot(fake)
  expect_equal(graphics::par("mfrow"), before)
}))---")


.step("Checking the classifier against the corrected expectations, directly")
.loaded <- tryCatch({ devtools::load_all(".", quiet = TRUE); TRUE },
                    error = function(e) conditionMessage(e))
if (!isTRUE(.loaded)) .die("load_all() failed:\n  ", .loaded)
.cl <- get(".bd_ttt_shape", envir = asNamespace("BetaDanish"))
.u  <- seq(0.02, 1, length.out = 50)
.cases <- list(
  list(lab = "flat, on the diagonal",      phi = .u,
       want = "constant (exponential)"),
  list(lab = "above throughout",           phi = pmin(.u + 0.15, 1),
       want = "increasing"),
  list(lab = "below throughout",           phi = pmax(.u - 0.15, 0),
       want = "decreasing"),
  list(lab = "above then below (concave-convex)",
       phi = .u + 0.2 * sin(2 * pi * .u),
       want = "unimodal (upside-down bathtub)"),
  list(lab = "below then above (convex-concave)",
       phi = .u - 0.2 * sin(2 * pi * .u),
       want = "bathtub"))
.bad <- character(0)
for (cs in .cases) {
  got <- .cl(.u, cs$phi)
  if (identical(got, cs$want)) {
    .ok(sprintf("%-36s -> %s", cs$lab, got))
  } else {
    .warn(sprintf("%-36s -> %s (expected %s)", cs$lab, got, cs$want))
    .bad <- c(.bad, cs$lab)
  }
}
if (length(.bad))
  .die("The TTT classifier disagrees with Barlow-Campo on: ",
       paste(.bad, collapse = "; "), "\nBackups: ", BACKUP_DIR)
.ok("classifier matches the theory on all five reference shapes")

.step("Parsing all R and test files")
.targets <- c(list.files("R", pattern = "[.]R$", full.names = TRUE),
              list.files("tests", pattern = "[.]R$", recursive = TRUE, full.names = TRUE))
.badp <- character(0)
for (f in .targets) {
  e <- tryCatch({ parse(f); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(e)) .badp <- c(.badp, paste0("  ", f, ": ", e))
}
if (length(.badp)) .die("These files do not parse:\n", paste(.badp, collapse = "\n"),
                        "\n\nBackups: ", BACKUP_DIR)
.ok(sprintf("%d file(s) parse cleanly", length(.targets)))

.step("Knitting the introduction vignette on its own, before the full check")
if (requireNamespace("rmarkdown", quietly = TRUE) && file.exists(.vp)) {
  .out <- file.path(tempdir(), "bd_intro_test.html")
  .kn <- tryCatch({
    rmarkdown::render(.vp, output_file = .out, quiet = TRUE,
                      envir = new.env())
    TRUE
  }, error = function(e) conditionMessage(e))
  if (isTRUE(.kn)) {
    .ok("the introduction vignette builds")
  } else {
    .warn(paste("vignette build failed:", .kn))
    .die("Fix the vignette before running check(), which would take minutes ",
         "to report the same thing.\nBackups: ", BACKUP_DIR)
  }
} else {
  .info("rmarkdown not available; the full check will exercise the vignette")
}

.step("devtools::test()")
.t <- tryCatch(devtools::test(), error = function(e) { .warn(conditionMessage(e)); NULL })

.step("devtools::check() -- several minutes, do not interrupt")
.chk <- tryCatch(devtools::check(document = TRUE, args = "--as-cran", error_on = "never"),
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
cat("  PATCH 3f-fix3 COMPLETE\n")
cat(strrep("=", 78), "\n\n")
cat("  vignette  the appended section now attaches the package\n")
cat("  tests     the two TTT expectations corrected to match Barlow-Campo\n\n")
cat("  If this reports 0 / 0 / 0, all 46 recommendations are implemented and\n")
cat("  the package is ready. NEXT:\n\n")
cat("      source(\"dev/BetaDanish_Patch3d_release.R\")\n\n")
cat("  Then: confirm the JAMSI citation, check_win_devel(), urlchecker,\n")
cat("  push and tag, devtools::build(), submit to CRAN.\n\n")
cat("  Backups: ", BACKUP_DIR, "\n\n", sep = "")
