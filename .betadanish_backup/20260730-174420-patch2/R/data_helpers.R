#' Read and Prepare Survival Data
#'
#' Reads survival data from a CSV or Excel file and returns a clean data frame
#' ready for [fit_betadanish()]. Columns are selected by name, and covariates
#' keep their original type, so factors are not silently coerced to numbers.
#'
#' @param file Path to a `.csv`, `.xls` or `.xlsx` file.
#' @param time_col Name of the time column.
#' @param status_col Name of the event indicator column (1 = event,
#'   0 = censored). If `NULL` (default), all observations are treated as
#'   uncensored.
#' @param covar_cols Character vector of covariate columns to retain, or `NULL`.
#'
#' @return A data frame with columns `time`, `status` and any retained
#'   covariates. A `"bd_data_report"` attribute records rows read, rows
#'   dropped and the censoring proportion.
#'
#' @details
#' A covariate that happens to be named `time` or `status` is renamed with a
#' `_cov` suffix and a warning, rather than colliding with the response.
#'
#' Excel input requires the `readxl` package.
#'
#' @export
#'
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' write.csv(data.frame(survival_time = c(5, 8, 12, 16),
#'                      status = c(1, 1, 0, 1)),
#'           tmp, row.names = FALSE)
#' dat <- read_survival_data(tmp, time_col = "survival_time",
#'                           status_col = "status")
#' attr(dat, "bd_data_report")
#' unlink(tmp)
read_survival_data <- function(file, time_col, status_col = NULL,
                               covar_cols = NULL) {

  if (!is.character(file) || length(file) != 1L)
    stop("'file' must be a single file path.", call. = FALSE)
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)

  ext <- tolower(tools::file_ext(file))
  dat <- if (ext == "csv") {
    utils::read.csv(file, stringsAsFactors = FALSE, check.names = TRUE)
  } else if (ext %in% c("xls", "xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Reading Excel files requires the 'readxl' package.", call. = FALSE)
    as.data.frame(readxl::read_excel(file))
  } else {
    stop("Unsupported file extension '", ext,
         "'. Supply a .csv, .xls or .xlsx file.", call. = FALSE)
  }

  if (!nrow(dat)) stop("The file contains no data rows.", call. = FALSE)

  if (!time_col %in% names(dat))
    stop("Time column '", time_col, "' not found. Available columns: ",
         paste(names(dat), collapse = ", "), call. = FALSE)

  ## ---- build the response by name, never by position -----------------------
  clean <- data.frame(time = as.numeric(dat[[time_col]]))

  if (is.null(status_col)) {
    message("No status column supplied; treating all observations as ",
            "uncensored (status = 1).")
    clean$status <- 1
  } else {
    if (!status_col %in% names(dat))
      stop("Status column '", status_col, "' not found.", call. = FALSE)
    clean$status <- as.numeric(dat[[status_col]])
  }

  if (!is.null(covar_cols)) {
    missing_cov <- setdiff(covar_cols, names(dat))
    if (length(missing_cov))
      stop("Covariate column(s) not found: ",
           paste(missing_cov, collapse = ", "), call. = FALSE)

    extra <- dat[, covar_cols, drop = FALSE]
    clash <- intersect(names(extra), c("time", "status"))
    if (length(clash)) {
      warning("Covariate(s) named ", paste(clash, collapse = ", "),
              " collide with the response; renamed with a '_cov' suffix.",
              call. = FALSE)
      names(extra)[names(extra) %in% clash] <-
        paste0(names(extra)[names(extra) %in% clash], "_cov")
    }
    clean <- cbind(clean, extra)
  }

  n_read <- nrow(clean)
  clean  <- clean[stats::complete.cases(clean), , drop = FALSE]
  n_kept <- nrow(clean)

  if (n_kept < n_read)
    warning("Dropped ", n_read - n_kept, " row(s) with missing values.",
            call. = FALSE)
  if (!n_kept) stop("No complete rows remain after removing missing values.",
                    call. = FALSE)

  if (any(clean$time <= 0))
    warning("Some times are <= 0; survival models require strictly positive ",
            "times.", call. = FALSE)
  if (!all(clean$status %in% c(0, 1)))
    stop("The status column contains values other than 0 and 1. Recode so ",
         "that 1 = event and 0 = censored.", call. = FALSE)

  row.names(clean) <- NULL
  attr(clean, "bd_data_report") <- list(
    file            = normalizePath(file, mustWork = FALSE),
    rows_read       = n_read,
    rows_kept       = n_kept,
    rows_dropped    = n_read - n_kept,
    n_events        = sum(clean$status == 1),
    censoring_prop  = mean(clean$status == 0),
    covariates      = setdiff(names(clean), c("time", "status"))
  )
  clean
}
