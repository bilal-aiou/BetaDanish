# Read and Prepare Survival Data

A helper function to read survival data from CSV or Excel files and
prepare it for analysis with the Beta-Danish package. It automatically
handles missing status columns by assuming all observations are complete
(uncensored).

## Usage

``` r
read_survival_data(file, time_col, status_col = NULL, covar_cols = NULL)
```

## Arguments

- file:

  Character string specifying the path to the file.

- time_col:

  Character string specifying the name of the time/survival column.

- status_col:

  Character string specifying the name of the event/censoring indicator
  column. If \`NULL\` (default), the function assumes all observations
  are uncensored and creates a status column filled with 1s.

- covar_cols:

  Character vector specifying the names of covariate columns to keep. If
  \`NULL\` (default), no covariates are kept.

## Value

A clean \`data.frame\` containing the \`time\`, \`status\`, and any
specified covariates, ready to be passed to \`fit_betadanish()\`.

## Details

The function checks the file extension to determine how to read the
data. For \`.xlsx\` or \`.xls\` files, the \`readxl\` package must be
installed. Missing values (\`NA\`) in the specified columns will cause
those rows to be dropped with a warning.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming a CSV file exists at "data.csv"
# Read uncensored data
dat <- read_survival_data("data.csv", time_col = "survival_time")

# Read censored data with covariates
dat <- read_survival_data("data.csv",
                          time_col = "time",
                          status_col = "event",
                          covar_cols = c("age", "treatment"))
} # }
```
