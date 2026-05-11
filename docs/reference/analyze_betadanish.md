# Comprehensive Beta-Danish Analysis Pipeline

Runs a complete end-to-end analysis: reads data, fits the 4-parameter
and 3-parameter models, compares them, benchmarks against standard
distributions, and generates diagnostic plots.

## Usage

``` r
analyze_betadanish(file, time_col, status_col = NULL)
```

## Arguments

- file:

  Path to the CSV or Excel file containing the data.

- time_col:

  Name of the time column.

- status_col:

  Name of the status column (optional).
