# Compare Nested Beta-Danish Models

Performs a Likelihood Ratio Test (LRT) between the 4-parameter full
model and the 3-parameter submodel.

## Usage

``` r
compare_models(full_model, sub_model)
```

## Arguments

- full_model:

  A fitted 4-parameter \`betadanish\` object.

- sub_model:

  A fitted 3-parameter \`betadanish\` object.

## Value

A data frame with the test statistic and p-value.
