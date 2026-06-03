# Plot Diagnostics for Beta-Danish Fit

Generates diagnostic plots for a fitted Beta-Danish model, including
survival, hazard, density, PP, and QQ plots.

## Usage

``` r
# S3 method for class 'betadanish'
plot(x, type = c("survival", "hazard", "density", "pp", "qq", "all"), ...)
```

## Arguments

- x:

  A fitted \`betadanish\` object.

- type:

  Character string specifying the plot type: \`"survival"\`,
  \`"hazard"\`, \`"density"\`, \`"pp"\`, \`"qq"\`, or \`"all"\`.

- ...:

  Additional arguments passed to the base \`plot\` function.

## Value

Invisibly returns the input betadanish object. Called mainly for its
side effect of producing diagnostic plots.
