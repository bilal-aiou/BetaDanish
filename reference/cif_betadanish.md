# Compute Cumulative Incidence Function (CIF)

Computes the CIF for a specific cause from a fitted Beta-Danish
competing risks model using numerical integration.

## Usage

``` r
cif_betadanish(fit, tvec, cause_idx)
```

## Arguments

- fit:

  An object of class \`bd_competing\`.

- tvec:

  Numeric vector of times at which to evaluate the CIF.

- cause_idx:

  Integer; the specific cause to evaluate (must match one of the causes
  in the fitted model).

## Value

A numeric vector of CIF probabilities corresponding to \`tvec\`.
