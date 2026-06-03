# Shannon Entropy of the Beta-Danish Distribution

Computes the differential Shannon entropy \\H(f) = -\int_0^\infty f(t)
\log f(t)\\ dt\\ for the four-parameter Beta-Danish distribution by
adaptive Gauss-Kronrod quadrature on the log-pdf.

## Usage

``` r
bd_entropy_shannon(a, b, c, k, subdivisions = 2000, rel.tol = 1e-08)
```

## Arguments

- a, b, c, k:

  Positive parameters of the Beta-Danish distribution.

- subdivisions, rel.tol:

  Passed to
  [`stats::integrate`](https://rdrr.io/r/stats/integrate.html).

## Value

Scalar Shannon entropy (in nats); `NA_real_` on integration failure.

## Examples

``` r
bd_entropy_shannon(a = 1.5, b = 2.5, c = 2, k = 1)
#> [1] 1.736391
```
