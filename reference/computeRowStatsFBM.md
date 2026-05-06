# Compute row-wise sum and sum of squares for a Filebacked Big Matrix

Efficiently computes row sums and row sum-of-squares for a
[`bigstatsr::FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
using column-wise chunking, suitable for large datasets that cannot be
loaded fully into memory.

## Usage

``` r
computeRowStatsFBM(fbm, ncores = 1)
```

## Arguments

- fbm:

  A
  [`bigstatsr::FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
  object.

- ncores:

  Integer; number of cores to use for parallel operations (default 1).

## Value

A list with two numeric vectors:

- row_sums:

  Sum of each row.

- row_sums_sq:

  Sum of squares of each row.
