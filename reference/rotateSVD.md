# Rotate SVD components to make dominant directions positive

Ensures consistency in SVD output by flipping signs so that each left
singular vector has a majority of positive entries.

## Usage

``` r
rotateSVD(svdres)
```

## Arguments

- svdres:

  A list as returned by [`svd()`](https://rdrr.io/r/base/svd.html), with
  components `u`, `d`, and `v`.

## Value

A modified `svd`-result list where each column of `$u` has been
sign-flipped so that its entries sum to a nonnegative value; `$v` is
flipped correspondingly.
