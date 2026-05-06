# Filter rows of a Filebacked Big Matrix based on mean and variance

Filters an FBM based on row-level mean and variance thresholds,
returning a new FBM with only the selected rows.

## Usage

``` r
filterFBM(
  fbm,
  rowStats,
  keep_samples_idx = NULL,
  mean_cutoff = NULL,
  var_cutoff = NULL,
  backingfile = "filtered_fbm"
)
```

## Arguments

- fbm:

  A
  [`bigstatsr::FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
  object.

- rowStats:

  A list with numeric vectors `row_means` and `row_variances`.

- keep_samples_idx:

  Optional integer vector of column indices to retain. Default is
  `filtered_fbm`.

- mean_cutoff:

  Optional minimum mean threshold; rows with means below this are
  removed.

- var_cutoff:

  Optional minimum variance threshold; rows with variances below this
  are removed.

- backingfile:

  A character string specifying the filename (without extension) for the
  new FBM.

## Value

A list with:

- fbm_filtered:

  A new FBM object containing only filtered rows.

- kept_rows:

  Indices of rows retained in the filtering step.

## Details

This function creates a new FBM and copies over only the rows that pass
the filtering criteria. The original FBM is unchanged.

## Examples

``` r
fbm <- bigstatsr::FBM(5, 3, init = matrix(rnorm(15), nrow = 5))
rs <- list(
    row_means = rowMeans(fbm[]),
    row_variances = apply(fbm[], 1, var)
)
out <- filterFBM(fbm, rs,
    mean_cutoff = -0.2, var_cutoff = 0.5,
    backingfile = tempfile()
)
```
