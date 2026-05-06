# Select default number of CLAMP latent variables from an SVD

Chooses the default `clamp_k` used by the CLAMP solvers when the user
does not provide one, and returns the scale used for downstream L1/L2
regularization. Multiple methods are available via `method`:

## Usage

``` r
select_clamp_k(
  svdres,
  n_samples,
  svd_k,
  method = c("elbow", "permutation", "gavish_donoho", "scaleSVs"),
  data = NULL,
  B = 20
)
```

## Arguments

- svdres:

  An SVD result with a `d` component (output of `compute_svd`).

- n_samples:

  Integer number of samples in the original matrix (i.e. `ncol(Y)`).
  Used by `"scaleSVs"` and `"gavish_donoho"`.

- svd_k:

  Integer upper bound on `clamp_k` (number of components actually
  computed in the SVD).

- method:

  One of `"elbow"`, `"permutation"`, `"gavish_donoho"`, `"scaleSVs"`.
  Defaults to `"elbow"`.

- data:

  Raw data matrix. Required for `"permutation"` (row-normalized
  internally) and `"gavish_donoho"` (used for `n_genes`).

- B:

  Number of permutations for `"permutation"`.

## Value

An integer: the selected number of latent variables.

## Details

- `"elbow"` (default):

  Elbow heuristic on the singular-value spectrum via
  `num.pc(svdres, method = "elbow")`. `scale = svdres$d[clamp_k]`.

- `"permutation"`:

  Permutation test via `num.pc(data, method = "permutation", B = B)`.
  Requires the raw row-normalized `data` matrix.
  `scale = svdres$d[clamp_k]`.

- `"gavish_donoho"`:

  Gavish-Donoho optimal singular-value threshold via
  [`PCAtools::chooseGavishDonoho()`](https://rdrr.io/pkg/PCAtools/man/chooseGavishDonoho.html).
  Requires the raw `data` matrix (used for `n_genes`).
  `scale = svdres$d[clamp_k]`.

- `"scaleSVs"`:

  Previous behavior:
  [`getScaleFromSVs()`](https://chikinalab.github.io/CLAMP/reference/getScaleFromSVs.md)
  linear-tail fit, `clamp_k <- min(floor(k * 1.5), svd_k)`, scale from
  the fit.

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
svdres <- compute_svd(Y, k = 25)
select_clamp_k(svdres, n_samples = ncol(Y), svd_k = 25)
#> [1] 22
```
