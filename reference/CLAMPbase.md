# CLAMP base matrix factorization

Runs the core matrix factorization procedure of CLAMP, decomposing the
gene expression matrix `Y` into latent variables `Z` and loadings `B`.
It supports sparse, dense, and Filebacked Big Matrices (FBM) as input
and includes options for adaptive sparsity, positive constraints, and
regularization.

## Usage

``` r
CLAMPbase(
  Y,
  clamp_k = NULL,
  svd_k = NULL,
  svdres = NULL,
  L1 = NULL,
  L2 = NULL,
  Zpos = TRUE,
  max.iter = 200,
  tol = 5e-04,
  trace = FALSE,
  rseed = NULL,
  B = NULL,
  scale = 1,
  pos.adj = 3,
  adaptive.p = 0.05,
  adaptive.iter = 20,
  cutoff = 0,
  ncores = 1,
  clamp_k_method = c("elbow", "permutation", "gavish_donoho", "scaleSVs")
)
```

## Arguments

- Y:

  Input gene expression matrix (genes x samples). Can be dense, sparse
  (`dgCMatrix`), or FBM.

- clamp_k:

  Number of latent variables for CLAMP (final model rank). If `NULL`, it
  is chosen automatically via
  [`select_clamp_k()`](https://chikinalab.org/CLAMP/reference/select_clamp_k.md).

- svd_k:

  Number of singular values/components to compute in the SVD. If `NULL`,
  defaults to `max(2, min(n_genes, n_samples) - 1)`.

- svdres:

  Optional precomputed SVD result. If not supplied, it is computed
  internally.

- L1:

  L1 regularization strength for Z. Defaults to scaled singular value.

- L2:

  L2 regularization strength for B. Defaults to scaled singular value.

- Zpos:

  Logical; if `TRUE`, negative entries in Z are zeroed. Default is
  `TRUE`.

- max.iter:

  Maximum number of optimization iterations. Default is 200.

- tol:

  Convergence tolerance for B update. Default is 5e-4.

- trace:

  Logical; if `TRUE`, prints progress. Default is `FALSE`.

- rseed:

  Optional integer for reproducible random initialization of B.

- B:

  Optional initial matrix for B. If not provided, initialized from SVD.

- scale:

  Scaling factor for L1 and L2 when not provided. Default is 1.

- pos.adj:

  Positive constraint adjustment divisor for L1. Default is 3.

- adaptive.p:

  Controls adaptive sparsity in `Z`. After each ALS update, negative
  entries in `Z` are assumed to reflect noise. The cutoff for
  thresholding is set according to the probability of positive values
  under a reflected negative distribution-effectively zeroing out small
  positive entries likely to be noise. Smaller values lead to more
  sparsity. Default is 0.05.

- adaptive.iter:

  Number of iterations before adaptive sparsity is applied. Default is
  20.

- cutoff:

  Scalar threshold to zero Z values when `Zpos = TRUE` and adaptive
  thresholding is not used. Default is 0.

- ncores:

  Number of cores to use for parallel computation (only used if Y is an
  FBM). Default is 1.

- clamp_k_method:

  Method for selecting `clamp_k` when not provided. One of `"elbow"`
  (default), `"permutation"`, `"gavish_donoho"`, or `"scaleSVs"`. Passed
  to
  [`select_clamp_k()`](https://chikinalab.org/CLAMP/reference/select_clamp_k.md).

## Value

A list with components:

- `B`:

  Latent variable loadings (LVs x genes)

- `Z`:

  Latent variable scores (LVs x samples)

- `Zraw`:

  Raw Z matrix before thresholding

- `L1`:

  Final value of L1 used

- `L2`:

  Final value of L2 used

## Details

This function is the low-level implementation of CLAMP. It alternates
between solving for `Z` given `B` and solving for `B` given `Z`, with
optional sparsity and non-negativity constraints on `Z`. Convergence is
assessed via relative change in `B`.

## Examples

``` r
# small toy dataset: 5 genes x 4 samples
Y <- matrix(rnorm(5 * 4), nrow = 5, ncol = 4)
# run a single iteration for speed
res <- CLAMPbase(Y, clamp_k = 2, max.iter = 1, trace = FALSE)
#> ****
#> Computing SVD
#> CLAMP k is set to 2
#> L1 is set to 0.926903271141505
#> L2 is set to 2.78070981342452
# inspect dimensions of B and Z
dim(res$B)
#> [1] 2 4
dim(res$Z)
#> [1] 5 2
```
