# Runs the streamlined full CLAMP model.

This version performs latent-variable decomposition of a gene expression
matrix `Y` guided by prior pathway annotations `priorMat`, with
simplified and lighter regularization compared to the original extended
CLAMP variant. The algorithm alternates updates of `Z`, `B`, and `U`,
where `U` captures pathway-latent variable associations inferred
directly from the data without ridge-regularized projections (`Chat` is
not used).

## Usage

``` r
CLAMPfull(
  Y,
  priorMat,
  Chat = NULL,
  svdres = NULL,
  clamp.base.result = NULL,
  clamp_k = NULL,
  svd_k = NULL,
  L1 = NULL,
  L2 = NULL,
  cvn = 5,
  max.iter = 30,
  trace = TRUE,
  maxPath = 10,
  doCrossval = TRUE,
  penalty.factor = rep(1, ncol(priorMat)),
  glm_alpha = 0.9,
  minGenes = 0,
  tol = 5e-04,
  seed = 123456,
  allGenes = FALSE,
  rseed = NULL,
  max.U.updates = Inf,
  pathwaySelection = c("fast", "complete"),
  multiplier = 5,
  adaptive.p = 0.05,
  useNNLS = TRUE,
  useRaw = TRUE,
  refitEvery = 3,
  useSE = FALSE,
  var.prior = TRUE,
  Uscale = FALSE,
  robust.vp = TRUE,
  use_cpp = FALSE,
  clamp_k_method = c("elbow", "permutation", "gavish_donoho", "scaleSVs")
)
```

## Arguments

- Y:

  Gene expression matrix (genes x samples). Can be dense, sparse
  (dgCMatrix), or FBM.

- priorMat:

  Binary or weighted prior matrix (genes x pathways) linking genes to
  pathways.

- Chat:

  Ignored in this version (kept for interface compatibility).

- svdres:

  Optional precomputed SVD result for initialization.

- clamp.base.result:

  Optional result from
  [`CLAMPbase()`](https://chikinalab.org/CLAMP/reference/CLAMPbase.md)
  providing initial values.

- clamp_k:

  Number of latent variables for CLAMP (final model rank). If `NULL`, it
  is chosen automatically via
  [`select_clamp_k()`](https://chikinalab.org/CLAMP/reference/select_clamp_k.md).

- svd_k:

  Number of singular values/components to compute in the SVD. If `NULL`,
  defaults to `max(2, min(n_genes, n_samples) - 1)`.

- L1, L2:

  Regularization parameters for `Z` and `B`. Defaults use values from
  `clamp.base.result`.

- cvn:

  Number of folds for pathway-level cross-validation. Default: 5.

- max.iter:

  Maximum number of outer iterations. Default: 30.

- trace:

  Logical; print iteration progress. Default: `TRUE`.

- maxPath:

  Maximum number of pathways per LV during U-fitting. Default: 10.

- doCrossval:

  Whether to mask prior entries for CV evaluation. Default: `TRUE`.

- penalty.factor:

  Optional vector of per-pathway penalties for glmnet. Default: 1.

- glm_alpha:

  Elastic net mixing parameter for U estimation. Default: 0.9.

- minGenes:

  Minimum number of genes per pathway. Default: 0.

- tol:

  Convergence tolerance for B updates. Default: 5e-4.

- seed:

  Random seed for CV masking. Default: 123456.

- allGenes:

  If `TRUE`, adds zero-filled rows for missing genes. Default: `FALSE`.

- rseed:

  Reproducibility, coordinate descent updates are done in random order.

- max.U.updates:

  Maximum number of U updates (capped by max.iter).

- pathwaySelection:

  Pathway selection mode for U fitting (`"fast"` or `"complete"`).

- multiplier:

  Variance-prior scaling factor. Default: 5.

- adaptive.p:

  Quantile of negative Z values used to define adaptive thresholding.
  Default: 0.05.

- useNNLS:

  Whether to use non-negative least squares for U estimation. Default:
  `TRUE`.

- useRaw:

  If `TRUE`, uses unthresholded Z in U updates. Default: `TRUE`.

- refitEvery:

  Frequency (in U updates) of full refits. Default: 3.

- useSE:

  Logical; whether to use the 1-standard-error rule for internal glmnet
  fitting. Default is FALSE.

- var.prior:

  Logical; if `TRUE`, enables adaptive variance prior updates for Z.
  Default: `TRUE`.

- Uscale:

  Logical; whether to scale U columns. Default: `FALSE`.

- robust.vp:

  Logical; winsorize prior-predicted Z2 values to reduce outlier
  effects. Default: `TRUE`.

- use_cpp:

  Logical; if TRUE, use C++ implementation for Z updates. Default is
  FALSE.

- clamp_k_method:

  Method for selecting `clamp_k` when not provided. One of `"elbow"`
  (default), `"permutation"`, `"gavish_donoho"`, or `"scaleSVs"`. Passed
  to
  [`select_clamp_k()`](https://chikinalab.org/CLAMP/reference/select_clamp_k.md).

## Value

A list with elements:

- `B`:

  LV loadings on samples (k × samples)

- `Z`:

  Gene loadings (genes × k)

- `U`:

  Pathway loadings (pathways × k)

- `C`:

  Prior matrix used during training (masked if CV)

- `Z2`:

  Predicted Z from pathway priors

- `heldOutGenes`:

  Held-out gene lists per pathway (CV mode)

- `Uauc, Up, summary`:

  CV evaluation metrics if CV is enabled

- `priorMat, priorMatCV`:

  Final and masked prior matrices

- `call`:

  Function call

## Details

Cross-validation can be used to evaluate pathway-LV specificity, and a
variance-based prior (`var.prior = TRUE`) introduces adaptive shrinkage
of `Z` based on how strongly each latent component aligns with prior
pathways. The scaling factor `multiplier` (default 5) controls the
strength of this adaptive shrinkage.

This implementation omits ridge-projected priors (`Chat`) and uses a
lighter variance prior with a lower default `multiplier = 5`, allowing
more flexible latent representations. Setting `var.prior = FALSE`
reproduces standard CLAMP-like updates. Cross-validation, if enabled,
masks 20% of gene-pathway associations per column to estimate pathway-LV
specificity (reported via AUC and p-values).

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
base <- CLAMPbase(mat, clamp_k = 5, trace = FALSE, max.iter = 5)
#> ****
#> Computing SVD
#> CLAMP k is set to 5
#> L1 is set to 0.736045659723086
#> L2 is set to 2.20813697916926
prior <- matrix(sample(0:1, 10 * 6, TRUE, prob = c(0.9, 0.1)),
    nrow = 10, ncol = 6
)
fit <- CLAMPfull(
    Y = mat,
    priorMat = prior,
    clamp.base.result = base,
    doCrossval = FALSE,
    adaptive.p = 0,
    max.U.updates = 0,
    max.iter = 1,
    trace = FALSE
)
#> ** CLAMPfull **
#> using provided CLAMPbase result
#> CLAMP k is set to 5
#> L1=0.736045659723086; L2=2.20813697916926
#> Not using cross-validation. No AUCs or p-values
```
