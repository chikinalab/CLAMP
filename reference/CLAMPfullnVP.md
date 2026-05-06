# Full CLAMP model with prior information and cross-validation

Runs the full CLAMP model using a gene expression matrix and prior
pathway annotation matrix. This function performs latent variable
decomposition guided by prior knowledge and includes optional
cross-validation to evaluate pathway associations.

## Usage

``` r
CLAMPfullnVP(
  Y,
  priorMat,
  svdres = NULL,
  clamp.base.result = NULL,
  clamp_k = NULL,
  svd_k = NULL,
  L1 = NULL,
  L2 = NULL,
  top = NULL,
  cvn = 5,
  max.iter = 350,
  trace = FALSE,
  Chat = NULL,
  maxPath = 10,
  doCrossval = TRUE,
  penalty.factor = rep(1, ncol(priorMat)),
  glm_alpha = 0.9,
  minGenes = 10,
  tol = 5e-04,
  seed = 123456,
  allGenes = FALSE,
  rseed = NULL,
  max.U.updates = 5,
  pathwaySelection = c("fast"),
  multiplier = 1,
  adaptive.p = 0.05,
  useNNLS = TRUE,
  useRaw = TRUE,
  refitAll = FALSE,
  useSE = FALSE,
  ncores = 1,
  clamp_k_method = "elbow"
)
```

## Arguments

- Y:

  Gene expression matrix (genes x samples). Can be dense, sparse
  (dgCMatrix), or FBM.

- priorMat:

  Binary matrix (genes x pathways) representing prior annotations.

- svdres:

  Optional SVD result used for initialization.

- clamp.base.result:

  Optional result from
  [`CLAMPbase()`](https://chikinalab.github.io/CLAMP/reference/CLAMPbase.md)
  to initialize B.

- clamp_k:

  Number of latent variables for CLAMP (final model rank). If `NULL`, it
  is chosen automatically via
  [`select_clamp_k()`](https://chikinalab.github.io/CLAMP/reference/select_clamp_k.md).

- svd_k:

  Number of singular values/components to compute in the SVD. If `NULL`,
  defaults to `max(2, min(n_genes, n_samples) - 1)`.

- L1:

  Regularization strength for Z. If `NULL`, initialized from SVD or
  `clamp.base.result`.

- L2:

  Regularization strength for B. If `NULL`, initialized from SVD or
  `clamp.base.result`.

- top:

  If set, keeps only top-n values per column in Z during U updates.

- cvn:

  Number of folds for cross-validation in U updates. Default is 5.

- max.iter:

  Maximum number of iterations. Default is 350.

- trace:

  Logical; if `TRUE`, prints iteration progress.

- Chat:

  Optional precomputed matrix for solving U.

- maxPath:

  Maximum number of pathways/features selected per LV. Default is 10.

- doCrossval:

  Whether to perform pathway-level cross-validation. Default is `TRUE`.

- penalty.factor:

  Vector of feature-specific penalties for glmnet. Default: all ones.

- glm_alpha:

  Elastic net mixing parameter for glmnet. Default is 0.9.

- minGenes:

  Minimum number of genes per pathway to retain. Default is 10.

- tol:

  Convergence tolerance on relative change in B. Default is 5e-4.

- seed:

  Seed for reproducibility of cross-validation masking. Default is
  123456.

- allGenes:

  If `TRUE`, zero-fills `priorMat` for genes not present. Default is
  `FALSE`.

- rseed:

  Optional seed for randomly reinitializing B and Z.

- max.U.updates:

  Maximum number of U updates. Default is 5.

- pathwaySelection:

  Pathway selection mode: `"fast"` or `"complete"`.

- multiplier:

  Scaling factor for adjusting L1 and L2.

- adaptive.p:

  Quantile threshold for adaptively zeroing small Z values. After each
  update, the `adaptive.p`-quantile of negative entries in Z is used
  (flipped positive) to threshold small positive values, assuming they
  reflect noise. Default is 0.05.

- useNNLS:

  If `TRUE`, uses non-negative least squares in U estimation. Default is
  `TRUE`.

- useRaw:

  If `TRUE`, uses unthresholded Z for solving U. Default is `TRUE`.

- refitAll:

  If `TRUE`, refits all U columns every update. Default is `FALSE`.

- useSE:

  Logical; passed to the internal
  [`solveU()`](https://chikinalab.github.io/CLAMP/reference/solveU.md)
  call. If `TRUE`, enables standard-error-aware selection when fitting U
  (pathway coefficients). Default is `FALSE`.

- ncores:

  Number of cores to use for parallel computation (only used if Y is an
  FBM). Default is 1.

- clamp_k_method:

  Method for selecting `clamp_k` when not provided. One of `"elbow"`
  (default), `"permutation"`, `"gavish_donoho"`, or `"scaleSVs"`. Passed
  to
  [`select_clamp_k()`](https://chikinalab.github.io/CLAMP/reference/select_clamp_k.md).

## Value

A list with the following components:

- `B`:

  Latent variable loadings (LVs x genes)

- `Z`:

  Latent variable matrix (LVs x samples)

- `U`:

  Pathway loadings matrix (pathways x LVs)

- `C`:

  Masked prior matrix used for training

- `L1`, `L2`:

  Regularization parameters

- `heldOutGenes`:

  List of held-out genes per pathway (if CV is enabled)

- `Uauc`:

  AUC matrix from CV evaluation (if enabled)

- `Up`:

  \-`log10(p)` values from CV evaluation (if enabled)

- `summary`:

  Data frame of AUC, p-values, and FDR per pathway x LV (if enabled)

- `priorMatCV`:

  Masked prior matrix used during CV

- `priorMat`:

  Final filtered prior matrix

- `withPrior`:

  Indices of LVs with non-zero pathway loadings

- `call`:

  Function call

## Details

The model alternates between solving `Z`, `B`, and `U`. Adaptive
sparsity is applied to `Z` using a dynamic threshold based on the
negative tail of its distribution. Cross-validation is used to hold out
gene annotations in `priorMat` and evaluate latent variable specificity.

## Examples

``` r
mat <- matrix(rnorm(100), 10, 10)
svdres <- rsvd::rsvd(mat, k = 5)
base <- CLAMPbase(Y = mat, clamp_k = 5, svdres = svdres, trace = FALSE)
#> ****
#> CLAMP k is set to 5
#> L1 is set to 0.903063161872916
#> L2 is set to 2.70918948561875
#> Converged at iteration= 31 | Bdiff=0.000000,  tol=0.000500     
priorMat <- matrix(1, nrow(mat), 5)
full <- CLAMPfullnVP(
    Y = mat, priorMat = priorMat, svdres = svdres,
    clamp.base.result = base, clamp_k = 5,
    doCrossval = FALSE, trace = FALSE, max.U.updates = 0
)
#> **CLAMPfullnVP v2 **
#> using provided CLAMPbase result
#> CLAMP k is set to 5
#> L1=0.903063161872916; L2=2.70918948561875
#> 
Converged at 8 / 350 | Bdiff=0.000000, minCor=1.000000
#> Not using cross-validation. No AUCs or p-values
```
