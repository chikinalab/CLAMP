# Fit the loading matrix Z using sparse regression of prior information U

For each column of a target matrix `Z` using pathway or prior annotation
`priorMat`. It performs regularization selection using cross-validation
and can apply either Supports continuous or binary response models.
Relaxed refitting is supported for final coefficient estimation.

## Usage

``` r
solveU(
  Z,
  Chat = NULL,
  priorMat,
  penalty.factor,
  pathwaySelection = "fast",
  alpha = 0.9,
  maxPath = 10,
  nfolds = 5,
  useSE = FALSE,
  top = NULL,
  binary = FALSE,
  nlambda = 20,
  scale = TRUE,
  refit = TRUE,
  Uprev = NULL,
  useAUC = TRUE,
  intercept = TRUE,
  ...
)
```

## Arguments

- Z:

  A numeric matrix with features (rows) and samples (columns).

- Chat:

  (Optional) Precomputed pseudo-inverse of `priorMat`; if `NULL`, it is
  calculated using ridge regularization.

- priorMat:

  A numeric matrix with prior information (features x pathways).

- penalty.factor:

  Optional penalty weights for features in `priorMat`.

- pathwaySelection:

  Method to select candidate pathways: `"fast"` (default) or
  `"complete"`.

- alpha:

  Elastic net mixing parameter (0 = ridge, 1 = lasso). Default is 0.9.

- maxPath:

  Maximum number of pathways/features selected per column. Default is
  10.

- nfolds:

  Number of cross-validation folds. Default is 5.

- useSE:

  Whether to use the 1-standard-error rule for lambda selection. Default
  is `FALSE`.

- top:

  If set, sets to 0 all but the top entries of `Z` per column before
  fitting.

- binary:

  If `TRUE`, fits a binomial model (e.g., classification) to `Z`\>0. Can
  be used in combination with `top`. Default is `FALSE`.

- nlambda:

  Number of lambda values for glmnet. Default is 20.

- scale:

  Whether to standardize predictors in glmnet. Default is `TRUE`.

- refit:

  Whether to perform relaxed refitting using selected predictors.
  Default is `TRUE`.

- Uprev:

  (Optional) Previous U matrix to reuse. In this mode only the columns
  of `U` that are all zero are estimated. Used internally in `CLAMP`.

- useAUC:

  Logical; whether to compute pathway-LV associations using AUC (default
  TRUE) instead of OLS.

- intercept:

  Logical; whether to include an intercept term in glmnet models.
  Default is TRUE.

- ...:

  Additional arguments passed to `glmnet()` or `cv.glmnet()`.

## Value

A list with one element:

- `U`:

  A matrix of loadings (features x components). Columns are named `LV1`,
  `LV2`, ...

## Examples

``` r
set.seed(123)
genes <- paste0("G", 1:200)
lvs <- paste0("LV", 1:4)
paths <- paste0("Path", 1:60)

Z <- matrix(rnorm(200 * 4), nrow = 200, dimnames = list(genes, lvs))

priorMat <- matrix(rbinom(200 * 60, 1, 0.07),
    nrow = 200, dimnames = list(genes, paths)
)

fit1 <- solveU(
    Z = Z,
    priorMat = priorMat,
    pathwaySelection = "fast",
    alpha = 0.9,
    maxPath = 10,
    nfolds = 5,
    binary = FALSE,
    refit = TRUE
)
#> Warning: Passing 'dfmax' to glmnet() is deprecated. Use control = list(dfmax = ...) instead.
```
