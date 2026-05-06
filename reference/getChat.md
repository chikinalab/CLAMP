# Compute Chat matrix from prior annotation

Computes the transformation matrix `Chat` used to map from observed data
to latent space, based on a pseudo-inverse of the prior annotation
matrix. Optionally standardizes the columns of `priorMat` before
computing.

## Usage

``` r
getChat(priorMat, scale = TRUE)
```

## Arguments

- priorMat:

  A numeric or sparse matrix (features x pathways) containing prior
  annotations.

- scale:

  Logical; if `TRUE` (default), standardizes the columns of `priorMat`
  before computing `Chat`.

## Value

A numeric matrix `Chat` of dimensions (pathways x features).

## Examples

``` r
# simple toy prior: 3 features x 2 pathways
priorMat <- matrix(
    c(
        1, 0, 1,
        0, 1, 0
    ),
    nrow = 3, ncol = 2,
    dimnames = list(
        paste0("gene", seq_len(3)),
        paste0("path", seq_len(2))
    )
)
# compute Chat (2 pathways x 3 features)
Chat <- getChat(priorMat)
#> Inverting...
#> done
```
