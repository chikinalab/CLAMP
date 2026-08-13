# Cross-validation AUC for CLAMP latent variables and pathways

Evaluates how well each latent variable in a CLAMP model captures
held-out pathway annotations, using cross-validation over the prior
matrix. For each latent variable and associated pathway, held-out genes
are selected and the AUC is computed using their scores in `clampRes$Z`.

## Usage

``` r
crossVal(clampRes, priorMat, priorMatcv)
```

## Arguments

- clampRes:

  A list containing `U` (loadings) and `Z` (scores) from a CLAMP model.

- priorMat:

  A binary matrix (genes x pathways) indicating original pathway
  annotations.

- priorMatcv:

  A version of `priorMat` used to mask held-out annotations for
  cross-validation.

## Value

A list with:

- `Uauc`:

  Matrix of AUC values (pathways x LVs)

- `Upval`:

  Matrix of `-log10(p)` values (pathways x LVs)

- `summary`:

  Data frame with pathway, LV index, AUC, p-value, and FDR
