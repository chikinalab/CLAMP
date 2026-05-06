# Changelog

## CLAMP 0.99.0

### New features

- Initial Bioconductor submission.
- Implements prior-informed latent variable decomposition for gene
  expression.
- Supports file-backed matrices
  ([`bigstatsr::FBM`](https://privefl.github.io/bigstatsr/reference/FBM-class.html))
  for large datasets.
- Added CPM, z-score, and filtering functions for preprocessing.
- Projection of new datasets into pre-trained models.
- Vignettes with example workflows and detailed documentation.

### Improvements

- Optimized in-place operations for large matrices.
- Added cross-validation utilities and adaptive sparsity support.
- Roxygen2 documentation for all exported functions.

### Bug fixes

- Fixed edge case in CPM normalization for zero-count columns.
- Corrected NA handling in preprocessing.
