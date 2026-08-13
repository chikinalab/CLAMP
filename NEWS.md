# CLAMP 0.99.5

## Improvements
- Added automated BiocCheck and cross-platform R CMD check workflows.
- Added automatic pkgdown deployment after updates to `devel`.

# CLAMP 0.99.4

## Bug fixes
- Made `preprocessCLAMP()` conditionally log2-transform and replace missing
  values before filtering, matching `preprocessCLAMPFBM()`.
- Made `CLAMPplotU()` handle missing FDR values when no pathways pass the
  requested thresholds.

# CLAMP 0.99.0

## New features
- Initial Bioconductor submission.
- Implements prior-informed latent variable decomposition for gene expression.
- Supports file-backed matrices (`bigstatsr::FBM`) for large datasets.
- Added CPM, z-score, and filtering functions for preprocessing.
- Projection of new datasets into pre-trained models.
- Vignettes with example workflows and detailed documentation.

## Improvements
- Optimized in-place operations for large matrices.
- Added cross-validation utilities and adaptive sparsity support.
- Roxygen2 documentation for all exported functions.

## Bug fixes
- Fixed edge case in CPM normalization for zero-count columns.
- Corrected NA handling in preprocessing.
