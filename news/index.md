# Changelog

## CLAMP 0.99.12

### Documentation

- Completed NEWS entries for earlier submission versions ahead of the
  Bioconductor 3.24 release.

## CLAMP 0.99.11

### Bug fixes

- Kept a single documentation URL in DESCRIPTION to avoid a merged
  hyperlink on the Bioconductor package page.
- Replaced unavailable Bioconductor release badges with “Not yet
  released” and clarified installation from Bioconductor devel in the
  README.

## CLAMP 0.99.10

### Improvements

- Updated the package hex logo to use a fully opaque hexagon background
  and refreshed `man/figures/clamp.png` and `man/figures/logo.png`
  accordingly.

## CLAMP 0.99.9

### Documentation

- Added Alexander McKim as a contributor.

## CLAMP 0.99.8

### Bug fixes

- Fixed whole-blood examples in the introductory vignette and README to
  use the supplied normalized expression values directly, avoiding
  redundant CPM normalization and log2 transformation that caused the
  pkgdown build to fail during
  [`CLAMPfull()`](https://chikinalab.org/CLAMP/reference/CLAMPfull.md).

## CLAMP 0.99.7

### Improvements

- Added `log2_transform = TRUE` to
  [`preprocessCLAMP()`](https://chikinalab.org/CLAMP/reference/preprocessCLAMP.md),
  [`preprocessCLAMPFBM()`](https://chikinalab.org/CLAMP/reference/preprocessCLAMPFBM.md),
  and
  [`cleanFBM()`](https://chikinalab.org/CLAMP/reference/cleanFBM.md).
  Set it to `FALSE` to skip log2 transformation while retaining
  missing-value handling and filtering.

### Bug fixes

- Changed
  [`preprocessCLAMP()`](https://chikinalab.org/CLAMP/reference/preprocessCLAMP.md)
  to use population variance, matching
  [`preprocessCLAMPFBM()`](https://chikinalab.org/CLAMP/reference/preprocessCLAMPFBM.md)
  for variance filtering and downstream scaling.

## CLAMP 0.99.6

### Bug fixes

- Corrected the pathway dot-plot heading and Visualization examples in
  the introductory vignette so all six plots are rendered.

## CLAMP 0.99.5

### Improvements

- Added automated BiocCheck and cross-platform R CMD check workflows.
- Added automatic pkgdown deployment after updates to `devel`.

## CLAMP 0.99.4

### Bug fixes

- Made
  [`preprocessCLAMP()`](https://chikinalab.org/CLAMP/reference/preprocessCLAMP.md)
  conditionally log2-transform and replace missing values before
  filtering, matching
  [`preprocessCLAMPFBM()`](https://chikinalab.org/CLAMP/reference/preprocessCLAMPFBM.md).
- Made
  [`CLAMPplotU()`](https://chikinalab.org/CLAMP/reference/CLAMPplotU.md)
  handle missing FDR values when no pathways pass the requested
  thresholds.

## CLAMP 0.99.3

### New features

- Added
  [`differentialLVActivity()`](https://chikinalab.org/CLAMP/reference/differentialLVActivity.md)
  to test latent variable activity between sample groups.

### Improvements

- Added input validation with informative error messages for argument
  types and ranges across the solvers and utility functions.
- Documented the source (NCBI GEO GSE130824) and processing of the
  `dataWholeBlood` example dataset.
- Replaced `hdf5r` with `rhdf5` for HDF5 access.
- Updated the introductory vignette to load data with the CLAMPData
  reader functions and to summarize the model results.

## CLAMP 0.99.2

### Other changes

- Build housekeeping only; no user-visible changes.

## CLAMP 0.99.1

### New features

- Renamed the package API from PLIER2 to CLAMP:
  [`CLAMPbase()`](https://chikinalab.org/CLAMP/reference/CLAMPbase.md),
  [`CLAMPfull()`](https://chikinalab.org/CLAMP/reference/CLAMPfull.md),
  [`CLAMPfullnVP()`](https://chikinalab.org/CLAMP/reference/CLAMPfullnVP.md),
  [`projectCLAMP()`](https://chikinalab.org/CLAMP/reference/projectCLAMP.md),
  [`preprocessCLAMP()`](https://chikinalab.org/CLAMP/reference/preprocessCLAMP.md),
  [`cpmCLAMP()`](https://chikinalab.org/CLAMP/reference/cpmCLAMP.md),
  [`zscoreCLAMP()`](https://chikinalab.org/CLAMP/reference/zscoreCLAMP.md)
  and their FBM variants.
- Added plotting functions
  [`CLAMPplotU()`](https://chikinalab.org/CLAMP/reference/CLAMPplotU.md),
  [`CLAMPdotplot()`](https://chikinalab.org/CLAMP/reference/CLAMPdotplot.md),
  [`CLAMPdotplotAll()`](https://chikinalab.org/CLAMP/reference/CLAMPdotplotAll.md),
  [`CLAMPplotTopZ()`](https://chikinalab.org/CLAMP/reference/CLAMPplotTopZ.md)
  and
  [`plotTopZ_Complex()`](https://chikinalab.org/CLAMP/reference/plotTopZ_Complex.md).
- Added SVD and rank selection helpers, including the Gavish-Donoho
  method:
  [`compute_svd()`](https://chikinalab.org/CLAMP/reference/compute_svd.md),
  [`select_svd_k()`](https://chikinalab.org/CLAMP/reference/select_svd_k.md)
  and
  [`select_clamp_k()`](https://chikinalab.org/CLAMP/reference/select_clamp_k.md).
- Exported
  [`read_gmt()`](https://chikinalab.org/CLAMP/reference/read_gmt.md) for
  reading pathway GMT files.
- Added `keep_samples_idx` to
  [`filterFBM()`](https://chikinalab.org/CLAMP/reference/filterFBM.md).
- Added a vignette comparing CLAMPbase and CLAMPfull.

### Improvements

- Moved the Z update to C++ via Rcpp.
- Added `src/Makevars` and `src/Makevars.win` for cross-platform
  BLAS/LAPACK linking.

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
