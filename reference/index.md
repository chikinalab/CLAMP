# Package index

## Main CLAMP functions

- [`CLAMPbase()`](https://chikinalab.github.io/CLAMP/reference/CLAMPbase.md)
  : CLAMP base matrix factorization
- [`CLAMPfull()`](https://chikinalab.github.io/CLAMP/reference/CLAMPfull.md)
  : Runs the streamlined full CLAMP model.
- [`CLAMPfullnVP()`](https://chikinalab.github.io/CLAMP/reference/CLAMPfullnVP.md)
  : Full CLAMP model with prior information and cross-validation
- [`projectCLAMP()`](https://chikinalab.github.io/CLAMP/reference/projectCLAMP.md)
  : Project new data into CLAMP latent space

## SVD and rank selection

- [`select_svd_k()`](https://chikinalab.github.io/CLAMP/reference/select_svd_k.md)
  : Select default number of components for a CLAMP solver SVD
- [`compute_svd()`](https://chikinalab.github.io/CLAMP/reference/compute_svd.md)
  : Compute a truncated SVD for a CLAMP input matrix
- [`select_clamp_k()`](https://chikinalab.github.io/CLAMP/reference/select_clamp_k.md)
  : Select default number of CLAMP latent variables from an SVD
- [`num.pc()`](https://chikinalab.github.io/CLAMP/reference/num.pc.md) :
  Estimate number of principal components via elbow or permutation
  method
- [`getScaleFromSVs()`](https://chikinalab.github.io/CLAMP/reference/getScaleFromSVs.md)
  : Estimate noise scale from singular values with linear tail
  extrapolation

## Preprocessing

- [`cpmCLAMP()`](https://chikinalab.github.io/CLAMP/reference/cpmCLAMP.md)
  : Compute counts-per-million (CPM) for CLAMP pipelines
- [`cpmCLAMPFBM()`](https://chikinalab.github.io/CLAMP/reference/cpmCLAMPFBM.md)
  : Compute CPM on a file-backed matrix for CLAMP (in-place)
- [`preprocessCLAMP()`](https://chikinalab.github.io/CLAMP/reference/preprocessCLAMP.md)
  : Preprocess an expression matrix for CLAMP
- [`preprocessCLAMPFBM()`](https://chikinalab.github.io/CLAMP/reference/preprocessCLAMPFBM.md)
  : Preprocess a bigstatsr FBM for CLAMP
- [`zscoreCLAMP()`](https://chikinalab.github.io/CLAMP/reference/zscoreCLAMP.md)
  : Z-score a filtered expression matrix for CLAMP
- [`zscoreCLAMPFBM()`](https://chikinalab.github.io/CLAMP/reference/zscoreCLAMPFBM.md)
  : Z-score a filtered FBM in-place
- [`tscale()`](https://chikinalab.github.io/CLAMP/reference/tscale.md) :
  Row-wise scaling (mean 0, sd 1)

## Prior knowledge

- [`getGMT()`](https://chikinalab.github.io/CLAMP/reference/getGMT.md) :
  Download and read a GMT file from a URL
- [`read_gmt()`](https://chikinalab.github.io/CLAMP/reference/read_gmt.md)
  : Read a GMT file into a list
- [`gmtListToSparseMat()`](https://chikinalab.github.io/CLAMP/reference/gmtListToSparseMat.md)
  : Convert a list of GMT gene sets to a sparse matrix
- [`getMatchedPathwayMat()`](https://chikinalab.github.io/CLAMP/reference/getMatchedPathwayMat.md)
  : Subset and filter pathway matrix to match target genes
- [`getMatchedPathwayMatList()`](https://chikinalab.github.io/CLAMP/reference/getMatchedPathwayMatList.md)
  : Subset and filter multiple pathway matrices to match target genes
- [`getChat()`](https://chikinalab.github.io/CLAMP/reference/getChat.md)
  : Compute Chat matrix from prior annotation

## Evaluation and visualization

- [`compareBs()`](https://chikinalab.github.io/CLAMP/reference/compareBs.md)
  : Compare two sets of factor loadings or embeddings
- [`allAgainstAllAUCs()`](https://chikinalab.github.io/CLAMP/reference/allAgainstAllAUCs.md)
  : Compute all-vs-all AUC matrix
- [`CLAMPplotU()`](https://chikinalab.github.io/CLAMP/reference/CLAMPplotU.md)
  : Plot the U matrix (pathway-LV associations) as a heatmap
- [`CLAMPplotTopZ()`](https://chikinalab.github.io/CLAMP/reference/CLAMPplotTopZ.md)
  : Plot top genes per LV by Z loading
- [`CLAMPdotplot()`](https://chikinalab.github.io/CLAMP/reference/CLAMPdotplot.md)
  : Dot plot of top pathways for a single latent variable
- [`CLAMPdotplotAll()`](https://chikinalab.github.io/CLAMP/reference/CLAMPdotplotAll.md)
  : Dot plot of pathway-LV associations across all latent variables
- [`plotTopZ_Complex()`](https://chikinalab.github.io/CLAMP/reference/plotTopZ_Complex.md)
  : ComplexHeatmap visualization of top genes by latent variable

## FBM utilities

- [`cleanFBM()`](https://chikinalab.github.io/CLAMP/reference/cleanFBM.md)
  : Clean a Filebacked Big Matrix (FBM) by log-transforming and handling
  NAs
- [`filterFBM()`](https://chikinalab.github.io/CLAMP/reference/filterFBM.md)
  : Filter rows of a Filebacked Big Matrix based on mean and variance

## Solvers and helpers

- [`cross_ZY()`](https://chikinalab.github.io/CLAMP/reference/cross_ZY.md)
  : Cross-product Z^T Y with FBM or dense matrices
- [`ridge_B()`](https://chikinalab.github.io/CLAMP/reference/ridge_B.md)
  : Ridge regression update for B
- [`solveU()`](https://chikinalab.github.io/CLAMP/reference/solveU.md) :
  Fit the loading matrix Z using sparse regression of prior information
  U
- [`mat_mult()`](https://chikinalab.github.io/CLAMP/reference/mat_mult.md)
  : Matrix multiplication with support for FBM objects
- [`findSplineMax()`](https://chikinalab.github.io/CLAMP/reference/findSplineMax.md)
  : Find the location of the maximum of a smoothing spline
- [`oneToOneMask()`](https://chikinalab.github.io/CLAMP/reference/oneToOneMask.md)
  : One-to-one masking of maximum associations
- [`squashZscore()`](https://chikinalab.github.io/CLAMP/reference/squashZscore.md)
  : Squash extreme z-scores
- [`winsor_topk()`](https://chikinalab.github.io/CLAMP/reference/winsor_topk.md)
  : Winsorize matrix columns by capping the top-k values

## Datasets

- [`dataWholeBlood`](https://chikinalab.github.io/CLAMP/reference/dataWholeBlood.md)
  : Whole-blood reference expression matrix
- [`celltypeTargets`](https://chikinalab.github.io/CLAMP/reference/celltypeTargets.md)
  : Cell-type deconvolution matrix
- [`majorCellTypes`](https://chikinalab.github.io/CLAMP/reference/majorCellTypes.md)
  : Major cell-type annotations
- [`panDB`](https://chikinalab.github.io/CLAMP/reference/panDB.md) :
  panDB gene-set database
- [`xCell`](https://chikinalab.github.io/CLAMP/reference/xCell.md) :
  xCell cell-signature matrix
