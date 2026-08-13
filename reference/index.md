# Package index

## All functions

- [`AUC()`](https://chikinalab.org/CLAMP/reference/AUC.md) : Compute AUC
  using Wilcoxon rank-sum test
- [`BH()`](https://chikinalab.org/CLAMP/reference/BH.md) : Adjust
  p-values using Benjamini-Hochberg method
- [`CLAMPbase()`](https://chikinalab.org/CLAMP/reference/CLAMPbase.md) :
  CLAMP base matrix factorization
- [`CLAMPdotplot()`](https://chikinalab.org/CLAMP/reference/CLAMPdotplot.md)
  : Dot plot of top pathways for a single latent variable
- [`CLAMPdotplotAll()`](https://chikinalab.org/CLAMP/reference/CLAMPdotplotAll.md)
  : Dot plot of pathway-LV associations across all latent variables
- [`CLAMPfull()`](https://chikinalab.org/CLAMP/reference/CLAMPfull.md) :
  Runs the streamlined full CLAMP model.
- [`CLAMPfullnVP()`](https://chikinalab.org/CLAMP/reference/CLAMPfullnVP.md)
  : Full CLAMP model with prior information and cross-validation
- [`CLAMPplotTopZ()`](https://chikinalab.org/CLAMP/reference/CLAMPplotTopZ.md)
  : Plot top genes per LV by Z loading
- [`CLAMPplotU()`](https://chikinalab.org/CLAMP/reference/CLAMPplotU.md)
  : Plot the U matrix (pathway-LV associations) as a heatmap
- [`allAgainstAllAUCs()`](https://chikinalab.org/CLAMP/reference/allAgainstAllAUCs.md)
  : Compute all-vs-all AUC matrix
- [`binarizeTop()`](https://chikinalab.org/CLAMP/reference/binarizeTop.md)
  : Binarize matrix by top-k values per column
- [`celltypeTargets`](https://chikinalab.org/CLAMP/reference/celltypeTargets.md)
  : Cell-type deconvolution matrix
- [`cleanFBM()`](https://chikinalab.org/CLAMP/reference/cleanFBM.md) :
  Clean a Filebacked Big Matrix (FBM) by log-transforming and handling
  NAs
- [`commonRows()`](https://chikinalab.org/CLAMP/reference/commonRows.md)
  : Find common row names between two matrices or data frames
- [`compareBs()`](https://chikinalab.org/CLAMP/reference/compareBs.md) :
  Compare two sets of factor loadings or embeddings
- [`computeRowStatsFBM()`](https://chikinalab.org/CLAMP/reference/computeRowStatsFBM.md)
  : Compute row-wise sum and sum of squares for a Filebacked Big Matrix
- [`compute_svd()`](https://chikinalab.org/CLAMP/reference/compute_svd.md)
  : Compute a truncated SVD for a CLAMP input matrix
- [`cpmCLAMP()`](https://chikinalab.org/CLAMP/reference/cpmCLAMP.md) :
  Compute counts-per-million (CPM) for CLAMP pipelines
- [`cpmCLAMPFBM()`](https://chikinalab.org/CLAMP/reference/cpmCLAMPFBM.md)
  : Compute CPM on a file-backed matrix for CLAMP (in-place)
- [`crossVal()`](https://chikinalab.org/CLAMP/reference/crossVal.md) :
  Cross-validation AUC for CLAMP latent variables and pathways
- [`cross_ZY()`](https://chikinalab.org/CLAMP/reference/cross_ZY.md) :
  Cross-product Z^T Y with FBM or dense matrices
- [`dataWholeBlood`](https://chikinalab.org/CLAMP/reference/dataWholeBlood.md)
  : Whole-blood reference expression matrix
- [`differentialLVActivity()`](https://chikinalab.org/CLAMP/reference/differentialLVActivity.md)
  : Differential latent-variable activity between sample groups
- [`filterFBM()`](https://chikinalab.org/CLAMP/reference/filterFBM.md) :
  Filter rows of a Filebacked Big Matrix based on mean and variance
- [`findSplineMax()`](https://chikinalab.org/CLAMP/reference/findSplineMax.md)
  : Find the location of the maximum of a smoothing spline
- [`getAUCstats()`](https://chikinalab.org/CLAMP/reference/getAUCstats.md)
  : Count number of latent variables exceeding AUC thresholds
- [`getChat()`](https://chikinalab.org/CLAMP/reference/getChat.md) :
  Compute Chat matrix from prior annotation
- [`getGMT()`](https://chikinalab.org/CLAMP/reference/getGMT.md) :
  Download and read a GMT file from a URL
- [`getMatchedPathwayMat()`](https://chikinalab.org/CLAMP/reference/getMatchedPathwayMat.md)
  : Subset and filter pathway matrix to match target genes
- [`getMatchedPathwayMat2()`](https://chikinalab.org/CLAMP/reference/getMatchedPathwayMat2.md)
  : Subset and filter multiple pathway matrices to match target genes
- [`getMatchedPathwayMatList()`](https://chikinalab.org/CLAMP/reference/getMatchedPathwayMatList.md)
  : Subset and filter multiple pathway matrices to match target genes
- [`getMatchedPathwayMatOld()`](https://chikinalab.org/CLAMP/reference/getMatchedPathwayMatOld.md)
  : Subset and filter pathway matrix to match target genes
- [`getMaxAUC()`](https://chikinalab.org/CLAMP/reference/getMaxAUC.md) :
  Get maximum AUC per latent variable
- [`getScaleFromSVs()`](https://chikinalab.org/CLAMP/reference/getScaleFromSVs.md)
  : Estimate noise scale from singular values with linear tail
  extrapolation
- [`gmtListToSparseMat()`](https://chikinalab.org/CLAMP/reference/gmtListToSparseMat.md)
  : Convert a list of GMT gene sets to a sparse matrix
- [`majorCellTypes`](https://chikinalab.org/CLAMP/reference/majorCellTypes.md)
  : Major cell-type annotations
- [`mat_mult()`](https://chikinalab.org/CLAMP/reference/mat_mult.md) :
  Matrix multiplication with support for FBM objects
- [`max_correspondence_greedy()`](https://chikinalab.org/CLAMP/reference/max_correspondence_greedy.md)
  : Greedy maximum correspondence from correlation matrix
- [`mymessage()`](https://chikinalab.org/CLAMP/reference/mymessage.md) :
  Print a concatenated message
- [`num.pc()`](https://chikinalab.org/CLAMP/reference/num.pc.md) :
  Estimate number of principal components via elbow or permutation
  method
- [`oneToOneMask()`](https://chikinalab.org/CLAMP/reference/oneToOneMask.md)
  : One-to-one masking of maximum associations
- [`panDB`](https://chikinalab.org/CLAMP/reference/panDB.md) : panDB
  gene-set database
- [`pinv.ridge()`](https://chikinalab.org/CLAMP/reference/pinv.ridge.md)
  : Ridge-regularized pseudoinverse via SVD
- [`plotTopZ_Complex()`](https://chikinalab.org/CLAMP/reference/plotTopZ_Complex.md)
  : ComplexHeatmap visualization of top genes by latent variable
- [`preprocessCLAMP()`](https://chikinalab.org/CLAMP/reference/preprocessCLAMP.md)
  : Preprocess an expression matrix for CLAMP
- [`preprocessCLAMPFBM()`](https://chikinalab.org/CLAMP/reference/preprocessCLAMPFBM.md)
  : Preprocess a bigstatsr FBM for CLAMP
- [`projectCLAMP()`](https://chikinalab.org/CLAMP/reference/projectCLAMP.md)
  : Project new data into CLAMP latent space
- [`read_gmt()`](https://chikinalab.org/CLAMP/reference/read_gmt.md) :
  Read a GMT file into a list
- [`ridge_B()`](https://chikinalab.org/CLAMP/reference/ridge_B.md) :
  Ridge regression update for B
- [`rotateSVD()`](https://chikinalab.org/CLAMP/reference/rotateSVD.md) :
  Rotate SVD components to make dominant directions positive
- [`row_cor()`](https://chikinalab.org/CLAMP/reference/row_cor.md) :
  Row-wise correlation between two matrices
- [`run_elbow()`](https://chikinalab.org/CLAMP/reference/run_elbow.md) :
  Run elbow method to estimate number of PCs
- [`run_permutation()`](https://chikinalab.org/CLAMP/reference/run_permutation.md)
  : Run permutation method to estimate number of PCs
- [`select_clamp_k()`](https://chikinalab.org/CLAMP/reference/select_clamp_k.md)
  : Select default number of CLAMP latent variables from an SVD
- [`select_svd_k()`](https://chikinalab.org/CLAMP/reference/select_svd_k.md)
  : Select default number of components for a CLAMP solver SVD
- [`solveU()`](https://chikinalab.org/CLAMP/reference/solveU.md) : Fit
  the loading matrix Z using sparse regression of prior information U
- [`squashZscore()`](https://chikinalab.org/CLAMP/reference/squashZscore.md)
  : Squash extreme z-scores
- [`tscale()`](https://chikinalab.org/CLAMP/reference/tscale.md) :
  Row-wise scaling (mean 0, sd 1)
- [`winsor_topk()`](https://chikinalab.org/CLAMP/reference/winsor_topk.md)
  : Winsorize matrix columns by capping the top-k values
- [`xCell`](https://chikinalab.org/CLAMP/reference/xCell.md) : xCell
  cell-signature matrix
- [`zscoreCLAMP()`](https://chikinalab.org/CLAMP/reference/zscoreCLAMP.md)
  : Z-score a filtered expression matrix for CLAMP
- [`zscoreCLAMPFBM()`](https://chikinalab.org/CLAMP/reference/zscoreCLAMPFBM.md)
  : Z-score a filtered FBM in-place
