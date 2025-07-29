#' @importFrom stats wilcox.test quantile median p.adjust
#' @importFrom utils download.file
#' @importFrom dplyr %>% group_by summarize ungroup
#' @importFrom bigstatsr big_cprodMat big_prodMat big_apply rows_along FBM
#' @importFrom stats coef
#' @importFrom utils data flush.console
#' @importFrom rsvd rsvd
#' @importFrom irlba irlba
#' @importFrom Matrix sparseMatrix
#' @importFrom rlang .data
#' @importFrom glmnet glmnet
#' @importFrom glmnet cv.glmnet
#' @importFrom Matrix t
# library(doParallel)
# library(foreach)
library(bigstatsr)
library(Matrix)
library(glmnet)
#' Adjust p-values using Benjamini-Hochberg method
#'
#' Applies the BH (Benjamini-Hochberg) correction for multiple hypothesis testing.
#'
#' @param p Numeric vector of p-values.
#' @return Adjusted p-values.
BH <- function(p) {
    p.adjust(p, method = "BH")
}

#' Row-wise correlation between two matrices
#'
#' Computes the Pearson correlation for each row between matrices \code{A} and \code{B}.
#'
#' @param A A numeric matrix.
#' @param B A numeric matrix of the same dimensions as \code{A}.
#' @return A numeric vector of correlations, one per row.
row_cor <- function(A, B) {
    A_mat <- as.matrix(A)
    B_mat <- as.matrix(B)
    if (!all(dim(A_mat) == dim(B_mat))) {
        stop("row_cor(): dimensions of A and B must match")
    }
    A_centered <- A_mat - rowMeans(A_mat)
    B_centered <- B_mat - rowMeans(B_mat)
    num <- rowSums(A_centered * B_centered)
    den <- sqrt(rowSums(A_centered^2) * rowSums(B_centered^2))
    num / den
}

#' Matrix multiplication with support for FBM objects
#'
#' Multiplies two matrices, using optimized multiplication if the first is a Filebacked Big Matrix (FBM).
#'
#' @param mat1 A matrix or an object of class \code{FBM}.
#' @param mat2 A numeric matrix.
#' @return Matrix product of \code{mat1} and \code{mat2}.
mat_mult <- function(mat1, mat2) {
    is_fbm <- inherits(mat1, "FBM")
    if (is_fbm) {
        # For FBM objects, use the specific multiplication method
        return(big_prodMat(mat1, as.matrix(mat2)))
    } else {
        # Regular matrix multiplication
        return(mat1 %*% mat2)
    }
}


#' Ridge-regularized pseudoinverse via SVD
#'
#' Computes a stable pseudoinverse of a symmetric positive semi-definite matrix
#' using singular value decomposition (SVD) and ridge regularization. This is
#' useful when the matrix is ill-conditioned or rank-deficient.
#'
#' @param m A symmetric numeric matrix (e.g., from \code{crossprod()}).
#' @param alpha Non-negative scalar specifying the ridge penalty. A small positive value
#'        stabilizes the inversion by shrinking large singular values.
#'
#' @return A numeric matrix representing the ridge-regularized pseudoinverse of \code{m}.
pinv.ridge <- function(m, alpha = 0) {
    msvd <- svd(m)
    d <- msvd$d

    if (length(d) == 0 || all(d == 0)) {
        return(matrix(0, nrow = ncol(m), ncol = nrow(m)))
    }

    d_inv <- if (alpha > 0) d / (d^2 + alpha^2) else 1 / d

    msvd$v %*% (d_inv * t(msvd$u))
}

#' Rotate SVD components to make dominant directions positive
#'
#' Ensures consistency in SVD output by flipping signs so that each left singular vector
#' has a majority of positive entries.
#'
rotateSVD <- function(svdres) {
    upos <- svdres$u
    uneg <- svdres$u
    upos[upos < 0] <- 0
    uneg[uneg >= 0] <- 0
    uneg <- -uneg
    sumposu <- colSums(upos)
    sumnegu <- colSums(uneg)


    for (i in 1:ncol(svdres$u)) {
        if (sumnegu[i] > sumposu[i]) {
            svdres$u[, i] <- -svdres$u[, i]
            svdres$v[, i] <- -svdres$v[, i]
        }
    }
    svdres
}


#' Binarize matrix by top-k values per column
#'
#' Keeps only the top \code{top} values in each column of a matrix, setting others to 0.
#'
#' @param Z A numeric matrix.
#' @param top Number of top entries to keep in each column.
#' @param keepVals If \code{TRUE}, retains original values above the cutoff; otherwise, sets them to 1.
#' @return A modified matrix with only top entries retained per column.
binarizeTop <- function(Z, top, keepVals = T) {
    for (i in 1:ncol(Z)) {
        cutoff <- sort(Z[, i], T)[top + 1]


        if (cutoff == 0) {
            cutoff <- min(Z[Z[, i] > 0, i])
        }
        Z[Z[, i] < cutoff, i] <- 0
        if (!keepVals) {
            Z[Z[, i] > 0, i] <- 1
        }
    }
    Z
}

#' Fit the loading matrix Z using sparse regression of prior information U
#'

#' For each column of a target matrix \code{Z} using pathway or prior annotation \code{priorMat}.
#' It performs regularization selection using cross-validation and can apply either
#' Supports continuous or binary response models. Relaxed refitting is supported for final coefficient estimation.
#'
#' @param Z A numeric matrix with features (rows) and samples (columns).
#' @param Chat (Optional) Precomputed pseudo-inverse of \code{priorMat}; if \code{NULL}, it is calculated using ridge regularization.
#' @param priorMat A numeric matrix with prior information (features × pathways).
#' @param penalty.factor Optional penalty weights for features in \code{priorMat}.
#' @param pathwaySelection Method to select candidate pathways: \code{"fast"} (default) or \code{"complete"}.
#' @param alpha Elastic net mixing parameter (0 = ridge, 1 = lasso). Default is 0.9.
#' @param maxPath Maximum number of pathways/features selected per column. Default is 10.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param useSE Whether to use the 1-standard-error rule for lambda selection. Default is \code{FALSE}.
#' @param top If set, sets to 0 all but the top entries of \code{Z} per column before fitting.
#' @param binary If \code{TRUE}, fits a binomial model (e.g., classification) to \code{Z}>0. Can be used incombination with \code{top}. Default is \code{FALSE}.
#' @param nlambda Number of lambda values for glmnet. Default is 20.
#' @param scale Whether to standardize predictors in glmnet. Default is \code{TRUE}.
#' @param refit Whether to perform relaxed refitting using selected predictors. Default is \code{TRUE}.
#' @param Uprev (Optional) Previous U matrix to reuse. In this mode only the columns of \code{U} that are all zero are estimated. Used internally in \code{PLIER2}.
#' @param ... Additional arguments passed to \code{glmnet()} or \code{cv.glmnet()}.
#'
#' @return A list with one element:
#' \describe{
#'   \item{\code{U}}{A matrix of loadings (features × components). Columns are named \code{LV1}, \code{LV2}, ...}
#' }
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom Matrix crossprod
#' @export

solveU <- function(Z, Chat = NULL, priorMat, penalty.factor, pathwaySelection = "fast", alpha = 0.9, maxPath = 10, nfolds = 5, useSE = F, top = NULL, binary = F, nlambda = 20, scale = T, refit = T, Uprev = NULL, ...) {
    if (nrow(Z) != nrow(priorMat)) {
        cm <- commonRows(Z, priorMat)
        Z <- Z[cm, ]
        iim <- match(cm, rownames(priorMat))
        priorMat <- priorMat[iim, ]

        message("matching rows")
    }
    if (scale) {
        col_means <- Matrix::colMeans(priorMat)
        col_sds <- sqrt(Matrix::colMeans(priorMat^2) - col_means^2)
        col_sds[col_sds == 0] <- 1 # avoid divide-by-zero
        priorMat <- Matrix::t((Matrix::t(priorMat) - col_means) / col_sds)
    }
    if (is.null(Chat)) {
        Chat <- pinv.ridge(crossprod(priorMat), 5) %*% (t(priorMat))
    }
    # printmessage("New solve U")
    Ur <- Chat %*% Z # get U by OLS

    Ur <- apply(-Ur, 2, rank) # rank
    Urm <- apply(Ur, 1, min)

    # Zhat <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
    if (is.null(Uprev)) {
        U <- Matrix::Matrix(0, nrow = ncol(priorMat), ncol = ncol(Z), sparse = TRUE)
        #  U=matrix(0,nrow=ncol(priorMat), ncol=ncol(Z))
    } else {
        U <- Uprev
    }
    pathwaySelection <- match.arg(pathwaySelection, c("fast", "complete"))

    if (!is.null(top)) {
        Znew <- binarizeTop(Z, top, keepVals = TRUE)

        Z <- Znew
    }
    if (pathwaySelection == "complete") {
        iip <- which(Urm <= maxPath * 5)
        message(paste("Picked", length(iip), "pathways"))
    }


    for (i in 1:ncol(Z)) {
        if (any(U[, i] > 0)) {
            next
        }
        if (pathwaySelection == "fast") {
            iip <- which(Ur[, i] <= maxPath)
        }

        if (!binary) { # not doing a binary prediction

            set.seed(1)
            gres <- cv.glmnet(y = Z[, i], x = priorMat[, iip], alpha = alpha, lower.limits = 0, foldid = ((1:nrow(Z)) %% nfolds) + 1, keep = T, nfolds = nfolds, standardize = scale, dfmax = maxPath, nlambda = nlambda, ...)


            # plot(gres)
        } else {
            set.seed(1)
            gres <- cv.glmnet(y = (Z[, i] > 0) + 1 - 1, x = priorMat[, iip], family = "binomial", alpha = alpha, lower.limits = 0, foldid = ((1:nrow(Z)) %% nfolds) + 1, keep = F, nfolds = nfolds, type.measure = "auc", dfmax = maxPath, nlambda = nlambda, standardize = scale, nlambda = nlambda, ...)
        }


        if (refit) {
            # Select lambda
            s_best <- if (useSE) gres$lambda.1se else gres$lambda.min
            active_coef <- coef(gres, s = s_best)
            active_ix <- which(active_coef[-1] != 0)
            selected_features <- iip[active_ix]
            if (length(selected_features) == 0) {
                next
            }
            X_sel <- priorMat[, selected_features, drop = FALSE]

            # Add dummy column if only one predictor
            add_dummy <- ncol(X_sel) == 1
            if (add_dummy) {
                X_sel <- cbind(X_sel, dummy = 0)
            }

            # Fit relaxed model
            if (!binary) {
                fit_relaxed <- glmnet(X_sel, Z[, i],
                    alpha = 0, lambda = 0,
                    lower.limits = 0, standardize = scale
                )
                coefs <- as.vector(coef(fit_relaxed))[-1]
                if (add_dummy) coefs <- coefs[-length(coefs)]
                U[selected_features, i] <- coefs
                #  Zhat[, i] <- predict(fit_relaxed, newx = X_sel, s = 0)[,1]
            } else {
                fit_relaxed <- glmnet(X_sel, (Z[, i] > 0) + 0,
                    family = "binomial", alpha = 0, lambda = 0,
                    lower.limits = 0, standardize = scale
                )
                coefs <- as.vector(coef(fit_relaxed))[-1]
                if (add_dummy) coefs <- coefs[-length(coefs)]
                U[selected_features, i] <- coefs
                # Zhat[, i] <- predict(fit_relaxed, newx = X_sel, s = 0, type = "response")[,1]
            }
        } else {
            # old code with extra functions
            #      betaI=getNonZeroBetas(gres,se = useSE, index=T)
            #     beta <- getNonZeroBetas(gres, se = useSE, index = FALSE)
            #    U[iip[betaI], i] <- as.vector(beta)

            s_best <- if (useSE) gres$lambda.1se else gres$lambda.min
            coef_vec <- as.vector(coef(gres, s = s_best))[-1] # drop intercept
            betaI <- which(coef_vec != 0)
            U[iip[betaI], i] <- coef_vec[betaI]
        }
        # end for i in Z
    }



    cat(paste(", Number of annotated columns is", sum(Matrix::colSums(U) > 0)))
    # rownames(U)=substr(colnames(priorMat),1,30)
    rownames(U) <- colnames(priorMat)
    colnames(U) <- paste("LV", 1:ncol(U))


    U <- as.matrix(U)

    return(list(U = U))



    message(paste("Number of annotated columns is", sum(Matrix::colSums(U) > 0)))

    rownames(U) <- colnames(priorMat)
    colnames(U) <- paste("LV", 1:ncol(U))
    return(U)
}
#' Compute Chat matrix from prior annotation
#'
#' Computes the transformation matrix \code{Chat} used to map from observed data to latent space,
#' based on a pseudo-inverse of the prior annotation matrix. Optionally standardizes the columns
#' of \code{priorMat} before computing.
#'
#' @param priorMat A numeric or sparse matrix (features × pathways) containing prior annotations.
#' @param scale Logical; if \code{TRUE} (default), standardizes the columns of \code{priorMat} before computing \code{Chat}.
#'
#' @return A numeric matrix \code{Chat} of dimensions (pathways × features).
#' @export
getChat <- function(priorMat, scale = TRUE) {
    if (scale) {
        col_means <- Matrix::colMeans(priorMat)
        col_sds <- sqrt(Matrix::colMeans(priorMat^2) - col_means^2)
        col_sds[col_sds == 0] <- 1
        priorMat <- Matrix::t((Matrix::t(priorMat) - col_means) / col_sds)
    }

    message("Inverting...")
    Chat <- pinv.ridge(crossprod(priorMat), 5) %*% t(priorMat)

    message("done")
    Chat
}












#' Subset and filter pathway matrix to match target genes
#'
#' Filters a gene-by-pathway annotation matrix to retain only pathways
#' with sufficient overlap with a given gene set. The result is a sparse matrix
#' aligned to \code{new.genes}, with columns (pathways) retained only if they
#' have at least \code{min.genes} matched genes.
#'
#' @param pathMat A sparse binary matrix of genes (rows) × pathways (columns).
#' @param new.genes Character vector of gene names to match.
#' @param min.genes Minimum number of overlapping genes required to keep a pathway.
#'
#' @return A sparse matrix of dimensions \code{length(new.genes)} × filtered pathways.
#' @export
getMatchedPathwayMat <- function(pathMat, new.genes, min.genes = 10) {
    cm <- intersect(rownames(pathMat), new.genes)
    mymessage("there are ", length(cm), " genes in the intersection between data and prior")

    matchPathMat <- Matrix::sparseMatrix(
        i = match(cm, new.genes),
        j = rep(1, length(cm)), # temporary, will be overwritten below
        dims = c(length(new.genes), ncol(pathMat)),
        x = 0, # fill with zeros for now
        dimnames = list(new.genes, colnames(pathMat))
    )
    matchPathMat[cm, ] <- pathMat[cm, ]

    genesInPath <- Matrix::colSums(matchPathMat)
    ii <- which(genesInPath >= min.genes)
    message(sprintf("Removing %d pathways", ncol(matchPathMat) - length(ii)))

    matchPathMat[, ii]
}







#' Compute AUC using Wilcoxon rank-sum test
#'
#' Computes the area under the ROC curve (AUC) by applying a Wilcoxon rank-sum test
#' between predicted values for positive and negative labels. This is equivalent to
#' computing the Mann–Whitney U statistic.
#'
#' @param labels A numeric or logical vector indicating class labels. Values > 0 are treated as positive.
#' @param values A numeric vector of prediction scores corresponding to \code{labels}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{auc}}{Estimated AUC, or 0.5 if one class is missing}
#'   \item{\code{pval}}{Wilcoxon test p-value, or \code{NA} if one class is missing}
#' }
#'
#' @export
AUC <- function(labels, values) {
    pos <- labels > 0
    neg <- !pos
    posn <- sum(pos)
    negn <- sum(neg)

    if (posn > 0 && negn > 0) {
        res <- suppressWarnings(wilcox.test(values[pos], values[neg], alternative = "greater"))
        auc <- unname(res$statistic) / (posn * negn)
        pval <- res$p.value
    } else {
        auc <- 0.5
        pval <- NA
    }

    list(auc = auc, pval = pval)
}







#' Cross-validation AUC for PLIER latent variables and pathways
#'
#' Evaluates how well each latent variable in a PLIER model captures held-out pathway annotations,
#' using cross-validation over the prior matrix. For each latent variable and associated pathway,
#' held-out genes are selected and the AUC is computed using their scores in \code{plierRes$Z}.
#'
#' @param plierRes A list containing \code{U} (loadings) and \code{Z} (scores) from a PLIER model.
#' @param priorMat A binary matrix (genes × pathways) indicating original pathway annotations.
#' @param priorMatcv A version of \code{priorMat} used to mask held-out annotations for cross-validation.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{Uauc}}{Matrix of AUC values (pathways × LVs)}
#'   \item{\code{Upval}}{Matrix of \code{-log10(p)} values (pathways × LVs)}
#'   \item{\code{summary}}{Data frame with pathway, LV index, AUC, p-value, and FDR}
#' }
#'
#' @export
crossVal <- function(plierRes, priorMat, priorMatcv) {
    out <- matrix(ncol = 4, nrow = 0)
    ii <- which(Matrix::colSums(plierRes$U) > 0)

    # Uauc <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = dim(plierRes$U))
    #  Up <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = dim(plierRes$U))
    Uauc <- Matrix(0, nrow = nrow(plierRes$U), ncol = ncol(plierRes$U), sparse = T)
    Up <- Matrix(0, nrow = nrow(plierRes$U), ncol = ncol(plierRes$U), sparse = T)
    for (i in ii) { # for each column in U

        iipath <- which(plierRes$U[, i] > 0) # get the pathways

        if (length(iipath) > 1) { # more than one pathway
            for (j in iipath) {
                iiheldout <- which((rowSums(priorMat[, iipath, drop = F]) == 0) | (priorMat[, j] > 0 & priorMatcv[, j] == 0))
                aucres <- AUC(priorMat[iiheldout, j], plierRes$Z[iiheldout, i])
                out <- rbind(out, c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
                Uauc[j, i] <- aucres$auc
                Up[j, i] <- -log10(aucres$pval)
            }
        } else {
            j <- iipath
            iiheldout <- which((rowSums(matrix(priorMat[, iipath], ncol = 1)) == 0) | (priorMat[, j] > 0 & priorMatcv[, j] == 0))
            aucres <- AUC(priorMat[iiheldout, j], plierRes$Z[iiheldout, i])
            out <- rbind(out, c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
            Uauc[j, i] <- aucres$auc
            Up[j, i] <- -log10(aucres$pval)
        } # else
    }
    out <- data.frame(out, stringsAsFactors = F)
    out[, 3] <- as.numeric(out[, 3])
    out[, 4] <- as.numeric(out[, 4])
    out[, 5] <- BH(out[, 4])
    colnames(out) <- c("pathway", "LV index", "AUC", "p-value", "FDR")
    return(list(Uauc = Uauc, Upval = Up, summary = out))
}




#'  PLIER base matrix factorization
#'
#' Runs the core matrix factorization procedure of PLIER (Pathway-Level Information ExtractoR),
#' decomposing the gene expression matrix \code{Y} into latent variables \code{Z} and loadings \code{B}.
#' It supports sparse, dense, and Filebacked Big Matrices (FBM) as input and includes options for
#' adaptive sparsity, positive constraints, and regularization.
#'
#' @param Y Input gene expression matrix (genes × samples). Can be dense, sparse (\code{dgCMatrix}), or FBM.
#' @param k Number of latent variables.
#' @param svdres Optional precomputed SVD result. If not supplied, it is computed internally.
#' @param L1 L1 regularization strength for Z. Defaults to scaled singular value.
#' @param L2 L2 regularization strength for B. Defaults to scaled singular value.
#' @param Zpos Logical; if \code{TRUE}, negative entries in Z are zeroed. Default is \code{TRUE}.
#' @param max.iter Maximum number of optimization iterations. Default is 200.
#' @param tol Convergence tolerance for B update. Default is 5e-4.
#' @param trace Logical; if \code{TRUE}, prints progress. Default is \code{FALSE}.
#' @param rseed Optional integer for reproducible random initialization of B.
#' @param B Optional initial matrix for B. If not provided, initialized from SVD.
#' @param scale Scaling factor for L1 and L2 when not provided. Default is 1.
#' @param pos.adj Positive constraint adjustment divisor for L1. Default is 3.
#' @param adaptive.p Controls adaptive sparsity in \code{Z}. After each ALS update, negative entries in \code{Z} are assumed to reflect noise. The cutoff for thresholding is set according to the probability of positive values under a reflected negative distribution—effectively zeroing out small positive entries likely to be noise. Smaller values lead to more sparsity. Default is 0.05.
#'
#' @param adaptive.iter Number of iterations before adaptive sparsity is applied. Default is 20.
#' @param cutoff Scalar threshold to zero Z values when \code{Zpos = TRUE} and adaptive thresholding is not used. Default is 0.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{B}}{Latent variable loadings (LVs × genes)}
#'   \item{\code{Z}}{Latent variable scores (LVs × samples)}
#'   \item{\code{Zraw}}{Raw Z matrix before thresholding}
#'   \item{\code{L1}}{Final value of L1 used}
#'   \item{\code{L2}}{Final value of L2 used}
#' }
#'
#' @details
#' This function is the low-level implementation of PLIER. It alternates between solving for \code{Z}
#' given \code{B} and solving for \code{B} given \code{Z}, with optional sparsity and non-negativity
#' constraints on \code{Z}. Convergence is assessed via relative change in \code{B}.
#'
#' @export
PLIERbase <- function(Y, k, svdres = NULL, L1 = NULL, L2 = NULL,
    Zpos = T, max.iter = 200, tol = 5e-4, trace = F,
    rseed = NULL, B = NULL, scale = 1, pos.adj = 3, adaptive.p = 0.05, adaptive.iter = 20, cutoff = 0) {
    # message("Checking type")
    # Detect matrix type
    is_fbm <- inherits(Y, "FBM")
    is_sparse <- inherits(Y, "dgCMatrix")

    ng <- nrow(Y)
    ns <- ncol(Y)

    Bdiff <- Inf
    BdiffTrace <- double()
    BdiffCount <- 0
    message("****")
    k <- min(k, min(dim(Y)) - 1)

    if (is.null(svdres) && is.null(B)) {
        message("Computing SVD")
        if (is_fbm) {
            # For FBM, we need special handling for SVD
            set.seed(123)
            if (requireNamespace("bigstatsr", quietly = TRUE)) {
                # Use big_SVD from bigstatsr if available
                svdres <- bigstatsr::big_SVD(X = Y, k = k)
            } else {
                # Fallback: convert to regular matrix for SVD
                # This might be memory-intensive for large matrices
                svdres <- rsvd(Y, k = k)
            }
        } else if (is_sparse) {
            # For sparse matrices, use irlba or other sparse SVD methods
            if (requireNamespace("irlba", quietly = TRUE)) {
                set.seed(123)
                svdres <- irlba::irlba(Y, nv = k)
            } else {
                set.seed(123)
                svdres <- rsvd(Y, k = k)
            }
        } else {
            # Regular matrix
            set.seed(123)
            svdres <- rsvd(Y, k = k)
        }

        svdres <- rotateSVD(svdres)
    }



    if (is.null(L1)) {
        L1 <- svdres$d[k] * scale
        if (!is.null(pos.adj)) {
            L1 <- L1 / pos.adj
        }
    }
    if (is.null(L2)) {
        L2 <- svdres$d[k] * scale
    }
    L2k <- L2 * diag(k)
    #    L1=svdres$d[k]/2*scale
    print(paste0("L1 is set to ", L1))
    print(paste0("L2 is set to ", L2))

    if (is.null(B)) {
        # initialize B with svd

        B <- t(svdres$v[, 1:k] %*% diag(sqrt(svdres$d[1:k])))
        # alternative initializations
        # seem to be not as good
        #   B=t(svdres$v[1:ncol(Y), 1:k]%*%diag((svdres$d[1:k])))
        #   B=t(svdres$v[1:ncol(Y), 1:k])
    } else {
        message("B given")
    }





    if (!is.null(rseed)) {
        message("using random start")
        set.seed(rseed)
        B <- t(apply(B, 1, sample))
    }


    round2 <- function(x) {
        signif(x, 4)
    }

    getT <- function(x) {
        -quantile(x[x < 0], adaptive.p)
    }



    for (i in 1:max.iter) {
        # main loop
        Zraw <- Z <- mat_mult(Y, t(B)) %*% solve(tcrossprod(B) + L1 * diag(k))

        if (i >= adaptive.iter && adaptive.p > 0) {
            cutoffs <- apply(Zraw, 2, getT)

            for (j in 1:ncol(Z)) {
                Z[Z[, j] < cutoffs[j], j] <- 0
            }
        } else if (Zpos) {
            Z[Z < cutoff] <- 0
        }

        oldB <- B


        if (is_fbm) {
            ZYt <- big_cprodMat(Y, as.matrix(Z))
            ZY <- t(ZYt)
            B <- solve(t(Z) %*% Z + L2k) %*% ZY
        } else {
            B <- solve(t(Z) %*% Z + L2k) %*% mat_mult(t(Z), Y)
        }

        # update error
        Bdiff <- sum((B - oldB)^2) / sum(B^2)

        # keeping track of this in case this can be useful for convergence
        minCor <- min(row_cor(B, oldB))


        BdiffTrace <- c(BdiffTrace, Bdiff)


        if (trace) {
            cat(sprintf("\rProgress %d / %d | Bdiff=%.6f, minCor=%.6f", i, max.iter, Bdiff, minCor))
            flush.console()
        }

        # check for convergence
        if (i > 52 && Bdiff > BdiffTrace[i - 50]) {
            BdiffCount <- BdiffCount + 1
        } else if (BdiffCount > 1) {
            BdiffCount <- BdiffCount - 1
        }

        if (Bdiff < tol && i > adaptive.iter + 10) {
            cat(sprintf("Converged at iteration= %d | Bdiff=%.6f,  tol=%.6f     ", i, Bdiff, tol))
            break
        }
        if (BdiffCount > 5 && i > adaptive.iter + 10) {
            message(paste0("stopped at  iteration ", i, " Bdiff is not decreasing"))
            break
        }
    }
    rownames(B) <- colnames(Z) <- paste("LV", 1:k)
    return(list(B = B, Z = Z, Zraw = Zraw, L1 = L1, L2 = L2))
}

#' Full PLIER model with prior information and cross-validation
#'
#' Runs the full PLIER (Pathway-Level Information ExtractoR) model using a gene expression matrix
#' and prior pathway annotation matrix. This function performs latent variable decomposition
#' guided by prior knowledge and includes optional cross-validation to evaluate pathway associations.
#'
#' @param Y Gene expression matrix (genes × samples). Can be dense, sparse (dgCMatrix), or FBM.
#' @param priorMat Binary matrix (genes × pathways) representing prior annotations.
#' @param svdres Optional SVD result used for initialization.
#' @param plier.base.result Optional result from \code{PLIERbase()} to initialize B.
#' @param k Number of latent variables. If \code{NULL}, estimated from SVD.
#' @param L1 Regularization strength for Z. If \code{NULL}, initialized from SVD or \code{plier.base.result}.
#' @param L2 Regularization strength for B. If \code{NULL}, initialized from SVD or \code{plier.base.result}.
#' @param top If set, keeps only top-n values per column in Z during U updates.
#' @param cvn Number of folds for cross-validation in U updates. Default is 5.
#' @param max.iter Maximum number of iterations. Default is 350.
#' @param trace Logical; if \code{TRUE}, prints iteration progress.
#' @param Chat Optional precomputed matrix for solving U.
#' @param maxPath Maximum number of pathways/features selected per LV. Default is 10.
#' @param doCrossval Whether to perform pathway-level cross-validation. Default is \code{TRUE}.
#' @param penalty.factor Vector of feature-specific penalties for glmnet. Default: all ones.
#' @param glm_alpha Elastic net mixing parameter for glmnet. Default is 0.9.
#' @param minGenes Minimum number of genes per pathway to retain. Default is 10.
#' @param tol Convergence tolerance on relative change in B. Default is 5e-4.
#' @param seed Seed for reproducibility of cross-validation masking. Default is 123456.
#' @param allGenes If \code{TRUE}, zero-fills \code{priorMat} for genes not present. Default is \code{FALSE}.
#' @param rseed Optional seed for randomly reinitializing B and Z.
#' @param max.U.updates Maximum number of U updates. Default is 5.
#' @param pathwaySelection Pathway selection mode: \code{"fast"} or \code{"complete"}.
#' @param multiplier Scaling factor for adjusting L1 and L2.
#' @param adaptive.p Quantile threshold for adaptively zeroing small Z values. After each update,
#'   the \code{adaptive.p}-quantile of negative entries in Z is used (flipped positive) to threshold
#'   small positive values, assuming they reflect noise. Default is 0.05.
#' @param useNNLS If \code{TRUE}, uses non-negative least squares in U estimation. Default is \code{TRUE}.
#' @param useRaw If \code{TRUE}, uses unthresholded Z for solving U. Default is \code{TRUE}.
#' @param refitAll If \code{TRUE}, refits all U columns every update. Default is \code{FALSE}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{B}}{Latent variable loadings (LVs × genes)}
#'   \item{\code{Z}}{Latent variable matrix (LVs × samples)}
#'   \item{\code{U}}{Pathway loadings matrix (pathways × LVs)}
#'   \item{\code{C}}{Masked prior matrix used for training}
#'   \item{\code{L1}, \code{L2}}{Regularization parameters}
#'   \item{\code{heldOutGenes}}{List of held-out genes per pathway (if CV is enabled)}
#'   \item{\code{Uauc}}{AUC matrix from CV evaluation (if enabled)}
#'   \item{\code{Up}}{-\code{log10(p)} values from CV evaluation (if enabled)}
#'   \item{\code{summary}}{Data frame of AUC, p-values, and FDR per pathway × LV (if enabled)}
#'   \item{\code{priorMatCV}}{Masked prior matrix used during CV}
#'   \item{\code{priorMat}}{Final filtered prior matrix}
#'   \item{\code{withPrior}}{Indices of LVs with non-zero pathway loadings}
#'   \item{\code{call}}{Function call}
#' }
#'
#' @details
#' The model alternates between solving \code{Z}, \code{B}, and \code{U}. Adaptive sparsity is applied
#' to \code{Z} using a dynamic threshold based on the negative tail of its distribution. Cross-validation
#' is used to hold out gene annotations in \code{priorMat} and evaluate latent variable specificity.
#'
#' @export
PLIERfull <- function(Y, priorMat, svdres = NULL, plier.base.result = NULL, k = NULL, L1 = NULL, L2 = NULL, top = NULL, cvn = 5, max.iter = 350, trace = F, Chat = NULL, maxPath = 10, doCrossval = T, penalty.factor = rep(1, ncol(priorMat)), glm_alpha = 0.9, minGenes = 10, tol = 5e-4, seed = 123456, allGenes = F, rseed = NULL, max.U.updates = 5, pathwaySelection = c("fast"), multiplier = 1, adaptive.p = 0.05, useNNLS = T, useRaw = T, refitAll = F) {
    getT <- function(x) {
        -quantile(x[x < 0], adaptive.p)
    }

    pathwaySelection <- match.arg(pathwaySelection, c("complete", "fast"))

    priorMat <- as.matrix(priorMat)

    message("**PLIER v2 **")

    # Detect matrix type
    is_fbm <- inherits(Y, "FBM")
    is_sparse <- inherits(Y, "dgCMatrix")

    if (nrow(priorMat) != nrow(Y) || !all(rownames(priorMat) == rownames(data))) {
        if (!allGenes) {
            cm <- commonRows(Y, priorMat)
            message(paste("Selecting common genes:", length(cm)))
            priorMat <- priorMat[cm, ]
            Y <- Y[cm, ]
        } else {
            extra.genes <- setdiff(rownames(Y), rownames(priorMat))
            eMat <- matrix(0, nrow = length(extra.genes), ncol = ncol(priorMat))
            rownames(eMat) <- extra.genes
            priorMat <- rbind(priorMat, eMat)
            priorMat <- priorMat[rownames(Y), ]
        }
    }

    numGenes <- Matrix::colSums(priorMat)

    heldOutGenes <- list()
    iibad <- which(numGenes < minGenes)
    if (length(iibad) > 0) {
        priorMat <- priorMat[, -iibad]
        message(paste("Removed", length(iibad), "pathways with too few genes"))
    }
    if (doCrossval) {
        priorMatCV <- as.matrix(priorMat)
        if (!is.null(seed)) {
            set.seed(seed)
        }
        for (j in 1:ncol(priorMatCV)) {
            iipos <- which(priorMatCV[, j] > 0)
            iiposs <- sample(iipos, length(iipos) / 5)
            priorMatCV[iiposs, j] <- 0
            heldOutGenes[[colnames(priorMat)[j]]] <- rownames(priorMat)[iiposs]
        }
        C <- priorMatCV
    } else {
        C <- priorMat
    }

    nc <- ncol(priorMat)
    ng <- nrow(Y)
    ns <- ncol(Y)

    Bdiff <- -1
    BdiffTrace <- double()
    BdiffCount <- 0
    if (is.null(Chat)) {
        Cp <- crossprod(C)
        Chat <- pinv.ridge(crossprod(C), 5) %*% (t(C))
    }
    # YsqSum=sum(Y^2)
    # compute svd and use that as the starting point

    if (!is.null(svdres) && nrow(svdres$v) != ncol(Y)) {
        message("SVD V has the wrong number of columns")
        svdres <- NULL
    }
    if (is.null(svdres) && is.null(plier.base.result)) {
        message("Computing SVD")
        if (ns > 500) {
            message("Using rsvd")
            set.seed(123456)
            svdres <- rsvd(Y, k = min(ns, max(200, ns / 4)), q = 3)
        } else {
            svdres <- svd(Y)
        }
        message("Done")
    }
    if (is.null(plier.base.result)) {
        if (is.null(k)) {
            k <- floor(sqrt(ncol(Y)))
            k <- min(k, floor(ncol(Y) * 0.9))
            message("k is set to ", k)
        }

        message("Running PLIERbase")
        if (is.null(plier.base.result)) {
            plier.base.result <- PLIERbase(Y, k = k)
        }
    } else {
        message("using provided PLIERbase result")

        if (nrow(Y) != nrow(plier.base.result$Z)) {
            if (is.null(rownames(Y)) | is.null(rownames(plier.base.result$Z))) {
                stop("Y and plier.base.result$Z must have equal row numbers or row names")
            }
            plier.base.result$Z <- plier.base.result$Z[rownames(Y), ]
        }
        u.iter <- 2
        k <- ncol(plier.base.result$Z)
    }

    Z <- plier.base.result$Z

    if (is.null(L1)) {
        L1 <- plier.base.result$L1
    }
    if (is.null(L2)) {
        L2 <- plier.base.result$L2
    }
    L1 <- L1 * multiplier
    L2 <- L2 / multiplier
    message(paste0("L1=", L1, "; L2=", L2))

    if (ncol(plier.base.result$B) == ncol(Y)) {
        B <- plier.base.result$B
    }

    oldB <- B

    if (!is.null(rseed)) {
        message("using random start")
        set.seed(rseed)
        B <- t(apply(B, 1, sample))
        Z <- apply(Z, 2, sample)
    }

    U <- matrix(0, nrow = ncol(C), ncol = k)


    round2 <- function(x) {
        signif(x, 4)
    }

    iter.full.start <- iter.full <- u.iter

    curfrac <- 0
    nposlast <- Inf
    npos <- -Inf
    num.U.updates <- 0
    L1k <- L1 * diag(k)
    L2k <- L2 * diag(k)

    if (is_fbm) {
        ZYt <- big_cprodMat(Y, as.matrix(Z))
        ZY <- t(ZYt)
        B <- solve(t(Z) %*% Z + L2k) %*% ZY
    } else {
        B <- solve(t(Z) %*% Z + L2k) %*% mat_mult(t(Z), Y)
    }
    Zraw <- Z
    Z2 <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))

    B <- as.matrix(B)

    for (iter in 1:max.iter) {
        if (iter >= iter.full.start) {
            if (iter >= iter.full && num.U.updates < max.U.updates & iter %% 2 == 1) {
                #    cat(paste(", Updating U"))
                if (any(Zraw < 0)) {
                    stop()
                }
                if (refitAll || num.U.updates %% 5 == 0) {
                    Uprev <- NULL
                    #  cat(", Refitting all U coefficinets")
                } else {
                    Uprev <- U
                }

                if (useRaw) {
                    res <- solveU(Zraw, Chat, C, penalty.factor, pathwaySelection, glm_alpha, maxPath, binary = F, nfolds = cvn, top = top, useNNLS = useNNLS, Uprev = Uprev)
                } else {
                    res <- solveU(Z, Chat, C, penalty.factor, pathwaySelection, glm_alpha, maxPath, binary = F, nfolds = cvn, top = top, useNNLS = useNNLS, Uprev = Uprev)
                }
                U <- res$U

                num.U.updates <- num.U.updates + 1

                iter.full <- iter.full + iter.full.start


                Z2 <- L1 * C %*% U
            }

            curfrac <- (npos <- sum(apply(U, 2, max) > 0)) / k
            # Z1=Y%*%t(B)
            Z1 <- mat_mult(Y, t(B))

            # ii=which(Z2>0)
            # ratio=median(Z2[ii]/abs(Z1[ii]))

            Z <- (Z1 + Z2) %*% solve(tcrossprod(B) + L1k)
        } else {
            Z <- mat_mult(Y, t(B)) %*% solve(tcrossprod(B) + L1k)
        }

        if (adaptive.p > 0) {
            Zraw <- Z
            Zraw[Zraw < 0] <- 0
            cutoffs <- apply(Z, 2, getT)

            for (j in 1:ncol(Z)) {
                Z[Z[, j] < cutoffs[j], j] <- 0
            }
        } else {
            Z[Z < 0] <- 0
            Zraw <- Z
        }

        oldB <- B

        if (is_fbm) {
            ZYt <- big_cprodMat(Y, as.matrix(Z))
            ZY <- t(ZYt)
            B <- solve(t(Z) %*% Z + L2k) %*% ZY
        } else {
            Z_mat <- if (inherits(Z, "matrix")) Z else as.matrix(Z)
            B <- solve(t(Z_mat) %*% Z + L2k) %*% mat_mult(t(Z), Y)
        }

        Bdiff <- sum((B - oldB)^2) / sum(B^2)
        minCor <- min(row_cor(B, oldB))

        BdiffTrace <- c(BdiffTrace, Bdiff)

        if (trace) {
            cat(sprintf("\rProgress %d / %d | Bdiff=%.6f", iter, max.iter, Bdiff))
            flush.console()
        }

        if (iter > 52 && Bdiff > BdiffTrace[iter - 50]) {
            BdiffCount <- BdiffCount + 1
            message("Bdiff is not decreasing")
        } else if (BdiffCount > 1) {
            BdiffCount <- BdiffCount - 1
        }

        if (Bdiff < tol & iter > u.iter + num.U.updates * 2 + 5) {
            cat(sprintf("\rConverged at %d / %d | Bdiff=%.6f, minCor=%.6f\n", iter, max.iter, Bdiff, minCor))
            break
        }
        if (BdiffCount > 5) {
            message(paste0("converged at  iteration ", iter, " Bdiff is not decreasing"))
            break
        }
    }
    rownames(U) <- colnames(priorMat)
    colnames(U) <- rownames(B) <- paste0("LV", 1:k)

    out <- list(B = B, Z = Z, U = U, C = C, L1 = L1, L2 = L2, heldOutGenes = heldOutGenes)

    if (doCrossval) {
        if (adaptive.p != 0) {
            message("Updating Z for CV")
            out$Z <- Zraw
            out$Z[out$Z < 0] <- 0
        }
        message("crossValidation")
        priorMat_m <- as.matrix(priorMat)
        priorMatCV_m <- as.matrix(priorMatCV)
        outAUC <- crossVal(out, priorMat, priorMatCV)
        out$Z <- Z
        out$Uauc <- outAUC$Uauc
        out$Up <- outAUC$Upval
        out$summary <- outAUC$summary
        out$priorMatCV <- priorMatCV
        out$priorMat <- priorMat
        out$withPrior <- which(colSums(out$U) > 0)

        tt <- apply(out$Uauc, 2, max)
        message(paste("There are", sum(tt > 0.70), " LVs with AUC>0.70"))
        message(paste("There are", sum(tt > 0.90), " LVs with AUC>0.90"))
    } else {
        message("Not using cross-validation. No AUCs or p-values")
    }

    # currently not working
    # rownames(out$B)=nameB(out)
    out$call <- call <- match.call()
    out
}

#' Project new data into PLIER latent space
#'
#' Computes the latent loadings \code{B} for new gene expression data using the latent variables
#' \code{Z} from a fitted PLIER model. This allows transfer of the learned latent structure to
#' new datasets with matched genes.
#'
#' @param PLIERres A result object from \code{PLIERfull()} or \code{PLIERbase()}, containing at least \code{Z} and \code{L2}.
#' @param newdata A gene expression matrix (genes × samples) to be projected. Must have the same genes (rows) as \code{PLIERres$Z}.
#'        Can be a standard matrix, sparse matrix, or FBM/big.matrix.
#' @param scale Optional numeric multiplier for the L2 regularization terms. Default is 1.
#'
#' @return A matrix \code{B} of projected latent loadings (LVs × samples) for the new dataset.
#'
#' @details
#' This function uses ridge-regularized least squares to compute \code{B = solve(ZᵗZ + L2·I) · ZᵗY}, where
#' \code{Z} is the latent matrix from the trained PLIER model and \code{Y} is the new dataset.
#' If \code{newdata} is a Filebacked Big Matrix (FBM), the computation is optimized using \code{bigstatsr::big_cprodMat()}.
#'
#' @export
projectPLIER <- function(PLIERres, newdata, scale = 1) {
    stopifnot(nrow(PLIERres$Z) == nrow(newdata))

    # Check if newdata is a FBM/big.matrix object
    is_fbm <- inherits(newdata, c("big.matrix", "FBM"))

    # Convert Matrix package matrices to standard R matrices for compatibility
    Z_matrix <- if (inherits(PLIERres$Z, "Matrix")) as.matrix(PLIERres$Z) else PLIERres$Z

    if (is_fbm) {
        # FBM implementation - use big_cprodMat for efficient computation
        # But ensure Z is a standard matrix, not a Matrix package object
        ZYt <- big_cprodMat(newdata, Z_matrix)
        ZY <- t(ZYt)
    } else {
        # Standard matrix implementation - use regular matrix multiplication
        ZY <- t(Z_matrix) %*% newdata
    }

    # Calculate the regularization matrix using standard matrix operations
    ZtZ <- t(Z_matrix) %*% Z_matrix
    L2k <- PLIERres$L2 * diag(ncol(Z_matrix)) * scale

    # Solve the regularized system
    B <- solve(ZtZ + L2k) %*% ZY

    return(B)
}

## Refactored PC estimation functions

#' Run elbow method to estimate number of PCs
#'
#' @param d Vector of singular values
#' @return Estimated number of PCs via elbow
run_elbow <- function(d) {
    # compute second differences
    x_raw <- abs(diff(diff(d)))
    # smoothing helper
    twiceit <- function(x, twiceit = TRUE) smooth(x, twiceit = TRUE)
    x_smooth <- twiceit(x_raw)
    cutoff <- quantile(x_smooth, 0.5)
    # first index below cutoff plus one for PC count
    which(x_smooth <= cutoff)[1] + 1
}

#' Run permutation method to estimate number of PCs
#'
#' @param data Raw data matrix (row-normalized)
#' @param d Vector of singular values
#' @param B Number of permutations
#' @return Estimated number of PCs via permutation test
run_permutation <- function(data, d, B = 20) {
    k <- length(d)
    # observed proportions
    obs_prop <- d^2 / sum(d^2)
    # permuted proportions matrix
    perm_mat <- matrix(0, nrow = B, ncol = k)
    for (i in seq_len(B)) {
        dat0 <- t(apply(data, 1, sample))
        if (ncol(dat0) == k) {
            uu0 <- svd(dat0)
        } else {
            set.seed(123456)
            uu0 <- rsvd(dat0, k = k, q = 3)
        }
        perm_mat[i, ] <- uu0$d[1:k]^2 / sum(uu0$d[1:k]^2)
    }
    p_vals <- apply(perm_mat >= obs_prop, 2, mean)
    p_vals <- cummax(p_vals)
    sum(p_vals <= 0.1)
}

#' Estimate number of principal components via elbow or permutation method
#'
#' @param data    Either a matrix (e.g. z-scored data) or an SVD result (list with $d).
#' @param method  One of "elbow" (fast) or "permutation" (slower).
#' @param B       Number of permutations (for method = "permutation").
#' @param seed    Seed for reproducibility.
#' @return        Estimated number of PCs.
#' @export
num.pc <- function(data, method = c("elbow", "permutation"), B = 20, seed = NULL) {
    method <- match.arg(method)
    if (!is.null(seed)) set.seed(seed)

    # Prepare SVD result
    if (!inherits(data, "list") || is.null(data$d)) {
        message("Computing SVD")
        # row-normalize
        row_sds <- apply(data, 1, sd)
        row_sds[row_sds == 0] <- 1
        data_norm <- sweep(data, 1, row_sds, "/")
        n <- ncol(data_norm)
        k <- if (n < 500) n else max(200, floor(n / 4))
        if (k == n) {
            uu <- svd(data_norm)
        } else {
            set.seed(123456)
            uu <- rsvd(data_norm, k = k, q = 3)
        }
    } else {
        uu <- data
        data_norm <- NULL # not needed for elbow
        if (method == "permutation") {
            message("Raw data required for permutation; switching to elbow")
            method <- "elbow"
        }
    }

    # Dispatch to specific method
    if (method == "permutation") {
        run_permutation(data_norm, uu$d, B)
    } else {
        run_elbow(uu$d)
    }
}
