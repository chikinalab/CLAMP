#' @useDynLib CLAMP, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Adjust p-values using Benjamini-Hochberg method
#'
#' Applies the BH (Benjamini-Hochberg) correction for
#' multiple hypothesis testing.
#'
#' @param p Numeric vector of p-values.
#' @return Adjusted p-values.
BH <- function(p) {
    p.adjust(p, method = "BH")
}

#' Row-wise correlation between two matrices
#'
#' Computes the Pearson correlation for each row between matrices
#' \code{A} and \code{B}.
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
#' Multiplies two matrices, using optimized multiplication if the first
#' is a Filebacked Big Matrix (FBM).
#'
#' @param mat1 A matrix or an object of class \code{FBM}.
#' @param mat2 A numeric matrix.
#' @param ncores Number of cores to use for parallel computation
#' (only used if mat1 is an FBM). Default is 1.
#' @return Matrix product of \code{mat1} and \code{mat2}.
#' @export
#' @examples
#' set.seed(123)
#' mat1 <- matrix(rnorm(20), nrow = 4)
#' mat2 <- matrix(rnorm(15), nrow = 5)
#' res1 <- mat_mult(mat1, mat2)
#' res1
mat_mult <- function(mat1, mat2, ncores = 1) {
    is_fbm <- inherits(mat1, "FBM")
    if (is_fbm) {
        # For FBM objects, use the specific multiplication method
        return(bigstatsr::big_prodMat(mat1, as.matrix(mat2), ncores = ncores))
    } else {
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
#' @param alpha Non-negative scalar specifying the ridge penalty.
#' A small positive value stabilizes the inversion by shrinking large
#' singular values.
#'
#' @return A numeric matrix representing the ridge-regularized
#' pseudoinverse of \code{m}.
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
#' Ensures consistency in SVD output by flipping signs so that each left
#' singular vector has a majority of positive entries.
#'
#' @param svdres A list as returned by \code{svd()}, with components
#'   \code{u}, \code{d}, and \code{v}.
#' @return A modified \code{svd}-result list where each column of \code{$u}
#'   has been sign-flipped so that its entries sum to a nonnegative value;
#'   \code{$v} is flipped correspondingly.
rotateSVD <- function(svdres) {
    upos <- svdres$u
    uneg <- svdres$u
    upos[upos < 0] <- 0
    uneg[uneg >= 0] <- 0
    uneg <- -uneg
    sumposu <- colSums(upos)
    sumnegu <- colSums(uneg)


    for (i in seq_len(ncol(svdres$u))) {
        if (sumnegu[i] > sumposu[i]) {
            svdres$u[, i] <- -svdres$u[, i]
            svdres$v[, i] <- -svdres$v[, i]
        }
    }
    svdres
}

#' Binarize matrix by top-k values per column
#'
#' Keeps only the top \code{top} values in each column of a matrix, setting
#' others to 0.
#'
#' @param Z A numeric matrix.
#' @param top Number of top entries to keep in each column.
#' @param keepVals If \code{TRUE}, retains original values above the cutoff;
#'   otherwise, sets them to 1.
#' @return A modified matrix with only top entries retained per column.
binarizeTop <- function(Z, top, keepVals = TRUE) {
    for (i in seq_len(ncol(Z))) {
        cutoff <- sort(Z[, i], decreasing = TRUE)[top + 1]
        if (cutoff == 0) {
            cutoff <- min(Z[Z[, i] > 0, i])
        }
        Z[Z[, i] < cutoff, i] <- 0
        if (!keepVals) Z[Z[, i] > 0, i] <- 1
    }
    Z
}

#' Fit the loading matrix Z using sparse regression of prior information U
#'
#' For each column of a target matrix \code{Z} using pathway or prior
#' annotation \code{priorMat}. It performs regularization selection using
#' cross-validation and can apply either Supports continuous or binary
#' response models. Relaxed refitting is supported for final coefficient
#' estimation.
#'
#' @param Z A numeric matrix with features (rows) and samples (columns).
#' @param Chat (Optional) Precomputed pseudo-inverse of \code{priorMat};
#'   if \code{NULL}, it is calculated using ridge regularization.
#' @param priorMat A numeric matrix with prior information
#'   (features x pathways).
#' @param penalty.factor Optional penalty weights for features in
#'   \code{priorMat}.
#' @param pathwaySelection Method to select candidate pathways: \code{"fast"}
#'   (default) or \code{"complete"}.
#' @param alpha Elastic net mixing parameter (0 = ridge, 1 = lasso).
#'   Default is 0.9.
#' @param maxPath Maximum number of pathways/features selected per column.
#'   Default is 10.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param useSE Whether to use the 1-standard-error rule for lambda selection.
#'   Default is \code{FALSE}.
#' @param top If set, sets to 0 all but the top entries of \code{Z} per
#'   column before fitting.
#' @param binary If \code{TRUE}, fits a binomial model (e.g., classification)
#'   to \code{Z}>0. Can be used in combination with \code{top}.
#'   Default is \code{FALSE}.
#' @param nlambda Number of lambda values for glmnet. Default is 20.
#' @param scale Whether to standardize predictors in glmnet.
#'   Default is \code{TRUE}.
#' @param refit Whether to perform relaxed refitting using selected
#'   predictors. Default is \code{TRUE}.
#' @param Uprev (Optional) Previous U matrix to reuse. In this mode only
#'   the columns of \code{U} that are all zero are estimated. Used
#'   internally in \code{CLAMP}.
#' @param ... Additional arguments passed to \code{glmnet()} or
#'   \code{cv.glmnet()}.
#' @param useAUC Logical; whether to compute pathway-LV associations using
#'   AUC (default TRUE) instead of OLS.
#' @param intercept Logical; whether to include an intercept term in glmnet
#'   models. Default is TRUE.
#' @return A list with one element:
#' \describe{
#'   \item{\code{U}}{A matrix of loadings (features x components).
#'    Columns are named \code{LV1}, \code{LV2}, ...}
#' }
#' @export
#' @examples
#' set.seed(123)
#' genes <- paste0("G", 1:200)
#' lvs <- paste0("LV", 1:4)
#' paths <- paste0("Path", 1:60)
#'
#' Z <- matrix(rnorm(200 * 4), nrow = 200, dimnames = list(genes, lvs))
#'
#' priorMat <- matrix(rbinom(200 * 60, 1, 0.07),
#'     nrow = 200, dimnames = list(genes, paths)
#' )
#'
#' fit1 <- solveU(
#'     Z = Z,
#'     priorMat = priorMat,
#'     pathwaySelection = "fast",
#'     alpha = 0.9,
#'     maxPath = 10,
#'     nfolds = 5,
#'     binary = FALSE,
#'     refit = TRUE
#' )
solveU <- function(
  Z, Chat = NULL, priorMat, penalty.factor, pathwaySelection = "fast",
  alpha = 0.9, maxPath = 10, nfolds = 5, useSE = FALSE, top = NULL,
    binary = FALSE, nlambda = 20, scale = TRUE, refit = TRUE,
    Uprev = NULL, useAUC = TRUE, intercept = TRUE, ...) {
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
    if (is.null(Chat) & !useAUC) {
        Chat <- pinv.ridge(crossprod(priorMat), 5) %*% (t(priorMat))
    }

    if (!useAUC) {
        Ur <- Chat %*% Z # get U by OLS
    } else {
        Ur <- t(allAgainstAllAUCs(Z, priorMat))
    }
    Ur <- apply(-Ur, 2, rank) # rank
    Urm <- apply(Ur, 1, min)

    # Zhat <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
    if (is.null(Uprev)) {
        U <- Matrix::Matrix(
            0,
            nrow = ncol(priorMat), ncol = ncol(Z), sparse = TRUE
        )
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
        message("Picked ", length(iip), " pathways")
    }

    for (i in seq_len(ncol(Z))) {
        if (all(U[, i] == 0)) {
            if (var(Z[, i]) < 1e-9) {
                # message("skipping")
                next
            }
            if (pathwaySelection == "fast") {
                iip <- which(Ur[, i] <= maxPath)
            }

            if (!binary) { # not doing a binary prediction
                gres <- cv.glmnet(
                    y = Z[, i], x = priorMat[, iip],
                    alpha = alpha, lower.limits = 0,
                    foldid = ((seq_len(nrow(Z))) %% nfolds) + 1,
                    keep = TRUE, nfolds = nfolds,
                    standardize = scale, dfmax = maxPath,
                    nlambda = nlambda, intercept = intercept, ...
                )

                # plot(gres)
            } else {
                gres <- cv.glmnet(
                    y = (Z[, i] > 0) + 1 - 1, x = priorMat[, iip],
                    family = "binomial", alpha = alpha, lower.limits = 0,
                    foldid = ((seq_len(nrow(Z))) %% nfolds) + 1,
                    keep = FALSE,
                    nfolds = nfolds, type.measure = "auc",
                    dfmax = maxPath, nlambda = nlambda,
                    standardize = scale, intercept = intercept, ...
                )
            }
        } # end if all U[,i]==0
        if (refit) {
            if (any(U[, i] != 0)) {
                selected_features <- which(U[, i] != 0)
            } else {
                s_best <- if (useSE) gres$lambda.1se else gres$lambda.min
                active_coef <- coef(gres, s = s_best)
                active_ix <- which(active_coef[-1] != 0)
                selected_features <- iip[active_ix]
            }

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
                # Zhat[, i] <- predict(fit_relaxed, newx = X_sel,
                #                      s = 0, type = "response")[,1]
            }
        } else {
            s_best <- if (useSE) gres$lambda.1se else gres$lambda.min
            coef_vec <- as.vector(coef(gres, s = s_best))[-1] # drop intercept
            betaI <- which(coef_vec != 0)
            U[iip[betaI], i] <- coef_vec[betaI]
        }
        # end for i in Z
    }
    # message(sprintf(", Number of annotated columns is %d",
    #                 sum(Matrix::colSums(U) > 0)))
    # rownames(U)=substr(colnames(priorMat),1,30)
    rownames(U) <- colnames(priorMat)
    colnames(U) <- colnames(Z)
    # flag this for dicussion
    U <- as.matrix(U)

    return(list(U = U))

    message("Number of annotated columns is ", sum(Matrix::colSums(U) > 0))

    rownames(U) <- colnames(priorMat)
    colnames(U) <- paste0("LV", seq_len(ncol(U)))
    return(U)
}

#' Compute Chat matrix from prior annotation
#'
#' Computes the transformation matrix \code{Chat} used to map from observed
#' data to latent space, based on a pseudo-inverse of the prior annotation
#' matrix. Optionally standardizes the columns
#' of \code{priorMat} before computing.
#'
#' @param priorMat A numeric or sparse matrix (features x pathways)
#'   containing prior annotations.
#' @param scale Logical; if \code{TRUE} (default), standardizes the columns
#'   of \code{priorMat}
#' before computing \code{Chat}.
#'
#' @return A numeric matrix \code{Chat} of dimensions (pathways x features).
#' @examples
#' # simple toy prior: 3 features x 2 pathways
#' priorMat <- matrix(
#'     c(
#'         1, 0, 1,
#'         0, 1, 0
#'     ),
#'     nrow = 3, ncol = 2,
#'     dimnames = list(
#'         paste0("gene", seq_len(3)),
#'         paste0("path", seq_len(2))
#'     )
#' )
#' # compute Chat (2 pathways x 3 features)
#' Chat <- getChat(priorMat)
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
#' @param pathMat A sparse binary matrix of genes (rows) x pathways (columns).
#' @param new.genes Character vector of gene names to match.
#' @param min.genes Minimum number of overlapping genes required to keep
#'   a pathway.
#'
#' @return A sparse matrix of dimensions \code{length(new.genes)} x
#'   filtered pathways.
#' @examples
#' library(Matrix)
#' # create a toy gene-by-pathway sparse matrix
#' genes <- paste0("g", seq_len(6))
#' pathways <- c("Path1", "Path2", "Path3")
#' # Path1: g1, g2; Path2: g2, g3, g4; Path3: g5
#' pathMat <- sparseMatrix(
#'     i = c(1, 2, 2, 3, 4, 5),
#'     j = c(1, 1, 2, 2, 2, 3),
#'     dims = c(length(genes), length(pathways)),
#'     dimnames = list(genes, pathways)
#' )
#' new.genes <- genes
#' filtered <- getMatchedPathwayMat(pathMat, new.genes, min.genes = 2)
#' @export
getMatchedPathwayMat <- function(pathMat, new.genes, min.genes = 10) {
    cm <- intersect(rownames(pathMat), new.genes)
    mymessage(
        "There are ", length(cm),
        " genes in the intersection between data and prior"
    )

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

#' Subset and filter multiple pathway matrices to match target genes
#'
#' Filters gene-by-pathway annotation matrices to retain only pathways
#' with sufficient overlap with a given gene set. The result is a sparse matrix
#' aligned to \code{new.genes}, combining all inputs column-wise.
#'
#' @param ... One or more sparse binary matrices (genes x pathways).
#' @param new.genes Character vector of gene names to match.
#' @param min.genes Minimum number of overlapping genes required to keep
#'   a pathway.
#'
#' @return A sparse matrix with rows = \code{new.genes} and columns =
#'   filtered pathways from all inputs.
getMatchedPathwayMat2 <- function(..., new.genes, min.genes = 10) {
    pathMats <- list(...)

    filtered <- lapply(pathMats, function(pathMat) {
        cm <- intersect(rownames(pathMat), new.genes)
        message(
            "There are ", length(cm),
            " genes in the intersection between data and prior"
        )

        matchPathMat <- Matrix::sparseMatrix(
            i = match(cm, new.genes),
            j = rep(1, length(cm)), # temporary
            dims = c(length(new.genes), ncol(pathMat)),
            x = 0,
            dimnames = list(new.genes, colnames(pathMat))
        )
        matchPathMat[cm, ] <- pathMat[cm, ]

        genesInPath <- Matrix::colSums(matchPathMat)

        ii <- which(genesInPath >= min.genes)

        message(sprintf(
            "Removing %d pathways", ncol(matchPathMat) - length(ii)
        ))

        matchPathMat[, ii, drop = FALSE]
    })

    if (length(filtered) == 1) {
        return(filtered[[1]])
    }
    do.call(Matrix::cBind, filtered)
}

#' Subset and filter pathway matrix to match target genes
#'
#' Filters a gene-by-pathway annotation matrix to retain only pathways
#' with sufficient overlap with a given gene set. The result is a sparse matrix
#' aligned to \code{new.genes}, with columns (pathways) retained only if they
#' have at least \code{min.genes} matched genes.
#'
#' @param pathMat A sparse binary matrix of genes (rows) x pathways (columns).
#' @param new.genes Character vector of gene names to match.
#' @param min.genes Minimum number of overlapping genes required to keep
#'   a pathway.
#'
#' @return A sparse matrix of dimensions \code{length(new.genes)} x
#'   filtered pathways.
getMatchedPathwayMatOld <- function(pathMat, new.genes, min.genes = 10) {
    cm <- intersect(rownames(pathMat), new.genes)
    mymessage(
        "There are ", length(cm),
        " genes in the intersection between data and prior"
    )

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
#' Computes the area under the ROC curve (AUC) by applying a Wilcoxon
#' rank-sum test between predicted values for positive and negative labels.
#' This is equivalent to
#' computing the Mann-Whitney U statistic.
#'
#' @param labels A numeric or logical vector indicating class labels.
#'   Values > 0 are treated as positive.
#' @param values A numeric vector of prediction scores corresponding to
#'   \code{labels}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{auc}}{Estimated AUC, or 0.5 if one class is missing}
#'   \item{\code{pval}}{Wilcoxon test p-value, or \code{NA} if one class
#'     is missing}
#' }
AUC <- function(labels, values) {
    pos <- labels > 0
    neg <- !pos
    posn <- sum(pos)
    negn <- sum(neg)

    if (posn > 0 && negn > 0) {
        res <- stats::wilcox.test(values[pos], values[neg],
            alternative = "greater", exact = FALSE
        )
        auc <- unname(res$statistic) / (posn * negn)
        pval <- res$p.value
    } else {
        auc <- 0.5
        pval <- NA
    }
    list(auc = auc, pval = pval, npos = posn, nneg = negn)
}

#' Cross-validation AUC for CLAMP latent variables and pathways
#'
#' Evaluates how well each latent variable in a CLAMP model captures
#' held-out pathway annotations, using cross-validation over the prior
#' matrix. For each latent variable and associated pathway, held-out genes
#' are selected and the AUC is computed using their scores in
#' \code{clampRes$Z}.
#'
#' @param clampRes A list containing \code{U} (loadings) and \code{Z}
#'   (scores) from a CLAMP model.
#' @param priorMat A binary matrix (genes x pathways) indicating original
#'   pathway annotations.
#' @param priorMatcv A version of \code{priorMat} used to mask held-out
#'   annotations for cross-validation.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{Uauc}}{Matrix of AUC values (pathways x LVs)}
#'   \item{\code{Upval}}{Matrix of \code{-log10(p)} values (pathways x LVs)}
#'   \item{\code{summary}}{Data frame with pathway, LV index, AUC, p-value,
#'     and FDR}
#' }
crossVal <- function(clampRes, priorMat, priorMatcv) {
    ii <- which(Matrix::colSums(clampRes$U) > 0)

    Uauc <- Matrix::Matrix(
        0,
        nrow = nrow(clampRes$U), ncol = ncol(clampRes$U), sparse = TRUE
    )
    Up <- Matrix::Matrix(
        0,
        nrow = nrow(clampRes$U), ncol = ncol(clampRes$U), sparse = TRUE
    )

    results <- list()

    for (i in ii) {
        iipath <- which(clampRes$U[, i] > 0)

        for (j in iipath) {
            iiheldout <- which(
                (rowSums(priorMat[, iipath, drop = FALSE]) == 0) |
                    (priorMat[, j] > 0 & priorMatcv[, j] == 0)
            )

            aucres <- AUC(priorMat[iiheldout, j], clampRes$Z[iiheldout, i])


            results[[length(results) + 1]] <- data.frame(
                pathway = colnames(priorMat)[j],
                LV = paste0("LV", i),
                AUC = as.numeric(aucres$auc),
                p_value = as.numeric(aucres$pval),
                FDR = NA_real_,
                npos = aucres$npos,
                nneg = aucres$nneg,
                stringsAsFactors = FALSE
            )

            Uauc[j, i] <- aucres$auc
            Up[j, i] <- -log10(aucres$pval)
        }
    }

    out <- do.call(rbind, results)
    out$FDR <- as.numeric(BH(out$p_value))
    out$AUC <- as.numeric(out$AUC)

    return(list(Uauc = Uauc, Upval = Up, summary = out))
}

#' CLAMP base matrix factorization
#'
#' Runs the core matrix factorization procedure of CLAMP,
#' decomposing the gene expression matrix \code{Y} into latent variables
#' \code{Z} and loadings \code{B}. It supports sparse, dense, and Filebacked
#' Big Matrices (FBM) as input and includes options for
#' adaptive sparsity, positive constraints, and regularization.
#'
#' @param Y Input gene expression matrix (genes x samples). Can be dense,
#'   sparse (\code{dgCMatrix}), or FBM.
#' @param clamp_k Number of latent variables for CLAMP (final model rank).
#'   If \code{NULL}, it is chosen automatically via \code{select_clamp_k()}.
#' @param svd_k Number of singular values/components to compute in the SVD.
#'   If \code{NULL}, defaults to \code{max(2, min(n_genes, n_samples) - 1)}.
#' @param svdres Optional precomputed SVD result. If not supplied, it is
#'   computed internally.
#' @param L1 L1 regularization strength for Z. Defaults to scaled singular
#'   value.
#' @param L2 L2 regularization strength for B. Defaults to scaled singular
#'   value.
#' @param Zpos Logical; if \code{TRUE}, negative entries in Z are zeroed.
#'   Default is \code{TRUE}.
#' @param max.iter Maximum number of optimization iterations. Default is 200.
#' @param tol Convergence tolerance for B update. Default is 5e-4.
#' @param trace Logical; if \code{TRUE}, prints progress. Default is
#'   \code{FALSE}.
#' @param rseed Optional integer for reproducible random initialization of B.
#' @param B Optional initial matrix for B. If not provided, initialized
#'   from SVD.
#' @param scale Scaling factor for L1 and L2 when not provided. Default is 1.
#' @param pos.adj Positive constraint adjustment divisor for L1. Default is 3.
#' @param adaptive.p Controls adaptive sparsity in \code{Z}. After each ALS
#'   update,
#' negative entries in \code{Z} are assumed to reflect noise.
#' The cutoff for thresholding is set according to the probability of
#' positive values under a reflected negative distribution-effectively
#' zeroing out small positive entries likely to be noise.
#' Smaller values lead to more sparsity. Default is 0.05.
#' @param adaptive.iter Number of iterations before adaptive sparsity is
#'   applied. Default is 20.
#' @param cutoff Scalar threshold to zero Z values when \code{Zpos = TRUE}
#'   and adaptive thresholding
#' is not used. Default is 0.
#' @param ncores Number of cores to use for parallel computation (only used
#'   if Y is an FBM). Default is 1.
#' @param clamp_k_method Method for selecting `clamp_k` when not provided.
#'   One of `"elbow"` (default), `"permutation"`, `"gavish_donoho"`, or
#'   `"scaleSVs"`. Passed to [select_clamp_k()].
#' @return A list with components:
#' \describe{
#'   \item{\code{B}}{Latent variable loadings (LVs x genes)}
#'   \item{\code{Z}}{Latent variable scores (LVs x samples)}
#'   \item{\code{Zraw}}{Raw Z matrix before thresholding}
#'   \item{\code{L1}}{Final value of L1 used}
#'   \item{\code{L2}}{Final value of L2 used}
#' }
#'
#' @details
#' This function is the low-level implementation of CLAMP. It alternates
#' between solving for \code{Z} given \code{B} and solving for \code{B}
#' given \code{Z}, with optional sparsity and non-negativity constraints on
#' \code{Z}. Convergence is assessed via relative change in \code{B}.
#'
#' @examples
#' # small toy dataset: 5 genes x 4 samples
#' Y <- matrix(rnorm(5 * 4), nrow = 5, ncol = 4)
#' # run a single iteration for speed
#' res <- CLAMPbase(Y, clamp_k = 2, max.iter = 1, trace = FALSE)
#' # inspect dimensions of B and Z
#' dim(res$B)
#' dim(res$Z)
#'
#' @export
CLAMPbase <- function(
  Y, clamp_k = NULL, svd_k = NULL, svdres = NULL, L1 = NULL, L2 = NULL,
  Zpos = TRUE, max.iter = 200, tol = 5e-4, trace = FALSE,
  rseed = NULL, B = NULL, scale = 1, pos.adj = 3,
  adaptive.p = 0.05, adaptive.iter = 20,
  cutoff = 0, ncores = 1, clamp_k_method = "elbow"
) {
    if (ncores > 1) {
        # if we are parallelizing, then disable BLAS parallelization
        options(bigstatsr.check.parallel.blas = FALSE)
        blas_nproc <- getOption("default.nproc.blas")
        options(default.nproc.blas = NULL)
    }

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

    if (is.null(svd_k)) {
        svd_k <- select_svd_k(Y)
        if (!is.null(clamp_k)) svd_k <- max(svd_k, clamp_k)
    }

    if (is.null(svdres) && is.null(B)) {
        message("Computing SVD")
        svdres <- compute_svd(Y, k = svd_k)
    }

    svdres <- rotateSVD(svdres)

    if (is.null(clamp_k)) {
        clamp_k <- select_clamp_k(svdres,
            n_samples = ncol(Y), svd_k = svd_k,
            method = clamp_k_method, data = Y
        )
        d <- svdres$d[clamp_k]
    } else {
        d <- svdres$d[clamp_k]
    }

    message("CLAMP k is set to ", clamp_k)

    if (is.null(L1)) {
        # L1 <- svdres$d[k] * scale
        L1 <- d * scale
        if (!is.null(pos.adj)) {
            L1 <- L1 / pos.adj
        }
    }

    if (is.null(L2)) {
        #   L2 <- svdres$d[k] * scale
        L2 <- d * scale
    }

    L2k <- L2 * diag(clamp_k)
    #    L1 <- svdres$d[k]/2*scale
    message("L1 is set to ", L1)
    message("L2 is set to ", L2)

    if (is.null(B)) {
        # initialize B with svd

        B <- t(
            svdres$v[, seq_len(clamp_k)] %*%
                diag(sqrt(svdres$d[seq_len(clamp_k)]))
        )

        # alternative initializations
        # seem to be not as good
        #   B <- t(svdres$v[seq_len(ncol(Y)), seq_len(k)] %*%
        #          diag(svdres$d[seq_len(k)]))
        #   B <- t(svdres$v[seq_len(ncol(Y)), seq_len(k)])
    } else {
        message("B given")
    }

    if (!is.null(rseed)) {
        message("Using random start")
        B <- t(apply(B, 1, sample))
    }

    round2 <- function(x) {
        signif(x, 4)
    }

    getT <- function(x) {
        -quantile(x[x < 0], adaptive.p)
    }

    for (i in seq_len(max.iter)) {
        # main loop
        Zraw <- Z <- mat_mult(Y, t(B), ncores = ncores) %*%
            solve(tcrossprod(B) + L1 * diag(clamp_k))

        if (i >= adaptive.iter && adaptive.p > 0) {
            cutoffs <- apply(Zraw, 2, getT)

            for (j in seq_len(ncol(Z))) {
                Z[Z[, j] < cutoffs[j], j] <- 0
            }
        } else if (Zpos) {
            Z[Z < cutoff] <- 0
        }

        oldB <- B

        if (is_fbm) {
            ZYt <- big_cprodMat(Y, as.matrix(Z), ncores = ncores)
            ZY <- Matrix::t(ZYt)
            B <- solve(Matrix::t(Z) %*% Z + L2k) %*% ZY
        } else {
            B <- solve(Matrix::t(Z) %*% Z + L2k) %*%
                mat_mult(Matrix::t(Z), Y, ncores = ncores)
        }

        # update error
        Bdiff <- sum((B - oldB)^2) / sum(B^2)

        # keeping track of this in case this can be useful for convergence
        minCor <- min(row_cor(B, oldB))

        BdiffTrace <- c(BdiffTrace, Bdiff)

        if (trace) {
            message(sprintf(
                "\rProgress %d / %d | Bdiff=%.6f, minCor=%.6f",
                i, max.iter, Bdiff, minCor
            ))
            flush.console()
        }

        # check for convergence
        if (i > 52 && Bdiff > BdiffTrace[i - 50]) {
            BdiffCount <- BdiffCount + 1
        } else if (BdiffCount > 1) {
            BdiffCount <- BdiffCount - 1
        }

        if (Bdiff < tol && i > adaptive.iter + 10) {
            message(sprintf(
                "Converged at iteration= %d | Bdiff=%.6f,  tol=%.6f     ",
                i, Bdiff, tol
            ))
            break
        }
        if (BdiffCount > 5 && i > adaptive.iter + 10) {
            message("stopped at iteration ", i, " Bdiff is not decreasing")
            break
        }
    }

    rownames(B) <- colnames(Z) <- paste0("LV", seq_len(clamp_k))
    if (!is.null(rownames(Y))) rownames(Z) <- rownames(Y)
    if (!is.null(colnames(Y))) colnames(B) <- colnames(Y)

    if (ncores > 1) {
        # restore previous state
        options(bigstatsr.check.parallel.blas = TRUE)
        options(default.nproc.blas = blas_nproc)
    }

    return(list(
        B = as.matrix(B), Z = as.matrix(Z),
        Zraw = Zraw, L1 = L1, L2 = L2
    ))
}

#' Full CLAMP model with prior information and cross-validation
#'
#' Runs the full CLAMP model using a gene expression matrix
#' and prior pathway annotation matrix. This function performs latent
#' variable decomposition guided by prior knowledge and includes optional
#' cross-validation to evaluate pathway associations.
#'
#' @param Y Gene expression matrix (genes x samples). Can be dense, sparse
#'   (dgCMatrix), or FBM.
#' @param priorMat Binary matrix (genes x pathways) representing prior
#'   annotations.
#' @param svdres Optional SVD result used for initialization.
#' @param clamp.base.result Optional result from \code{CLAMPbase()} to
#'   initialize B.
#' @param clamp_k Number of latent variables for CLAMP (final model rank).
#'   If \code{NULL}, it is chosen automatically via \code{select_clamp_k()}.
#' @param svd_k Number of singular values/components to compute in the SVD.
#'   If \code{NULL}, defaults to \code{max(2, min(n_genes, n_samples) - 1)}.
#' @param L1 Regularization strength for Z. If \code{NULL}, initialized from
#'   SVD or \code{clamp.base.result}.
#' @param L2 Regularization strength for B. If \code{NULL}, initialized from
#'   SVD or \code{clamp.base.result}.
#' @param top If set, keeps only top-n values per column in Z during U updates.
#' @param cvn Number of folds for cross-validation in U updates. Default is 5.
#' @param max.iter Maximum number of iterations. Default is 350.
#' @param trace Logical; if \code{TRUE}, prints iteration progress.
#' @param Chat Optional precomputed matrix for solving U.
#' @param maxPath Maximum number of pathways/features selected per LV.
#'   Default is 10.
#' @param doCrossval Whether to perform pathway-level cross-validation.
#'   Default is \code{TRUE}.
#' @param penalty.factor Vector of feature-specific penalties for glmnet.
#'   Default: all ones.
#' @param glm_alpha Elastic net mixing parameter for glmnet. Default is 0.9.
#' @param minGenes Minimum number of genes per pathway to retain. Default is 10.
#' @param tol Convergence tolerance on relative change in B. Default is 5e-4.
#' @param seed Seed for reproducibility of cross-validation masking.
#'   Default is 123456.
#' @param allGenes If \code{TRUE}, zero-fills \code{priorMat} for genes not
#'   present. Default is \code{FALSE}.
#' @param rseed Optional seed for randomly reinitializing B and Z.
#' @param max.U.updates Maximum number of U updates. Default is 5.
#' @param pathwaySelection Pathway selection mode: \code{"fast"} or
#'   \code{"complete"}.
#' @param multiplier Scaling factor for adjusting L1 and L2.
#' @param adaptive.p Quantile threshold for adaptively zeroing small Z
#'   values. After each update, the \code{adaptive.p}-quantile of negative
#'   entries in Z is used (flipped positive) to threshold
#'   small positive values, assuming they reflect noise. Default is 0.05.
#' @param useNNLS If \code{TRUE}, uses non-negative least squares in U
#'   estimation. Default is \code{TRUE}.
#' @param useRaw If \code{TRUE}, uses unthresholded Z for solving U.
#'   Default is \code{TRUE}.
#' @param refitAll If \code{TRUE}, refits all U columns every update.
#'   Default is \code{FALSE}.
#' @param useSE Logical; passed to the internal \code{solveU()} call. If
#'   \code{TRUE}, enables standard-error-aware selection when fitting U
#'   (pathway coefficients). Default is \code{FALSE}.
#' @param ncores Number of cores to use for parallel computation (only used
#'   if Y is an FBM). Default is 1.
#' @param clamp_k_method Method for selecting `clamp_k` when not provided.
#'   One of `"elbow"` (default), `"permutation"`, `"gavish_donoho"`, or
#'   `"scaleSVs"`. Passed to [select_clamp_k()].
#' @return A list with the following components:
#' \describe{
#'   \item{\code{B}}{Latent variable loadings (LVs x genes)}
#'   \item{\code{Z}}{Latent variable matrix (LVs x samples)}
#'   \item{\code{U}}{Pathway loadings matrix (pathways x LVs)}
#'   \item{\code{C}}{Masked prior matrix used for training}
#'   \item{\code{L1}, \code{L2}}{Regularization parameters}
#'   \item{\code{heldOutGenes}}{List of held-out genes per pathway (if CV
#'     is enabled)}
#'   \item{\code{Uauc}}{AUC matrix from CV evaluation (if enabled)}
#'   \item{\code{Up}}{-\code{log10(p)} values from CV evaluation (if enabled)}
#'   \item{\code{summary}}{Data frame of AUC, p-values, and FDR per pathway
#'     x LV (if enabled)}
#'   \item{\code{priorMatCV}}{Masked prior matrix used during CV}
#'   \item{\code{priorMat}}{Final filtered prior matrix}
#'   \item{\code{withPrior}}{Indices of LVs with non-zero pathway loadings}
#'   \item{\code{call}}{Function call}
#' }
#'
#' @details
#' The model alternates between solving \code{Z}, \code{B}, and \code{U}.
#' Adaptive sparsity is applied to \code{Z} using a dynamic threshold based
#' on the negative tail of its distribution. Cross-validation is used to
#' hold out gene annotations in \code{priorMat} and evaluate latent variable
#' specificity.
#'
#' @examples
#' mat <- matrix(rnorm(100), 10, 10)
#' svdres <- rsvd::rsvd(mat, k = 5)
#' base <- CLAMPbase(Y = mat, clamp_k = 5, svdres = svdres, trace = FALSE)
#' priorMat <- matrix(1, nrow(mat), 5)
#' full <- CLAMPfullnVP(
#'     Y = mat, priorMat = priorMat, svdres = svdres,
#'     clamp.base.result = base, clamp_k = 5,
#'     doCrossval = FALSE, trace = FALSE, max.U.updates = 0
#' )
#' @export
CLAMPfullnVP <- function(
  Y, priorMat, svdres = NULL, clamp.base.result = NULL, clamp_k = NULL,
  svd_k = NULL, L1 = NULL, L2 = NULL, top = NULL,
  cvn = 5, max.iter = 350, trace = FALSE, Chat = NULL,
  maxPath = 10, doCrossval = TRUE,
  penalty.factor = rep(1, ncol(priorMat)), glm_alpha = 0.9,
  minGenes = 10, tol = 5e-4, seed = 123456,
  allGenes = FALSE, rseed = NULL,
  max.U.updates = 5, pathwaySelection = c("fast"), multiplier = 1,
  adaptive.p = 0.05, useNNLS = TRUE, useRaw = TRUE,
  refitAll = FALSE, useSE = FALSE, ncores = 1,
  clamp_k_method = "elbow") {
    if (ncores > 1) {
        # if we are parallelizing, then disable BLAS parallelization
        options(bigstatsr.check.parallel.blas = FALSE)
        blas_nproc <- getOption("default.nproc.blas")
        options(default.nproc.blas = NULL)
    }

    getT <- function(x) {
        -quantile(x[x < 0], adaptive.p)
    }

    pathwaySelection <- match.arg(pathwaySelection, c("complete", "fast"))

    priorMat <- as.matrix(priorMat)

    message("**CLAMPfullnVP v2 **")

    # Detect matrix type
    is_fbm <- inherits(Y, "FBM")
    is_sparse <- inherits(Y, "dgCMatrix")

    if (nrow(priorMat) != nrow(Y) || !all(rownames(priorMat) == rownames(Y))) {
        if (!allGenes) {
            cm <- commonRows(Y, priorMat)
            message("Selecting common genes: ", length(cm))
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
        message("Removed ", length(iibad), " pathways with too few genes")
    }
    if (doCrossval) {
        priorMatCV <- as.matrix(priorMat)
        if (!is.null(seed)) {
            warning(
                "`seed` is deprecated and ignored. Use set.seed(seed) ",
                "before calling this function.",
                call. = FALSE
            )
        }
        for (j in seq_len(ncol(priorMatCV))) {
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

    # YsqSum=sum(Y^2)
    # compute svd and use that as the starting point

    if (!is.null(svdres) && nrow(svdres$v) != ncol(Y)) {
        message("SVD V has the wrong number of columns")
        svdres <- NULL
    }

    if (is.null(svd_k) && is.null(clamp.base.result)) {
        svd_k <- select_svd_k(Y)
        if (!is.null(clamp_k)) svd_k <- max(svd_k, clamp_k)
    }

    if (is.null(svdres) && is.null(clamp.base.result)) {
        message("Computing SVD")
        svdres <- compute_svd(Y, k = svd_k)
    }

    if (is.null(svdres) && is.null(clamp.base.result)) {
        svdres <- rotateSVD(svdres)
    }

    if (is.null(clamp.base.result)) {
        clamp_k <- select_clamp_k(svdres,
            n_samples = ncol(Y), svd_k = svd_k,
            method = clamp_k_method, data = Y
        )
        d <- svdres$d[clamp_k]
    } else {
        d <- svdres$d[clamp_k]
    }

    if (is.null(clamp.base.result)) {
        message("Running CLAMPbase")
        clamp.base.result <- CLAMPbase(Y, clamp_k = clamp_k, svdres = svdres)
    } else {
        message("using provided CLAMPbase result")
        if (nrow(Y) != nrow(clamp.base.result$Z)) {
            if (is.null(rownames(Y)) | is.null(rownames(clamp.base.result$Z))) {
                stop(
                    "Y and clamp.base.result$Z must have equal row ",
                    "numbers or row names"
                )
            }
            clamp.base.result$Z <- clamp.base.result$Z[rownames(Y), ]
        }
        clamp_k <- ncol(clamp.base.result$Z)
    }

    message("CLAMP k is set to ", clamp_k)

    Z <- clamp.base.result$Z

    if (is.null(L1)) {
        L1 <- clamp.base.result$L1
    }
    if (is.null(L2)) {
        L2 <- clamp.base.result$L2
    }
    L1 <- L1 * multiplier
    L2 <- L2 / multiplier
    message("L1=", L1, "; L2=", L2)

    if (ncol(clamp.base.result$B) == ncol(Y)) {
        B <- clamp.base.result$B
    }

    oldB <- B

    if (!is.null(rseed)) {
        message("Using random start")
        # reproducibility is controlled by user calling set.seed before
        # this function
        B <- t(apply(B, 1, sample))
        Z <- apply(Z, 2, sample)
    }

    U <- matrix(0, nrow = ncol(C), ncol = clamp_k)

    round2 <- function(x) {
        signif(x, 4)
    }

    u.iter <- 2
    curfrac <- 0
    nposlast <- Inf
    npos <- -Inf
    num.U.updates <- 0
    L1k <- L1 * diag(clamp_k)
    L2k <- L2 * diag(clamp_k)

    if (is_fbm) {
        ZYt <- big_cprodMat(Y, as.matrix(Z), ncores = ncores)
        ZY <- Matrix::t(ZYt)
        B <- solve(Matrix::t(Z) %*% Z + L2k) %*% ZY
    } else {
        B <- solve(Matrix::t(Z) %*% Z + L2k) %*%
            mat_mult(Matrix::t(Z), Y, ncores = ncores)
    }

    Zraw <- Z
    Z2 <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))

    B <- as.matrix(B)

    for (iter in seq_len(max.iter)) {
        if (iter >= u.iter) { # actually do a U iteration
            if (num.U.updates < max.U.updates & iter %% 2 == 1) {
                #  Updating U
                if (any(Zraw < 0)) {
                    stop()
                }
                if (refitAll || num.U.updates %% 5 == 0) {
                    Uprev <- NULL
                    #  mRefitting all U coefficinets
                } else {
                    Uprev <- U
                }

                if (useRaw) {
                    res <- solveU(
                        Zraw, Chat, C,
                        penalty.factor, pathwaySelection,
                        glm_alpha, maxPath,
                        binary = FALSE, nfolds = cvn, top = top,
                        useNNLS = useNNLS, Uprev = Uprev, useSE = useSE
                    )
                } else {
                    res <- solveU(Z, Chat, C, penalty.factor, pathwaySelection,
                        glm_alpha, maxPath,
                        binary = FALSE, nfolds = cvn, top = top,
                        useNNLS = useNNLS, Uprev = Uprev, useSE = useSE
                    )
                }
                U <- res$U

                num.U.updates <- num.U.updates + 1

                iter.full <- iter.full + iter.full.start


                Z2 <- L1 * C %*% U
            }

            curfrac <- (npos <- sum(apply(U, 2, max) > 0)) / clamp_k
            # Z1=Y%*%t(B)
            Z1 <- mat_mult(Y, t(B), ncores = ncores)

            # ii=which(Z2>0)
            # ratio=median(Z2[ii]/abs(Z1[ii]))

            Z <- (Z1 + Z2) %*% solve(tcrossprod(B) + L1k)
        } else {
            Z <- mat_mult(Y, t(B), ncores = ncores) %*%
                solve(tcrossprod(B) + L1k)
        }

        if (adaptive.p > 0) {
            Zraw <- Z
            Zraw[Zraw < 0] <- 0
            cutoffs <- apply(Z, 2, getT)

            for (j in seq_len(ncol(Z))) {
                Z[Z[, j] < cutoffs[j], j] <- 0
            }
        } else {
            Z[Z < 0] <- 0
            Zraw <- Z
        }

        oldB <- B

        if (is_fbm) {
            ZYt <- big_cprodMat(Y, as.matrix(Z), ncores = ncores)
            ZY <- Matrix::t(ZYt)
            B <- solve(Matrix::t(Z) %*% Z + L2k) %*% ZY
        } else {
            Z_mat <- if (inherits(Z, "matrix")) Z else as.matrix(Z)
            B <- solve(Matrix::t(Z_mat) %*% Z + L2k) %*%
                mat_mult(Matrix::t(Z), Y, ncores = ncores)
        }

        Bdiff <- sum((B - oldB)^2) / sum(B^2)
        minCor <- min(row_cor(B, oldB))

        BdiffTrace <- c(BdiffTrace, Bdiff)

        if (trace) {
            message(sprintf(
                "\rProgress %d / %d | Bdiff=%.6f",
                iter, max.iter, Bdiff
            ))
            flush.console()
        }

        if (iter > 52 && Bdiff > BdiffTrace[iter - 50]) {
            BdiffCount <- BdiffCount + 1
            # message("Bdiff is not decreasing")
        } else if (BdiffCount > 1) {
            BdiffCount <- BdiffCount - 1
        }

        if (Bdiff < tol & iter > u.iter + num.U.updates * 2 + 5) {
            message(sprintf(
                "\rConverged at %d / %d | Bdiff=%.6f, minCor=%.6f\n",
                iter, max.iter, Bdiff, minCor
            ))
            break
        }
        if (BdiffCount > 5) {
            message("converged at iteration", iter, "Bdiff is not decreasing")
            break
        }
    }

    rownames(U) <- colnames(priorMat)
    colnames(U) <- rownames(B) <- colnames(Z) <- paste0("LV", seq_len(clamp_k))
    if (!is.null(rownames(Y))) rownames(Z) <- rownames(Y)
    if (!is.null(colnames(Y))) colnames(B) <- colnames(Y)

    out <- list(
        B = B, Z = Z, U = U, C = C,
        L1 = L1, L2 = L2, heldOutGenes = heldOutGenes
    )

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
        message("There are ", sum(tt > 0.70), " LVs with AUC>0.70")
        message("There are ", sum(tt > 0.90), " LVs with AUC>0.90")
    } else {
        message("Not using cross-validation. No AUCs or p-values")
    }

    out$call <- call <- match.call()
    out$Z <- as.matrix(out$Z)
    out$B <- as.matrix(out$B)

    if (ncores > 1) {
        # restore previous state
        options(bigstatsr.check.parallel.blas = TRUE)
        options(default.nproc.blas = blas_nproc)
    }

    return(out)
}

#' Project new data into CLAMP latent space
#'
#' Computes the latent loadings \code{B} for new gene expression data using
#' the latent variables \code{Z} from a fitted CLAMP model. This allows
#' transfer of the learned latent structure to new datasets with matched genes.
#'
#' @param CLAMPres A result object from \code{CLAMPfull()} or
#'   \code{CLAMPbase()}, containing at least \code{Z} and \code{L2}.
#' @param newdata A gene expression matrix (genes x samples) to be projected.
#'   Can be a standard matrix, sparse matrix, or FBM/big.matrix.
#' @param scale Optional numeric multiplier for the L2 regularization terms.
#'   Default is 1.
#' @param ncores Number of cores to use for parallel computation (only used
#'   if newdata is an FBM).
#'  Default is 1.
#' @param align Logical; if \code{TRUE} (default), row names shared by
#'   \code{CLAMPres$Z} and \code{newdata} are used to align both matrices to
#'   the same genes in the same order before projection. If row names are not
#'   available, dimensions must already match.
#' @param verbose Logical; if \code{TRUE} (default), report how many common
#'   rows are used for projection.
#' @return A matrix \code{B} of projected latent loadings (LVs x samples)
#'   for the new dataset.
#'
#' @details
#' This function uses ridge-regularized least squares to compute
#' \code{B = solve(Z'Z + L2 * I) * Z'Y}, where \code{Z} is the latent matrix
#' from the trained CLAMP model and \code{Y} is the new dataset. If
#' \code{newdata} is a Filebacked Big Matrix (FBM) and does not need row-name
#' subsetting, the computation is optimized using
#' \code{bigstatsr::big_cprodMat()}.
#'
#' @examples
#' # fit a tiny CLAMP model for projection
#' Y0 <- matrix(rnorm(5 * 3),
#'     nrow = 5,
#'     dimnames = list(paste0("Gene", 1:5), paste0("S", 1:3))
#' )
#' base <- CLAMPbase(Y0, clamp_k = 2, max.iter = 1, trace = FALSE)
#' # new data can be provided in a different row order
#' newY <- matrix(rnorm(5 * 2),
#'     nrow = 5,
#'     dimnames = list(rev(rownames(Y0)), paste0("N", 1:2))
#' )
#' projB <- projectCLAMP(base, newdata = newY)
#' # check dimensions: 2 latent vars x 2 samples
#' dim(projB)
#'
#' @export
projectCLAMP <- function(CLAMPres, newdata, scale = 1, ncores = 1,
                         align = TRUE, verbose = TRUE) {
    if (is.null(CLAMPres$Z)) stop("'CLAMPres' must contain a 'Z' matrix.")
    if (is.null(CLAMPres$L2)) stop("'CLAMPres' must contain an 'L2' value.")

    Z_matrix <- if (inherits(CLAMPres$Z, "Matrix")) {
        as.matrix(CLAMPres$Z)
    } else {
        CLAMPres$Z
    }

    z_genes <- rownames(Z_matrix)
    new_genes <- rownames(newdata)

    if (align && !is.null(z_genes) && !is.null(new_genes)) {
        if (anyDuplicated(z_genes) > 0 || anyDuplicated(new_genes) > 0) {
            stop("Projection row-name alignment requires unique gene names.")
        }

        cm <- commonRows(Z_matrix, newdata)
        if (length(cm) == 0) {
            stop("No common row names found between CLAMPres$Z and newdata.")
        }
        if (verbose) message(length(cm), " common rows found")

        Z_matrix <- Z_matrix[cm, , drop = FALSE]
        newdata <- newdata[cm, , drop = FALSE]
    } else {
        if (nrow(Z_matrix) != nrow(newdata)) {
            stop(
                "CLAMPres$Z and newdata have different numbers of rows ",
                "and cannot be aligned because row names are missing."
            )
        }

        if (!is.null(z_genes) && !is.null(new_genes) &&
            !identical(z_genes, new_genes)) {
            stop(
                "CLAMPres$Z and newdata row names are not in the same order. ",
                "Use align = TRUE to align common genes automatically."
            )
        }

        if (verbose && (is.null(z_genes) || is.null(new_genes))) {
            message(
                "Row names unavailable; assuming CLAMPres$Z and newdata ",
                "are already aligned."
            )
        }
    }

    if (ncores > 1) {
        # if we are parallelizing, then disable BLAS parallelization
        options(bigstatsr.check.parallel.blas = FALSE)
        blas_nproc <- getOption("default.nproc.blas")
        options(default.nproc.blas = NULL)
        on.exit(
            {
                options(bigstatsr.check.parallel.blas = TRUE)
                options(default.nproc.blas = blas_nproc)
            },
            add = TRUE
        )
    }

    # Check if newdata is a FBM/big.matrix object
    is_fbm <- inherits(newdata, c("big.matrix", "FBM"))

    if (is_fbm) {
        # FBM implementation - use big_cprodMat for efficient computation
        # But ensure Z is a standard matrix, not a Matrix package object
        ZYt <- big_cprodMat(newdata, Z_matrix, ncores = ncores)
        ZY <- t(ZYt)
    } else {
        # Standard matrix implementation - use regular matrix multiplication
        ZY <- t(Z_matrix) %*% newdata
    }

    # Calculate the regularization matrix using standard matrix operations
    ZtZ <- t(Z_matrix) %*% Z_matrix
    L2k <- CLAMPres$L2 * diag(ncol(Z_matrix)) * scale

    # Solve the regularized system
    B <- solve(ZtZ + L2k) %*% ZY

    rownames(B) <- colnames(Z_matrix)
    colnames(B) <- colnames(newdata)

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
            uu0 <- rsvd(dat0, k = k, q = 3)
        }
        perm_mat[i, ] <- uu0$d[seq_len(k)]^2 / sum(uu0$d[seq_len(k)]^2)
    }
    p_vals <- apply(perm_mat >= obs_prop, 2, mean)
    p_vals <- cummax(p_vals)
    sum(p_vals <= 0.1)
}

#' Estimate number of principal components via elbow or permutation method
#'
#' @param data    Either a matrix (e.g. z-scored data) or an SVD result
#'   (list with $d).
#' @param method  One of "elbow" (fast) or "permutation" (slower).
#' @param B       Number of permutations (for method = "permutation").
#' @param seed    Seed for reproducibility.
#' @return        Estimated number of PCs.
#' @examples
#' # generate a small random matrix: 5 features x 10 samples
#' mat <- matrix(rnorm(5 * 10), nrow = 5)
#' # fast elbow estimate
#' num.pc(mat, method = "elbow")
#' # slower permutation estimate (use fewer perms for example speed)
#' num.pc(mat, method = "permutation", B = 5)
#' @export
num.pc <- function(data, method = c("elbow", "permutation"),
                   B = 20, seed = NULL) {
    method <- match.arg(method)
    if (!is.null(seed)) {
        warning(
            "`seed` is deprecated and ignored. ",
            "For reproducibility, call set.seed(seed) *before* this function.",
            call. = FALSE
        )
    }
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

#' Winsorize matrix columns by capping the top-k values
#'
#' For each column, replaces values greater than the k-th largest
#' entry with that threshold.
#'
#' @param M A numeric matrix.
#' @param k Integer; number of top elements to cap. Must be >= 1 and <= nrow(M).
#' @return A numeric matrix of the same dimensions as \code{M}, winsorized
#'   per column.
#' @export
#' @examples
#' set.seed(123)
#' M <- matrix(rnorm(20 * 5, mean = 0, sd = 2),
#'     nrow = 20,
#'     dimnames = list(paste0("Gene", 1:20), paste0("S", 1:5))
#' )
#'
#' # Display column maxima before winsorization
#' apply(M, 2, max)
#'
#' # Winsorize each column by capping top 3 values
#' M_winsor <- winsor_topk(M, k = 3)
winsor_topk <- function(M, k) {
    if (nrow(M) < 10 * k) {
        return(M)
    }
    stopifnot(is.matrix(M), is.numeric(M), k >= 1L, k <= nrow(M))
    n <- nrow(M)
    thr <- apply(M, 2L, function(x) {
        sort.int(x, partial = n - k + 1L)[n - k + 1L]
    })
    sweep(M, 2L, thr, pmin)
}

#' Cross-product Z^T Y with FBM or dense matrices
#'
#' Computes \eqn{Z^T Y}, handling FBM objects from bigstatsr
#' as well as base R matrices.
#'
#' @param Y Gene expression matrix (genes x samples), dense or FBM.
#' @param Z Latent variable matrix (genes x k).
#' @return A numeric matrix giving Z^T Y.
#' @importFrom Matrix t
#' @export
#' @examples
#' set.seed(123)
#'
#' genes <- 40
#' samples <- 10
#' k <- 4
#'
#' Y <- matrix(rnorm(genes * samples), nrow = genes)
#' Z <- matrix(rnorm(genes * k), nrow = genes)
#'
#' # Compute Z^T Y
#' res1 <- cross_ZY(Y, Z)
#' dim(res1) # k × samples
cross_ZY <- function(Y, Z) {
    if (inherits(Y, "FBM")) {
        Matrix::t(bigstatsr::big_cprodMat(Y, as.matrix(Z)))
    } else {
        base::crossprod(as.matrix(Z), as.matrix(Y))
    }
}

#' Ridge regression update for B
#'
#' Solves \eqn{B = (Z^T Z + L2)^{-1} Z^T Y} with ridge penalty matrix L2.
#'
#' @param Y Gene expression matrix (genes x samples), dense or FBM.
#' @param Z Latent variable matrix (genes x k).
#' @param L2k Ridge penalty matrix (k x k).
#' @return A numeric matrix of size k x samples.
#' @importFrom Matrix crossprod
#' @export
#' @examples
#' set.seed(123)
#'
#' genes <- paste0("Gene", 1:50)
#' samples <- paste0("S", 1:20)
#' k <- 5
#'
#' Y <- matrix(rnorm(50 * 20), nrow = 50, dimnames = list(genes, samples))
#' Z <- matrix(rnorm(50 * k),
#'     nrow = 50,
#'     dimnames = list(genes, paste0("LV", 1:k))
#' )
#'
#' lambda <- 0.1
#' L2k <- diag(lambda, k)
#'
#' # Solve for B = (Z'Z + L2)^(-1) Z'Y
#' B <- ridge_B(Y, Z, L2k)
ridge_B <- function(Y, Z, L2k) {
    Zm <- as.matrix(Z)
    ZtZ <- Matrix::crossprod(Zm) # Z^T Z
    ZY <- cross_ZY(Y, Zm) # Z^T Y
    out <- base::solve(ZtZ + L2k, ZY)
    as.matrix(out)
}


#' Runs the streamlined full CLAMP model.
#'
#' This version performs latent-variable decomposition of a gene expression
#' matrix \code{Y} guided by prior pathway annotations \code{priorMat}, with
#' simplified and lighter regularization compared to the original extended
#' CLAMP variant. The algorithm alternates updates of \code{Z}, \code{B},
#' and \code{U}, where \code{U} captures pathway-latent variable
#' associations inferred directly from the data without ridge-regularized
#' projections (\code{Chat} is not used).
#'
#' Cross-validation can be used to evaluate pathway-LV specificity, and a
#' variance-based prior (\code{var.prior = TRUE}) introduces adaptive
#' shrinkage of \code{Z} based on how strongly each latent component aligns
#' with prior pathways. The scaling factor \code{multiplier} (default 5)
#' controls the strength of this adaptive shrinkage.
#'
#' @param Y Gene expression matrix (genes x samples). Can be dense, sparse
#'   (dgCMatrix), or FBM.
#' @param priorMat Binary or weighted prior matrix (genes x pathways)
#'   linking genes to pathways.
#' @param Chat Ignored in this version (kept for interface compatibility).
#' @param svdres Optional precomputed SVD result for initialization.
#' @param clamp.base.result Optional result from \code{CLAMPbase()}
#'   providing initial values.
#' @param clamp_k Number of latent variables for CLAMP (final model rank).
#'   If \code{NULL}, it is chosen automatically via \code{select_clamp_k()}.
#' @param svd_k Number of singular values/components to compute in the SVD.
#'   If \code{NULL}, defaults to \code{max(2, min(n_genes, n_samples) - 1)}.
#' @param L1,L2 Regularization parameters for \code{Z} and \code{B}.
#'   Defaults use values from
#'   \code{clamp.base.result}.
#' @param cvn Number of folds for pathway-level cross-validation. Default: 5.
#' @param max.iter Maximum number of outer iterations. Default: 30.
#' @param trace Logical; print iteration progress. Default: \code{TRUE}.
#' @param maxPath Maximum number of pathways per LV during U-fitting.
#'   Default: 10.
#' @param doCrossval Whether to mask prior entries for CV evaluation.
#'   Default: \code{TRUE}.
#' @param penalty.factor Optional vector of per-pathway penalties for glmnet.
#'   Default: 1.
#' @param glm_alpha Elastic net mixing parameter for U estimation. Default: 0.9.
#' @param minGenes Minimum number of genes per pathway. Default: 0.
#' @param tol Convergence tolerance for B updates. Default: 5e-4.
#' @param seed Random seed for CV masking. Default: 123456.
#' @param allGenes If \code{TRUE}, adds zero-filled rows for missing genes.
#'   Default: \code{FALSE}.
#' @param rseed Reproducibility, coordinate descent updates are done in
#'   random order.
#' @param max.U.updates Maximum number of U updates (capped by max.iter).
#' @param pathwaySelection Pathway selection mode for U fitting
#'   (\code{"fast"} or \code{"complete"}).
#' @param multiplier Variance-prior scaling factor. Default: 5.
#' @param adaptive.p Quantile of negative Z values used to define adaptive
#'   thresholding. Default: 0.05.
#' @param useNNLS Whether to use non-negative least squares for U estimation.
#'   Default: \code{TRUE}.
#' @param useRaw If \code{TRUE}, uses unthresholded Z in U updates.
#'   Default: \code{TRUE}.
#' @param refitEvery Frequency (in U updates) of full refits. Default: 3.
#' @param var.prior Logical; if \code{TRUE}, enables adaptive variance prior
#'   updates for Z.
#'   Default: \code{TRUE}.
#' @param Uscale Logical; whether to scale U columns. Default: \code{FALSE}.
#' @param robust.vp Logical; winsorize prior-predicted Z2 values to reduce
#'   outlier effects. Default: \code{TRUE}.
#' @param useSE Logical; whether to use the 1-standard-error rule for
#'   internal glmnet fitting. Default is FALSE.
#' @param use_cpp Logical; if TRUE, use C++ implementation for Z updates.
#'   Default is FALSE.
#' @param clamp_k_method Method for selecting `clamp_k` when not provided.
#'   One of `"elbow"` (default), `"permutation"`, `"gavish_donoho"`, or
#'   `"scaleSVs"`. Passed to [select_clamp_k()].
#' @return A list with elements:
#' \describe{
#'   \item{\code{B}}{LV loadings on samples (k × samples)}
#'   \item{\code{Z}}{Gene loadings (genes × k)}
#'   \item{\code{U}}{Pathway loadings (pathways × k)}
#'   \item{\code{C}}{Prior matrix used during training (masked if CV)}
#'   \item{\code{Z2}}{Predicted Z from pathway priors}
#'   \item{\code{heldOutGenes}}{Held-out gene lists per pathway (CV mode)}
#'   \item{\code{Uauc, Up, summary}}{CV evaluation metrics if CV is enabled}
#'   \item{\code{priorMat, priorMatCV}}{Final and masked prior matrices}
#'   \item{\code{call}}{Function call}
#' }
#'
#' @details
#' This implementation omits ridge-projected priors (\code{Chat}) and uses
#' a lighter variance prior with a lower default \code{multiplier = 5},
#' allowing more flexible latent representations. Setting
#' \code{var.prior = FALSE} reproduces standard CLAMP-like updates.
#' Cross-validation, if
#' enabled, masks 20% of gene-pathway associations per column to estimate
#' pathway-LV specificity
#' (reported via AUC and p-values).
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' base <- CLAMPbase(mat, clamp_k = 5, trace = FALSE, max.iter = 5)
#' prior <- matrix(sample(0:1, 10 * 6, TRUE, prob = c(0.9, 0.1)),
#'     nrow = 10, ncol = 6
#' )
#' fit <- CLAMPfull(
#'     Y = mat,
#'     priorMat = prior,
#'     clamp.base.result = base,
#'     doCrossval = FALSE,
#'     adaptive.p = 0,
#'     max.U.updates = 0,
#'     max.iter = 1,
#'     trace = FALSE
#' )
#' @export
CLAMPfull <- function(
  Y, priorMat, Chat = NULL, svdres = NULL,
  clamp.base.result = NULL, clamp_k = NULL,
  svd_k = NULL, L1 = NULL, L2 = NULL, cvn = 5,
  max.iter = 30, trace = TRUE, maxPath = 10, doCrossval = TRUE,
  penalty.factor = rep(1, ncol(priorMat)), glm_alpha = 0.9,
  minGenes = 0, tol = 5e-4, seed = 123456,
  allGenes = FALSE, rseed = NULL,
  max.U.updates = Inf,
  pathwaySelection = c("fast", "complete"), multiplier = 5,
  adaptive.p = 0.05, useNNLS = TRUE, useRaw = TRUE, refitEvery = 3,
  useSE = FALSE, var.prior = TRUE, Uscale = FALSE,
  robust.vp = TRUE, use_cpp = FALSE,
  clamp_k_method = "elbow") {
    if (is.infinite(max.U.updates)) max.U.updates <- max.iter

    getT <- function(x) -stats::quantile(x[x < 0], adaptive.p)

    pathwaySelection <- match.arg(pathwaySelection, c("fast", "complete"))

    priorMat <- as.matrix(priorMat)

    message("** CLAMPfull **")

    ## Detect matrix type (FBM block preserved for compatibility)
    is_fbm <- inherits(Y, "FBM")
    is_sparse <- inherits(Y, "dgCMatrix")

    getVarMultiplier <- function(Zinput, Z2) {
        Zmultiplier <- matrix(0, nrow = nrow(Zinput), ncol = ncol(Zinput))
        for (zi in seq_len(ncol(Zinput))) {
            if (all(Z2[, zi] == 0)) next
            iiback <- Z2[, zi] == 0 & Zinput[, zi] > 0
            iiiforward <- Z2[, zi] > 0 & Zinput[, zi] > 0
            # define proxives for variance
            # background in Z
            mvarBack <- mean(Zinput[iiback, zi]^2)
            # foreground in Z
            mvarForward <- mean(Zinput[iiiforward, zi]^2)
            # foreground predicted in Z2
            mvarPredicted <- mean(Z2[iiiforward, zi]^2)
            # this offset will use the implied foreground variance scale
            # from actual Z
            # as opposed to the predcicted on in Z2
            offset <- mvarForward / mvarPredicted
            # compute the predicted variance
            Z2var <- Z2[, zi]^2 * offset
            # compute a fold change - 1 and clip at 0
            # 1 will be added later ensuring that the background variance
            # scale is always L1
            Zmultiplier[, zi] <- pmax((Z2var / mvarBack) - 1, 0)
        }
        Zmultiplier
    }

    ## Align genes
    if (nrow(priorMat) != nrow(Y) || !all(rownames(priorMat) == rownames(Y))) {
        if (!allGenes) {
            cm <- commonRows(Y, priorMat)
            message("Selecting common genes: ", length(cm))
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

    ## Filter small pathways
    numGenes <- Matrix::colSums(priorMat)
    heldOutGenes <- list()
    iibad <- which(numGenes < minGenes)
    if (length(iibad) > 0) {
        priorMat <- priorMat[, -iibad, drop = FALSE]
        message("Removed ", length(iibad), " pathways with too few genes")
    }
    ## Cross-validation masking
    if (doCrossval) {
        priorMatCV <- as.matrix(priorMat)

        for (j in seq_len(ncol(priorMatCV))) {
            iipos <- which(priorMatCV[, j] > 0)
            iiposs <- sample(iipos, length(iipos) / 5)
            priorMatCV[iiposs, j] <- 0
            heldOutGenes[[colnames(priorMat)[j]]] <- rownames(priorMat)[iiposs]
        }
        C <- priorMatCV
    } else {
        C <- priorMat
    }

    ## Initialize Chat
    if (is.null(Chat)) {
        # Chat <- pinv.ridge(crossprod(C), 5) %*% t(C)
    }

    ## SVD / base init
    ns <- ncol(Y)

    if (!is.null(svdres) && nrow(svdres$v) != ncol(Y)) {
        message("SVD V has the wrong number of columns")
        svdres <- NULL
    }

    if (is.null(svd_k) && is.null(clamp.base.result)) {
        svd_k <- select_svd_k(Y)
        if (!is.null(clamp_k)) svd_k <- max(svd_k, clamp_k)
    }

    if (is.null(svdres) && is.null(clamp.base.result)) {
        message("Computing SVD")
        svdres <- compute_svd(Y, k = svd_k)
    }

    if (is.null(svdres) && is.null(clamp.base.result)) {
        svdres <- rotateSVD(svdres)
    }

    if (is.null(clamp.base.result)) {
        clamp_k <- select_clamp_k(svdres,
            n_samples = ncol(Y), svd_k = svd_k,
            method = clamp_k_method, data = Y
        )
        d <- svdres$d[clamp_k]
    } else {
        d <- svdres$d[clamp_k]
    }

    if (is.null(clamp.base.result)) {
        message("Running CLAMPbase")
        clamp.base.result <- CLAMPbase(Y, clamp_k = clamp_k, svdres = svdres)
    } else {
        message("using provided CLAMPbase result")
        if (nrow(Y) != nrow(clamp.base.result$Z)) {
            if (is.null(rownames(Y)) | is.null(rownames(clamp.base.result$Z))) {
                stop(
                    "Y and clamp.base.result$Z must have equal row ",
                    "numbers or row names"
                )
            }
            clamp.base.result$Z <- clamp.base.result$Z[rownames(Y), ]
        }
        clamp_k <- ncol(clamp.base.result$Z)
    }

    message("CLAMP k is set to ", clamp_k)

    u.iter <- 2

    Z <- clamp.base.result$Z
    if (is.null(L1)) L1 <- clamp.base.result$L1
    if (is.null(L2)) L2 <- clamp.base.result$L2
    message("L1=", L1, "; L2=", L2)

    if (ncol(clamp.base.result$B) == ncol(Y)) {
        B <- clamp.base.result$B
        L1scalevec <- rowSums(B^2)
    } else {
        stop("clamp.base.result$B does not match number of samples in Y")
    }

    if (!is.null(rseed)) {
        message("Using random start")
        B <- t(apply(B, 1, sample))
        Z <- apply(Z, 2, sample)
    }

    U <- matrix(0, nrow = ncol(C), ncol = clamp_k)
    L1k <- L1 * diag(clamp_k)
    L2k <- L2 * diag(clamp_k)

    ## Optional FBM path (kept in comments for compatibility)
    # if (is_fbm) {
    #   ZYt <- bigstatsr::big_cprodMat(Y, as.matrix(Z))
    #   ZY  <- Matrix::t(ZYt)
    #   B   <- solve(Matrix::t(Z) %*% Z + L2k) %*% ZY
    # } else {
    #   B   <- solve(Matrix::t(Z) %*% Z + L2k) %*% mat_mult(Matrix::t(Z), Y)
    # }

    Z_mat <- if (inherits(Z, "matrix")) Z else as.matrix(Z)
    B <- ridge_B(Y, Z_mat, L2k)

    Zraw <- Z
    Z2 <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))

    B <- as.matrix(B)
    n <- nrow(Z)


    BdiffTrace <- numeric(0)
    BdiffCount <- 0
    oldB <- B

    num.U.updates <- 0
    start_time <- Sys.time()
    for (iter in seq_len(max.iter)) {
        if (iter == 4) {
            iter_time <- as.numeric(
                difftime(Sys.time(), start_time, units = "secs")
            )
            est_total <- iter_time * max.iter / 4
            message(sprintf(
                "Estimated total runtime: ~%.1f min", est_total / 60
            ))
        }
        if (iter >= u.iter) { # do U
            if (num.U.updates < max.U.updates) { # actually update U
                if (any(Zraw < 0)) {
                    stop("Zraw has negative entries before U update")
                }
                Uprev <- if (
                    num.U.updates %% refitEvery == 0 || iter == u.iter
                ) {
                    NULL
                } else {
                    U
                }
                if (!is.null(Uprev)) {
                    #  print("Reusing previous")
                } else {
                    #  print("Fitting new")
                }

                Zinput <- if (useRaw) Zraw else Z
                res <- solveU(Zinput, NULL, C, penalty.factor, pathwaySelection,
                    glm_alpha, maxPath,
                    binary = FALSE, nfolds = cvn,
                    useNNLS = useNNLS, Uprev = Uprev, useSE = useSE,
                    scale = Uscale, useAUC = TRUE
                )

                U <- res$U
                num.U.updates <- num.U.updates + 1
            } # end

            Z1 <- mat_mult(Y, t(B))

            if (!var.prior) {
                Z2 <- C %*% U
                Z <- (Z1 + L1 * Z2) %*% solve(Matrix::tcrossprod(B) + L1k)
            } else {
                Z2 <- as.matrix(C %*% U)
                if (robust.vp) Z2 <- winsor_topk(Z2, 20)
                Zmultiplier <- getVarMultiplier(Zinput, Z2)

                B2 <- B %*% t(B) # k×k
                bk2_all <- diag(B2) # length k
                YBt <- Y %*% t(B) # n×k
                Q <- Z %*% B2 # n×k cache (ZB projected via B)

                if (use_cpp) {
                    if (!is_fbm) {
                        updateZcpp(
                            Z, Y, B, B2, bk2_all, YBt, Q,
                            Zmultiplier, L1, multiplier, 3, iter
                        )
                    } else {
                        updateZcpp(
                            Z, Y$address, B, B2, bk2_all, YBt, Q,
                            Zmultiplier, L1, multiplier, 3, iter,
                            Y_is_fbm = TRUE, n = nrow(Y), p = ncol(Y)
                        )
                    }
                } else {
                    for (inner.iter in seq_len(3)) {
                        for (k_index in sample.int(clamp_k)) {
                            gene_var <- L1 *
                                (1 / (multiplier * Zmultiplier[, k_index] + 1))

                            bk2 <- bk2_all[k_index]
                            denom <- bk2 + gene_var
                            num <- YBt[, k_index] - Q[, k_index] +
                                Z[, k_index] * bk2

                            newZk <- num / denom
                            delta <- newZk - Z[, k_index]
                            Z[, k_index] <- newZk

                            # rank-1 n*1 * 1*k
                            Q <- Q + tcrossprod(delta, B2[k_index, ])

                            # if (anyNA(Z)) stop()
                        } # end for k
                    } # end inner iter
                } # end no CPP
            } # end var prior
        } # end if full update
        else {
            Z <- mat_mult(Y, t(B)) %*% solve(Matrix::tcrossprod(B) + L1k)
        }

        ## Adaptive thresholding
        if (adaptive.p > 0) {
            Zraw <- Z
            Zraw[Zraw < 0] <- 0
            cutoffs <- apply(Z, 2, getT)
            for (j in seq_len(ncol(Z))) Z[Z[, j] < cutoffs[j], j] <- 0
        } else {
            Z[Z < 0] <- 0
            Zraw <- Z
        }

        ## Update B
        oldB <- B

        Z_mat <- if (inherits(Z, "matrix")) Z else as.matrix(Z)
        B <- ridge_B(Y, Z_mat, L2k)

        ## Convergence diagnostics
        Bdiff <- sum((B - oldB)^2) / sum(B^2)
        BdiffTrace <- c(BdiffTrace, Bdiff)
        if (trace) {
            message(sprintf(
                "\rProgress %d / %d | Bdiff=%.6f",
                iter, max.iter, Bdiff
            ))
            utils::flush.console()
        }
        if (iter > 52 && Bdiff > BdiffTrace[iter - 50]) {
            BdiffCount <- BdiffCount + 1
            message("Bdiff is not decreasing")
        } else if (BdiffCount > 1) {
            BdiffCount <- BdiffCount - 1
        }
        if (Bdiff < tol && iter > u.iter + num.U.updates * 2 + 5) {
            message(sprintf(
                "\rConverged at %d / %d | Bdiff=%.6f",
                iter, max.iter, Bdiff
            ))
            break
        }
        if (BdiffCount > 5) {
            message("Converged early: Bdiff not decreasing")
            break
        }
    }

    rownames(U) <- colnames(priorMat)
    colnames(U) <- rownames(B) <- colnames(Z) <- paste0("LV", seq_len(clamp_k))
    rownames(Z) <- rownames(Y)
    if (!is.null(colnames(Y))) colnames(B) <- colnames(Y)

    out <- list(
        B = B, Z = Z, U = U, C = C,
        L1 = L1, L2 = L2, Z2 = Z2,
        heldOutGenes = heldOutGenes
    )

    if (doCrossval) {
        if (adaptive.p != 0) {
            message("Updating Z for CV")
            out$Z <- Zraw
            out$Z[out$Z < 0] <- 0
        }
        message("crossValidation")
        priorMat_m <- as.matrix(priorMat)
        priorMatCV_m <- if (exists("priorMatCV")) {
            as.matrix(priorMatCV)
        } else {
            as.matrix(priorMat)
        }
        outAUC <- crossVal(
            out, priorMat,
            if (exists("priorMatCV")) priorMatCV else priorMat
        )
        out$Z <- Z
        out$Uauc <- outAUC$Uauc
        out$Up <- outAUC$Upval
        out$summary <- outAUC$summary
        if (exists("priorMatCV")) out$priorMatCV <- priorMatCV
        out$priorMat <- priorMat
        out$withPrior <- which(colSums(out$U) > 0)
        tt <- apply(out$Uauc, 2, max)
        message("There are ", sum(tt > 0.70), " LVs with AUC>0.70")
        message("There are ", sum(tt > 0.90), " LVs with AUC>0.90")
    } else {
        message("Not using cross-validation. No AUCs or p-values")
    }

    out$call <- match.call()
    out$Z <- as.matrix(out$Z)
    out$B <- as.matrix(out$B)
    out
}

#' Subset and filter multiple pathway matrices to match target genes
#'
#' Filters one or more gene-by-pathway annotation matrices to retain only
#' pathways
#' with sufficient overlap with a given target gene set. Each input matrix is
#' restricted to \code{new.genes}, and pathways with fewer than \code{min.genes}
#' overlapping genes are removed. The resulting matrices are column-bound into a
#' single sparse matrix aligned to \code{new.genes}.
#'
#' @param ... One or more binary matrices (genes x pathways), either base
#'   \code{matrix}
#'   or sparse \code{Matrix} objects. Row names must be gene identifiers.
#' @param new.genes Character vector of gene names to align all pathway
#'   matrices to.
#' @param min.genes Integer; minimum number of overlapping genes required
#'   for a pathway
#'   to be retained. Default is 10.
#'
#' @return A sparse binary matrix with rows equal to \code{new.genes} and
#'   columns equal to
#'   the union of filtered pathways from all input matrices.
#' @export
#' @examples
#' set.seed(123)
#' library(Matrix)
#'
#' # Simulate two small pathway matrices (genes × pathways)
#' genes <- paste0("Gene", 1:100)
#' pathways1 <- paste0("Path", 1:5)
#' pathways2 <- paste0("Path", 6:10)
#'
#' mat1 <- matrix(sample(c(0, 1), 100 * 5, replace = TRUE, prob = c(0.9, 0.1)),
#'     nrow = 100, ncol = 5,
#'     dimnames = list(genes, pathways1)
#' )
#' mat2 <- matrix(sample(c(0, 1), 100 * 5, replace = TRUE, prob = c(0.9, 0.1)),
#'     nrow = 100, ncol = 5,
#'     dimnames = list(genes, pathways2)
#' )
#'
#' # Define target genes (subset of total)
#' new.genes <- sample(genes, 50)
#'
#' # Match and filter pathways with at least 5 genes
#' matched <- getMatchedPathwayMatList(mat1, mat2,
#'     new.genes = new.genes, min.genes = 5
#' )
getMatchedPathwayMatList <- function(..., new.genes, min.genes = 10) {
    pathMats <- list(...)

    stopifnot(is.character(new.genes), length(new.genes) > 0)

    filtered <- lapply(pathMats, function(pathMat) {
        # coerce to sparse
        if (!inherits(pathMat, "Matrix")) {
            pathMat <- Matrix::Matrix(pathMat, sparse = TRUE)
        }

        cm <- intersect(rownames(pathMat), new.genes)
        message(
            "There are ", length(cm),
            " genes in the intersection between data and prior"
        )

        if (length(cm) == 0L) {
            return(Matrix::Matrix(0,
                nrow = length(new.genes), ncol = 0,
                sparse = TRUE,
                dimnames = list(new.genes, character())
            ))
        }

        # create aligned matrix
        matchPathMat <- Matrix::Matrix(0,
            nrow = length(new.genes),
            ncol = ncol(pathMat), sparse = TRUE,
            dimnames = list(new.genes, colnames(pathMat))
        )

        matchPathMat[cm, ] <- pathMat[cm, , drop = FALSE]

        genesInPath <- Matrix::colSums(matchPathMat)
        ii <- which(genesInPath >= min.genes)

        message(sprintf(
            "Removing %d pathways",
            ncol(matchPathMat) - length(ii)
        ))

        matchPathMat[, ii, drop = FALSE]
    })

    if (length(filtered) == 1L) {
        return(filtered[[1L]])
    }
    do.call(cbind, filtered)
}
