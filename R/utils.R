library(dplyr)


#' Row-wise scaling (mean 0, sd 1)
#'
#' Standardize each row of a numeric matrix to have mean 0 and standard deviation 1.
#'
#' @param x A numeric matrix.
#' @return A matrix of the same shape, with each row scaled independently.
#' @export
tscale <- function(x) {
    row_means <- rowMeans(x)
    row_sds <- sqrt(rowMeans((x - row_means)^2))
    row_sds[row_sds == 0] <- 1 # avoid division by zero
    sweep(sweep(x, 1, row_means), 1, row_sds, "/")
}


#' Compare latent variable loadings against a target using correlation or other statistics
#'
#' This function compares two sets of latent variable loadings (`res1`, `res2`) with respect to a binary or continuous
#' target matrix. It computes the maximal association (e.g., correlation, t-statistic, or AUC) between each latent variable
#' and each column in the target, and returns a paired comparison plot and summary.
#' @export
compareBs <- function(res1, res2, target, method = "p", xlab = "1", ylab = "2", stat.method = "t") {
    extract_B <- function(res) {
        if (class(res)[1] == "list") {
            return(as.matrix(res$B))
        } else if (class(res)[1] == "rsvd") {
            return(t(res$v))
        } else {
            return(as.matrix(res))
        }
    }

    B1 <- extract_B(res1)
    B2 <- extract_B(res2)
    noNA <- !apply(target, 1, function(x) {
        any(is.na(x))
    })
    B1 <- B1[, noNA]
    B2 <- B2[, noNA]
    target <- target[noNA, ]

    if (method %in% c("s", "p")) {
        mat1 <- cor(t(B1), target, method = method)
        mat2 <- cor(t(B2), target, method = method)
        cor1 <- apply(mat1, 2, max)
        cor2 <- apply(mat2, 2, max)
        idx1 <- apply(mat1, 2, which.max)
        idx2 <- apply(mat2, 2, which.max)
    } else if (method == "a") {
        mat1 <- allAgainstAllAUCs(t(B1), target)
        mat2 <- allAgainstAllAUCs(t(B2), target)
        cor1 <- apply(mat1, 2, max)
        cor2 <- apply(mat2, 2, max)
        idx1 <- apply(mat1, 2, which.max)
        idx2 <- apply(mat2, 2, which.max)
    } else if (method == "t") {
        mat1 <- allAgainstAllTstats(t(B1), target)
        mat2 <- allAgainstAllTstats(t(B2), target)
        cor1 <- apply(mat1, 1, max)
        cor2 <- apply(mat2, 1, max)
        idx1 <- apply(mat1, 1, which.max)
        idx2 <- apply(mat2, 1, which.max)
    }
    Labels <- colnames(target)

    df <- data.frame(Cor1 = cor1, Cor2 = cor2, BestIdxRes1 = idx1, BestIdxRes2 = idx2, Label = Labels)
    df$corMean <- (df$Cor1 + df$Cor2) / 2
    if (stat.method == "t") {
        pval <- t.test(df$Cor2, df$Cor1, paired = TRUE, alternative = "greater")$p.value
    } else {
        pval <- wilcox.test(df$Cor2, df$Cor1, paired = TRUE, alternative = "greater")$p.value
    }
    method_label <- switch(method,
        p = "Pearson correlation",
        s = "Spearman correlation",
        a = "AUC",
        t = "T-statistic"
    )
    pl <- ggplot(df, aes(x = Cor1, y = Cor2, label = Label)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
        geom_text_repel() +
        labs(
            x = paste("Max", method_label, xlab),
            y = paste("Max", method_label, ylab)
        ) +
        annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = paste("p =", signif(pval, 3))) +
        #    annotate("text", x = Inf, y = -Inf, label = paste("p =", signif(pval, 3)),
        #            hjust = 1.1, vjust = -0.5, size = 4)+
        theme_minimal()
    return(list(plot = pl, df = df))
}


#' Print a concatenated message
#'
#' Wrapper around \code{message()} that pastes arguments together into a single string.
#'
#' @param ... Character strings to concatenate and print.
#'
#' @export
mymessage <- function(...) {
    message(paste(...))
}

#' Get maximum AUC per latent variable
#'
#' Summarizes a cross-validation results data frame to extract the highest AUC value
#' associated with each latent variable (LV).
#'
#' @param summary A data frame (e.g., from \code{crossVal()}) .
#' @param verbose Logical; if \code{TRUE}, prints counts of LVs exceeding AUC thresholds. Default is \code{FALSE}.
#'
#' @return A data frame with columns \code{LV index} and \code{max_AUC}.
#'
#' @export
getMaxAUC <- function(summary, verbose = F) {
    max_auc_per_lv <- summary %>%
        group_by(.data$`LV index`) %>%
        summarize(max_AUC = max(.data$AUC, na.rm = TRUE)) %>%
        ungroup()

    if (verbose) {
        message(paste("There are", sum(max_auc_per_lv$max_AUC > 0.70), " LVs with AUC>0.70"))
        message(paste("There are", sum(max_auc_per_lv$max_AUC > 0.90), " LVs with AUC>0.90"))
    }
    max_auc_per_lv
}


#' Count number of latent variables exceeding AUC thresholds
#'
#' Given a summary data frame from cross-validation, reports the number of latent variables
#' with maximum AUC exceeding 0.7, 0.8, and 0.9.
#'
#' @param summary A data frame with \code{LV index} and \code{AUC} columns.
#'
#' @return A named numeric vector with counts for thresholds 0.7, 0.8, and 0.9.
#'
#' @export
getAUCstats <- function(summary) {
    out <- getMaxAUC(summary)
    unlist(lapply(c(0.7, 0.8, 0.9), function(x) {
        sum(out$max_AUC > x)
    }))
}


#' Greedy maximum correspondence from correlation matrix
#'
#' Finds a one-to-one assignment (permutation) between rows and columns of a square correlation matrix
#' that maximizes the total correlation score, using a greedy algorithm.
#'
max_correspondence_greedy <- function(cor_mat) {
    n <- nrow(cor_mat)
    used_rows <- rep(FALSE, n)
    used_cols <- rep(FALSE, n)
    assignment <- integer(n)

    corr_entries <- as.data.frame(which(!is.na(cor_mat), arr.ind = TRUE))
    corr_entries$val <- cor_mat[cbind(corr_entries$row, corr_entries$col)]
    corr_entries <- corr_entries[order(-corr_entries$val), ]

    for (i in seq_len(nrow(corr_entries))) {
        r <- corr_entries$row[i]
        c <- corr_entries$col[i]
        if (!used_rows[r] && !used_cols[c]) {
            assignment[r] <- c
            used_rows[r] <- TRUE
            used_cols[c] <- TRUE
        }
    }

    list(permutation = assignment, sum = sum(cor_mat[cbind(1:n, assignment)]))
}

#' Download and read a GMT file from a URL
#'
#' Downloads a Gene Matrix Transposed (GMT) file from a specified URL, reads it into R as a list,
#' and removes the temporary file afterward.
#'
#' @param url A character string specifying the URL to a GMT file.
#'
#' @return A named list where each element is a character vector of gene names for a given gene set.
#' @export
getGMT <- function(url) {
    tmp_file <- tempfile(fileext = ".gmt")

    download.file(url, tmp_file)

    result <- read_gmt(tmp_file)

    unlink(tmp_file)

    return(result)
}

#' Read a GMT file into a list
#'
#' Parses a local GMT file and returns a list of gene sets, with each gene set represented as
#' a character vector of unique gene names.
#'
#' @param filename A character string giving the path to a .gmt file.
#'
#' @return A named list where each element is a character vector of gene names.
#' @export
read_gmt <- function(filename) {
    gmt <- list()
    lines <- readLines(filename)
    for (line in lines) {
        line <- gsub("\"", "", trimws(line))
        sp <- unlist(strsplit(line, "\t"))
        sp[3:length(sp)] <- gsub(",.*$", "", sp[3:length(sp)])
        gmt[[sp[1]]] <- sort(unique(sp[3:length(sp)]))
    }
    return(gmt)
}

#' Convert a list of GMT gene sets to a sparse matrix
#'
#' Converts a list of named gene sets (e.g., from \code{getGMT()}) into a sparse binary matrix
#' where rows are genes, columns are gene sets, and entries are 1 if the gene is in the set.
#'
#' @param gmtList A nested list of gene sets. Outer names are gene set names; each entry is a character vector of gene names.
#'
#' @return A sparse binary matrix with genes as rows and gene sets as columns.
#' @export
gmtListToSparseMat <- function(gmtList) {
    allnames <- unlist(lapply(gmtList, names))
    # there are usually no duplicates
    stopifnot(all(table(allnames) == 1))
    allGenes <- unique(unlist(lapply(gmtList, unlist)))

    row_indices <- integer(0)
    col_indices <- integer(0)
    values <- integer(0)
    for (gmt in seq_along(gmtList)) {
        for (path in names(gmtList[[gmt]])) {
            thisPathGenes <- gmtList[[gmt]][[path]]
            iiGenes <- match(thisPathGenes, allGenes)
            iPath <- match(path, allnames)

            # Store indices and values
            row_indices <- c(row_indices, iiGenes)
            col_indices <- c(col_indices, rep(iPath, length(iiGenes)))
            values <- c(values, rep(1, length(iiGenes)))
        }
    }

    # Use sparseMatrix to create the matrix in one go

    pathMat <- sparseMatrix(
        i = row_indices, j = col_indices, x = values,
        dims = c(length(allGenes), length(allnames))
    )
    rownames(pathMat) <- allGenes
    colnames(pathMat) <- allnames
    pathMat
}

#' Find common row names between two matrices or data frames
#'
#' Returns the intersection of row names shared by two input objects.
#'
#' @param data1 A matrix, data frame, or similar object with row names.
#' @param data2 A matrix, data frame, or similar object with row names.
#'
#' @return A character vector of row names common to both inputs.
#'
#' @export
commonRows <- function(data1, data2) {
    intersect(rownames(data1), rownames(data2))
}


#' Clean a Filebacked Big Matrix (FBM) by log-transforming and handling NAs
#'
#' This function inspects an FBM to determine if log-transformation is needed (based on value range)
#' and whether NA values are present. If the maximum value is ≥ 100, it applies a log2(x + 1)
#' transformation in-place. If any NA values are detected, they are replaced with 0.
#'
#' @param fbm A \code{bigmemory::FBM} or \code{bigstatsr::FBM} object.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{max_value}}{The maximum value encountered in the FBM (after log transformation if applied).}
#'   \item{\code{had_na}}{Logical indicating whether any NA values were found and filled.}
#' }
#'
#' @details
#' Modifies the FBM in place. Uses \code{bigstatsr::big_apply()} to process in parallel-safe chunks.
#'
#' @export
cleanFBM <- function(fbm) {
    # Step 1: Check for NA and max value
    max_value <- -Inf
    has_na <- FALSE

    big_apply(fbm, a.FUN = function(X, ind) {
        max_value <<- max(max_value, max(X[ind, ], na.rm = TRUE))
        if (anyNA(X[ind, ])) {
            has_na <<- TRUE
        }
        NULL # No return, just updating global values
    }, ind = rows_along(fbm))

    # Step 2: Log2 transform if necessary
    if (max_value >= 100) {
        message("Applying log2 transformation")
        big_apply(fbm, a.FUN = function(X, ind) {
            X[ind, ] <- log2(X[ind, ] + 1)
        }, ind = rows_along(fbm))
    } else {
        message("Already on log scale")
    }

    # Step 3: Fill NAs with 0 if necessary
    if (has_na) {
        message("Filling NAs with 0")
        big_apply(fbm, a.FUN = function(X, ind) {
            X[ind, ][is.na(X[ind, ])] <- 0
            NULL
        }, ind = rows_along(fbm))
    } else {
        message("No NA values found")
    }

    return(list(max_value = max_value, had_na = has_na))
}


#' Compute row-wise sum and sum of squares for a Filebacked Big Matrix
#'
#' Efficiently computes row sums and row sum-of-squares for a \code{bigstatsr::FBM} using
#' column-wise chunking, suitable for large datasets that cannot be loaded fully into memory.
#'
#' @param fbm A \code{bigstatsr::FBM} object.
#' @param chunk_size Number of columns to process at a time. Default is 1000.
#'
#' @return A list with two numeric vectors:
#' \describe{
#'   \item{\code{row_sums}}{Sum of each row.}
#'   \item{\code{row_sums_sq}}{Sum of squares of each row.}
#' }
#'
#' @export
computeRowStatsFBM <- function(fbm, chunk_size = 1000) {
    n_rows <- nrow(fbm)
    n_cols <- ncol(fbm)

    # Initialize vectors to store row sums and row sums of squares
    row_sums <- numeric(n_rows)
    row_sums_sq <- numeric(n_rows)

    # Iterate over columns in chunks
    for (start_col in seq(1, n_cols, by = chunk_size)) {
        end_col <- min(start_col + chunk_size - 1, n_cols)

        # Extract the current chunk of columns
        col_chunk <- fbm[, start_col:end_col]

        # Update row sums and sums of squares
        row_sums <- row_sums + rowSums(col_chunk)
        row_sums_sq <- row_sums_sq + rowSums(col_chunk^2)
    }

    # Compute row means and variances
    row_means <- row_sums / n_cols
    row_variances <- (row_sums_sq / n_cols) - (row_means^2)

    return(list(row_means = row_means, row_variances = row_variances))
}

#' Filter rows of a Filebacked Big Matrix based on mean and variance
#'
#' Filters an FBM based on row-level mean and variance thresholds, returning a new FBM
#' with only the selected rows.
#'
#' @param fbm A \code{bigstatsr::FBM} object.
#' @param rowStats A list with numeric vectors \code{row_means} and \code{row_variances}.
#' @param mean_cutoff Optional minimum mean threshold; rows with means below this are removed.
#' @param var_cutoff Optional minimum variance threshold; rows with variances below this are removed.
#' @param backingfile A character string specifying the filename (without extension) for the new FBM. Default is \code{"filtered_fbm"}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{fbm_filtered}}{A new FBM object containing only filtered rows.}
#'   \item{\code{kept_rows}}{Indices of rows retained in the filtering step.}
#' }
#' @details
#' This function creates a new FBM and copies over only the rows that pass the filtering criteria.
#' The original FBM is unchanged.
#' @export
filterFBM <- function(fbm, rowStats, mean_cutoff = NULL, var_cutoff = NULL, backingfile = "filtered_fbm") {
    row_means <- rowStats$row_means
    row_variances <- rowStats$row_variances

    # Determine rows to keep based on cutoffs
    keep_rows <- rep(TRUE, length(row_means)) # Default: keep all rows

    if (!is.null(mean_cutoff)) {
        keep_rows <- keep_rows & (row_means >= mean_cutoff)
    }

    if (!is.null(var_cutoff)) {
        keep_rows <- keep_rows & (row_variances >= var_cutoff)
    }

    # Number of rows to keep
    n_kept <- sum(keep_rows)

    if (n_kept == 0) {
        stop("No rows meet the filtering criteria.")
    }

    # Create a new FBM with the filtered data
    fbm_filtered <- FBM(n_kept, ncol(fbm), backingfile = backingfile)
    fbm_filtered[] <- fbm[keep_rows, ]

    return(list(fbm_filtered = fbm_filtered, kept_rows = which(keep_rows)))
}

#' Z-score a filtered expression matrix for PLIER2
#'
#' Centers each gene to mean 0 and scales to unit variance.
#'
#' @param Y_filtered Numeric matrix (genes × samples) returned by preprocessPLIER2
#' @param rowStats   Data frame with numeric columns `mean` and `variance`,
#'                   row-named to match `rownames(Y_filtered)`
#'
#' @return Numeric matrix of the same dimensions as `Y_filtered`, with each
#'         row centered and scaled.
#' @export
zscorePLIER2 <- function(Y_filtered, rowStats) {
    # 1) Input validation
    if (!is.matrix(Y_filtered) || !is.numeric(Y_filtered)) {
        stop("`Y_filtered` must be a numeric matrix (genes × samples).")
    }
    if (!is.data.frame(rowStats) ||
        !all(c("mean", "variance") %in% colnames(rowStats))) {
        stop("`rowStats` must be a data.frame with columns 'mean' and 'variance'.")
    }
    # 2) Align rowStats to Y_filtered
    if (!all(rownames(Y_filtered) %in% rownames(rowStats))) {
        stop("Row names of `Y_filtered` and `rowStats` do not match.")
    }
    rowStats <- rowStats[rownames(Y_filtered), , drop = FALSE]
    # 3) Ensure numeric
    mu <- as.numeric(rowStats$mean)
    var <- as.numeric(rowStats$variance)
    if (any(is.na(mu)) || any(is.na(var))) {
        stop("Missing values detected in 'mean' or 'variance'.")
    }
    if (any(var <= 0)) {
        stop("All variances must be positive; zero or negative found.")
    }
    # 4) Compute standard deviation
    sd <- sqrt(var)
    # 5) Center and scale
    #    subtract mu from each row, then divide by sd
    Y_centered <- sweep(Y_filtered, 1L, mu, "-")
    Y_scaled <- sweep(Y_centered, 1L, sd, "/")
    # 6) Return
    return(Y_scaled)
}


#' Preprocess a bigstatsr FBM for PLIER2
#'
#' Makes a writable copy of the input FBM, cleans it (log2 transform if needed, fill NAs),
#' filters rows by mean/variance, and returns the filtered FBM plus stats and indices.
#'
#' @param fbm A bigstatsr::FBM (genes × samples), possibly read-only.
#' @param mean_cutoff Numeric or NULL. Minimum row mean to keep (NULL = no mean filter).
#' @param var_cutoff  Numeric or NULL. Minimum row variance to keep (NULL = no var filter).
#' @param backingfile Character or NULL. Base name for the *copy* FBM and filtered FBM on disk.
#'                    If NULL, defaults to paste0(fbm$backingfile, "_preproc") and "_filtered".
#' @return A list with:
#'   \item{fbm_filtered}{The filtered FBM (writable).}
#'   \item{rowStats}{List with row_means & row_variances for fbm_filtered.}
#'   \item{kept_rows}{Integer vector of original row indices that were retained.}
#' @export
preprocessPLIER2FBM <- function(
        fbm,
        mean_cutoff = NULL,
        var_cutoff = NULL,
        backingfile = NULL) {
    # 1. Choose base names
    base_bk <- if (is.null(backingfile)) paste0(fbm$backingfile, "_preproc") else backingfile

    # 2. Make a writable copy
    fbm_copy <- FBM(
        nrow        = nrow(fbm),
        ncol        = ncol(fbm),
        backingfile = base_bk,
        create_bk   = TRUE
    )
    # copy all data
    fbm_copy[] <- fbm[]

    # 3. Clean in-place (log2 if needed, fill NAs)
    cleanFBM(fbm_copy)

    # 4. Compute row stats on cleaned copy
    rs_all <- computeRowStatsFBM(fbm_copy)

    # 5. Filter rows, writing to a new filtered FBM
    filt_bk <- paste0(base_bk, "_filtered")
    filter_res <- filterFBM(
        fbm_copy,
        rowStats    = rs_all,
        mean_cutoff = mean_cutoff,
        var_cutoff  = var_cutoff,
        backingfile = filt_bk
    )
    fbm_filtered <- filter_res$fbm_filtered
    kept_rows <- filter_res$kept_rows

    # 6. Subset stats to kept rows
    stats_filt <- list(
        row_means     = rs_all$row_means[kept_rows],
        row_variances = rs_all$row_variances[kept_rows]
    )

    list(
        fbm_filtered = fbm_filtered,
        rowStats     = stats_filt,
        kept_rows    = kept_rows
    )
}

#' Z-score a filtered FBM in-place
#'
#' Standardizes each row of an FBM using provided row means and variances.
#'
#' @param fbm_filtered A bigstatsr::FBM produced by preprocessPLIER2FBM().
#' @param rowStats A list with row_means and row_variances from that FBM.
#' @param chunk_size Columns per block (default 1000).
#' @export
zscorePLIER2FBM <- function(fbm_filtered, rowStats, chunk_size = 1000) {
    message("Applying Z-score transformation")
    means <- rowStats$row_means
    sds <- sqrt(rowStats$row_variances)
    sds[sds == 0] <- 1

    for (start in seq(1, ncol(fbm_filtered), by = chunk_size)) {
        end <- min(start + chunk_size - 1, ncol(fbm_filtered))
        mat <- fbm_filtered[, start:end]
        mat <- sweep(mat, 1, means, FUN = "-")
        mat <- sweep(mat, 1, sds, FUN = "/")
        fbm_filtered[, start:end] <- mat
    }

    invisible(NULL)
}

#' Preprocess an expression matrix for PLIER2
#'
#' Filters genes by mean expression and variance, returning the filtered matrix
#' and per-gene statistics.
#'
#' @param Y Numeric matrix of gene expression (rows = genes, cols = samples)
#' @param mean_cutoff Numeric. Minimum row-mean required to keep a gene (default 0).
#' @param var_cutoff  Numeric. Minimum row-variance required to keep a gene (default 0).
#'
#' @return A list with components:
#'   - **fbm_filtered**: filtered matrix (genes × samples)
#'   - **rowStats**: data.frame with columns `mean` and `variance` for each kept gene
#'   - **kept_rows**: integer vector of the original row indices that were kept
#'
#' @export
preprocessPLIER2 <- function(Y, mean_cutoff = 0, var_cutoff = 0) {
    if (!is.matrix(Y) || !is.numeric(Y)) {
        stop("`Y` must be a numeric matrix (genes × samples).")
    }
    # Compute per‐gene statistics
    row_mean <- rowMeans(Y, na.rm = TRUE)
    row_var <- apply(Y, 1, stats::var, na.rm = TRUE)

    rowStats <- data.frame(
        mean = row_mean,
        variance = row_var,
        stringsAsFactors = FALSE
    )
    rownames(rowStats) <- rownames(Y)

    # Identify genes passing both thresholds
    keep <- which(rowStats$mean >= mean_cutoff & rowStats$variance >= var_cutoff)
    if (length(keep) == 0) {
        stop("No genes passed the mean/variance filters.")
    }

    # Subset matrix and stats
    Y_filtered <- Y[keep, , drop = FALSE]
    rowStats_filtered <- rowStats[keep, , drop = FALSE]

    return(list(
        Y_filtered = Y_filtered,
        rowStats   = rowStats_filtered,
        kept_rows  = keep
    ))
}
