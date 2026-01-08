#' Compute all-vs-all AUC matrix
#'
#' Calculates the area under the ROC curve (AUC) for all pairs of columns
#' between a prediction matrix `B` and binary targets in `target`.
#' Each column of `B` is ranked, and AUC is computed based on how well
#' the ranks separate positive vs. negative samples in each target column.
#'
#' @param B A numeric matrix of predictions (samples × features).
#' @param target A binary matrix of the same number of rows as `B`
#'   (samples × targets), where 1 indicates positive and 0 indicates negative.
#'
#' @return A numeric matrix of AUC values (features × targets).
#'
#' @examples
#' set.seed(1)
#' B <- matrix(rnorm(100), nrow = 20)
#' target <- matrix(sample(0:1, 40, replace = TRUE), nrow = 20)
#' allAgainstAllAUCs(B, target)
#'
#' @importFrom matrixStats colRanks
#' @export
allAgainstAllAUCs <- function(B, target) {
  B <- as.matrix(B)
  target <- as.matrix(target)
  if (!all(dim(B)[1] == dim(target)[1]))
    stop("B and target must have the same number of rows")

  ranks <- matrixStats::colRanks(B, ties.method = "average")

  n_pos <- colSums(target == 1, na.rm = TRUE)
  n_neg <- colSums(target == 0, na.rm = TRUE)

  pos_mean_rank <- ranks %*% (target == 1)
  pos_mean_rank <- sweep(pos_mean_rank, 2, n_pos * (n_pos + 1) / 2, "-")

  auc_matrix <- sweep(pos_mean_rank, 2, n_pos * n_neg, "/")
  return(auc_matrix)
}


#' Row-wise scaling (mean 0, sd 1)
#'
#' Standardize each row of a numeric matrix to have mean 0 and
#' standard deviation 1.
#'
#' @param x A numeric matrix.
#' @return A matrix of the same shape, with each row scaled independently.
#' @export
tscale <- function(x) {
  row_means <- rowMeans(x)
  row_sds <- sqrt(rowMeans((x - row_means)^2))
  row_sds[row_sds == 0] <- 1  # avoid division by zero
  sweep(sweep(x, 1, row_means), 1, row_sds, "/")
}



#' Print a concatenated message
#'
#' Wrapper around `message()` that pastes arguments together into a single string.
#' @param ... Character strings to concatenate and print.
#' @return Invisibly returns NULL. Called for side effects (messages).
mymessage <- function(...) {
  message(...)
}

#' Get maximum AUC per latent variable
#'
#' Summarizes a cross-validation results data frame to extract the highest AUC value
#' associated with each latent variable (LV).
#'
#' @param summary A data frame (e.g., from `crossVal()`) .
#' @param verbose Logical; if `TRUE`, prints counts of LVs exceeding AUC thresholds. Default is `FALSE`.
#'
#' @return A data frame with columns LV index and max_AUC.
getMaxAUC <- function(summary, verbose = FALSE) {
  max_auc_per_lv <- summary %>%
    group_by(.data$`LV index`) %>%
    summarize(max_AUC = max(.data$AUC, na.rm = TRUE)) %>%
    ungroup()

  if (verbose) {
    message("There are ", sum(max_auc_per_lv$max_AUC > 0.7), " LVs with AUC>0.70")
    message("There are ", sum(max_auc_per_lv$max_AUC > 0.9), " LVs with AUC>0.90")
  }
  max_auc_per_lv
}

#' Count number of latent variables exceeding AUC thresholds
#'
#' Given a summary data frame from cross-validation, reports the number of latent variables
#' with maximum AUC exceeding 0.7, 0.8, and 0.9.
#'
#' @param summary A data frame with `LV index` and `AUC` columns.
#'
#' @return A named numeric vector with counts for thresholds 0.7, 0.8, and 0.9.
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
#' @param cor_mat A square numeric matrix of pairwise correlations (rows = items, cols = items).
#' @return A vector of assignments (integer indices).
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

  list(permutation = assignment, sum = sum(cor_mat[cbind(seq_len(n), assignment)]))
}

#' Download and read a GMT file from a URL
#'
#' Downloads a Gene Matrix Transposed (GMT) file from a specified URL, reads it into R as a list,
#' and removes the temporary file afterward.
#'
#' @param url A character string specifying the URL to a GMT file.
#' @param name Optional name for the GMT file (defaults to the portion after the last '=' in the URL).
#' @param cache_dir Optional directory in which to cache/download the GMT file.
#' @param redownload Logical; if `TRUE`, forces re-download even if cached.
#'
#' @return A named list where each element is a character vector of gene names for a given gene set.
#'
#' @examples
#' url <- 'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human'
#' gmt_list <- getGMT(url)
#' # list available gene sets
#' names(gmt_list)
#' # inspect the first few genes in the first gene set
#' head(gmt_list[[1]])
#' @export
getGMT <- function(url, name = NULL, cache_dir = NULL, redownload = FALSE) {
  if (is.null(name)) {
    name <- sub(".*[=]", "", url)
    message("Auto-detected name: ", name)
  }
  if (is.null(cache_dir)) {
    cache_dir <- system.file("extdata", package = "PLIER2")
  }
  cache_file <- file.path(cache_dir, paste0(name, ".gmt"))
  if (!file.exists(cache_file) || redownload) {
    message("Downloading ", name, " from Enrichr...")
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
    download.file(url, cache_file)
  } else {
    message("Using cached file for ", name)
  }
  read_gmt(cache_file)
}

#' Read a GMT file into a list
#'
#' Parses a local GMT file and returns a list of gene sets, with each gene
#' set represented as
#' a character vector of unique gene names.
#'
#' @param filename A character string giving the path to a .gmt file.
#'
#' @return A named list where each element is a character vector of gene names.
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
#' Converts a list of named gene sets (e.g., from `getGMT()`) into a sparse binary matrix
#' where rows are genes, columns are gene sets, and entries are 1 if the gene is in the set.
#'
#' @param gmtList A nested list of gene sets. Outer names are gene set names; each entry is a character
#' vector of gene names.
#'
#' @return A sparse binary matrix with genes as rows and gene sets as columns.
#' @examples
#' # define a simple nested GMT list
#' gmt1 <- list(
#'     PathwayA = c('Gene1', 'Gene2', 'Gene3'),
#'     PathwayB = c('Gene2', 'Gene4')
#' )
#' gmt2 <- list(
#'     PathwayC = c('Gene1', 'Gene4'),
#'     PathwayD = c('Gene3', 'Gene5')
#' )
#' # combine into a nested list
#' nestedList <- list(gmt1 = gmt1, gmt2 = gmt2)
#' # convert to sparse matrix
#' sparseMat <- gmtListToSparseMat(nestedList)
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

  pathMat <- Matrix::sparseMatrix(i = row_indices, j = col_indices, x = values, dims = c(length(allGenes),
                                                                                         length(allnames)))
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
commonRows <- function(data1, data2) {
  intersect(rownames(data1), rownames(data2))
}

#' Clean a Filebacked Big Matrix (FBM) by log-transforming and handling NAs
#'
#' This function inspects an FBM to determine if log-transformation is needed (based on value range)
#' and whether NA values are present. If the maximum value is >= 100, it applies a log2(x + 1)
#' transformation in-place. If any NA values are detected, they are replaced with 0.
#'
#' @param fbm A `bigmemory::FBM` or `bigstatsr::FBM` object.
#' @param ncores Integer; number of cores to use for parallel operations (default 1).
#' @return A list with:
#'   \describe{
#'     \item{max_value}{The maximum value encountered in the FBM (after log transformation if applied).}
#'     \item{had_na}{Logical indicating whether any NA values were found and filled.}
#'   }
#'
#' @details
#' Modifies the FBM in place. Uses `bigstatsr::big_apply()` to process in parallel-safe chunks.
#' @importFrom bigstatsr big_apply rows_along FBM
cleanFBM <- function(fbm, ncores = 1) {
  # Block‐wise scan for max and NA
  stats <- big_apply(fbm, a.FUN = function(X, ind) {
    vals <- X[, ind, drop = FALSE]
    list(max = if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE), na = anyNA(vals))
  }, a.combine = function(...) {
    Reduce(function(a, b) {
      list(max = max(a$max, b$max, na.rm = TRUE), na = a$na || b$na)
    }, list(...))
  }, ind = bigstatsr::cols_along(fbm), ncores = ncores, )

  max_value <- stats$max
  has_na <- stats$na

  # Log2 transform if necessary
  if (!is.na(max_value) && max_value >= 100) {
    message("Applying log2 transformation")
    big_apply(fbm, a.FUN = function(X, ind) {
      X[, ind] <- log2(X[, ind] + 1)
      NULL
    }, ind = bigstatsr::cols_along(fbm), ncores = ncores, )
  } else {
    message("Already on log scale or all NA")
  }

  # Fill NAs if present
  if (has_na) {
    message("Filling NAs with 0")
    big_apply(fbm, a.FUN = function(X, ind) {
      X[, ind][is.na(X[, ind])] <- 0
      NULL
    }, ind = bigstatsr::cols_along(fbm), ncores = ncores, )
  } else {
    message("No NA values found")
  }

  return(list(max_value = max_value, had_na = has_na))
}

#' Compute row-wise sum and sum of squares for a Filebacked Big Matrix
#'
#' Efficiently computes row sums and row sum-of-squares for a `bigstatsr::FBM` using
#' column-wise chunking, suitable for large datasets that cannot be loaded fully into memory.
#'
#' @param fbm A `bigstatsr::FBM` object.
#' @param ncores Integer; number of cores to use for parallel operations (default 1).
#' @return A list with two numeric vectors:
#' \describe{
#'   \item{row_sums}{Sum of each row.}
#'   \item{row_sums_sq}{Sum of squares of each row.}
#' }
#'
computeRowStatsFBM <- function(fbm, ncores = 1) {
  # Compute row sums in blocks
  row_sums <- big_apply(fbm, a.FUN = function(X, ind) rowSums(X[, ind]), a.combine = "plus",
                        ncores = ncores)

  # Compute row sums of squares in blocks
  row_sums_sq <- big_apply(fbm, a.FUN = function(X, ind) rowSums(X[, ind]^2), a.combine = "plus",
                           ncores = ncores)

  n_cols <- ncol(fbm)
  # Final means and variances
  row_means <- row_sums/n_cols
  row_variances <- (row_sums_sq/n_cols) - (row_means^2)

  list(row_means = row_means, row_variances = row_variances)
}

#' Filter rows of a Filebacked Big Matrix based on mean and variance
#'
#' Filters an FBM based on row-level mean and variance thresholds, returning a new FBM
#' with only the selected rows.
#'
#' @param fbm A `bigstatsr::FBM` object.
#' @param rowStats A list with numeric vectors `row_means` and `row_variances`.
#' @param mean_cutoff Optional minimum mean threshold; rows with means below this are removed.
#' @param var_cutoff Optional minimum variance threshold; rows with variances below this are removed.
#' @param backingfile A character string specifying the filename (without extension) for the new FBM.
#' Default is `filtered_fbm`.
#'
#' @return A list with:
#' \describe{
#'   \item{fbm_filtered}{A new FBM object containing only filtered rows.}
#'   \item{kept_rows}{Indices of rows retained in the filtering step.}
#' }
#' @details
#' This function creates a new FBM and copies over only the rows that pass the filtering criteria.
#' The original FBM is unchanged.
filterFBM <- function(fbm, rowStats, mean_cutoff = NULL, var_cutoff = NULL, backingfile = "filtered_fbm") {
  row_means <- rowStats$row_means
  row_variances <- rowStats$row_variances

  # Determine rows to keep based on cutoffs
  keep_rows <- rep(TRUE, length(row_means))  # Default: keep all rows

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
#' @param Y_filtered Numeric matrix (genes x samples) returned by preprocessPLIER2
#' @param rowStats   Data frame with numeric columns `mean` and `variance`,
#'                   row-named to match `rownames(Y_filtered)`
#'
#' @return Numeric matrix of the same dimensions as Y_filtered, with each row centered and scaled.
#' @examples
#' # simple 2 genes x 3 samples matrix
#' Y <- matrix(
#'     c(
#'         2, 4, 6, # gene1 counts
#'         8, 10, 12 # gene2 counts
#'     ),
#'     nrow = 2, byrow = TRUE,
#'     dimnames = list(c('gene1', 'gene2'), paste0('sample', seq_len(3)))
#' )
#'
#' # compute per‐gene mean and variance
#' rowStats <- data.frame(
#'     mean = rowMeans(Y),
#'     variance = apply(Y, 1, var),
#'     row.names = rownames(Y)
#' )
#'
#' # z‐score each row
#' Y_z <- zscorePLIER2(Y, rowStats)
#' @export
zscorePLIER2 <- function(Y_filtered, rowStats) {
  # Input validation
  if (!is.matrix(Y_filtered) || !is.numeric(Y_filtered)) {
    stop("`Y_filtered` must be a numeric matrix (genes x samples).")
  }
  if (!is.data.frame(rowStats) || !all(c("mean", "variance") %in% colnames(rowStats))) {
    stop("`rowStats` must be a data.frame with columns 'mean' and 'variance'.")
  }
  # Align rowStats to Y_filtered
  if (!all(rownames(Y_filtered) %in% rownames(rowStats))) {
    stop("Row names of `Y_filtered` and `rowStats` do not match.")
  }
  rowStats <- rowStats[rownames(Y_filtered), , drop = FALSE]
  # Ensure numeric
  mu <- as.numeric(rowStats$mean)
  var <- as.numeric(rowStats$variance)
  if (any(is.na(mu)) || any(is.na(var))) {
    stop("Missing values detected in 'mean' or 'variance'.")
  }
  if (any(var <= 0)) {
    stop("All variances must be positive; zero or negative found.")
  }
  # Compute standard deviation
  sd <- sqrt(var)
  # Center and scale subtract mu from each row, then divide by sd
  Y_centered <- sweep(Y_filtered, 1L, mu, "-")
  Y_scaled <- sweep(Y_centered, 1L, sd, "/")
  # Return
  return(Y_scaled)
}

#' Preprocess a bigstatsr FBM for PLIER2
#'
#' Makes a writable copy of the input FBM, cleans it (log2 transform if needed, fill NAs),
#' filters rows by mean/variance, and returns the filtered FBM plus stats and indices.
#'
#' @param fbm A bigstatsr::FBM (genes x samples), possibly read-only.
#' @param mean_cutoff Numeric or NULL. Minimum row mean to keep (NULL = no mean filter).
#' @param var_cutoff  Numeric or NULL. Minimum row variance to keep (NULL = no var filter).
#' @param backingfile Character or NULL. Base name for the *copy* FBM and filtered FBM on disk.
#'                    If NULL, defaults to paste0(fbm$backingfile, '_preproc') and '_filtered'.
#' @param ncores Integer; number of cores to use for parallel operations (default 1).
#' @param block_size Number of rows to process at a time when copying data. Default is 1000.
#' @return A list with:
#'   \item{fbm_filtered}{The filtered FBM (writable).}
#'   \item{rowStats}{List with row_means & row_variances for fbm_filtered.}
#'   \item{kept_rows}{Integer vector of original row indices that were retained.}
#' @examples
#' library(bigstatsr)
#' # create a toy matrix and back it with an FBM
#' mat <- matrix(
#'     c(
#'         1, 2, 3, # geneA
#'         10, 20, 30, # geneB
#'         100, 200, 300 # geneC
#'     ),
#'     nrow = 3, byrow = TRUE,
#'     dimnames = list(c('geneA', 'geneB', 'geneC'), paste0('s', seq_len(3)))
#' )
#' fbm <- FBM(nrow(mat), ncol(mat), init = mat)
#'
#' # preprocess without filtering (all genes kept)
#' res_all <- preprocessPLIER2FBM(fbm)
#' @export
preprocessPLIER2FBM <- function(fbm, mean_cutoff = NULL, var_cutoff = NULL, backingfile = NULL,
                                block_size = 1000, ncores = 1) {
  n_r <- nrow(fbm)
  n_c <- ncol(fbm)

  # Choose base names
  base_bk <- if (is.null(backingfile)) {
    paste0(fbm$backingfile, "_preproc")
  } else {
    backingfile
  }

  # Make a writable copy
  fbm_copy <- FBM(nrow = n_r, ncol = n_c, backingfile = base_bk, create_bk = TRUE)

  # copy all data
  for (rs in seq(1, n_r, by = block_size)) {
    rows <- rs:min(rs + block_size - 1L, n_r)
    fbm_copy[rows, ] <- fbm[rows, ]
  }

  if (ncores > 1) {
    # if we are parallelizing, then disable BLAS parallelization
    options(bigstatsr.check.parallel.blas = FALSE)
    blas_nproc <- getOption("default.nproc.blas")
    options(default.nproc.blas = NULL)
  }

  # Clean in-place (log2 if needed, fill NAs)
  cleanFBM(fbm_copy, ncores)

  # Compute row stats on cleaned copy
  rs_all <- computeRowStatsFBM(fbm_copy, ncores)

  if (ncores > 1) {
    options(bigstatsr.check.parallel.blas = TRUE)
    options(default.nproc.blas = blas_nproc)
  }

  # Filter rows, writing to a new filtered FBM
  filt_bk <- paste0(base_bk, "_filtered")
  filter_res <- filterFBM(fbm_copy, rowStats = rs_all, mean_cutoff = mean_cutoff,
                          var_cutoff = var_cutoff, backingfile = filt_bk)

  fbm_filtered <- filter_res$fbm_filtered
  kept_rows <- filter_res$kept_rows

  # Subset stats to kept rows
  stats_filt <- list(row_means = rs_all$row_means[kept_rows], row_variances = rs_all$row_variances[kept_rows])

  list(fbm_filtered = fbm_filtered, rowStats = stats_filt, kept_rows = kept_rows)
}

#' Z-score a filtered FBM in-place
#'
#' Standardizes each row of an FBM using provided row means and variances.
#'
#' @param fbm_filtered A bigstatsr::FBM produced by preprocessPLIER2FBM().
#' @param rowStats A list with row_means and row_variances from that FBM.
#' @param chunk_size Columns per block (default 1000).
#' @param ncores Integer; number of cores to use for parallel operations (default 1).
#' @examples
#' library(bigstatsr)
#' fbm <- FBM(
#'     nrow = 2, ncol = 4,
#'     init = matrix(seq_len(8), nrow = 2)
#' )
#' stats <- list(
#'     row_means     = rowMeans(fbm[]),
#'     row_variances = apply(fbm[], 1, var)
#' )
#' zscorePLIER2FBM(fbm, stats, chunk_size = 2)
#' @return A normalized FBM with z-scored rows.
#' @export
zscorePLIER2FBM <- function(fbm_filtered, rowStats, chunk_size = 1000, ncores = 1) {
  message("Applying Z-score transformation")
  means <- rowStats$row_means
  sds <- sqrt(rowStats$row_variances)
  sds[!is.finite(sds) | sds == 0] <- 1

  if (ncores > 1) {
    options(bigstatsr.check.parallel.blas = FALSE)
    old_blas <- getOption("default.nproc.blas")
    options(default.nproc.blas = NULL)
    on.exit({
      options(bigstatsr.check.parallel.blas = TRUE)
      options(default.nproc.blas = old_blas)
    }, add = TRUE)
  }

  bigstatsr::big_apply(fbm_filtered, a.FUN = function(X, ind, means, sds) {
    block <- X[, ind, drop = FALSE]
    block <- sweep(block, 1, means, "-")
    block <- sweep(block, 1, sds, "/")
    X[, ind] <- block
    integer(0)
  }, a.combine = "c", ind = bigstatsr::cols_along(fbm_filtered), block.size = chunk_size,
  ncores = ncores, means = means, sds = sds)

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
#'   - Y_filtered: filtered matrix (genes x samples)
#'   - rowStats: data.frame with columns mean and variance for each kept gene
#'   - kept_rows: integer vector of the original row indices that were kept
#'
#' @examples
#' # construct a small example matrix
#' mat <- matrix(
#'     c(
#'         1, 5, 10,
#'         2, 6, 11,
#'         3, 7, 12,
#'         4, 8, 13
#'     ),
#'     nrow = 4, byrow = FALSE,
#'     dimnames = list(paste0('gene', seq_len(4)), paste0('sample', seq_len(3)))
#' )
#'
#' # keep genes with mean >= 6 and variance >= 2
#' res <- preprocessPLIER2(mat, mean_cutoff = 6, var_cutoff = 2)
#' @export
preprocessPLIER2 <- function(Y, mean_cutoff = 0, var_cutoff = 0) {
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("`Y` must be a numeric matrix (genes x samples).")
  }
  # Compute per‐gene statistics
  row_mean <- rowMeans(Y, na.rm = TRUE)
  row_var <- apply(Y, 1, stats::var, na.rm = TRUE)

  rowStats <- data.frame(mean = row_mean, variance = row_var, stringsAsFactors = FALSE)
  rownames(rowStats) <- rownames(Y)

  # Identify genes passing both thresholds
  keep <- which(rowStats$mean >= mean_cutoff & rowStats$variance >= var_cutoff)
  if (length(keep) == 0) {
    stop("No genes passed the mean/variance filters.")
  }

  # Subset matrix and stats
  Y_filtered <- Y[keep, , drop = FALSE]
  rowStats_filtered <- rowStats[keep, , drop = FALSE]

  return(list(Y_filtered = Y_filtered, rowStats = rowStats_filtered, kept_rows = keep))
}

#' Compute counts-per-million (CPM) for PLIER2 pipelines
#'
#' @param counts A numeric matrix or data.frame of raw counts (genes x samples).
#' @return A numeric matrix of CPM values (same dimensions), ready for PLIER2 input.
#' @examples
#' mat <- matrix(seq_len(12), nrow = 3)
#' cpmPLIER2(mat)
#' @export
cpmPLIER2 <- function(counts) {
  mat <- if (is.data.frame(counts)) {
    as.matrix(counts)
  } else {
    counts
  }
  stopifnot(is.numeric(mat), length(dim(mat)) == 2)
  lib_sizes <- colSums(mat, na.rm = TRUE)
  if (any(lib_sizes == 0)) {
    warning("Some samples have zero total counts - CPM will be Inf/NaN.")
  }
  sweep(mat, 2, lib_sizes, "/") * 1e+06
}


#' Compute CPM on a file-backed matrix for PLIER2 (in-place)
#'
#' @param fbm_counts A bigstatsr::FBM of raw counts (genes x samples).
#' @param block_size Integer; columns per block (default 1000).
#' @param ncores Integer; number of cores to use for parallel operations (default 1).
#' @return Invisibly returns the modified FBM (now holding CPM values).
#' @examples
#' library(bigstatsr)
#' mat <- matrix(c(10, 20, 30, 40, 50, 60),
#'     nrow = 2,
#'     dimnames = list(c('gene1', 'gene2'), paste0('sample', seq_len(3)))
#' )
#' fbm <- FBM(nrow(mat), ncol(mat), init = mat)
#' cpmPLIER2FBM(fbm, block_size = 1)
#' @export
cpmPLIER2FBM <- function(fbm_counts, block_size = 1000, ncores = 1) {
  if (!inherits(fbm_counts, "FBM")) {
    stop("`fbm_counts` must be a bigstatsr::FBM object.")
  }
  block_size <- as.integer(block_size)
  if (block_size <= 0) {
    stop("`block_size` must be a positive integer.")
  }

  # Avoid BLAS oversubscription when parallelizing
  if (ncores > 1) {
    options(bigstatsr.check.parallel.blas = FALSE)
    old_blas <- getOption("default.nproc.blas")
    options(default.nproc.blas = NULL)
    on.exit({
      options(bigstatsr.check.parallel.blas = TRUE)
      options(default.nproc.blas = old_blas)
    }, add = TRUE)
  }

  # Library sizes (sum per column), processed in column chunks
  lib_sizes <- bigstatsr::big_apply(fbm_counts, a.FUN = function(X, ind) {
    colSums(X[, ind, drop = FALSE])
  }, a.combine = "c", ind = bigstatsr::cols_along(fbm_counts), block.size = block_size,
  ncores = ncores)
  lib_sizes[lib_sizes == 0] <- 1

  # Divide each column by its library size and scale to CPM, in-place
  bigstatsr::big_apply(fbm_counts, a.FUN = function(X, ind, libs) {
    blk <- X[, ind, drop = FALSE]
    blk <- sweep(blk, 2, libs[ind], "/") * 1e+06
    X[, ind] <- blk
    integer(0)
  }, a.combine = "c", ind = bigstatsr::cols_along(fbm_counts), block.size = block_size,
  ncores = ncores, libs = lib_sizes)

  invisible(fbm_counts)
}
#' Find the location of the maximum of a smoothing spline
#'
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values (same length as \code{x}).
#' @param n Integer, number of grid points to evaluate (default 1000).
#' @param spar Smoothing parameter passed to \code{stats::smooth.spline}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{x}{The x coordinate at which the spline reaches its maximum.}
#'   \item{y}{The corresponding maximum fitted y value.}
#' }
#' @examples
#' x <- 1:10
#' y <- sin(x) + rnorm(10, 0, 0.1)
#' findSplineMax(x, y)
#'
#' @importFrom stats smooth.spline predict
#' @export
findSplineMax <- function(x, y, n = 1000, spar = NULL) {
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
  fit <- stats::smooth.spline(x, y, spar = spar)
  grid <- seq(min(x), max(x), length.out = n)
  pred <- stats::predict(fit, grid)
  max_idx <- which.max(pred$y)
  list(x = pred$x[max_idx], y = pred$y[max_idx])
}

#' Squash extreme z-scores
#'
#' @param zdata Numeric vector or matrix of z-scores.
#' @param maxScore Numeric scalar, maximum absolute score (default 2).
#'
#' @return A numeric object of same dimensions as \code{zdata}, with values
#' shrunk by a hyperbolic tangent transformation.
#'
#' @examples
#' z <- rnorm(10, 0, 5)
#' squashZscore(z)
#'
#' @export
squashZscore <- function(zdata, maxScore = 2) {
  stopifnot(is.numeric(zdata), is.numeric(maxScore), length(maxScore) == 1)
  maxScore * tanh(zdata / maxScore)
}


#' Estimate noise scale from singular values with linear tail extrapolation
#'
#' This function estimates a characteristic scale from a vector of singular values
#' by fitting a linear model to the tail and extrapolating to length \code{n}.
#' The median of the extrapolated values is returned as the estimate.
#' If no sufficiently linear tail is detected (based on \code{min_r2}) or
#' the extrapolation produces negative values, a fallback estimate is returned
#' using the 75\% quantile singular value.
#'
#' @param sv Numeric vector of singular values sorted in decreasing order.
#' @param n Integer, total length to which the linear tail is extrapolated.
#' @param min_r2 Numeric, minimum R-squared value required for accepting
#'   the linear tail fit. Defaults to \code{0.95}.
#'
#' @return A numeric scalar giving the estimated scale. If linear extrapolation
#'   is unreliable, returns \code{sv[ceiling(0.75 * length(sv))]}.
#'
#' @examples
#' sv <- exp(-seq(0, 5, length.out = 50)) + rnorm(50, 0, 0.01)
#' getScaleFromSVs(sv, n = 100)
#'
#' @export
getScaleFromSVs <- function(sv, n, min_r2 = 0.95) {
  k <- length(sv)
  #works best on log scale
sv=log(sv)
  if(k<20){

    stop("Need at least 20 singular values ")
  }
  if (k ==n){
    #drop the last few
    sv=sv[1:(k-5)]
    k <- length(sv)
  }
  fallback <- sv[ceiling(0.75 * k)]

  best_r2 <- -Inf
  best_drop <- 0

  for (frac_drop in seq(0, 0.8, by = 0.01)) {
    drop <- ceiling(k * frac_drop)
    if (k - drop < 4) next
    x <- (drop + 1):k
    fit <- lm(sv[x] ~ x)
    r2 <- summary(fit)$r.squared
    if (r2 > best_r2) {
      best_r2 <- r2
      best_drop <- drop

    }
  }

  if (best_r2 < min_r2){
    message(paste("Minimum r2 not reached, best_r2 is", best_r2))
    scale=exp(median(y_pred))
    k=NULL
    return(list(scale=scale))
  }


  x <- (best_drop + 1):k
  fit <- lm(sv[x] ~ x)
  y_pred <- predict(fit, newdata = data.frame(x = 1:n))

  if (mean(y_pred < 0)>0.1){

    message("y_pred is <0, consider supplying more SVs")
    scale=exp(median(y_pred))
    k=best_drop
    return(list(scale=scale,k=k))
  }
  scale=exp(median(y_pred))
  k=best_drop
  return(list(scale=scale, k=k))
}


#' Scale matrix rows by inverse SVD row residuals
#'
#' Computes a low-rank reconstruction from an `rsvd` result, calculates the
#' per-row root-mean-squared residuals, and rescales each row of the input
#' matrix by the inverse of its residual. This is equivalent to applying
#' \eqn{\sqrt{w_i} = 1 / r_i} row-weights for weighted least squares.
#'
#' @param X Numeric matrix. Input data matrix.
#' @param svdres List returned by `rsvd::rsvd()`, containing components
#'   `u`, `d`, and `v`.
#'
#' @return A matrix of the same dimensions as `X` with rowwise scaling applied.
#'
#' @examples
#' # svdres <- rsvd::rsvd(X, k = 5)
#' # X_scaled <- scale_rows(X, svdres)
#'
#' @export
scale_rows <- function(X, svdres, return.scale=F) {
  Xhat <- (svdres$u %*% diag(svdres$d)) %*% t(svdres$v)
  r <- sqrt(rowSums((X - Xhat)^2) / ncol(X))  # RMSE per row
  if(!return.scale){
    X * (1 / r)
  }
  else{
    1/r
  }
}


