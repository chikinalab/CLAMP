#' Compare two sets of factor loadings or embeddings
#'
#' Compares correspondence between two result matrices (e.g., factor loadings)
#' across multiple targets using correlation, AUC, or t-statistics.
#' Produces a paired comparison plot and summary statistics.
#'
#' @param res1, res2 Result objects or matrices containing factor loadings.
#'   If a list, must contain element `B`; if of class `"rsvd"`,
#'   `t(res$v)` is used.
#' @param target Numeric matrix of target variables (samples × targets).
#' @param method Character, one of `"p"`, `"s"`, `"a"`, or `"t"`, indicating
#'   Pearson, Spearman, AUC, or t-statistic comparison.
#' @param xlab, ylab Labels for x and y axes in the plot.
#' @param stat.method Statistical test to compare correlations
#'   (`"t"` or `"wilcox"`).
#' @param oneToOne Logical, whether to apply one-to-one masking of associations.
#' @param res2 Second result object or matrix to compare against.
#' @param ylab Label for y-axis (used in plot).
#' @return A list with:
#' \describe{
#'   \item{plot}{A ggplot object comparing maximal correlations across targets.}
#'   \item{df}{A data.frame with per-target metrics and (if available)
#'     top genes.}
#' }
#'
#' @importFrom matrixStats colRanks
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs annotate
#' @importFrom ggplot2 theme_minimal
#' @export
#' @examples
#' set.seed(123)
#' # Simulated example with 50 genes × 20 samples
#' Y <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)
#'
#' # Run two CLAMP-like decompositions (here using simple SVD)
#' svd1 <- rsvd::rsvd(Y, k = 5)
#' svd2 <- rsvd::rsvd(Y + matrix(rnorm(50 * 20, 0, 0.1), 50, 20), k = 5)
#'
#' # Define a target variable (e.g., binary or continuous trait)
#' target <- matrix(rnorm(20 * 3), ncol = 3)
#' colnames(target) <- c("Trait1", "Trait2", "Trait3")
#'
#' # Compare the two sets of embeddings
#' res <- compareBs(svd1, svd2, target,
#'     method = "p", xlab = "SVD1", ylab = "SVD2"
#' )
compareBs <- function(res1, res2, target, method = "p", xlab = "1", ylab = "2",
                      stat.method = "t", oneToOne = TRUE) {
    extract_B <- function(res) {
        if (is.list(res) & !is.null(res$B)) {
            return(as.matrix(res$B))
        } else if (inherits(res, "rsvd")) {
            return(t(res$v))
        } else {
            return(as.matrix(res))
        }
    }

    B1 <- extract_B(res1)
    B2 <- extract_B(res2)
    show(dim(B1))
    show(dim(B2))
    target <- as.matrix(target)
    noNA <- !apply(target, 1, function(x) any(is.na(x)))
    B1 <- B1[, noNA, drop = FALSE]
    B2 <- B2[, noNA, drop = FALSE]
    target <- target[noNA, , drop = FALSE]

    # Only correlation-based comparison supported
    if (method %in% c("s", "p")) {
        mat1 <- cor(t(B1), target, method = method)
        mat2 <- cor(t(B2), target, method = method)
    }

    # if (method %in% c("s", "p")) {
    #   mat1 <- cor(t(B1), target, method = method)
    #   mat2 <- cor(t(B2), target, method = method)
    # } else if (method == "a") {
    #   mat1 <- allAgainstAllAUCs(t(B1), target)
    #   mat2 <- allAgainstAllAUCs(t(B2), target)
    # } else if (method == "t") {
    #   mat1 <- allAgainstAllTstats(t(B1), target)
    #   mat2 <- allAgainstAllTstats(t(B2), target)
    # }

    if (inherits(res1, "rsvd")) {
        mat1 <- abs(mat1)
    }
    if (inherits(res2, "rsvd")) {
        mat1 <- abs(mat1)
    }
    mat1[is.na(mat1)] <- 0
    mat2[is.na(mat2)] <- 0
    if (oneToOne) {
        mat1 <- oneToOneMask(mat1)
        mat2 <- oneToOneMask(mat2)
    }

    cor1 <- apply(mat1, 2, max)
    cor2 <- apply(mat2, 2, max)
    cor1[cor1 < 0.01] <- NA
    cor2[cor2 < 0.01] <- NA
    idx1 <- apply(mat1, 2, which.max)
    idx2 <- apply(mat2, 2, which.max)
    Labels <- colnames(target)

    df <- data.frame(
        Cor1 = cor1, Cor2 = cor2,
        BestIdxRes1 = idx1, BestIdxRes2 = idx2,
        Label = Labels
    )
    df$corMean <- (df$Cor1 + df$Cor2) / 2

    pval <- if (stat.method == "t") {
        t.test(df$Cor2, df$Cor1, paired = TRUE, alternative = "greater")$p.value
    } else {
        wilcox.test(
            df$Cor2, df$Cor1,
            paired = TRUE, alternative = "greater"
        )$p.value
    }

    method_label <- switch(method,
        p = "Pearson correlation",
        s = "Spearman correlation",
        a = "AUC",
        t = "T-statistic"
    )

    pl <- ggplot2::ggplot(df, ggplot2::aes(x = Cor1, y = Cor2, label = Label)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(
            slope = 1, intercept = 0,
            linetype = "dashed", color = "red", linewidth = 1
        ) +
        ggrepel::geom_text_repel() +
        ggplot2::labs(
            x = sprintf("Max %s %s", method_label, xlab),
            y = sprintf("Max %s %s", method_label, ylab)
        ) +
        ggplot2::annotate(
            "text",
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
            label = sprintf("p = %.3g", pval)
        ) +
        ggplot2::theme_minimal()

    if (is.list(res1) && !is.null(res1$Z)) {
        get_top_genes <- function(Z, k) {
            apply(Z, 2, function(col) {
                top_n <- seq_len(min(k, length(col)))
                names(sort(col, decreasing = TRUE))[top_n]
            })
        }
        topk1 <- get_top_genes(res1$Z, 7)
        topk2 <- get_top_genes(res2$Z, 7)
        df$TopGenesRes1 <- vapply(
            idx1, function(i) paste(topk1[, i], collapse = ","), ""
        )
        df$TopGenesRes2 <- vapply(
            idx2, function(i) paste(topk2[, i], collapse = ","), ""
        )
    }

    return(list(plot = pl, df = df))
}

#' ComplexHeatmap visualization of top genes by latent variable
#'
#' Produces a combined heatmap showing expression, pathway membership,
#' and optionally Z-loadings for top genes per LV using ComplexHeatmap.
#'
#' @param clampRes A CLAMP result list containing matrices \code{Z},
#'   \code{U}, etc.
#' @param data Expression matrix with genes as rows.
#' @param priorMat Binary gene × pathway matrix.
#' @param top Integer, number of top genes per LV.
#' @param top.pathway Integer, number of top pathways per LV to annotate.
#' @param index Optional vector of LVs to include.
#' @param allLVs Logical; include all LVs.
#' @param Zheat Logical; whether to include the Z matrix as an additional
#'   heatmap.
#' @param LV.names Optional character vector for LV names.
#' @param max.genes Maximum number of genes allowed in the plot.
#' @param max.col Maximum number of columns (samples).
#' @param seed Random seed for column subsampling.
#'
#' @return Invisibly returns the drawn ComplexHeatmap object.
#' @import ComplexHeatmap
#' @import circlize
#' @importFrom stats setNames
#' @importFrom Matrix rowSums
#' @export
#' @examples
#' library(ComplexHeatmap)
#' library(Matrix)
#'
#' # Simulate small CLAMP-like results
#' set.seed(123)
#' genes <- paste0("Gene", 1:100)
#' samples <- paste0("S", 1:20)
#' lvs <- paste0("LV", 1:3)
#'
#' # Simulated Z (gene loadings) and U (pathway loadings)
#' Z <- matrix(rnorm(100 * 3), nrow = 100, dimnames = list(genes, lvs))
#' U <- matrix(abs(rnorm(50 * 3)),
#'     nrow = 50,
#'     dimnames = list(paste0("Path", 1:50), lvs)
#' )
#'
#' # Expression data
#' expr_data <- matrix(rnorm(100 * 20),
#'     nrow = 100,
#'     dimnames = list(genes, samples)
#' )
#'
#' # Binary gene × pathway matrix
#' priorMat <- matrix(sample(0:1, 100 * 50, replace = TRUE, prob = c(0.9, 0.1)),
#'     nrow = 100, ncol = 50,
#'     dimnames = list(genes, paste0("Path", 1:50))
#' )
#'
#' # Create a CLAMP-like result list
#' clampRes <- list(Z = Z, U = U)
#'
#' # Plot top genes and pathway memberships
#' plotTopZ_Complex(clampRes, expr_data, priorMat,
#'     top = 5, top.pathway = 3, index = 1:2, Zheat = TRUE
#' )
plotTopZ_Complex <- function(clampRes, data, priorMat,
                             top = 10, top.pathway = 5,
                             index = NULL, allLVs = FALSE, Zheat = FALSE,
                             LV.names = NULL, max.genes = 100, max.col = 50,
                             seed = 1234) {
    data <- data[rownames(clampRes$Z), , drop = FALSE]
    if (top * length(index) > max.genes) {
        stop(
            "Too many genes. Reduce number of LVs or 'top', ",
            "or increase 'max.genes'."
        )
    }

    if (ncol(data) > max.col) {
        data <- data[, sample(ncol(data), max.col)]
    }

    priorMat <- priorMat[rownames(clampRes$Z), , drop = FALSE]
    ii <- which(colSums(clampRes$U) > 0)
    if (!allLVs) {
        if (!is.null(index)) ii <- intersect(ii, index)
    } else {
        ii <- index
    }

    tmp <- apply(-clampRes$Z[, ii, drop = FALSE], 2, rank)
    nn <- unique(unlist(apply(tmp, 2, function(x) names(which(x <= top)))))
    nn <- sort(unique(nn))
    data_sub <- t(scale(t(data[nn, , drop = FALSE])))

    nnpath <- vapply(
        seq_along(ii),
        function(i) {
            gene_idx <- match(nn, rownames(priorMat))
            col_idx <- which(clampRes$U[, ii[i]] > 0)
            Matrix::rowSums(priorMat[gene_idx, col_idx, drop = FALSE]) > 0
        },
        logical(length(nn)) # expected output type
    )

    nnpath <- rowSums(nnpath) > 0

    gene_annot <- ComplexHeatmap::rowAnnotation(
        present = nnpath,
        col = list(present = c("TRUE" = "black", "FALSE" = "beige"))
    )

    top_pathways <- unique(unlist(lapply(ii, function(i) {
        names(sort(clampRes$U[, i], decreasing = TRUE))[seq_len(top.pathway)]
    })))
    gene_idx <- match(nn, rownames(priorMat))
    path_idx <- match(top_pathways, colnames(priorMat))
    pathway_mat <- priorMat[gene_idx, path_idx, drop = FALSE]
    pathway_mat <- pathway_mat[, colSums(pathway_mat) > 0, drop = FALSE]

    ht1 <- ComplexHeatmap::Heatmap(
        data_sub,
        name = "expression",
        show_row_names = TRUE, show_column_names = FALSE,
        cluster_rows = TRUE, cluster_columns = TRUE,
        width = grid::unit(7, "cm"),
        row_dend_width = grid::unit(0, "mm"),
        use_raster = FALSE
    )

    col_fun <- circlize::colorRamp2(c(0, 1), c("white", "black"))
    ht2 <- ComplexHeatmap::Heatmap(
        as.matrix(pathway_mat > 0),
        name = "in pathway", col = col_fun,
        show_row_names = TRUE, show_column_names = TRUE,
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_names_rot = 45,
        row_names_gp = grid::gpar(fontsize = 8),
        column_names_gp = grid::gpar(fontsize = 8),
        width = grid::unit(6, "cm"),
        use_raster = FALSE
    )

    if (Zheat) {
        gene_idx <- match(nn, rownames(clampRes$Z))
        z_sub <- scale(
            clampRes$Z[gene_idx, index, drop = FALSE],
            center = FALSE
        )
        if (!is.null(LV.names)) colnames(z_sub) <- LV.names[index]
        z_col_fun <- circlize::colorRamp2(
            c(0, max(z_sub, na.rm = TRUE)),
            c("#e5f5e0", "#31a354")
        )
        ht_z <- ComplexHeatmap::Heatmap(
            z_sub,
            name = "Z", col = z_col_fun,
            cluster_rows = TRUE, cluster_columns = TRUE,
            show_row_names = TRUE, show_column_names = TRUE,
            column_names_rot = 90,
            width = grid::unit(1, "cm"),
            use_raster = FALSE
        )
    }

    lv_labels <- colnames(clampRes$Z)[
        max.col(scale(clampRes$Z, center = FALSE), ties.method = "first")
    ]
    names(lv_labels) <- rownames(clampRes$Z)
    gene_lv <- lv_labels[rownames(data_sub)]
    row_annot <- ComplexHeatmap::rowAnnotation(
        LV = gene_lv, show_annotation_name = FALSE
    )

    if (!Zheat) {
        ComplexHeatmap::draw(row_annot + ht1 + ht2, row_dend_side = "left")
    } else {
        ComplexHeatmap::draw(
            row_annot + ht_z + ht1 + ht2,
            row_dend_side = "left"
        )
    }
}

# Plotting functions (adapted from PLIER, Mao et al.)
#' Plot the U matrix (pathway-LV associations) as a heatmap
#'
#' Displays the pathway loading matrix `U` after filtering by AUC and FDR
#' thresholds.  Only the top-`top` pathways per LV are shown.
#'
#' @param clampRes A CLAMP result list containing at least `U`, `Uauc`, `Up`,
#'   and `summary`.
#' @param auc.cutoff Minimum AUC threshold; entries below this are set to zero.
#'   Default `0.6`.
#' @param fdr.cutoff Maximum FDR threshold for pathway significance filtering.
#'   Default `0.05`.
#' @param indexCol Integer vector of LV column indices to include. `NULL` uses
#'   all LVs.
#' @param indexRow Integer vector of pathway row indices to include. `NULL` uses
#'   all pathways.
#' @param top Number of top pathways to retain per LV. Default `3`.
#' @param sort.row Logical; if `TRUE`, rows are sorted by the dominant LV.
#'   Default `FALSE`.
#'
#' @param cluster.rows Logical; if `TRUE` (default), rows are reordered by
#'   hierarchical clustering (overridden when `sort.row = TRUE`).
#'
#' @return Invisibly returns a [ggplot2::ggplot()] object.
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient theme_minimal
#'   theme element_text element_blank labs unit
#' @importFrom rlang .data
#' @importFrom stats hclust dist
#' @export
#' @examples
#' set.seed(42)
#' pathways <- paste0("Path", 1:30)
#' lvs <- paste0("LV", 1:5)
#'
#' U <- matrix(abs(rnorm(30 * 5)),
#'     nrow = 30,
#'     dimnames = list(pathways, lvs)
#' )
#' Uauc <- matrix(runif(30 * 5, 0.5, 1.0),
#'     nrow = 30,
#'     dimnames = list(pathways, lvs)
#' )
#' Up <- matrix(runif(30 * 5, 0, 3),
#'     nrow = 30,
#'     dimnames = list(pathways, lvs)
#' )
#'
#' # Build a minimal summary table
#' nr <- 30 * 5
#' summ <- data.frame(
#'     pathway = rep(pathways, 5),
#'     LV = rep(lvs, each = 30),
#'     AUC = as.vector(Uauc),
#'     p_value = runif(nr, 0, 0.1),
#'     FDR = runif(nr, 0, 0.1),
#'     stringsAsFactors = FALSE
#' )
#'
#' clampRes <- list(U = U, Uauc = Uauc, Up = Up, summary = summ)
#' CLAMPplotU(clampRes, auc.cutoff = 0.6, fdr.cutoff = 0.1, top = 3)
CLAMPplotU <- function(clampRes, auc.cutoff = 0.6, fdr.cutoff = 0.05,
                       indexCol = NULL, indexRow = NULL, top = 3,
                       sort.row = FALSE, cluster.rows = TRUE) {
    indexCol <- if (is.null(indexCol)) seq_len(ncol(clampRes$U)) else indexCol
    indexRow <- if (is.null(indexRow)) seq_len(nrow(clampRes$U)) else indexRow

    U <- as.matrix(clampRes$U)

    pval.cutoff <- if (!is.null(clampRes$summary) &&
        any(clampRes$summary[, 5] < fdr.cutoff)) {
        max(clampRes$summary[clampRes$summary[, 5] < fdr.cutoff, 4])
    } else {
        Inf
    }

    U[as.matrix(clampRes$Uauc) < auc.cutoff] <- 0
    U[as.matrix(clampRes$Up) > pval.cutoff] <- 0

    U <- U[indexRow, indexCol, drop = FALSE]

    for (i in seq_len(ncol(U))) {
        ct <- sort(U[, i], decreasing = TRUE)[top]
        if (!is.na(ct)) U[U[, i] < ct, i] <- 0
    }

    rownames(U) <- make.unique(strtrim(rownames(U), 30))
    colnames(U) <- make.unique(strtrim(colnames(U), 30))

    keep_row <- rowSums(abs(U)) > 0
    keep_col <- colSums(abs(U)) > 0
    U <- U[keep_row, keep_col, drop = FALSE]

    if (nrow(U) == 0 || ncol(U) == 0) {
        message("No entries pass the AUC/FDR thresholds.")
        return(invisible(NULL))
    }

    if (sort.row) {
        Utmp <- sweep(sign(U), 2, seq_len(ncol(U)) * 100, "*")
        Um <- apply(Utmp, 1, max)
        U <- U[order(-Um), , drop = FALSE]
    } else if (cluster.rows && nrow(U) > 1) {
        U <- U[hclust(dist(U))$order, , drop = FALSE]
    }

    df <- data.frame(
        pathway = factor(rep(rownames(U), ncol(U)), levels = rev(rownames(U))),
        LV = factor(rep(colnames(U), each = nrow(U)), levels = colnames(U)),
        value = as.vector(U),
        stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(
        x    = .data$LV,
        y    = .data$pathway,
        fill = .data$value
    )) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
        ggplot2::scale_fill_gradient(
            name   = "U",
            low    = "#F7FBFF",
            high   = "#08519C",
            limits = c(0, NA)
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            axis.text.x      = ggplot2::element_text(
                angle = 45, hjust = 1, size = 9
            ),
            axis.text.y      = ggplot2::element_text(size = 7),
            panel.grid       = ggplot2::element_blank(),
            legend.key.width = ggplot2::unit(0.4, "cm")
        ) +
        ggplot2::labs(x = "Latent Variable", y = NULL)

    p
}

#' Plot top genes per LV by Z loading
#'
#' For each selected latent variable, ranks genes by their Z loading and plots
#' the top genes as a loading-versus-rank scatter plot. The highest ranking
#' genes are labelled with `ggrepel`.
#'
#' @param clampRes A CLAMP result list containing at least `Z`, and optionally
#'   `U` for selecting LVs with pathway support.
#' @param data Deprecated; retained for backward compatibility and ignored.
#' @param priorMat Deprecated; retained for backward compatibility and ignored.
#' @param top Number of top genes to plot per LV. Default `50`.
#' @param index Integer or character vector of LV columns to include. `NULL`
#'   keeps LVs with non-zero `U` entries when `U` is present.
#' @param allLVs Logical; if `TRUE`, all LVs are eligible when `index = NULL`.
#'   Default `FALSE`.
#' @param label.top Number of top genes to label per LV. Default `min(10, top)`.
#' @param max.name.len Maximum characters for displayed gene labels.
#'
#' @return A [ggplot2::ggplot()] object for one LV, or a `patchwork` object for
#'   multiple LVs.
#' @importFrom ggplot2 ggplot aes geom_point theme labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang .data
#' @export
#' @examples
#' set.seed(1)
#' genes <- paste0("Gene", 1:80)
#' lvs <- paste0("LV", 1:4)
#' paths <- paste0("Path", 1:10)
#' Z <- matrix(abs(rnorm(80 * 4)),
#'     nrow = 80,
#'     dimnames = list(genes, lvs)
#' )
#' U <- matrix(abs(rnorm(10 * 4)),
#'     nrow = 10,
#'     dimnames = list(paths, lvs)
#' )
#' clampRes <- list(Z = Z, U = U)
#' CLAMPplotTopZ(clampRes, top = 20, index = 1:2)
CLAMPplotTopZ <- function(clampRes, data = NULL, priorMat = NULL, top = 50,
                          index = NULL, allLVs = FALSE,
                          label.top = min(10, top), max.name.len = 50) {
    if (is.null(clampRes$Z)) {
        stop("'clampRes' must contain a 'Z' matrix.")
    }

    Z <- as.matrix(clampRes$Z)
    if (is.null(rownames(Z))) rownames(Z) <- paste0("Gene", seq_len(nrow(Z)))
    if (is.null(colnames(Z))) colnames(Z) <- paste0("LV", seq_len(ncol(Z)))

    if (is.null(index)) {
        if (!allLVs && !is.null(clampRes$U)) {
            ii <- which(colSums(as.matrix(clampRes$U), na.rm = TRUE) > 0)
            if (length(ii) == 0) ii <- seq_len(ncol(Z))
        } else {
            ii <- seq_len(ncol(Z))
        }
    } else if (is.character(index)) {
        ii <- match(index, colnames(Z))
        if (anyNA(ii)) {
            stop(
                "Unknown latent variable(s): ",
                paste(index[is.na(ii)], collapse = ", ")
            )
        }
    } else {
        ii <- as.integer(index)
    }

    ii <- unique(ii[!is.na(ii) & ii >= 1 & ii <= ncol(Z)])
    if (length(ii) == 0) stop("No latent variables selected.")

    top <- as.integer(top)[1]
    if (!is.finite(top) || top < 1) stop("'top' must be a positive integer.")

    label.top <- as.integer(label.top)[1]
    if (!is.finite(label.top) || label.top < 0) {
        stop("'label.top' must be a non-negative integer.")
    }
    label.top <- min(label.top, top)

    make_lv_df <- function(i) {
        loadings <- Z[, i]
        keep <- which(!is.na(loadings))
        if (length(keep) == 0) {
            return(data.frame())
        }

        loadings <- loadings[keep]
        genes <- rownames(Z)[keep]
        ord <- order(loadings, decreasing = TRUE)
        n_top <- min(top, length(ord))
        ord <- ord[seq_len(n_top)]

        data.frame(
            LV = colnames(Z)[i],
            rank = seq_len(n_top),
            gene = genes[ord],
            loading = loadings[ord],
            label = ifelse(seq_len(n_top) <= label.top,
                strtrim(genes[ord], max.name.len),
                NA_character_
            ),
            stringsAsFactors = FALSE
        )
    }

    plot_df <- do.call(rbind, lapply(ii, make_lv_df))
    if (nrow(plot_df) == 0) {
        stop("Selected latent variables have no finite Z loadings.")
    }

    plot_df$LV <- factor(plot_df$LV, levels = colnames(Z)[ii])

    make_plot <- function(lv_name) {
        df <- plot_df[plot_df$LV == lv_name, , drop = FALSE]
        label_df <- df[!is.na(df$label), , drop = FALSE]

        ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$loading)) +
            ggplot2::geom_point(size = 1.2, colour = "grey75") +
            ggplot2::geom_point(
                data = label_df, colour = "#C23B22", size = 1.6
            ) +
            ggrepel::geom_text_repel(
                data = label_df,
                ggplot2::aes(label = .data$label),
                size = 3,
                fontface = "italic",
                direction = "y",
                hjust = 0,
                seed = 42,
                max.overlaps = Inf,
                min.segment.length = 0,
                box.padding = 0.3,
                point.padding = 0.2,
                nudge_x = max(1, top * 0.04),
                segment.color = "grey60",
                segment.size = 0.25
            ) +
            ggplot2::coord_cartesian(clip = "off") +
            ggplot2::scale_x_continuous(
                expand = ggplot2::expansion(mult = c(0.02, 0.18))
            ) +
            ggplot2::labs(
                title = as.character(lv_name), x = NULL, y = "Loadings"
            ) +
            ggplot2::theme_classic(base_size = 11) +
            ggplot2::theme(
                axis.line    = ggplot2::element_line(
                    colour = "black", linewidth = 0.4
                ),
                axis.ticks   = ggplot2::element_line(
                    colour = "black", linewidth = 0.35
                ),
                axis.text.x  = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                plot.title   = ggplot2::element_text(
                    face = "bold", hjust = 0.5, size = 12
                ),
                plot.margin  = ggplot2::margin(5.5, 28, 5.5, 5.5)
            )
    }

    plots <- lapply(levels(plot_df$LV), make_plot)
    if (length(plots) == 1) {
        return(plots[[1]])
    }

    patchwork::wrap_plots(plots, ncol = min(2, length(plots)))
}

#' Dot plot of top pathways for a single latent variable
#'
#' Lollipop-style dot plot showing the top pathways associated with one
#' selected LV. Dot size encodes AUC; dot colour encodes `-log10(FDR)`.
#' The x-axis and pathway ordering can use either AUC or `-log10(FDR)`.
#'
#' @param clampRes A CLAMP result list containing a `summary` data frame, or
#'   the summary data frame itself.
#' @param lv LV to plot: either a numeric index (e.g. `1` -> `"LV1"`) or a
#'   character name (e.g. `"LV3"`). Default `1`.
#' @param top Maximum number of pathways to display, chosen by `order.by`.
#'   Default `20`.
#' @param auc.cutoff Minimum AUC to display. Default `0.6`.
#' @param fdr.cutoff Maximum FDR to display. Default `0.05`.
#' @param max.name.len Maximum characters for pathway label trimming.
#'   Default `50`.
#' @param x.axis Metric to place on the x-axis. Use `"AUC"` or
#'   `"-log10(FDR)"`. `"log10FDR"` is also accepted. Default `"AUC"`.
#' @param order.by Metric used to choose the top pathways and order the y-axis.
#'   Use `"AUC"` or `"-log10(FDR)"`. `"log10FDR"` is also accepted. Defaults
#'   to `x.axis`.
#'
#' @return Invisibly returns a [ggplot2::ggplot()] object.
#' @importFrom ggplot2 ggplot aes geom_segment geom_point scale_size_continuous
#'   scale_colour_gradient scale_x_continuous theme_minimal theme element_text
#'   element_blank element_line labs unit
#' @importFrom rlang .data
#' @export
#' @examples
#' set.seed(9)
#' pathways <- paste0("Pathway_", 1:20)
#' lvs <- paste0("LV", 1:5)
#' nr <- length(pathways) * length(lvs)
#'
#' summ <- data.frame(
#'     pathway = rep(pathways, length(lvs)),
#'     LV = rep(lvs, each = length(pathways)),
#'     AUC = runif(nr, 0.5, 1.0),
#'     FDR = runif(nr, 0, 0.05),
#'     stringsAsFactors = FALSE
#' )
#'
#' CLAMPdotplot(list(summary = summ), lv = 1, top = 10)
#' CLAMPdotplot(list(summary = summ),
#'     lv = "LV2", x.axis = "-log10(FDR)",
#'     order.by = "-log10(FDR)"
#' )
CLAMPdotplot <- function(clampRes, lv = 1, top = 20, auc.cutoff = 0.6,
                         fdr.cutoff = 0.05, max.name.len = 50,
                         x.axis = c("AUC", "-log10(FDR)", "log10FDR"),
                         order.by = x.axis) {
    summ <- if (is.data.frame(clampRes)) clampRes else clampRes$summary
    if (is.null(summ)) {
        stop("'clampRes' must contain a 'summary' data frame or be one itself.")
    }

    if ("LV_index" %in% colnames(summ) && !"LV" %in% colnames(summ)) {
        summ$LV <- paste0("LV", summ$LV_index)
    }

    required <- c("pathway", "LV", "AUC", "FDR")
    missing <- setdiff(required, colnames(summ))
    if (length(missing) > 0) {
        stop("summary is missing columns: ", paste(missing, collapse = ", "))
    }

    normalize_dotplot_metric <- function(x) {
        x <- match.arg(x, c("AUC", "-log10(FDR)", "log10FDR"))
        if (x == "log10FDR") "-log10(FDR)" else x
    }

    x.axis <- normalize_dotplot_metric(x.axis)
    order.by <- normalize_dotplot_metric(order.by)

    lv_name <- if (is.numeric(lv)) paste0("LV", lv) else as.character(lv)
    df <- summ[summ$LV == lv_name &
        summ$AUC >= auc.cutoff &
        summ$FDR <= fdr.cutoff, , drop = FALSE]

    if (nrow(df) == 0) {
        message("No associations pass the thresholds for ", lv_name, ".")
        return(invisible(NULL))
    }

    df$log10_fdr <- -log10(pmax(df$FDR, .Machine$double.eps))
    x_col <- if (x.axis == "AUC") "AUC" else "log10_fdr"
    order_col <- if (order.by == "AUC") "AUC" else "log10_fdr"

    df <- df[order(-df[[order_col]], -df$AUC), , drop = FALSE]
    df <- df[seq_len(min(top, nrow(df))), , drop = FALSE]
    df$x_value <- df[[x_col]]
    df$pathway_label <- make.unique(strtrim(df$pathway, max.name.len))
    df$pathway <- factor(
        df$pathway_label,
        levels = df$pathway_label[order(df[[order_col]], df$AUC)]
    )

    log10_fdr_min <- -log10(max(fdr.cutoff, .Machine$double.eps))
    log10_fdr_max <- max(df$log10_fdr, na.rm = TRUE)
    if (log10_fdr_max <= log10_fdr_min) {
        log10_fdr_max <- log10_fdr_min + 1e-6
    }

    x_baseline <- if (x.axis == "AUC") auc.cutoff else log10_fdr_min
    if (x.axis == "AUC") {
        x_limits <- c(max(0, auc.cutoff - 0.02), 1)
        x_label <- "AUC"
    } else {
        x_range <- range(c(x_baseline, df$x_value), na.rm = TRUE)
        x_pad <- diff(x_range) * 0.05
        if (!is.finite(x_pad) || x_pad == 0) x_pad <- 0.1
        x_limits <- c(max(0, x_range[1] - x_pad), x_range[2] + x_pad)
        x_label <- "-log10(FDR)"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(
        x      = .data$x_value,
        y      = .data$pathway,
        colour = .data$log10_fdr
    )) +
        ggplot2::geom_segment(
            ggplot2::aes(xend = .data$x_value, yend = .data$pathway),
            x = x_baseline,
            colour = "#A1D99B",
            linewidth = 1
        ) +
        ggplot2::geom_point(ggplot2::aes(size = .data$AUC), alpha = 0.9) +
        ggplot2::scale_size_continuous(
            name   = "AUC",
            range  = c(3, 8),
            limits = c(auc.cutoff, 1)
        ) +
        ggplot2::scale_colour_gradient(
            name   = "-log10(FDR)",
            low    = "#DCEFD8",
            high   = "#1B7837",
            limits = c(log10_fdr_min, log10_fdr_max)
        ) +
        ggplot2::scale_x_continuous(
            limits = x_limits,
            expand = c(0, 0.01)
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            plot.title         = ggplot2::element_text(hjust = 0.5),
            axis.text.y        = ggplot2::element_text(size = 9),
            panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
            panel.grid.major.x = ggplot2::element_line(colour = "grey92"),
            panel.grid.minor   = ggplot2::element_blank()
        ) +
        ggplot2::labs(title = lv_name, x = x_label, y = NULL)

    p
}

#' Dot plot of pathway-LV associations across all latent variables
#'
#' Produces a dot plot where each point represents one pathway × LV pair.
#' Dot size encodes the AUC value and dot colour encodes `-log10(FDR)`. Only
#' associations passing `auc.cutoff` and `fdr.cutoff` are shown.
#'
#' @param clampRes A CLAMP result list containing a `summary` data frame with
#'   columns `pathway`, `LV`, `AUC`, and `FDR`.  Alternatively, `clampRes`
#'   may be the summary data frame itself.
#' @param auc.cutoff Minimum AUC to display. Default `0.6`.
#' @param fdr.cutoff Maximum FDR to display. Default `0.05`.
#' @param top.per.lv Maximum number of pathways to display per LV, chosen by
#'   highest AUC. `NULL` shows all. Default `NULL`.
#' @param max.name.len Maximum characters for pathway label trimming.
#'   Default `40`.
#'
#' @return Invisibly returns a [ggplot2::ggplot()] object.
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_gradient
#'   scale_size_continuous theme_minimal theme element_text element_line labs
#' @importFrom dplyr group_by slice_max ungroup
#' @importFrom rlang .data
#' @export
#' @examples
#' set.seed(9)
#' pathways <- paste0("Pathway_", 1:20)
#' lvs <- paste0("LV", 1:5)
#' nr <- length(pathways) * length(lvs)
#'
#' summ <- data.frame(
#'     pathway = rep(pathways, length(lvs)),
#'     LV = rep(lvs, each = length(pathways)),
#'     AUC = runif(nr, 0.5, 1.0),
#'     FDR = runif(nr, 0, 0.2),
#'     stringsAsFactors = FALSE
#' )
#'
#' CLAMPdotplotAll(list(summary = summ), auc.cutoff = 0.6, fdr.cutoff = 0.15)
CLAMPdotplotAll <- function(clampRes, auc.cutoff = 0.6, fdr.cutoff = 0.05,
                            top.per.lv = NULL, max.name.len = 40) {
    summ <- if (is.data.frame(clampRes)) clampRes else clampRes$summary

    if (is.null(summ)) {
        stop("'clampRes' must contain a 'summary' data frame or be one itself.")
    }

    if ("LV_index" %in% colnames(summ) && !"LV" %in% colnames(summ)) {
        summ$LV <- paste0("LV", summ$LV_index)
    }

    required <- c("pathway", "LV", "AUC", "FDR")
    missing <- setdiff(required, colnames(summ))
    if (length(missing) > 0) {
        stop("summary is missing columns: ", paste(missing, collapse = ", "))
    }

    df <- summ[summ$AUC >= auc.cutoff & summ$FDR <= fdr.cutoff, , drop = FALSE]

    if (nrow(df) == 0) {
        message("No associations pass the AUC/FDR thresholds.")
        return(invisible(NULL))
    }

    if (!is.null(top.per.lv)) {
        df <- dplyr::group_by(df, .data$LV)
        df <- dplyr::slice_max(df,
            order_by = .data$AUC, n = top.per.lv,
            with_ties = FALSE
        )
        df <- dplyr::ungroup(df)
    }

    df$pathway <- strtrim(df$pathway, max.name.len)
    df$pathway <- factor(df$pathway,
        levels = rev(unique(df$pathway[order(df$AUC)]))
    )
    df$log10_fdr <- -log10(pmax(df$FDR, .Machine$double.eps))

    log10_fdr_min <- -log10(max(fdr.cutoff, .Machine$double.eps))
    log10_fdr_max <- max(df$log10_fdr, na.rm = TRUE)
    if (log10_fdr_max <= log10_fdr_min) {
        log10_fdr_max <- log10_fdr_min + 1e-6
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(
        x      = .data$LV,
        y      = .data$pathway,
        size   = .data$AUC,
        colour = .data$log10_fdr
    )) +
        ggplot2::geom_point(alpha = 0.85) +
        ggplot2::scale_size_continuous(
            name   = "AUC",
            range  = c(2, 8),
            limits = c(auc.cutoff, 1)
        ) +
        ggplot2::scale_colour_gradient(
            name   = "-log10(FDR)",
            low    = "#DCEFD8",
            high   = "#1B7837",
            limits = c(log10_fdr_min, log10_fdr_max)
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
            axis.text.y      = ggplot2::element_text(size = 8),
            panel.grid.major = ggplot2::element_line(colour = "grey92")
        ) +
        ggplot2::labs(x = "Latent Variable", y = "Pathway")

    p
}
