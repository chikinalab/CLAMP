#' Compare two sets of factor loadings or embeddings
#'
#' Compares correspondence between two result matrices (e.g., factor loadings)
#' across multiple targets using correlation, AUC, or t-statistics.
#' Produces a paired comparison plot and summary statistics.
#'
#' @param res1, res2 Result objects or matrices containing factor loadings.
#'   If a list, must contain element `B`; if of class `"rsvd"`, `t(res$v)` is used.
#' @param target Numeric matrix of target variables (samples × targets).
#' @param method Character, one of `"p"`, `"s"`, `"a"`, or `"t"`, indicating
#'   Pearson, Spearman, AUC, or t-statistic comparison.
#' @param xlab, ylab Labels for x and y axes in the plot.
#' @param stat.method Statistical test to compare correlations (`"t"` or `"wilcox"`).
#' @param oneToOne Logical, whether to apply one-to-one masking of associations.
#'
#' @return A list with:
#' \describe{
#'   \item{plot}{A ggplot object comparing maximal correlations across targets.}
#'   \item{df}{A data.frame with per-target metrics and (if available) top genes.}
#' }
#'
#' @importFrom matrixStats colRanks
#' @importFrom ggplot2 ggplot aes geom_point geom_abline  labs annotate theme_minimal
#' @export
compareBs <- function(res1, res2, target, method = "p", xlab = "1", ylab = "2",
                      stat.method = "t", oneToOne = TRUE) {

  extract_B <- function(res) {
    if (is.list(res) && !is.null(res$B)) {

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

  if (method %in% c("s", "p")) {
    mat1 <- cor(t(B1), target, method = method)
    mat2 <- cor(t(B2), target, method = method)
  } else if (method == "a") {
    mat1 <- allAgainstAllAUCs(t(B1), target)
    mat2 <- allAgainstAllAUCs(t(B2), target)
  } else if (method == "t") {
    mat1 <- allAgainstAllTstats(t(B1), target)
    mat2 <- allAgainstAllTstats(t(B2), target)
  }
  if (inherits(res1, "rsvd"))
    mat1=abs(mat1)
  if (inherits(res2, "rsvd"))
    mat1=abs(mat1)
  mat1[is.na(mat1)] <- 0
  mat2[is.na(mat2)] <- 0
  if (oneToOne) {
    mat1 <- oneToOneMask(mat1)
    mat2 <- oneToOneMask(mat2)
  }

  cor1 <- apply(mat1, 2, max)
  cor2 <- apply(mat2, 2, max)
  cor1[cor1 < -50] <- NA
  cor2[cor2 < -50] <- NA
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
    wilcox.test(df$Cor2, df$Cor1, paired = TRUE, alternative = "greater")$p.value
  }

  method_label <- switch(method,
                         p = "Pearson correlation",
                         s = "Spearman correlation",
                         a = "AUC",
                         t = "T-statistic"
  )

  pl <- ggplot2::ggplot(df, ggplot2::aes(x = Cor1, y = Cor2, label = Label)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    ggrepel::geom_text_repel() +
    ggplot2::labs(
      x = paste("Max", method_label, xlab),
      y = paste("Max", method_label, ylab)
    ) +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
                      label = paste("p =", signif(pval, 3))) +
    ggplot2::theme_minimal()

  if (is.list(res1) && !is.null(res1$Z)&&is.list(res2) && !is.null(res2$Z)) {
    get_top_genes <- function(Z, k) {
      apply(Z, 2, function(col) {
        names(sort(col, decreasing = TRUE))[seq_len(min(k, length(col)))]
      })
    }
    topk1 <- get_top_genes(res1$Z, 7)
    topk2 <- get_top_genes(res2$Z, 7)
    df$TopGenesRes1 <- vapply(idx1, function(i) paste(topk1[, i], collapse = ","), "")
    df$TopGenesRes2 <- vapply(idx2, function(i) paste(topk2[, i], collapse = ","), "")
  }

  return(list(plot = pl, df = df))
}


#' One-to-one masking of maximum associations
#'
#' Selects the highest-scoring one-to-one pairs between rows and columns
#' of a matrix, similar to a greedy bipartite matching. All other entries
#' are set to a sentinel value (default \code{-100}).
#'
#' @param cc A numeric matrix of association scores.
#'
#' @return A numeric matrix of the same dimensions as \code{cc},
#'   where only the selected one-to-one maxima are retained
#'   and all other entries are set to \code{-100}.
#'
#' @examples
#' set.seed(1)
#' m <- matrix(runif(16), 4, 4)
#' oneToOneMask(m)
#'
#' @export
oneToOneMask <- function(cc) {
  cc <- as.matrix(cc)
  cc_out <- matrix(-100, nrow(cc), ncol(cc))
  tmp <- cc
  for (i in seq_len(ncol(cc))) {
    imax <- which(tmp == max(tmp, na.rm = TRUE), arr.ind = TRUE)[1, ]
    cc_out[imax[1], imax[2]] <- cc[imax[1], imax[2]]
    tmp[imax[1], ] <- -100
    tmp[, imax[2]] <- -100
  }
  return(cc_out)
}



#' ComplexHeatmap visualization of top genes by latent variable
#'
#' Produces a combined heatmap showing expression, pathway membership,
#' and optionally Z-loadings for top genes per LV using ComplexHeatmap.
#'
#' @param plierRes A PLIER result list containing matrices \code{Z}, \code{U}, etc.
#' @param data Expression matrix with genes as rows.
#' @param priorMat Binary gene × pathway matrix.
#' @param top Integer, number of top genes per LV.
#' @param top.pathway Integer, number of top pathways per LV to annotate.
#' @param index Optional vector of LVs to include.
#' @param allLVs Logical; include all LVs.
#' @param Zheat Logical; whether to include the Z matrix as an additional heatmap.
#' @param LV.names Optional character vector for LV names.
#' @param max.genes Maximum number of genes allowed in the plot.
#' @param max.col Maximum number of columns (samples).
#' @param seed Random seed for column subsampling.
#'
#' @return Invisibly returns the drawn ComplexHeatmap object.
#'
#' @examples
#' \dontrun{
#' plotTopZ_Complex(plierRes, expr_data, priorMat, top = 15, index = 1:5, Zheat = TRUE)
#' }
#'
#' @import ComplexHeatmap
#' @import circlize
#' @importFrom stats setNames
#' @importFrom Matrix rowSums
#' @export
plotTopZ_Complex <- function(plierRes, data, priorMat, top = 10, top.pathway = 5,
                             index = NULL, allLVs = FALSE, Zheat = FALSE,
                             LV.names = NULL, max.genes = 100, max.col = 50, seed = 1234) {

  data <- data[rownames(plierRes$Z), , drop = FALSE]
  if (top * length(index) > max.genes)
    stop("Too many genes. Reduce number of LVs or 'top', or increase 'max.genes'.")

  if (ncol(data) > max.col) {
    set.seed(seed)
    data <- data[, sample(ncol(data), max.col)]
  }

  priorMat <- priorMat[rownames(plierRes$Z), , drop = FALSE]
  ii <- which(colSums(plierRes$U) > 0)
  if (!allLVs) {
    if (!is.null(index)) ii <- intersect(ii, index)
  } else {
    ii <- index
  }

  tmp <- apply(-plierRes$Z[, ii, drop = FALSE], 2, rank)
  nn <- unique(unlist(apply(tmp, 2, function(x) names(which(x <= top)))))
  nn <- sort(unique(nn))
  data_sub <- t(scale(t(data[nn, , drop = FALSE])))

  nnpath <- sapply(seq_along(ii), function(i) {
    gene_idx <- match(nn, rownames(priorMat))
    col_idx <- which(plierRes$U[, ii[i]] > 0)
    Matrix::rowSums(priorMat[gene_idx, col_idx, drop = FALSE]) > 0
  })
  nnpath <- rowSums(nnpath) > 0
  gene_annot <- ComplexHeatmap::rowAnnotation(
    present = nnpath,
    col = list(present = c("TRUE" = "black", "FALSE" = "beige"))
  )

  top_pathways <- unique(unlist(lapply(ii, function(i)
    names(sort(plierRes$U[, i], decreasing = TRUE))[seq_len(top.pathway)])))
  gene_idx <- match(nn, rownames(priorMat))
  path_idx <- match(top_pathways, colnames(priorMat))
  pathway_mat <- priorMat[gene_idx, path_idx, drop = FALSE]
  pathway_mat <- pathway_mat[, colSums(pathway_mat) > 0, drop = FALSE]

  ht1 <- ComplexHeatmap::Heatmap(
    data_sub, name = "expression",
    show_row_names = TRUE, show_column_names = FALSE,
    cluster_rows = TRUE, cluster_columns = TRUE,
    width = grid::unit(7, "cm"),
    row_dend_width = grid::unit(0, "mm")
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
    width = grid::unit(6, "cm")
  )

  if (Zheat) {
    gene_idx <- match(nn, rownames(plierRes$Z))
    z_sub <- scale(plierRes$Z[gene_idx, index, drop = FALSE], center = FALSE)
    if (!is.null(LV.names)) colnames(z_sub) <- LV.names[index]
    z_col_fun <- circlize::colorRamp2(c(0, max(z_sub, na.rm = TRUE)),
                                      c("#e5f5e0", "#31a354"))
    ht_z <- ComplexHeatmap::Heatmap(
      z_sub, name = "Z", col = z_col_fun,
      cluster_rows = TRUE, cluster_columns = TRUE,
      show_row_names = TRUE, show_column_names = TRUE,
      column_names_rot = 90,
      width = grid::unit(1, "cm")
    )
  }

  lv_labels <- colnames(plierRes$Z)[max.col(scale(plierRes$Z, center = FALSE), ties.method = "first")]
  names(lv_labels) <- rownames(plierRes$Z)
  gene_lv <- lv_labels[rownames(data_sub)]
  row_annot <- ComplexHeatmap::rowAnnotation(LV = gene_lv, show_annotation_name = FALSE)

  if (!Zheat) {
    ComplexHeatmap::draw(row_annot + ht1 + ht2, row_dend_side = "left")
  } else {
    ComplexHeatmap::draw(row_annot + ht_z + ht1 + ht2, row_dend_side = "left")
  }
}



library(ggplot2)
library(ggrepel)

library(ggplot2)
library(ggrepel)

plot_xy<- function(x, y, n_labels = 10, title = NULL, color = NULL) {
  stopifnot(length(x) == length(y))
  if (!is.null(color)) stopifnot(length(color) == length(x))

  # default labels = names(x) or fallback to row numbers
  lbl <- names(x)
  if (is.null(lbl)) lbl <- seq_along(x)

  # rescale to the same range
  r  <- range(c(x, y), na.rm = TRUE)
  xs <- scales::rescale(x, to = r)
  ys <- scales::rescale(y, to = r)

  df <- data.frame(x = xs, y = ys, label = lbl)
  if (!is.null(color)) df$color <- color

  # allocate labels: x-high, y-high, |x−y|-high
  k <- floor(n_labels / 3)
  idx_x  <- order(df$x, decreasing = TRUE)[1:k]
  idx_y  <- order(df$y, decreasing = TRUE)[1:k]
  idx_xy <- order(abs(df$x - df$y), decreasing = TRUE)[1:(n_labels - 2*k)]

  idx <- unique(c(idx_x, idx_y, idx_xy))

  df$lab_show <- NA
  df$lab_show[idx] <- df$label[idx]

  ggplot(df, aes(x, y)) +
    geom_point(aes(color = if (!is.null(color)) color else NULL), alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray60") +
    geom_text_repel(aes(label = lab_show), na.rm = TRUE, size = 3) +
    coord_equal() +
    labs(title = title, x = "x (scaled)", y = "y (scaled)", color = NULL) +
    scale_color_viridis_c(option = "A")+
    theme_bw()
}

