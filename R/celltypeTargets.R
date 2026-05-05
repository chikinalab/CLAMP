#' Cell-type deconvolution matrix
#'
#' A numeric matrix of estimated cell-type proportions for whole-blood samples.
#' Rows correspond to sample IDs and columns to major immune cell types.
#' This dataset can be used for validation or illustrative purposes in
#' CLAMP analyses.
#'
#' @format A numeric matrix with samples as rows and cell types as columns.
#'   Row names are sample identifiers.
#' @usage data(celltypeTargets)
#' @keywords datasets
#' @docType data
#' @name celltypeTargets
#' @aliases celltypeTargets
#' @examples
#' data(celltypeTargets)
#' @return A numeric matrix of cell-type proportions.
"celltypeTargets"
