#' Whole-blood reference expression matrix
#'
#' A numeric matrix of whole-blood gene expression where rows correspond to
#' genes and columns to samples.
#'
#' @format A numeric matrix with G genes (rows) and N samples (columns).
#'   Row names are gene symbols, and column names are sample IDs.
#'
#' @details
#' This object is a whole-blood RNA-seq expression matrix. The source data are
#' publicly available from NCBI GEO under accession
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130824}{GSE130824}
#' (Homo sapiens whole-blood RNA-seq, 36 samples). Raw sequencing data are
#' deposited in SRA, and the processed normalized expression file is provided
#' as \code{GSE130824_dataNormedFiltered.txt.gz}. The matrix bundled with
#' CLAMP is derived from that processed file and serves as a compact example
#' dataset for package demonstrations and unit tests.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130824}
#'
#' @usage data(dataWholeBlood)
#'
#' @keywords datasets
#' @docType data
#' @name dataWholeBlood
#' @aliases dataWholeBlood
#'
#' @return A numeric matrix of expression values.
#'
#' @examples
#' data(dataWholeBlood)
"dataWholeBlood"
