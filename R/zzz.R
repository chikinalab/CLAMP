#' @keywords internal
#' @noRd
#'
#' @importFrom stats cor sd var t.test smooth
#' @importFrom stats quantile median p.adjust coef setNames na.omit
#' @importFrom utils download.file
#' @importFrom utils data flush.console globalVariables
#' @importFrom Matrix Matrix sparseMatrix colSums
#' @importFrom Matrix colMeans t crossprod tcrossprod
#' @importFrom bigstatsr big_cprodMat big_prodMat big_apply rows_along FBM
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom rsvd rsvd
#' @importFrom irlba irlba
#' @importFrom ggplot2 ggplot aes annotate
#' @importFrom ggplot2 geom_point geom_abline labs theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate select arrange group_by ungroup summarize %>%

NULL

# silence NSE notes from ggplot/data.frame column names
utils::globalVariables(c("Cor1", "Cor2", "Label"))
