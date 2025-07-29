context("prior-related functions")

test_that("gmtListToSparseMat builds a sparse matrix", {
    gmtList <- list(
        lib1 = list(A = letters[1:5], B = letters[6:10]),
        lib2 = list(C = letters[5:8], D = letters[9:10])
    )
    mat <- gmtListToSparseMat(gmtList)
    expect_s4_class(mat, "dgCMatrix")

    total_sets <- sum(vapply(gmtList, length, integer(1)))
    expect_equal(ncol(mat), total_sets)

    all_genes <- sort(unique(unlist(unname(unlist(gmtList)))))
    expect_equal(sort(rownames(mat)), all_genes)
})

test_that("getMatchedPathwayMat returns zero-column matrix if no overlap", {
    gmtList <- list(lib = list(A = letters[1:5], B = letters[3:7]))
    mat <- gmtListToSparseMat(gmtList)
    genes <- LETTERS[1:5] # no overlap

    matched <- getMatchedPathwayMat(mat, genes, min.genes = 1)
    expect_s4_class(matched, "dgCMatrix")
    expect_equal(nrow(matched), length(genes))
    expect_equal(ncol(matched), 0)
})
