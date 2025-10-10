test_that("getMatchedPathwayMat returns zero‐column matrix if no overlap", {
    simple <- sparseMatrix(
        i = c(1, 2), j = c(1, 1), x = 1,
        dimnames = list(c("g1", "g2"), "p1")
    )
    out <- getMatchedPathwayMat(simple, new.genes = c("x", "y"), min.genes = 1)
    expect_s4_class(out, "dgCMatrix")
    expect_equal(ncol(out), 0)
})

test_that("getMatchedPathwayMat filters pathways by gene overlap", {
    m1 <- sparseMatrix(
        i = c(1, 2, 3), j = c(1, 1, 2), x = 1,
        dimnames = list(c("g1", "g2", "g3"), c("p1", "p2"))
    )
    genes <- c("g1", "g3", "g4")
    matched <- getMatchedPathwayMat(m1, genes, min.genes = 1)

    expect_s4_class(matched, "dgCMatrix")
    expect_equal(rownames(matched), genes)
    expect_equal(sort(colnames(matched)), c("p1", "p2"))
})
