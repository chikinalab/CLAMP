test_that("CLAMPbase and returns B, Z, U", {
    set.seed(1)

    data("dataWholeBlood", package = "CLAMP")
    data("xCell", package = "CLAMP")

    matchedPaths <- getMatchedPathwayMatList(
        xCell,
        new.genes = rownames(dataWholeBlood),
        min.genes = 3
    )

    base <- CLAMPbase(
        Y = dataWholeBlood,
        trace = FALSE,
        adaptive.p = 0.05
    )

    expect_type(base, "list")
    expect_true(all(c("B", "Z") %in% names(base)))
    expect_true(is.matrix(base$B))
    expect_true(is.matrix(base$Z))
    expect_equal(rownames(base$Z), rownames(dataWholeBlood))
    expect_equal(colnames(base$B), colnames(dataWholeBlood))
    expect_equal(rownames(base$B), colnames(base$Z))
    expect_equal(rownames(base$B), paste0("LV", seq_len(nrow(base$B))))
})
