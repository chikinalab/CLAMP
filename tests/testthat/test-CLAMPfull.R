test_that("CLAMPfull and returns B, Z, U", {
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

    full <- CLAMPfull(
        Y = dataWholeBlood,
        priorMat = matchedPaths,
        clamp.base.result = base,
        trace = FALSE,
        max.iter = 1,
        max.U.updates = 0
    )

    expect_type(full, "list")
    expect_true(all(c("B", "Z", "U") %in% names(full)))
    expect_true(is.matrix(full$B))
    expect_true(is.matrix(full$Z))
    expect_true(is.matrix(full$U))
    expect_equal(rownames(full$Z), rownames(dataWholeBlood))
    expect_equal(colnames(full$B), colnames(dataWholeBlood))
    expect_equal(rownames(full$B), colnames(full$Z))
    expect_equal(colnames(full$U), rownames(full$B))
    expect_equal(rownames(full$U), colnames(matchedPaths))
})

test_that("CLAMPfullnVP returns matrix B and Z with labels", {
    set.seed(2)

    mat <- matrix(
        rnorm(80),
        nrow = 10,
        dimnames = list(paste0("Gene", 1:10), paste0("Sample", 1:8))
    )
    prior <- matrix(
        sample(0:1, 30, replace = TRUE),
        nrow = 10,
        dimnames = list(rownames(mat), paste0("Path", 1:3))
    )
    prior[, 1] <- 1
    svdres <- rsvd::rsvd(mat, k = 4)

    base <- CLAMPbase(
        Y = mat,
        clamp_k = 3,
        svdres = svdres,
        trace = FALSE,
        max.iter = 1
    )

    full <- CLAMPfullnVP(
        Y = mat,
        priorMat = prior,
        svdres = svdres,
        clamp.base.result = base,
        clamp_k = 3,
        doCrossval = FALSE,
        trace = FALSE,
        max.iter = 1,
        max.U.updates = 0,
        minGenes = 0
    )

    expect_true(is.matrix(full$B))
    expect_true(is.matrix(full$Z))
    expect_equal(rownames(full$Z), rownames(mat))
    expect_equal(colnames(full$B), colnames(mat))
    expect_equal(rownames(full$B), colnames(full$Z))
})
