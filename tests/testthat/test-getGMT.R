test_that("getGMT reads a GMT file into a named list", {
    url <- paste0(
        "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=",
        "GTEx_Tissues_V8_2023"
    )
    gmt <- getGMT(url)

    expect_type(gmt, "list")
    expect_true(length(gmt) > 0)
    expect_true(is.character(gmt[[1]]))
    expect_true(length(gmt[[1]]) > 0)
})
