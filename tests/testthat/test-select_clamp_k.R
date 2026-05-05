test_that("select_clamp_k returns a list with clamp_k and scale", {
    set.seed(1)
    Y <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    svdres <- compute_svd(Y, k = 25)
    res <- select_clamp_k(svdres, n_samples = 30, svd_k = 25)
    expect_true(all(c("clamp_k", "scale") %in% names(res)))
    expect_true(is.numeric(res$clamp_k))
    expect_true(is.numeric(res$scale))
})

test_that("select_clamp_k caps clamp_k at svd_k (default method)", {
    set.seed(1)
    Y <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    svdres <- compute_svd(Y, k = 25)
    res <- select_clamp_k(svdres, n_samples = 30, svd_k = 3)
    expect_lte(res$clamp_k, 3)
})

test_that("select_clamp_k caps clamp_k at svd_k (scaleSVs method)", {
    set.seed(1)
    d <- sort(runif(30, min = 1, max = 100), decreasing = TRUE)
    svdres <- list(d = d)
    res <- select_clamp_k(svdres,
        n_samples = 50, svd_k = 5, method = "scaleSVs"
    )
    expect_lte(res$clamp_k, 5)
})

test_that("select_clamp_k scaleSVs matches inline behavior", {
    # regression check: scaleSVs method should reproduce the former inline logic
    set.seed(1)
    d <- sort(runif(30, min = 1, max = 100), decreasing = TRUE)
    svdres <- list(d = d)
    svd_k <- 30
    n_samples <- 50

    scale.res <- getScaleFromSVs(svdres$d, n_samples)
    expected_k <- min(floor(scale.res$k * 1.5), svd_k)
    expected_scale <- scale.res$scale

    res <- select_clamp_k(svdres,
        n_samples = n_samples, svd_k = svd_k,
        method = "scaleSVs"
    )
    expect_equal(res$clamp_k, expected_k)
    expect_equal(res$scale, expected_scale)
})

test_that("select_clamp_k default method is 'elbow'", {
    set.seed(1)
    Y <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    svdres <- compute_svd(Y, k = 25)
    res_default <- select_clamp_k(svdres, n_samples = ncol(Y), svd_k = 25)
    res_elbow <- select_clamp_k(svdres,
        n_samples = ncol(Y), svd_k = 25,
        method = "elbow"
    )
    expect_equal(res_default, res_elbow)
})

test_that(paste(
    "select_clamp_k non-scaleSVs methods return",
    "scale = svdres$d[clamp_k]"
), {
    set.seed(1)
    Y <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    svdres <- compute_svd(Y, k = 25)
    res <- select_clamp_k(svdres,
        n_samples = ncol(Y), svd_k = 25,
        method = "elbow"
    )
    expect_equal(res$scale, svdres$d[res$clamp_k])
})

test_that("select_clamp_k permutation method works with raw data", {
    set.seed(1)
    Y <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    svdres <- compute_svd(Y, k = 25)
    res <- select_clamp_k(svdres,
        n_samples = ncol(Y), svd_k = 25,
        method = "permutation", data = Y, B = 3
    )
    expect_true(is.numeric(res$clamp_k))
    expect_lte(res$clamp_k, 25)
    expect_equal(res$scale, svdres$d[res$clamp_k])
})

test_that("select_clamp_k permutation method errors without data", {
    set.seed(1)
    Y <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    svdres <- compute_svd(Y, k = 25)
    expect_error(
        select_clamp_k(svdres,
            n_samples = ncol(Y), svd_k = 25,
            method = "permutation"
        ),
        "data"
    )
})

test_that("select_clamp_k gavish_donoho method works with raw data", {
    skip_if_not_installed("PCAtools")
    set.seed(1)
    Y <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    svdres <- compute_svd(Y, k = 25)
    res <- select_clamp_k(svdres,
        n_samples = ncol(Y), svd_k = 25,
        method = "gavish_donoho", data = Y
    )
    expect_true(is.numeric(res$clamp_k))
    expect_lte(res$clamp_k, 25)
    expect_equal(res$scale, svdres$d[res$clamp_k])
})

test_that("select_clamp_k gavish_donoho method errors without data", {
    set.seed(1)
    Y <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    svdres <- compute_svd(Y, k = 25)
    expect_error(
        select_clamp_k(svdres,
            n_samples = ncol(Y), svd_k = 25,
            method = "gavish_donoho"
        ),
        "data"
    )
})
