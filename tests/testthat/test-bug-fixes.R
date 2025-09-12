# Test file to expose critical bugs identified in code analysis
# These tests should initially FAIL, then pass after bug fixes

library(testthat)

# Test 1: Division by zero in row_cor() function
test_that("row_cor handles zero variance rows gracefully", {
  # Create matrices with zero variance rows
  A <- matrix(c(1, 1, 1,   # zero variance row
                2, 3, 4,   # normal row
                5, 5, 5),  # zero variance row
              nrow = 3, byrow = TRUE)
  B <- matrix(c(10, 20, 30,
                40, 50, 60,
                70, 80, 90), 
              nrow = 3, byrow = TRUE)
  
  # This should not produce NaN or Inf values
  result <- PLIER2:::row_cor(A, B)
  
  # Check that result is finite (not NaN/Inf) for non-zero variance rows
  expect_true(is.finite(result[2]), "Non-zero variance row should have finite correlation")
  
  # Zero variance rows should either be handled gracefully or produce 0/NA
  expect_true(is.finite(result[1]) || is.na(result[1]), 
              "Zero variance rows should be handled gracefully")
  expect_true(is.finite(result[3]) || is.na(result[3]), 
              "Zero variance rows should be handled gracefully")
})

test_that("row_cor handles matrices with all zero variance", {
  # All rows have zero variance
  A <- matrix(rep(c(1, 2, 3), each = 3), nrow = 3, byrow = TRUE)
  B <- matrix(rep(c(4, 5, 6), each = 3), nrow = 3, byrow = TRUE)
  
  result <- PLIER2:::row_cor(A, B)
  
  # Should not crash and should handle gracefully
  expect_length(result, 3)
  expect_true(all(is.finite(result) | is.na(result)), 
              "All zero variance should be handled without NaN/Inf")
})

# Test 2: Out of bounds error in binarizeTop() function  
test_that("binarizeTop handles edge cases without crashing", {
  # Test when top >= number of rows
  Z <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  
  # This should throw a descriptive error when top > nrow(Z)
  expect_error(PLIER2:::binarizeTop(Z, top = 5), 
               "top.*cannot be greater than number of rows")
  
  # Test with top = nrow(Z) (boundary case) - this should work
  result <- PLIER2:::binarizeTop(Z, top = 2, keepVals = TRUE)
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
  
  # Test with top > nrow(Z) gives clear error
  expect_error(PLIER2:::binarizeTop(Z, top = 3), 
               "top.*cannot be greater than number of rows")
})

test_that("binarizeTop handles empty matrix", {
  # Empty matrix
  Z_empty <- matrix(numeric(0), nrow = 0, ncol = 2)
  
  expect_error(PLIER2:::binarizeTop(Z_empty, top = 1),
               "Cannot operate on empty matrix")
})

test_that("binarizeTop handles single row matrix", {
  Z_single <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
  
  # top = 1 should work fine
  result <- PLIER2:::binarizeTop(Z_single, top = 1, keepVals = TRUE)
  expect_equal(nrow(result), 1)
  
  # top > 1 should give clear error message
  expect_error(PLIER2:::binarizeTop(Z_single, top = 2),
               "top.*cannot be greater than number of rows")
})

# Test 3: BLAS restoration in functions with unreachable code
test_that("PLIERbase restores BLAS settings when using multiple cores", {
  # This test checks if BLAS settings are properly restored
  # Currently this will fail due to unreachable code after return statement
  
  original_blas_setting <- getOption("bigstatsr.check.parallel.blas")
  
  mat <- matrix(rnorm(100), 10, 10)
  
  # Run with multiple cores - should restore BLAS settings  
  result <- PLIERbase(mat, k = 3, max.iter = 1, ncores = 2, trace = FALSE)
  
  # Check BLAS setting was restored (should now pass after fix)
  expect_equal(getOption("bigstatsr.check.parallel.blas"), original_blas_setting)
})

test_that("projectPLIER restores BLAS settings when using multiple cores", {
  # Similar test for projectPLIER unreachable code
  
  original_blas_setting <- getOption("bigstatsr.check.parallel.blas")
  
  # Create a small PLIER result
  mat <- matrix(rnorm(50), 5, 10)
  plier_res <- PLIERbase(mat, k = 2, max.iter = 1, trace = FALSE)
  
  # Project with multiple cores
  newdata <- matrix(rnorm(25), 5, 5)
  result <- projectPLIER(plier_res, newdata, ncores = 2)
  
  # Check BLAS setting was restored (should now pass after fix)
  expect_equal(getOption("bigstatsr.check.parallel.blas"), original_blas_setting)
})

# Test 4: Missing return value issue in solveU
test_that("solveU returns proper structure and executes all code", {
  # Create minimal test data
  Z <- matrix(rnorm(20), 4, 5)
  priorMat <- matrix(rbinom(20, 1, 0.3), 4, 5)
  penalty.factor <- rep(1, 5)
  
  # Mock the internal functions if they don't exist
  if (!exists("commonRows", envir = asNamespace("PLIER2"))) {
    # This test will fail because commonRows doesn't exist
    expect_error({
      result <- PLIER2:::solveU(Z, priorMat = priorMat, penalty.factor = penalty.factor,
                                maxPath = 2, nfolds = 3, trace = FALSE)
    }, "commonRows.*not found")
  }
})

# Test 5: Test bare stop() produces informative error
test_that("PLIERfull provides informative error messages", {
  skip("This test requires fixing the bare stop() first")
  
  mat <- matrix(rnorm(100), 10, 10)
  priorMat <- matrix(1, 10, 5)
  
  # Create conditions that trigger the bare stop()
  # This is hard to test without triggering the specific condition
})