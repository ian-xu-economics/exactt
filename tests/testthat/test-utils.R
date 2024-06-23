# tests/testthat/test-utils.R

library(testthat)

# Does not test case where nBlocks >= 10
test_that("block_permute() works correctly", {
  
  blockIndexMatrix <- matrix(c(1:30), ncol = 5)
  possibleBlockPermutations <- do.call(rbind, combinat::permn(1:5))
  
  # Create X vector
  X <- matrix(sample(1:10, 30, replace = TRUE))
  
  # First Test
  result <- block_permute(X, blockIndexMatrix, possibleBlockPermutations)
  expect_equal(result$shuffledData, 
               X[c(blockIndexMatrix[, possibleBlockPermutations[1,]])])
  
  expect_equal(result$possibleBlockPermutations,
               possibleBlockPermutations[-1,])
  
  possibleBlockPermutations <- possibleBlockPermutations[-1,]
  
  # Second Test
  result <- block_permute(X, blockIndexMatrix, possibleBlockPermutations)
  expect_equal(result$shuffledData, 
               X[c(blockIndexMatrix[, possibleBlockPermutations[1,]])])
  
  expect_equal(result$possibleBlockPermutations,
               possibleBlockPermutations[-1,])
  
})

# Does not test case where nBlocks >= 10
test_that("build_GX() works correctly", {
  
  blockIndexMatrix <- matrix(c(1:30), ncol = 5)
  
  # Create X vector
  set.seed(10)
  X <- matrix(sample(1:10, 30, replace = TRUE))
  
  # First Test
  result <- build_GX(X, blockIndexMatrix)
  
  expect_equal(unname(result), 
               matrix(scan(file = "tests/testthat/expected_values_build_GX.txt", 
                           quiet = TRUE), 
                      ncol = 17))
})

test_that("build_GX2() and build_QGX2() works correctly", {
  
  blockIndexMatrix <- matrix(c(1:30), ncol = 5)
  
  # Create X vector
  set.seed(10)
  X1 <- matrix(sample(x = 1:10, size = 30, replace = TRUE))
  X2 <- matrix(sample(x = 11:20, size = 30, replace = TRUE))
  
  X <- cbind(X1, X2) 
  
  # Test build_GX2()
  result <- unname(build_GX2(X, blockIndexMatrix))
  
  expect_equal(unname(result), 
               matrix(scan(file = "tests/testthat/expected_values_build_GX2.txt", 
                           quiet = TRUE), 
                      ncol = 34))
  
  # Test build_QGX2()
  result <- build_QGX2(result)
  
  expect_equal(unname(result), 
               matrix(scan(file = "tests/testthat/expected_values_build_QGX2.txt", 
                           quiet = TRUE), 
                      ncol = 30),
               tolerance = 1e-6)
  
})


