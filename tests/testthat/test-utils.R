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
               matrix(scan(text = "9 10 7 8 6 7 3 8 10 7 10 2 8 8 7 6 7 6 2 5 
               9 2 10 5 10 1 7 10 2 9 9 10 7 8 6 7 3 8 10 7 10 2 8 8 7 6 7 6 
               10 1 7 10 2 9 2 5 9 2 10 5 9 10 7 8 6 7 3 8 10 7 10 2 10 1 7 10 
               2 9 8 8 7 6 7 6 2 5 9 2 10 5 9 10 7 8 6 7 10 1 7 10 2 9 3 8 10 
               7 10 2 8 8 7 6 7 6 2 5 9 2 10 5 10 1 7 10 2 9 9 10 7 8 6 7 3 8 
               10 7 10 2 8 8 7 6 7 6 2 5 9 2 10 5 10 1 7 10 2 9 9 10 7 8 6 7 3 
               8 10 7 10 2 2 5 9 2 10 5 8 8 7 6 7 6 9 10 7 8 6 7 3 8 10 7 10 2 
               2 5 9 2 10 5 10 1 7 10 2 9 8 8 7 6 7 6 9 10 7 8 6 7 2 5 9 2 10 5
               3 8 10 7 10 2 8 8 7 6 7 6 10 1 7 10 2 9 9 10 7 8 6 7 2 5 9 2 10 
               5 10 1 7 10 2 9 3 8 10 7 10 2 8 8 7 6 7 6 10 1 7 10 2 9 2 5 9 2 
               10 5 9 10 7 8 6 7 3 8 10 7 10 2 8 8 7 6 7 6 2 5 9 2 10 5 10 1 7 
               10 2 9 9 10 7 8 6 7 3 8 10 7 10 2 8 8 7 6 7 6 2 5 9 2 10 5 9 10 
               7 8 6 7 8 8 7 6 7 6 10 1 7 10 2 9 3 8 10 7 10 2 9 10 7 8 6 7 8 8 
               7 6 7 6 2 5 9 2 10 5 3 8 10 7 10 2 10 1 7 10 2 9 8 8 7 6 7 6 9 10 
               7 8 6 7 3 8 10 7 10 2 2 5 9 2 10 5 10 1 7 10 2 9 8 8 7 6 7 6 2 5 
               9 2 10 5 10 1 7 10 2 9 9 10 7 8 6 7 3 8 10 7 10 2 2 5 9 2 10 5 8 
               8 7 6 7 6 3 8 10 7 10 2 10 1 7 10 2 9 9 10 7 8 6 7 3 8 10 7 10 2 
               8 8 7 6 7 6 9 10 7 8 6 7 2 5 9 2 10 5 10 1 7 10 2 9", quiet = TRUE), 
                      ncol = 17))
})
