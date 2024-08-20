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
  
  # First Test
  result <- build_GX(blockIndexMatrix)
  
  #saveRDS(result, "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/_BFI Predoc/Pouliot/exactt/tests/testthat/expected_values_build_GX.rds")
  expect_equal(result, 
               readRDS(test_path("expected_values_build_GX.rds")))
})

test_that("build_GX2(), build_QGX2(), and build_QGX1GX2() work correctly", {
  
  n <- 50
  
  blockIndexMatrix <- matrix(c(1:n), ncol = 5)
  
  GX.indices <- build_GX(blockIndexMatrix)
  
  # Create X vector
  set.seed(10)
  X1 <- matrix(rexp(n))
  X2 <- matrix(rnorm(n))
  
  X <- cbind(X1, X2) 
  
  # Test build_GX2()
  result.GX2 <- build_GX2(X, GX.indices)
  
  #saveRDS(result.GX2, "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/_BFI Predoc/Pouliot/exactt/tests/testthat/expected_values_build_GX2.rds")
  expect_equal(result.GX2, 
               readRDS(test_path("expected_values_build_GX2.rds")))
  
  # Test build_QGX2()
  result.QGX2 <- build_QGX2(result.GX2)
  
  expect_equal(result.QGX2, 
               readRDS(test_path("expected_values_build_QGX2.rds")))
  
  # Test build_QGX1GX2()
  result.QGX1GX2 <- build_QGX1GX2(X1, result.GX2, GX.indices, GX1 = TRUE)
  
  expect_equal(result.QGX1GX2, 
               readRDS(test_path("expected_values_build_QGX1GX2.rds")))
})

test_that("exactt_pval() works correctly", {
  
  set.seed(31740)
  n = 50
  X1 <- matrix(rbinom(n, size = 8, prob = 0.5))
  X2 <- matrix(rexp(n))
  eps <- matrix(rnorm(n, sd = 2))
  Y <- matrix(3*X1 + 2*X2 + eps)
  
  nBlocks = 5
  blockIndexMatrix = matrix(1:n, ncol = nBlocks)
  blockPermutations <- do.call(rbind, combinat::permn(1:nBlocks))
  permIndices <- apply(blockPermutations, MARGIN = 1, function (x) c(blockIndexMatrix[, x]))
  GX.indices <- build_GX(blockIndexMatrix)
  
  betaNullVec = seq(2, 3, 0.1)
  
  # First test (GX1)
  result.1 <- exactt_pval(betaNullVec, 
                          Y.temp = Y, 
                          X1.temp = X1, 
                          X2.temp = X2, 
                          nBlocks, 
                          permIndices, 
                          GX.indices,
                          studentize = TRUE, 
                          GX1 = TRUE)
  # saveRDS(result.1, "/Users/ianxu/Library/Mobile Documents/com~apple~CloudDocs/Documents/_BFI Predoc/Pouliot/exactt/tests/testthat/expected_values_exactt_pval_GX1.rds")
  expect_equal(result.1, readRDS(test_path("expected_values_exactt_pval_GX1.rds")))
  
  # Second test (X1)
  result.2 <- exactt_pval(betaNullVec, 
                          Y.temp = Y, 
                          X1.temp = X1, 
                          X2.temp = X2, 
                          nBlocks,
                          permIndices, 
                          GX.indices,
                          studentize = TRUE, 
                          GX1 = FALSE)
  # saveRDS(result.2, "/Users/ianxu/Library/Mobile Documents/com~apple~CloudDocs/Documents/_BFI Predoc/Pouliot/exactt/tests/testthat/expected_values_exactt_pval_X1.rds")
  expect_equal(result.2, readRDS(test_path("expected_values_exactt_pval_X1.rds")))
  
  # Third test (unstudentized)
  result.3 <- exactt_pval(betaNullVec, 
                          Y.temp = Y, 
                          X1.temp = X1, 
                          X2.temp = X2, 
                          nBlocks,
                          permIndices,
                          GX.indices,
                          studentize = FALSE)
  # saveRDS(result.3, "/Users/ianxu/Library/Mobile Documents/com~apple~CloudDocs/Documents/_BFI Predoc/Pouliot/exactt/tests/testthat/expected_values_exactt_pval_unstudentized.rds")
  expect_equal(result.3, readRDS(test_path("expected_values_exactt_pval_unstudentized.rds")))
})

# test_that("getBetaNull() works correctly", {
#   
#   set.seed(31740)
#   n = 50
#   X1 <- matrix(rbinom(n, size = 8, prob = 0.5))
#   X2 <- matrix(rexp(n))
#   eps <- matrix(rnorm(n, sd = 2))
#   Y <- matrix(3*X1 + 2*X2 + eps)
#   data <- data.frame(Y, X1, X2)
#   
#   reg <- ivreg(Y ~ X1 + X2, data = data)
#   
#   nBlocks = 5
#   blockIndexMatrix = matrix(1:n, ncol = nBlocks)
#   blockPermutations <- do.call(rbind, combinat::permn(1:nBlocks))
#   permIndices <- apply(blockPermutations, MARGIN = 1, function (x) c(blockIndexMatrix[, x]))
#   
#   se = reg$cov.unscaled[2,2]
#   precisionToUse = floor(log(se, base = 10)) - 1
#   
#   # First Test
#   result.1 <- getBetaNull(Y.temp = Y, 
#                           X1.temp = X1, 
#                           X2.temp = X2, 
#                           Z.temp = NULL, 
#                           alpha = 0.1, 
#                           nBlocks, 
#                           permIndices, 
#                           beta_hat = reg$coefficients[2], 
#                           se = se,
#                           studentize = TRUE,
#                           precisionToUse = precisionToUse, 
#                           GX1 = TRUE)
# 
#   expect_equal(result.1, 
#                scan(file = test_path("expected_values_getBetaNull_GX1.txt"), 
#                     quiet = TRUE))
#   
#   # Second Test
#   result.2 <- getBetaNull(Y.temp = Y, 
#                           X1.temp = X1, 
#                           X2.temp = X2, 
#                           Z.temp = NULL, 
#                           alpha = 0.1, 
#                           nBlocks, 
#                           permIndices, 
#                           beta_hat = reg$coefficients[2], 
#                           se = se,
#                           studentize = TRUE,
#                           precisionToUse = precisionToUse, 
#                           GX1 = FALSE)
#   
#   expect_equal(result.2, 
#                scan(file = test_path("expected_values_getBetaNull_X1.txt"), 
#                     quiet = TRUE))
#   
#   # Third Test
#   result.3 <- getBetaNull(Y.temp = Y, 
#                           X1.temp = X1, 
#                           X2.temp = X2, 
#                           Z.temp = NULL, 
#                           alpha = 0.1, 
#                           nBlocks, 
#                           permIndices, 
#                           beta_hat = reg$coefficients[2], 
#                           se = se,
#                           studentize = FALSE,
#                           precisionToUse = precisionToUse)
#   
#   expect_equal(result.3, 
#                scan(file = test_path("expected_values_getBetaNull_unstudentized.txt"), 
#                     quiet = TRUE))
# })
