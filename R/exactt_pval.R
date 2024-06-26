#' Compute p-values for Exact t-tests
#'
#' This function calculates p-values based on permutation tests using exact t statistics.
#' It supports studentized and non-studentized approaches depending on the presence of
#' secondary predictors (X2.temp). The function handles different statistical models
#' based on the dimensionality and characteristics of the input matrices.
#'
#' @param betaNullVec A numeric vector of null hypothesis values to be tested.
#' @param Y.temp The response vector for which the test is being performed.
#' @param X1.temp A numeric column vector of the primary variable.
#' @param X2.temp A matrix of secondary variables.
#' @param nBlocks The number of blocks in the permutation test.
#' @param permIndices A matrix of permutation indices used in the test.
#' @param studentize Logical indicating whether to use studentized test statistics.
#'        Default is TRUE.
#' @param GX1 Logical indicating whether to use GX1 or X1 when constructing eps_hat.
#'
#' @return A list containing:
#'   - pval: A vector of p-values computed for each null hypothesis value.
#'   - randomizationDist: The matrix of randomization distributions used in computing
#'     the p-values.
#'
#' @details
#' The function uses a block matrix approach to generate permutations of the test statistics
#' based on the input matrices and the permutation indices. If `X2.temp` is non-empty,
#' the function builds several matrices including GX2 (generated from `X2.temp` and
#' block indices), QGX2 (projection matrix from GX2), and QGX1GX2 (a combined projection
#' matrix). These matrices are used to calculate test statistics and ultimately p-values.
#' If `X2.temp` is empty, the function simplifies the computations by directly using `X1.temp`.
#'
#' @noRd
exactt_pval <- function(betaNullVec, Y.temp, X1.temp, X2.temp, nBlocks, permIndices, studentize = TRUE, GX1){
  
  n <- nrow(X1.temp)
  
  betaNullVec <- matrix(betaNullVec, nrow = 1)
  
  blockIndexMatrix <- matrix(1:nrow(Y.temp), ncol = nBlocks)
  
  if(ncol(X2.temp) > 0){
    
    GX2.temp <- build_GX2(X2.temp, blockIndexMatrix)
    QGX2.temp <- build_QGX2(GX2.temp)
    QGX1GX2.temp <- build_QGX1GX2(X1.temp, GX2.temp, blockIndexMatrix, GX1)
    
    Q.X1.temp <- QGX2.temp%*%X1.temp
    eps_hat <- QGX1GX2.temp%*%Y.temp
    
    Q.X1.temp.permuted <- matrix(Q.X1.temp[permIndices], nrow = n)
    
    E <- replicate(length(betaNullVec), Y.temp, simplify = TRUE) - X1.temp%*%betaNullVec
    t_num <- t(Q.X1.temp.permuted) %*% E
    
    # n x nBlocks! matrix
    eps_hat.permuted <- matrix(eps_hat[permIndices], nrow = n)
    
    # nBlocks! x 1 matrix
    t_denom <- t(t(Q.X1.temp^2) %*% eps_hat.permuted^2/n)
    
    # t is nPerms x length(betaNullVec)
    # Each column of t is the randomization distribution of the studentized test statistics
    t <- sweep(t_num, 1, t_denom, "/")
    #t <- apply(t_num, 2, function(x) x/t_denom[,1]) OLD incorrect code.
    
    p_val_seq_beta_index <- apply(t, MARGIN = 2, function(x) mean(abs(x[1]) <= abs(x)))
    
  } else{
    X1.temp.permuted <- matrix(X1.temp[permIndices], nrow = n)
    
    E <- replicate(length(betaNullVec), Y.temp, simplify = TRUE) - X1.temp%*%t(betaNullVec)
    
    T_numerator_hat <- t(X1.temp)%*%E #Don't need to multiply by 1/sqrt(n)
    
    t_num <- t(X1.temp.permuted)%*%E #Don't need to multiply by 1/sqrt(n)
    
    p_val_seq_beta_index <- apply(t, MARGIN = 2, function(x) mean(abs(x[1]) <= abs(x)))
  }
  
  return(list(pval = p_val_seq_beta_index,
              randomizationDist = t(t)))
}
