#' @title Exact t-Test P-Value with Instrumental Variables
#' @description Computes the exact p-values for the t-test with instrumental variables, optionally studentized.
#' @param betaNullVec A numeric vector of null hypothesis values for the regression coefficients.
#' @param Y.temp A matrix of the dependent variable values.
#' @param X1.temp A matrix of the independent variable(s) values.
#' @param X2.temp A matrix of additional independent variable(s) values.
#' @param Z.temp A matrix of instrumental variable(s) values.
#' @param nBlocks The number of blocks to use for block permutations.
#' @param permIndices A matrix of permutation indices.
#' @param studentize A logical value indicating whether to studentize the test statistics.
#' @return A list containing the p-values and the randomization distribution.
#' @importFrom combinat permn
#' @importFrom parallel makeCluster detectCores stopCluster clusterExport clusterCall
#' @importFrom doParallel registerDoParallel
#' @importFrom tibble tibble
#' @importFrom dplyr filter
#' @keywords internal
exactt_pval_IV <- function(betaNullVec, Y.temp, X1.temp, X2.temp, Z.temp, nBlocks, permIndices, Q.Z.temp, studentize = TRUE, GX1){
  
  n <- nrow(X1.temp)
  
  nInstruments <- ncol(Z.temp)
  
  betaNullVec <- matrix(betaNullVec, nrow = 1)

  if(ncol(X2.temp) > 0){
    
    eps_hat <- QGX1GX2.temp%*%Y.temp
    
    # Q.Z.temp.permuted is an array with dimension n x M_G x ncol(Z.temp) 
    Q.Z.temp.permuted <- lapply(1:nInstruments, 
                                function(i){
                                  matrix(Q.Z.temp[permIndices, i, drop = FALSE], 
                                         nrow = n,
                                         byrow = FALSE) 
                                }) %>%
      simplify2array() 
    
    E <- replicate(length(betaNullVec), Y.temp, simplify = TRUE) - X1.temp%*%betaNullVec
    
    TauArray <- apply(Q.Z.temp.permuted, 
                      MARGIN = 3, 
                      FUN = function(x){
                        t(x) %*% E
                      },
                      simplify = FALSE) %>%
      simplify2array()
    
    # Then we need to build combinations of Tau (e.g. T1sq, T1T2, T2sq)
    TauCombosArray <- getHadamardCombinations(TauArray)
    
    if(studentize){
      # Then build GQZ combos (e.g. GQZ1sq, GQZ1*GQZ2, GQZ2sq)
      GQZComboArray <- getHadamardCombinations(Q.Z.temp.permuted)
      
      # Then build GSigma. For each slice in GQZComboArray, multiply each
      # column with eps_hat^2 element wise. Then for each column, sum. 
      GSigma <- apply(GQZComboArray, 
                      MARGIN = 3,
                      function(x){
                        apply(x * c(eps_hat^2),
                              MARGIN = 2,
                              sum)
                      })
      
      # Each column is the componets of GSigma (upper triangular).
      # There are nPerms columns.
      GSigmaInvMatrix = apply(GSigma,
                              MARGIN = 1,
                              function(x){
                                GSigmaMatLower <- matrix(0, nrow = nInstruments, ncol = nInstruments)
                                GSigmaMatLower[lower.tri(GSigmaMatLower,diag = TRUE)] <- x
                                
                                diagonalMat <- GSigmaMatLower
                                diagonalMat[row(diagonalMat) != col(diagonalMat)] <- 0
                                
                                GSigmaMat = GSigmaMatLower + t(GSigmaMatLower) - diagonalMat
                                
                                GSigmaMatInv = solve(GSigmaMat)
                                return(GSigmaMatInv[upper.tri(GSigmaMatInv, diag = TRUE)])
                              },
                              simplify = TRUE)
      
      if(!is.matrix(GSigmaInvMatrix)){
        GSigmaInvMatrix <- matrix(GSigmaInvMatrix)
      }
      
      # We know the order that is returned by getHadamardCombinations()
      # We will use this order to help determine which to multiply by two
      # e.g. T1sq*GSigmaInv[,1] + 2T1T2GismgaInv[,2] + T2sq*GSigmaInv[,3]
      trueFalseVectorForDouble = upper.tri(matrix(nrow = nInstruments, ncol = nInstruments))[upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = TRUE)]
      
      t <- lapply(1:dim(TauCombosArray)[3],
                  function(i){
                    if(trueFalseVectorForDouble[i] == TRUE){
                      return(2*TauCombosArray[,,i]*GSigmaInvMatrix[i,])
                    } else{
                      return(TauCombosArray[,,i]*GSigmaInvMatrix[i,])
                    }
                  }) %>%
        Reduce(f = '+') %>%
        matrix(nrow = factorial(nBlocks))
    } else{
      # If we don't studentize, then only need T1sq + T2sq + ...
      TausToAdd = !upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = FALSE)[upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = TRUE)]
      
      t <- apply(TauCombosArray[,,TausToAdd, drop = FALSE],
                 MARGIN = c(1,2),
                 sum)
    }
    
    p_val_seq_beta_index <- apply(t, MARGIN = 2, function(x) mean(abs(x[1]) <= abs(x) + 1e-9))
    
  } else{
    
    QGX1.temp <- build_QGX1GX2(X1.temp, GX2 = NULL, GX.indices)
    
    eps_hat <- QGX1.temp%*%Y.temp
    
    # Z.temp.permuted is an array with dimension n x M_G x ncol(Z.temp) 
    Z.temp.permuted <- lapply(1:nInstruments, 
                              function(i){
                                matrix(Z.temp[permIndices, i, drop = FALSE], 
                                       nrow = n,
                                       byrow = TRUE) 
                              }) %>%
      simplify2array()
    
    E <- replicate(length(betaNullVec), Y.temp, simplify = TRUE) - X1.temp%*%betaNullVec
    
    TauArray <- apply(Z.temp.permuted, 
                      MARGIN = 3, 
                      FUN = function(x){
                        t(x) %*% E
                      },
                      simplify = FALSE) %>%
      simplify2array()
    
    # Then we need to build combinations of Tau (e.g. T1sq, T1T2, T2sq)
    TauCombosArray <- getHadamardCombinations(TauArray)
    
    if(studentize){
      # Then build GQZ combos (e.g. GQZ1sq, GQZ1*GQZ2, GQZ2sq)
      GZComboArray <- getHadamardCombinations(Z.temp.permuted)
      
      # Then build GSigma. For each slice in GQZComboArray, multiply each
      # column with eps_hat^2 element wise. Then for each column, sum. 
      GSigma <- apply(GZComboArray, 
                      MARGIN = 3,
                      function(x){
                        apply(x * c(eps_hat^2),
                              MARGIN = 2,
                              sum)
                      })
      
      # Each column is the componets of GSigma (upper triangular).
      # There are nPerms columns.
      GSigmaInvMatrix = apply(GSigma,
                              MARGIN = 1,
                              function(x){
                                GSigmaMatLower <- matrix(0, nrow = nInstruments, ncol = nInstruments)
                                GSigmaMatLower[lower.tri(GSigmaMatLower,diag = TRUE)] <- x
                                GSigmaMat = GSigmaMatLower + t(GSigmaMatLower) - diag(diag(GSigmaMatLower))
                                
                                GSigmaMatInv = solve(GSigmaMat)
                                return(GSigmaMatInv[upper.tri(GSigmaMatInv, diag = TRUE)])
                              },
                              simplify = TRUE)
      
      # We know the order that is returned by getHadamardCombinations()
      # We will use this order to help determine which to multiply by two
      # e.g. T1sq*GSigmaInv[,1] + 2T1T2GismgaInv[,2] + T2sq*GSigmaInv[,3]
      trueFalseVectorForDouble = upper.tri(matrix(nrow = nInstruments, ncol = nInstruments))[upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = TRUE)]
      
      t <- lapply(1:dim(TauCombosArray)[3],
                  function(i){
                    if(trueFalseVectorForDouble[i] == TRUE){
                      return(2*TauCombosArray[,,i]*GSigmaInvMatrix[i,])
                    } else{
                      return(TauCombosArray[,,i]*GSigmaInvMatrix[i,])
                    }
                  }) %>%
        Reduce(f = '+')
    } else{
      # If we don't studentize, then only need T1sq + T2sq + ...
      TausToAdd = !upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = FALSE)[upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = TRUE)]
      
      t <- apply(TauCombosArray[,,TausToAdd],
                 MARGIN = c(1,2),
                 sum)
    }
    
    p_val_seq_beta_index <- apply(t, MARGIN = 2, function(x) mean(abs(x[1]) <= abs(x) + 1e-9))
  }
  
  return(list(pval = p_val_seq_beta_index,
              randomizationDist = t(t)))
}
