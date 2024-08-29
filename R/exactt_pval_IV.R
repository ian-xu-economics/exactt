#' @title Exact t-Test P-Value with Instrumental Variables
#' @description Computes the exact p-values for the t-test with instrumental variables, optionally studentized.
#' @param betaNullVec A numeric vector of null hypothesis values for the regression coefficients.
#' @param Y.temp A matrix of the dependent variable values.
#' @param X1.temp A matrix of the independent variable(s) values.
#' @param X2.temp A matrix of additional independent variable(s) values.
#' @param Z.temp A matrix of instrumental variable(s) values.
#' @param nBlocks The number of blocks to use for block permutations.
#' @param permIndices A matrix of permutation indices.
#' @return A list containing the p-values and the randomization distribution.
#' @importFrom combinat permn
#' @importFrom parallel makeCluster detectCores stopCluster clusterExport clusterCall
#' @importFrom doParallel registerDoParallel
#' @importFrom tibble tibble
#' @importFrom dplyr filter
#' @keywords internal
exactt_pval_IV <- function(beta0.df, Y.temp, X1.temp, X2.temp, Z.temp, nBlocks, permIndices, Q.Z.temp, QGX1GX2.temp, GX1){
  
  n <- nrow(X1.temp)
  nInstruments <- ncol(Z.temp)
  
  beta0.pval_na_indices <- which(is.na(beta0.df$beta0.pval))
  beta0.vec <- matrix(beta0.df$beta0[beta0.pval_na_indices], nrow = 1)
  
  blockIndexMatrix <- matrix(1:nrow(Y.temp), ncol = nBlocks)
  
  E <- replicate(length(beta0.vec), Y.temp, simplify = TRUE) - X1.temp%*%beta0.vec
  E.permuted <- apply(E,
                      MARGIN = 2,
                      function(x){
                        matrix(x[permIndices], nrow = n)
                      },
                      simplify = FALSE) %>%
    simplify2array()
  
  if(ncol(X2.temp) > 0){
    
    TauArray <- apply(E.permuted, 
                      MARGIN = 3, 
                      FUN = function(x){
                        t(Q.Z.temp) %*% x
                      },
                      simplify = FALSE) %>%
      simplify2array() %>%
      aperm(c(2,3,1)) # orients from nInstruments x numPerms x length(beta0.vec) to numPerms x length(beta0.vec) x nInstruments
    
    # Then we need to build combinations of Tau (e.g. T1sq, T1T2, T2sq)
    TauCombosArray <- getHadamardCombinations(TauArray)
    
    if(!is.null(QGX1GX2.temp)){ # If studentize
      
      QZComboArray <- getHadamardCombinations(Q.Z.temp)
      
      if(!GX1){
        
        eps_hat.permuted <- apply(E.permuted,
                                  MARGIN = 3,
                                  function(x){
                                    QGX1GX2.temp %*% x
                                  },
                                  simplify = FALSE) %>%
          simplify2array()
        
        GSigma <- apply(eps_hat.permuted,
                        MARGIN = 3,
                        function(x){
                          t(QZComboArray) %*% x^2
                        },
                        simplify = FALSE) %>%
          simplify2array() %>%
          aperm(c(2,1,3))
        
        GSigmaInvMatrix = apply(GSigma,
                                MARGIN = 3,
                                function(slice){
                                  
                                  apply(slice,
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
                                  
                                },
                                simplify = FALSE) %>%
          simplify2array()
        
        # We know the order that is returned by getHadamardCombinations()
        # We will use this order to help determine which to multiply by two
        # e.g. T1sq*GSigmaInv[,1] + 2T1T2GismgaInv[,2] + T2sq*GSigmaInv[,3]
        trueFalseVectorForDouble = upper.tri(matrix(nrow = nInstruments, ncol = nInstruments))[upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = TRUE)]
        
        t <- lapply(1:dim(TauCombosArray)[3],
                    function(i){
                      if(trueFalseVectorForDouble[i] == TRUE){
                        return(2*TauCombosArray[,,i]*GSigmaInvMatrix[i,,])
                      } else{
                        return(TauCombosArray[,,i]*GSigmaInvMatrix[i,,])
                      }
                    }) %>%
          Reduce(f = '+') %>%
          matrix(nrow = factorial(nBlocks))
      } else{ # If using GX1
        Y.temp.permuted <- matrix(Y.temp[permIndices], nrow = n)
        eps_hat.permuted <- QGX1GX2.temp %*% Y.temp.permuted
        
        GSigma <- t(QZComboArray) %*% eps_hat.permuted
        
        GSigmaInvMatrix = apply(GSigma,
                                MARGIN = 2,
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
      }
    } else{
      # If we don't studentize, then only need T1sq + T2sq + ...
      TausToAdd = !upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = FALSE)[upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = TRUE)]
      
      t <- apply(TauCombosArray[,,TausToAdd, drop = FALSE],
                 MARGIN = c(1,2),
                 sum)
    }
  } else{ # Like above, but don't partial out GX2
    
    TauArray <- apply(E.permuted, 
                      MARGIN = 3, 
                      FUN = function(x){
                        t(Z.temp) %*% x
                      },
                      simplify = FALSE) %>%
      simplify2array() %>%
      aperm(c(2,3,1)) 
    
    # Then we need to build combinations of Tau (e.g. T1sq, T1T2, T2sq)
    TauCombosArray <- getHadamardCombinations(TauArray)
    
    if(!is.null(QGX1GX2.temp)){
      
      ZComboArray <- getHadamardCombinations(Q.Z.temp)
      
      if(!GX1){
        Y.temp.minus.X1.temp.beta0.permuted <- lapply(beta0.vec,
                                                      function(x){
                                                        matrix((Y.temp - X1.temp*x)[permIndices], nrow = n)
                                                      }) %>%
          simplify2array()  
        
        eps_hat.permuted <- apply(Y.temp.minus.X1.temp.beta0.permuted,
                                  MARGIN = 3,
                                  function(x){
                                    QGX1GX2.temp %*% x
                                  },
                                  simplify = FALSE) %>%
          simplify2array()
        
        GSigma <- apply(eps_hat.permuted,
                        MARGIN = 3,
                        function(x){
                          t(ZComboArray) %*% x^2
                        },
                        simplify = FALSE) %>%
          simplify2array() %>%
          aperm(c(2,1,3))
        
        GSigmaInvMatrix = apply(GSigma,
                                MARGIN = 3,
                                function(slice){
                                  
                                  apply(slice,
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
                                  
                                },
                                simplify = FALSE) %>%
          simplify2array()
        
        # We know the order that is returned by getHadamardCombinations()
        # We will use this order to help determine which to multiply by two
        # e.g. T1sq*GSigmaInv[,1] + 2T1T2GismgaInv[,2] + T2sq*GSigmaInv[,3]
        trueFalseVectorForDouble = upper.tri(matrix(nrow = nInstruments, ncol = nInstruments))[upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = TRUE)]
        
        t <- lapply(1:dim(TauCombosArray)[3],
                    function(i){
                      if(trueFalseVectorForDouble[i] == TRUE){
                        return(2*TauCombosArray[,,i]*GSigmaInvMatrix[i,,])
                      } else{
                        return(TauCombosArray[,,i]*GSigmaInvMatrix[i,,])
                      }
                    }) %>%
          Reduce(f = '+') %>%
          matrix(nrow = factorial(nBlocks))
      } else{ # If using GX1
        Y.temp.permuted <- matrix(Y.temp[permIndices], nrow = n)
        eps_hat.permuted <- QGX1GX2.temp %*% Y.temp.permuted
        
        GSigma <- t(ZComboArray) %*% eps_hat.permuted
        
        GSigmaInvMatrix = apply(GSigma,
                                MARGIN = 2,
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
      }
    } else{
      # If we don't studentize, then only need T1sq + T2sq + ...
      TausToAdd = !upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = FALSE)[upper.tri(matrix(nrow = nInstruments, ncol = nInstruments), diag = TRUE)]
      
      t <- apply(TauCombosArray[,,TausToAdd, drop = FALSE],
                 MARGIN = c(1,2),
                 sum)
    }
  }
  
  p_val_seq_beta_index <- apply(t, MARGIN = 2, function(x) mean(abs(x[1]) <= abs(x) + 1e-9))
  
  beta0.df$beta0.pval[beta0.pval_na_indices] <- p_val_seq_beta_index
  
  return(list(beta0.df = beta0.df,
              randomizationDist = t(t)))
}
