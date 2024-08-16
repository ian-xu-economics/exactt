exactt_pval_IV_V2 <- function(betaNullVec, Y.temp, X1.temp, X2.temp, Z.temp, nBlocks, permIndices, studentize = TRUE, GX1){
  
  n <- nrow(X1.temp)
  
  nInstruments <- ncol(Z.temp)
  
  betaNullVec <- matrix(betaNullVec, nrow = 1)
  
  blockIndexMatrix <- matrix(1:nrow(Y.temp), ncol = nBlocks)
  
  E <- replicate(length(betaNullVec), Y.temp, simplify = TRUE) - X1.temp%*%betaNullVec
  E.permuted <- apply(E,
                      MARGIN = 2,
                      function(x){
                        matrix(x[permIndices], nrow = n)
                      },
                      simplify = FALSE) %>%
    simplify2array()
  
  if(ncol(X2.temp) > 0){
    
    GX2.temp <- build_GX2(X2.temp, blockIndexMatrix)
    QGX2.temp <- build_QGX2(GX2.temp)
    
    Q.Z.temp <- QGX2.temp%*%Z.temp
    
    TauArray <- apply(E.permuted, 
                      MARGIN = 3, 
                      FUN = function(x){
                        t(Q.Z.temp) %*% x
                      },
                      simplify = FALSE) %>%
      simplify2array() %>%
      aperm(c(2,3,1)) # orients from nInstruments x numPerms x length(betaNullVec) to numPerms x length(betaNullVec) x nInstruments
    
    # Then we need to build combinations of Tau (e.g. T1sq, T1T2, T2sq)
    TauCombosArray <- getHadamardCombinations(TauArray)
    
    if(studentize){
      
      QZComboArray <- getHadamardCombinations(Q.Z.temp)
      
      QGX1GX2.temp <- build_QGX1GX2(X1.temp, GX2.temp, blockIndexMatrix, GX1)
      
      if(!GX1){
        Y.temp.minus.X1.temp.beta0.permuted <- lapply(betaNullVec,
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
    
    if(studentize){
      
      ZComboArray <- getHadamardCombinations(Q.Z.temp)
      
      QGX1GX2.temp <- build_QGX1GX2(X1.temp, GX2.temp = NULL, blockIndexMatrix, GX1)
      
      if(!GX1){
        Y.temp.minus.X1.temp.beta0.permuted <- lapply(betaNullVec,
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
  
  return(list(pval = p_val_seq_beta_index,
              randomizationDist = t(t)))
}
