#' Internal Fitness Function for Genetic Algorithm
#'
#' This function computes the fitness of a permutation in the context of a genetic algorithm.
#' It is specifically designed to work with data partitioned into blocks and uses a
#' permutation matrix to evaluate block effects. This function is intended for internal use
#' and is not exported from the package.
#'
#' @param permutation A numeric vector representing the permutation of indices.
#' @param X1.temp The matrix of the primary variable.
#' @param X2.temp The matrix of non-fixed effect secondary variables.
#' @param Z.temp The matrix of instrumental variables.
#' @param indep.X2.index The indices of variables in X2 that are independent to X1.
#' @param blockIndexMatrix A matrix of indices specifying the blocks.
#' @param GX.indices A matrix specifying the indices to construct a GX Matrix with attempted maximum rank.
#' @param permIndices A matrix of permutation indices.
#'
#' @return Returns a numeric value representing the fitness of the permutation.
#' 
#' @importFrom stats median formula model.matrix lm
#' 
#' @noRd
fitness_function <- function(permutation, X1.temp, X2.temp, Z.temp = NULL, indep.X2.index, blockIndexMatrix, GX.indices, permIndices){
  
  n <- max(blockIndexMatrix)
  
  GX.indices.use <- GX.indices[1:n,, drop = FALSE]
  permIndices.use <- permIndices[1:n,, drop = FALSE]

  X1.temp.permuted <- X1.temp[permutation,, drop = FALSE][1:n,, drop = FALSE]
  
  if(is.null(Z.temp)){
    Z.temp.permuted <- X1.temp.permuted
  } else{
    Z.temp.permuted <- Z.temp[permutation,, drop = FALSE][1:n,, drop = FALSE]
  }
  
  if(ncol(X1.temp.permuted) == 1){
    if(ncol(X2.temp) == 0){
      gFF <- matrix(Z.temp.permuted, nrow = 1) %*%
        matrix(X1.temp.permuted[permIndices.use,], nrow = n)
    } else{
      gFF <- stats::lm(Z.temp.permuted ~ 0 + 
                         build_GX(X2.temp[permutation,, drop = FALSE][1:n,, drop = FALSE], 
                                  GX.indices.use,
                                  indep.X2.index),
                       model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals |>
        matrix(ncol = nrow(Z.temp.permuted), byrow = TRUE) %*%
        matrix(X1.temp.permuted[permIndices.use,], nrow = n)
    }
    
    if(nrow(gFF) != 1){
      gFF <- apply(gFF,
                   MARGIN = 2,
                   function(x){
                     crossprod(x)
                   })
    } else{
      gFF <- as.numeric(gFF)
    }
    
    return(gFF[1] - mean(gFF[-1]))
  } else{
    if(ncol(X2.temp) == 0){
      gFF <- apply(permIndices.use,
                   MARGIN = 2,
                   FUN = function(x){
                     t(Z.temp.permuted) %*%
                       X1.temp.permuted[x,, drop = FALSE]
                   },
                   simplify = FALSE) |>
        simplify2array()
    } else{
      gFF <- apply(permIndices.use,
                   MARGIN = 2,
                   FUN = function(x){
                     stats::lm(Z.temp.permuted ~ 0 + 
                                 build_GX(X2.temp[permutation,, drop = FALSE][1:n,, drop = FALSE], 
                                          GX.indices.use,
                                          indep.X2.index),
                               model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals |>
                       t() %*%
                       X1.temp.permuted[x,, drop = FALSE]
                   },
                   simplify = FALSE) |>
        simplify2array()
    }
    
    gFF.eigenvalues <- apply(gFF,
                             MARGIN = 3,
                             function(x){
                               eigen(crossprod(x))$values
                             }) 
    
    if(!is.matrix(gFF.eigenvalues)){
      gFF.eigenvalues <- matrix(gFF.eigenvalues, nrow = 1)
    }
    
    # L2 norm of [Xi - mean(Xi)] across i
    return(sqrt(sum((gFF.eigenvalues[,1] - apply(gFF.eigenvalues, MARGIN = 1, mean))^2)))
  }
}
