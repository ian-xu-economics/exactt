#' Internal Fitness Function for Genetic Algorithm
#'
#' This function computes the fitness of a permutation in the context of a genetic algorithm.
#' It is specifically designed to work with data partitioned into blocks and uses a
#' permutation matrix to evaluate block effects. This function is intended for internal use
#' and is not exported from the package.
#'
#' @param permutation A numeric vector representing the permutation of indices.
#' @param X2.temp The matrix of secondary variables, affected by permutation.
#' @param X1.temp The matrix of the primary variable.
#' @param blockIndexMatrix A matrix of indices specifying the blocks.
#' @param GX.indices A matrix specifying the indices to construct a GX Matrix with attempted maximum rank.
#' @param permIndices A matrix of permutation indices.
#' @param blockPermutations A matrix of indices plugged into blockIndexMatrix resulting in block permutations.
#'
#' @return Returns a numeric value representing the fitness of the permutation.
#' @noRd
fitness_function <- function(permutation, X1.temp, X2.temp, Z.temp = NULL, blockIndexMatrix, GX.indices, permIndices, blockPermutations){
  
  n <- nrow(X1.temp)
  
  GX2 <- build_GX2(X2.temp[permutation,, drop = FALSE], GX.indices)
  Q.GX2 <- build_QGX2(GX2)
  
  X1.temp.permuted <- X1.temp[permutation,, drop = FALSE]
  
  if(is.null(Z.temp)){
    gFF <- t(X1.temp.permuted) %*% 
      Q.GX2 %*% 
      matrix(X1.temp.permuted[permIndices,], nrow = n) |>
      as.numeric()
  } else{
    Z.temp.permuted = Z.temp[permutation,, drop = FALSE]
    
    gFF <- apply(blockPermutations,
                 MARGIN = 1,
                 function(x) {
                     t(Z.temp.permuted) %*% 
                     Q.GX2 %*% 
                     X1.temp.permuted[c(blockIndexMatrix[,x]), drop = FALSE] |>
                     as.numeric() |>
                     (\(y) y^2)() |>
                     sum()
                 })
  }
  
  return(gFF[1] - mean(gFF[-1]))
}
