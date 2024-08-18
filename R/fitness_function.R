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
#' @param blocks A matrix of indices specifying the blocks.
#' @param blockPermutations A matrix of block permutations used for calculations.
#'
#' @return Returns a numeric value representing the fitness of the permutation.
#' @importFrom MASS ginv
#' @importFrom methods as
#' @noRd
fitness_function <- function(permutation, X1.temp, X2.temp, Z.temp = NULL, blockIndexMatrix, blockPermutations, GX.indices){
  
  GX2 <- build_GX2(X2.temp[permutation,, drop = FALSE], GX.indices)
  
  Q <- build_QGX2(GX2)
  
  if(is.null(Z.temp)){
    gFF <- apply(blockPermutations,
                 MARGIN = 1,
                 function(x) {
                   t(X1.temp) %*% 
                     Matrix::t(as(permutation, "pMatrix")) %*% 
                     Matrix::t(as(c(blockIndexMatrix[,x]), "pMatrix")) %*%
                     Q %*% 
                     as(permutation, "pMatrix") %*% 
                     X1.temp %>%
                     as.numeric()
                 })
  } else{
    gFF <- apply(blockPermutations,
                 MARGIN = 1,
                 function(x) {
                     t(Z.temp) %*% 
                     Matrix::t(as(permutation, "pMatrix")) %*% 
                     Matrix::t(as(c(blockIndexMatrix[,x]), "pMatrix")) %*%
                     Q %*% 
                     as(permutation, "pMatrix") %*% 
                     X1.temp %>%
                     as.numeric() %>%
                     `^`(2) %>%
                     sum()
                 })
  }
  
  return(gFF[1] - mean(gFF[-1]))
}
