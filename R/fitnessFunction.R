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
#' 
#' @importFrom stats median formula model.matrix lm
#' 
#' @noRd
fitness_function <- function(permutation, X1.temp, X2.temp, Z.temp = NULL, blockIndexMatrix, GX.indices, permIndices, blockPermutations){
  
  n <- max(blockIndexMatrix)
  
  GX.indices <- GX.indices[1:n,, drop = FALSE]
  permIndices <- permIndices[1:n,, drop = FALSE]

  X1.temp.permuted <- X1.temp[permutation,, drop = FALSE][1:n,, drop = FALSE]
  X2.temp.permuted <- X2.temp[permutation,, drop = FALSE][1:n,, drop = FALSE]
  
  if(is.null(Z.temp)){
    gFF <- t(X1.temp.permuted) %*% 
      matrix(stats::lm(matrix(X1.temp.permuted[permIndices,], nrow = n) ~ 
                         0 + build_GX(X2.temp.permuted, GX.indices),
                       model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
             nrow = n) |>
      as.numeric()
  } else{
    Z.temp.permuted <- Z.temp[permutation,, drop = FALSE][1:n,, drop = FALSE]
    
    gFF <- t(Z.temp.permuted) %*% 
      matrix(stats::lm(matrix(X1.temp.permuted[permIndices], nrow = n) ~ 
                         0 + build_GX(X2.temp.permuted, GX.indices),
                       model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
             nrow = n) |>
      apply(MARGIN = 2,
            function(x){
              sum(x^2)
            })
  }
  
  return(gFF[1] - mean(gFF[-1]))
}
