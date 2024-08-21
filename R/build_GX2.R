#' Construct G Matrix for X2 Using Block Permutations
#'
#' This function constructs a matrix `GX2` where each column represents a block permutation 
#' of the columns of `X2`, which contains secondary regressors. Each column of `X2` is permuted 
#' within the block structure defined by `blocks`, including the identity permutation.
#'
#' @param X2 Matrix of secondary regressors.
#' @param blocks Vector or factor indicating the block structure for each observation in `X2`.
#'
#' @return A matrix where each column is a block permutation of a column in `X2`.
#' @noRd
build_GX2 <- function(X2.temp, GX.indices){
  
  n <- nrow(X2.temp)
  
  GX2.list <- apply(X2.temp,
                    MARGIN = 2,
                    function(x){
                      remove_dependent_columns(matrix(x[GX.indices], nrow = n))
                    },
                    simplify = FALSE)
  
  return(do.call('cbind', GX2.list))
}
