#' Build Combined Q Matrix from GX1 and GX2 (Internal Function)
#'
#' Computes the combined projection matrix Q using matrices GX1 and GX2. The matrix
#' GX1 is constructed from X1 and blockIndexMatrix, and then combined with GX2.
#' The resulting matrix, GX1X2, is used to compute Q as:
#' \dqen{Q = I - GX1X2 * Ginv(GX1X2' * GX1X2) * GX1X2'},
#' where Ginv denotes the generalized inverse and I is the identity matrix. This
#' type of computation is used internally for regression analysis or other
#' multivariate techniques where projections are needed.
#'
#' @param X1 A numeric column vector of the primary variable.
#' @param GX2 A numeric matrix, where each column is a block permutation of X2. 
#'        It is expected that GX2 is a non-empty square or rectangular matrix.
#' @param blockIndexMatrix Integer matrix indicating block indices, where each
#'        column represents a block and each element is an index in X.
#' @param GX1 Logical indicating whether to use GX1 or X1 when constructing eps_hat.
#' 
#' @return Returns the projection matrix Q, which is symmetric and has the same number
#' of rows and columns as the number of rows in GX2.
#'
#' @importFrom MASS ginv
#'
#' @noRd
build_QGX1GX2 <- function(X1, GX2, GX.indices, GX1 = TRUE){
  
  if(GX1){
    GX1 <- matrix(X1[GX.indices], nrow = nrow(GX.indices))
    GX1X2 <- cbind(1, GX1, GX2)
  } else{
    GX1X2 <- cbind(1, X1, GX2)
  }
  
  return(diag(nrow(GX2)) - GX1X2 %*% MASS::ginv(t(GX1X2) %*% GX1X2, max(dim(GX1X2)) * .Machine$double.eps) %*% t(GX1X2))
  
}
