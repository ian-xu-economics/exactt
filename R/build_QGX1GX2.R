#' Build Combined Q Matrix from GX1 and GX2 (Internal Function)
#'
#' Computes the combined projection matrix Q using matrices GX1 and GX2. The matrix
#' GX1 is constructed from X1 and blockIndexMatrix, and then combined with GX2.
#' The resulting matrix, GX1X2, is used to compute Q as:
#' Q = I - GX1X2 * Ginv(GX1X2' * GX1X2) * GX1X2',
#' where Ginv denotes the generalized inverse and I is the identity matrix. This
#' type of computation is used internally for regression analysis or other
#' multivariate techniques where projections are needed.
#'
#' @param X1 A numeric column vector of the primary variable.
#' @param GX2 A numeric matrix, where each column is a block permutation of X2. 
#'        It is expected that GX2 is a non-empty square or rectangular matrix.
#' @param blockIndexMatrix Integer matrix indicating block indices, where each
#'        column represents a block and each element is an index in X.
#'
#' @return Returns the projection matrix Q, which is symmetric and has the same number
#' of rows and columns as the number of rows in GX2.
#'
#' @importFrom MASS ginv
#'
#' @noRd
build_QGX1GX2 <- function(X1, GX2, blockIndexMatrix){
  
  GX1 <- build_GX(X1, blockIndexMatrix)
  
  GX1X2 <- cbind(GX1, GX2)
  
  return(diag(nrow(GX2)) - GX1X2 %*% MASS::ginv(t(GX1X2) %*% GX1X2, tol = max(dim(GX2)) * .Machine$double.eps) %*% t(GX1X2))
  
}
