#' Build Q Matrix for GX2
#'
#' Computes the orthogonal projection matrix Q from the matrix GX2. This matrix Q is defined as:
#' Q = I - GX2 * Ginv(GX2' * GX2) * GX2', where Ginv represents the generalized inverse
#' and I is the identity matrix. 
#'
#' @param GX2 A numeric matrix, where each column is a block permutation of X2. 
#' It is expected that GX2 is a non-empty square or rectangular matrix.
#'
#' @return Returns the Q matrix computed as described above. The Q matrix has the same
#' dimensions as the input matrix GX2 and is symmetric.
#'
#' @importFrom MASS ginv
#' 
#' @noRd
build_QGX2 <- function(GX2){
  
  return(diag(nrow(GX2)) - GX2 %*% MASS::ginv(crossprod(GX2), tol = max(dim(GX2)) * .Machine$double.eps) %*% t(GX2))
  
}
