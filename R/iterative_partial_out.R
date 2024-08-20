iterative_partial_out <- function(partial.out, GX) {
  
  n <- nrow(GX)
  k <- ncol(GX)
  
  # Start with the original vector to partial out
  residual <- partial.out
  
  # Create a matrix to hold the orthogonalized columns of GX
  Q <- matrix(0, nrow = n, ncol = k)
  
  # Iteratively orthogonalize and project out each column of GX
  for (j in 1:k) {
    
    # Start with the j-th column of GX
    v <- GX[, j]
    
    # Orthogonalize v against all previous orthogonal vectors in Q
    if (j > 1) {
      
      test <- apply(Q[,1:(j-1), drop = FALSE],
                    MARGIN = 2,
                    function(x){
                      if(sum(x^2) > 1e-10){
                        sum(x*GX[,j]) / sum(x^2) * x
                      } else{
                        return(rep(0, n))
                      }
                    })
      
      v <- v - apply(test, MARGIN = 1, sum)
      
      #for (i in 1:(j-1)) {
      #  if (sum(Q[, i]^2) > 1e-10) {  # Check if Q[, i] is non-zero
      #    v <- v - as.vector(crossprod(Q[, i], GX[, j]) / crossprod(Q[, i])) * Q[, i]
      #  }
      #}
    }

    # Check if the vector v is non-zero (i.e., check for linear dependence)
    v_norm <- sqrt(sum(v^2))
    if (v_norm > 1e-10) {  # Use a small threshold to avoid numerical issues
      # Normalize the orthogonalized vector
      Q[, j] <- v / v_norm
      
      # Project residual onto the orthogonalized vector Q[, j]
      proj_j <- outer(Q[, j], Q[, j]) / as.vector(crossprod(Q[, j]))
      
      # Update residual by subtracting the projection
      residual <- residual - proj_j %*% residual
    }
    
  }
  
  # The residual is now partial.out with the effects of GX partialed out using orthogonalized vectors
  return(residual)
}
