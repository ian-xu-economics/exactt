#' Construct a GX Matrix with Attempted Maximum Rank
#'
#' This function attempts to construct a GX matrix from the input matrix X,
#' using permutations of block indices provided in blockIndexMatrix to achieve
#' a matrix of the highest possible rank. It tries different permutations of
#' the block indices to add new columns to the GX matrix until it either reaches
#' the calculated maximum possible rank or exhausts all permutations. For blocks
#' of size 9 or less, it uses all permutations; for larger blocks, it iterates
#' over a set number of random permutations.
#'
#' @param X Numeric matrix, the original data from which GX is constructed.
#' @param blockIndexMatrix Integer matrix indicating block indices, where each
#'        column represents a block and each element is an index in X.
#'
#' @return Numeric matrix GX, potentially augmented from X with additional columns
#'         derived from permuting blocks to increase its rank.
#'
#' @importFrom Matrix rankMatrix
#' @noRd
build_GX <- function(blockIndexMatrix){
  
  nBlocks <- ncol(blockIndexMatrix)
  
  GX <- generate_block_permutations(nBlocks)
  
  GX.indices <- apply(GX,
                      MARGIN = 2,
                      function(x) blockIndexMatrix[,x])
  
  return(GX.indices)
}


#' Create block permutations from "An Exact t-Test".
#'
#' @param nBlocks The number of blocks to use for block permutations.
#'
#' @return A matrix containing (nBlocks - 1) * (nBlocks - 3) columns from S1,
#' nBlocks - 1 from S2 and nBlocks - 2 from S3+.
generate_block_permutations <- function(nBlocks){
  
  if(nBlocks == 2){
    s1.2 <- matrix(1:2)
    s2.2 <- matrix(2:1)
    
    return(cbind(s1.2, s2.2))
  }
  
  # Create s1
  s1 <- rbind(1, generate_block_permutations(nBlocks - 1) + 1)
  
  # Create s2
  s2.bot <- matrix(NA, nrow = nBlocks - 1, ncol = nBlocks - 1)
  diag(s2.bot) <- 1
  s2.bot[which(is.na(s2.bot))] <- 3:nBlocks
  s2 <- rbind(2, s2.bot)
  
  # Create s3.plus  
  s3.plus <- sapply(3:nBlocks,
                    function(x){
                      matrix(c(x, (1:nBlocks)[-x]))
                    })
  
  return(cbind(s1, s2, s3.plus))
  
}

#' Remove linearly dependent columns
#'
#' @param X A matrix to remove linearly dependent columns (if they exist).
#'
#' @return X matrix with linearly independent columns
remove_dependent_columns <- function(X) {
  qr_decomp <- qr(X, tol = max(dim(X)) * .Machine$double.eps)
  X_independent <- X[, qr_decomp$pivot[1:qr_decomp$rank], drop = FALSE]
  return(X_independent)
}

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
#' @param X1.temp A numeric column vector of the primary variable.
#' @param GX2 A numeric matrix, where each column is a block permutation of X2. 
#'        It is expected that GX2 is a non-empty square or rectangular matrix.
#' @param blockIndexMatrix Integer matrix indicating block indices, where each
#'        column represents a block and each element is an index in X.
#' @param denominator String argument indicating how to calculate epsilon hat.
#' 
#' @return Returns the projection matrix Q, which is symmetric and has the same number
#' of rows and columns as the number of rows in GX2.
#'
#' @importFrom MASS ginv
#'
#' @noRd
build_QGX1GX2 <- function(X1.temp, GX2.temp, GX.indices, denominator){
  
  if(denominator == "GX1"){
    GX1.temp <- apply(X1.temp,
                      MARGIN = 2,
                      function(x){
                        matrix(x[GX.indices], nrow = nrow(GX.indices)) |>
                          remove_dependent_columns()
                      },
                      simplify = FALSE)
    
    GX1.temp <- do.call('cbind', GX1.temp)
    
    GX1X2.temp <- cbind(1, GX1.temp, GX2.temp)
  } else{
    GX1X2.temp <- cbind(1, X1.temp, GX2.temp)
  }
  
  return(diag(nrow(GX2.temp)) - GX1X2.temp %*% MASS::ginv(crossprod(GX1X2.temp), max(dim(GX1X2.temp)) * .Machine$double.eps) %*% t(GX1X2.temp))
  
}

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

