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
