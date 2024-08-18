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
#' @importFrom combinat permn
#' @noRd
build_GX <- function(blockIndexMatrix){
  
  nBlocks <- ncol(blockIndexMatrix)
  
  # Create s1
  # IX 2024-08-18: need to find quicker way to construct s1. Don't need all columns
  permutation <- expand.grid(rep(list(2:nBlocks), nBlocks-1))
  
  s1 <- rbind(1,
              t(permutation[apply(permutation, MARGIN = 1, function(x) length(unique(x)) == (nBlocks - 1)),]))
  
  # Create s2
  s2.bot <- matrix(NA, nrow = nBlocks - 1, ncol = nBlocks - 1)
  diag(s2.bot) <- 1
  
  if(nBlocks >= 3){
    s2.bot[which(is.na(s2.bot))] <- 3:nBlocks
  }
  
  s2 <- rbind(2, s2.bot)
  
  # Create s3 and so on...
  if(nBlocks >= 3){
    s3.plus <- sapply(3:nBlocks,
                      function(x){
                        matrix(c(x, (1:nBlocks)[-x]))
                      })
  } else{
    s3.plus <- NULL
  }
  
  GX <- unname(cbind(s1, s2, s3.plus))
  
  GX.indices <- apply(GX,
                      MARGIN = 2,
                      function(x) blockIndexMatrix[,x])
  
  return(GX.indices)
}



