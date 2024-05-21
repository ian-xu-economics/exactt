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
build_GX = function(X, blockIndexMatrix){
  
  nBlocks <- ncol(blockIndexMatrix)
  
  max_rank <- nBlocks*(nBlocks - 2) + 2
  
  GX <- X
  rank <- 1
  
  if(nBlocks <= 9){
    possibleBlockPermutations = do.call(rbind, combinat::permn(1:nBlocks))
  
    triedAllBlockPerms <- FALSE
    
    while(rank < max_rank && triedAllBlockPerms == FALSE){
      
      temp <- block_permute(X, blockIndexMatrix, possibleBlockPermutations)
      
      shuffledData <- temp[[1]]
      possibleBlockPermutations <- temp[[2]]
      
      GX <- cbind(GX, shuffledData)
      
      rank <- Matrix::rankMatrix(GX)
      
      if(rank < max_rank && nrow(possibleBlockPermutations) == 0){
        warning(paste0("Unable to construct maximum rank GX. All possible block permutations attempted. \n Current GX rank: ", rank,
                       "\n Max GX rank: ", max_rank))
        
        triedAllBlockPerms <- TRUE
      }
    }
  } else{
    # Need to find a way to switch out of while loop.
    while(rank < max_rank){ 
      for(j in 1:(nBlocks^2)){
        GX <- cbind(GX, block_permute(X, blockIndexMatrix, possibleBlockPermutations = NULL)[[1]])
      }
      rank <- Matrix::rankMatrix(GX)
    }
  }
  return(GX)
}
