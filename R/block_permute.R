#' Internal Function to Block Permute a Vector
#'
#' Permute the elements of a vector `X` based on block indices specified in `blockIndexMatrix`.
#' This function either randomly permutes the blocks or uses a specified permutation from
#' `possibleBlockPermutations`, updating it by removing the used permutation. This function is
#' intended for internal use and is not exposed to package users.
#'
#' @param X A vector whose elements are grouped into blocks to be permuted.
#' @param blockIndexMatrix A matrix where each column corresponds to indices in `X` that form a block.
#' @param possibleBlockPermutations Optional matrix where each row is a permutation of the block indices.
#'
#' @return A list containing:
#' \itemize{
#'   \item{shuffledData}{The permuted vector `X`.}
#'   \item{possibleBlockPermutations}{The updated matrix of possible block permutations, with the used permutation removed.}
#' }
#'
#' @noRd
block_permute <- function(X, blockIndexMatrix, possibleBlockPermutations = NULL){
  
  nBlocks <- ncol(blockIndexMatrix)
  
  if(is.null(possibleBlockPermutations) == TRUE){
    block_order <- sample(1:nBlocks, nBlocks)
  } else{
    possibleBlockPermutationsLeft <- nrow(possibleBlockPermutations)
    
    row <- sample(possibleBlockPermutationsLeft, 1)
    
    block_order <- possibleBlockPermutations[row, ]
    
    possibleBlockPermutations <- possibleBlockPermutations[-row,, drop = FALSE]
  }
  
  return(list("shuffledData" = X[c(blockIndexMatrix[, block_order])],
              "possibleBlockPermutations" = possibleBlockPermutations))
  
}
