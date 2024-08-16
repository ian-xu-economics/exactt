#' Get Hadamard Combinations of Tensor Slices
#'
#' This function computes the Hadamard (element-wise) product for all unique combinations of 
#' the slices of a given 3-dimensional tensor. The order of the resulting combinations is 
#' such that combinations are treated as identical if they differ only by the order of slices.
#'
#' @param tensor A 3-dimensional array where the third dimension represents different slices.
#'
#' @return A 3-dimensional array where each slice represents the Hadamard product of a unique 
#' combination of the input tensor's slices. The resulting order is 11, 12, 22, 13, 23, 33, etc.
#'
#' @noRd
getHadamardCombinations <- function(tensor){
  
  num_dims <- length(dim(tensor))
  
  numSlices <- dim(tensor)[num_dims]
  n <- dim(tensor)[1]
  
  # Generate combinations of matrix indices
  allCombos = as.matrix(expand.grid(1:numSlices, 1:numSlices))
  
  # We treat 12 and 21 the same, so we only include those where the first column is <= the second
  allCombos = allCombos[allCombos[,1] <= allCombos[,2],, drop = FALSE]
  
  # Compute Hadamard products for each combination
  hadamardProducts <- apply(allCombos,
                            MARGIN = 1,
                            function(indices){
                              lastInd(tensor, num_dims, indices[1]) * lastInd(tensor, num_dims,  indices[2])
                            },
                            simplify = FALSE) %>%
    simplify2array()
  
  # Order is 11, 12, 22, 13, 23, 33, 14, 24, 34, 44,...
  return(hadamardProducts)
  
}


lastInd <- function(tensor, num_dims, index){
  
  # Create index lists for the two subsetting operations
  index_list <- rep(list(TRUE), num_dims)
  
  # Replace the last dimension with the desired indices
  index_list[[num_dims]] <- index
  
  return(do.call(`[`, c(list(tensor), index_list)))
}
