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
  
  numSlices <- dim(tensor)[3]
  n <- dim(tensor)[1]
  
  print(numSlices)
  
  # Generate combinations of matrix indices
  allCombos = as.matrix(expand.grid(1:numSlices, 1:numSlices))
  
  # We treat 12 and 21 the same, so we only include those where the first column is <= the second
  allCombos = allCombos[allCombos[,1] <= allCombos[,2],, drop = FALSE]
  
  # Compute Hadamard products for each combination
  hadamardProducts <- apply(allCombos,
                            MARGIN = 1,
                            function(indices){
                              matrix(tensor[,,indices[1]] * tensor[,,indices[2]],
                                     nrow = n)
                            },
                            simplify = FALSE) %>%
    simplify2array()
  
  # Order is 11, 12, 22, 13, 23, 33, 14, 24, 34, 44,...
  return(hadamardProducts)
  
}
