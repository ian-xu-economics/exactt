#' Internal Fitness Function for Genetic Algorithm
#'
#' This function computes the fitness of a permutation in the context of a genetic algorithm.
#' It is specifically designed to work with data partitioned into blocks and uses a
#' permutation matrix to evaluate block effects. This function is intended for internal use
#' and is not exported from the package.
#'
#' @param permutation A numeric vector representing the permutation of indices.
#' @param X1.temp The matrix of the primary variable.
#' @param X2.temp The matrix of non-fixed effect secondary variables.
#' @param blockIndexMatrix A matrix of indices specifying the blocks.
#' @param GX.indices A matrix specifying the indices to construct a GX Matrix with attempted maximum rank.
#' @param permIndices A matrix of permutation indices.
#'
#' @return Returns a numeric value representing the fitness of the permutation.
#' 
#' @importFrom stats median formula model.matrix lm
#' 
#' @noRd
fitness_function <- function(permutation, X1.temp, X2.temp, Z.temp = NULL, blockIndexMatrix, GX.indices, permIndices){
  
  n <- max(blockIndexMatrix)
  
  GX.indices.use <- GX.indices[1:n,, drop = FALSE]
  permIndices.use <- permIndices[1:n,, drop = FALSE]

  X1.temp.permuted <- X1.temp[permutation,, drop = FALSE][1:n,, drop = FALSE]
  # Don't want to store X2.temp.permuted for RAM
  # X2.temp.permuted <- X2.temp[permutation,, drop = FALSE][1:n,, drop = FALSE]
  
  if(is.null(Z.temp) && ncol(X1.temp) == 1){
    gFF <- stats::lm(X1.temp.permuted ~ 
                       0 + build_GX(X2.temp[permutation,, drop = FALSE][1:n,, drop = FALSE], 
                                    GX.indices.use),
                            model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals |>
      matrix(nrow = 1) %*%
      matrix(X1.temp.permuted[permIndices.use,], nrow = n) |>
      as.numeric()

  } else{
    #if(!is.null(Z.temp)){
      # Don't want to store Z.temp.permuted for RAM
      # Z.temp.permuted <- Z.temp[permutation,, drop = FALSE][1:n,, drop = FALSE]
      
      gFF <- lm(Z.temp[permutation,, drop = FALSE][1:n,, drop = FALSE] ~ 
                  0 + build_GX(X2.temp[permutation,, drop = FALSE][1:n,, drop = FALSE], 
                               GX.indices.use),
                model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals |>
        matrix(ncol = ncol(Z.temp)) |>
        t() %*%
        matrix(X1.temp.permuted[permIndices.use,], nrow = n) |>
        apply(MARGIN = 2,
              function(x){
                sum(x^2)
              })
    #} else{
      #gFF <- lm(X1.temp[permutation,, drop = FALSE][1:n,, drop = FALSE] ~ 
      #            0 + build_GX(X2.temp[permutation,, drop = FALSE][1:n,, drop = FALSE], 
      #                         GX.indices.use),
      #          model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals |>
      #  matrix(ncol = ncol(X1.temp)) |>
      #  t() %*%
      #  matrix(X1.temp.permuted[permIndices.use,], nrow = n) |>
      #  apply(MARGIN = 2,
      #        function(x){
      #          sum(x^2)
      #        })
      
      
      # conics <- lm(X1.temp[permutation,, drop = FALSE][1:n,, drop = FALSE] ~ 
      #      0 + build_GX(X2.temp[permutation,, drop = FALSE][1:n,, drop = FALSE], 
      #                   GX.indices.use),
      #    model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals |>
      #   matrix(ncol = ncol(X1.temp)) |>
      #   t() %*%
      #   matrix(X1.temp.permuted[permIndices.use,], nrow = n) |>
      #   array(dim = c(ncol(X1.temp), ncol(X1.temp), 120)) |>
      #   apply(MARGIN = 3,
      #         function(x){
      #           t(x) %*% x
      #         },
      #         simplify = FALSE) |>
      #   simplify2array()
      # 
      # return(apply(conics[,,-1],
      #              MARGIN = 3,
      #              function(x){
      #                eigen(conics[,,1] - x)$values > 0
      #              }) |>
      #          sum())
    
    #}
  }
  
  return(gFF[1] - mean(gFF[-1]))
}
