# Function: extreme_vals_quadric
# A: an n x n symmetric matrix representing the quadric hypersurface in projective (n-1) dimensions
#    (with homogeneous coordinates: [x; 1])
extremums.quadric <- function(A) {
  n <- nrow(A)
  if(n < 3) {
    stop("The input matrix A must be at least 3x3 to represent a quadric in projective (n-1) dimensions.")
  }
  
  # Partition A into A_aff (top-left (n-1)x(n-1)), b (first n-1 entries of last column), and c (bottom-right scalar)
  A_aff <- A[1:(n-1), 1:(n-1)]
  b <- A[1:(n-1), n]
  c_scalar <- A[n, n]
  
  # Compute the center in affine coordinates: c_aff = -A_aff^{-1} b
  center <- as.numeric(-solve(A_aff, b))
  
  # Compute the constant R = c - b^T A_aff^{-1} b, converting to numeric
  R <- as.numeric(c_scalar - t(b) %*% solve(A_aff, b))
  
  # For an ellipsoid, we expect F to be negative (so that after dividing by -F we get a positive-definite Q).
  if(R >= 0) {
    stop("The given quadric does not represent a bounded ellipsoid (F must be negative).")
  }
  
  # Normalize the quadratic form: Q = -A_aff / R,
  # so that the affine equation becomes: (x - center)^T Q (x - center) = 1.
  Q <- -A_aff / R
  
  # To get extreme values along the x-dimension (first coordinate), compute:
  #   alpha = e1^T Q^{-1} e1, where e1 = (1,0,...,0)^T in R^(n-1)
  Q_inv <- solve(Q)
  
  extremum.values <- sapply(1:ncol(Q_inv),
                            function(x){
                              alpha <- diag(Q_inv)[x]
                              
                              u_star <- as.numeric(Q_inv[,x]) / sqrt(alpha)
                              
                              return(matrix(c(center + u_star, center - u_star), 
                                            ncol = length(center), 
                                            byrow = TRUE))
                            },
                            simplify = FALSE) 
  
  return(do.call('rbind', extremum.values))
}

#' Calculate Confidence Interval Bounds
#'
#' @param omega.g The response vector for which the test is being performed.
#' @param gurobi.params A list of parameters to pass to Gurobi().
#' @param X1.names The column names of the `X1.temp` variable.
#' @param alpha The level of significance.
#'
#' @importFrom gurobi gurobi
exact.wald.projected.confidence.intervals <- function(omega.g, gurobi.params, X1.names, alpha){
  
  dimX1 <- dim(omega.g)[1] - 1
  nPerms <- dim(omega.g)[3]
  
  A.g <- lapply(2:nPerms,
                function(x){
                  omega.g[,,1] - omega.g[,,x]
                }) |>
    simplify2array()
  
  combos <- expand.grid("modelsense" = c("max", "min"),
                        "dimension" = 1:dimX1, 
                        stringsAsFactors = FALSE) |>
    cbind(name = rep(X1.names, each = 2))
  
  # Can change this to pblapply to speed it up.
  extremum.points <- lapply(1:dim(A.g)[3],
                            function(ellipse.index){
                              cbind(ellipse.index, 
                                    combos, 
                                    extremums.quadric(A.g[,,ellipse.index]))
                            })
  
  extremum.points <- do.call('rbind', extremum.points)
  
  orders <- apply(combos, 
                  MARGIN = 1, 
                  function(row) {
                    # Extract the dimension and modelsense values from the row
                    d <- row[["dimension"]]
                    ms <- row[["modelsense"]]
                    # The column name is the same as the dimension number (as a character)
                    colname <- as.character(d)
                    
                    # Filter the dataframe based on modelsense and dimension
                    filtered <- extremum.points[extremum.points$modelsense == ms & extremum.points$dimension == d, ]
                    
                    # Arrange based on modelsense: descending for "max", ascending for "min"
                    if (ms == "max") {
                      return(filtered[order(filtered[[colname]], decreasing = TRUE), ])
                    } else {
                      return(filtered[order(filtered[[colname]]), ])
                    }
                  },
                  simplify = FALSE)
  
  M <- 1e8
  
  gurobi.results <- lapply(orders,
                           function(order){
                     
                             modelsense <- order[1, "modelsense"]
                             dimension <- order[1, "dimension"]
                             
                             start.col <- which(colnames(order) == "1")
                             
                             for(i in floor(nPerms*alpha):nrow(order)){
                               count <- 1
                               if(i > 1){
                                 for(j in 1:(i-1)){
                                   quadratic.vector <- c(order[i, start.col:(start.col + dimX1 - 1)], 1) |>
                                     matrix() |>
                                     as.numeric()
                                   
                                   if(t(quadratic.vector) %*% A.g[,,order$ellipse.index[j]] %*% quadratic.vector <= 0){
                                     count <- count + 1
                                   } 
                                 }
                               }
                               
                               if(count >= floor(nPerms*alpha)){
                                 lower.bound.index <- i
                                 break
                               } else if(i == nrow(order)){
                                 lower.bound.index <- i
                               }
                             }
                             
                             if(lower.bound.index == floor(nPerms*alpha)){
                               return(list(point = as.numeric(order[nPerms*alpha, start.col:ncol(order)])))
                             }
                             
                             gurobi.model <- list()
                             
                             gurobi.model$vtype <- c(rep("C", dimX1),
                                                     rep("B", lower.bound.index))
                             
                             gurobi.model$ub <- c(rep(Inf, dimX1),
                                                  rep(1, lower.bound.index))
                             gurobi.model$lb <- c(rep(-Inf, dimX1),
                                                  rep(0, lower.bound.index))
                             
                             if(modelsense == "max"){
                               gurobi.model$lb[dimension] <- order[lower.bound.index, as.character(dimension)]
                               gurobi.model$ub[dimension] <- order[nPerms*alpha, as.character(dimension)]
                             } else{
                               gurobi.model$ub[dimension] <- order[lower.bound.index, as.character(dimension)]
                               gurobi.model$lb[dimension] <- order[floor(nPerms*alpha), as.character(dimension)]
                             }
                             
                             gurobi.model$A <- c(rep(0, dimX1),
                                                 rep(1, lower.bound.index)) |>
                               matrix(nrow = 1, byrow = TRUE)
                             
                             gurobi.model$sense <- ">="
                             gurobi.model$rhs <- floor(nPerms*alpha) # Reject when greater than alpha. If equal, reject.
                             
                             gurobi.model$quadcon <- lapply(1:lower.bound.index,
                                                            function(i){
                                                              
                                                              x <- A.g[,,order$ellipse.index[i]]
                                                              
                                                              quad.constraint <- list()
                                                              
                                                              quad.constraint$Qc <- Matrix::sparseMatrix(i = rep(1:dimX1, times = dimX1),
                                                                                                         j = rep(1:dimX1, each = dimX1), 
                                                                                                         x = as.vector(x[1:dimX1, 1:dimX1]),
                                                                                                         dims = c(dimX1 + lower.bound.index, dimX1 + lower.bound.index))
                                                              
                                                              quad.constraint$q <- c(2*x[1:dimX1, dimX1 + 1], 
                                                                                     rep(0, lower.bound.index))
                                                              quad.constraint$q[dimX1 + i] <- M
                                                              
                                                              quad.constraint$sense <- "<="
                                                              
                                                              quad.constraint$rhs <- M - x[dimX1 + 1, dimX1 + 1]
                                                              
                                                              return(quad.constraint)
                                                            })
                             
                             # Set the model to maximize
                             gurobi.model$obj <- rep(0, dimX1 + lower.bound.index)
                             gurobi.model$obj[dimension] <- 1
                             gurobi.model$modelsense <- modelsense
                             
                             gurobi.result <- gurobi::gurobi(gurobi.model, params = gurobi.params)
                             
                             return(list(point = gurobi.result$x[1:dimX1],
                                         gurobi.result = gurobi.result))
                           })
          
  conf.int.points <- cbind(combos,
                           sapply(gurobi.results,
                                  function(x){
                                    x$point
                                    }) |>
                             t())
  
  conf.ints <- sapply(1:dimX1,
                      function(i){
                        conf.int.points[(2+(i-1)*2):(1+(i-1)*2),3+i]
                      }) |>
    t()
  
  names(gurobi.results) <- apply(combos[, c("modelsense", "name")], 
                                 MARGIN = 1, 
                                 function(x) paste0(x, collapse = "."))
  
  final.result <- list(confidence.intervals = conf.ints,
                       confidence.interval.points = conf.int.points,
                       gurobi.results = gurobi.results)
  
  return(final.result)
  
}


