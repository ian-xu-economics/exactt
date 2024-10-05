# Calculate Sigma

compute_sigma_coefficients = function(e.g.Y, e.g.X1, Z.j){
  # Format is a + bx + cx^2
  
  a <- sum(e.g.Y^2 * Z.j)             # same as t(eY) %*% diag(Zj) %*% eY
  b <- -2 * sum(e.g.Y * e.g.X1 * Z.j) # same as -2 * t(eY) %*% diag(Zj) %*% eX1
  c <- sum(e.g.X1^2 * Z.j)            # same as t(eX1) %*% diag(Zj) %*% eX1
  
  return(polynom::polynomial(c(a, b, c)))
}

calculate_determinant <- function(matrix){

  # Base for 1x1
  if(all(dim(matrix) == c(1,1))){
    return(matrix[[1,1]])
  }

  # Base for 2x2 
  if(all(dim(matrix) == c(2,2))){
    sigma11 <- matrix[[1,1]]
    sigma22 <- matrix[[2,2]]
    sigma12 <- matrix[[1,2]]
    
    return(sigma11 * sigma22 - sigma12^2)
  }

  # Base for 3x3
  if(all(dim(matrix) == c(3,3))){
    sigma11 <- matrix[[1,1]]
    sigma12 <- matrix[[1,2]]
    sigma22 <- matrix[[2,2]]
    sigma13 <- matrix[[1,3]]
    sigma23 <- matrix[[2,3]]
    sigma33 <- matrix[[3,3]]

    return(sigma11 * sigma22 * sigma33 + 2 * sigma12 * sigma13 * sigma23 - sigma13^2 * sigma22 - sigma11 * sigma23^2 - sigma12^2 * sigma33)
  }

  # Recursive case
  determinant <- 0
  for (col in 1:dim(matrix)[2]) {
    sub_matrix <- matrix[-1, -col, drop = FALSE]
    cofactor <- (-1)^(1 + col) * matrix[[1,col]] * calculate_determinant(sub_matrix)
    determinant <- determinant + cofactor
  }

  return(determinant)
}

calculate_cofactor_matrix <- function(matrix) {
  
  cofactor_matrix <- sapply(1:nInstruments, function(i){
    sapply(1:nInstruments, function(j){
      minor <- matrix[-i, -j, drop = FALSE]
      list((-1)^(i + j) * calculate_determinant(minor))
    })
  })

  return(cofactor_matrix)
}

Y.X1.temp <- cbind(Y.temp, X1.temp)

# result is nInstruments x nPerms
T.g.Y.X1 <- apply(permIndices,
                  MARGIN = 2,
                  function(x){
                   apply(t(Q.Z.temp) %*% Y.X1.temp[x,],
                         MARGIN = 1,
                         function(y){
                          polynom::polynomial(y)
                         },
                         simplify = FALSE)
                  },
                  simplify = FALSE) |>
  simplify2array()

e.g.Y.X1 <- apply(permIndices,
                  MARGIN = 2,
                  function(x){
                    Q.X1.GX2 %*% Y.X1.temp[x,]
                  },
                  simplify = FALSE) |>
  simplify2array()

Z.temp <- array(Z.temp, dim = c(nrow(Z.temp), 1, ncol(Z.temp)))

Z.hadamard <- getHadamardCombinations(Z.temp)

# Compute sigma coefficients for each combination of e.g.Y and e.g.X1 with each Z.hadamard 
# Resulting array is coefficients (a, b, c) x Z.hadamard slices x permIndices
sigma_coefficients <- apply(e.g.Y.X1,
                            MARGIN = 3,
                            function(x){
                              matLower <- matrix(0, nrow = nInstruments, ncol = nInstruments)
                              polynomials <- apply(Z.hadamard,
                                    MARGIN = 2,
                                    function(Z){
                                      compute_sigma_coefficients(e.g.Y = x[,1], 
                                                                 e.g.X1 = x[,2], 
                                                                 Z)
                                    },
                                    simplify = FALSE)

                              # Initialize matLower as a list to hold polynomials
                              matLowerList <- as.list(matLower)
                              matLowerList[lower.tri(matLower, diag = TRUE)] <- polynomials

                              # Convert the list back to a matrix
                              matLower <- matrix(matLowerList, nrow = nInstruments, ncol = nInstruments)

                              matUpper.noDiag <- matrix(0, nrow = nrow(matLower), ncol = ncol(matLower))
                              matUpper.noDiagList <- as.list(matUpper.noDiag)
                              matUpper.noDiagList[upper.tri(matUpper.noDiag, diag = FALSE)] <- matLower[lower.tri(matLower, diag = FALSE)]
                              matUpper.noDiag <- matrix(matUpper.noDiagList, nrow = nInstruments, ncol = nInstruments)

                              matrix <- add_polynomial_matrices(matLower, matUpper.noDiag)
                            },
                            simplify = FALSE) |>
  simplify2array() 

cofactor_matrices <- apply(sigma_coefficients,
                           MARGIN = 3,
                           function(g){
                             calculate_cofactor_matrix(g)
                           },
                           simplify = FALSE) |>
                           simplify2array()

determinants <- sapply(1:nPerms,
                       function(g){
                        Reduce('+', sapply(1:nInstruments, function(i){
                          sigma_coefficients[[1,i,g]]*cofactor_matrices[[1,i,g]]
                        },
                        simplify = FALSE))
                       },
                       simplify = FALSE) |>
                       matrix()

determinantRoots <- apply(determinants,
MARGIN = 1,
function(x){
  solve(x[[1]])
})

t.g.numerator <- sapply(1:nPerms,
                        function(x){
                          multiply_polynomial_matrices(t(T.g.Y.X1[,x, drop = FALSE]), t(cofactor_matrices[,,x])) |>
                          multiply_polynomial_matrices(T.g.Y.X1[,x, drop = FALSE])
                        }) |>
                        matrix()

final.polynomials <- sapply(2:nPerms,
                           function(x){
                            t.g.numerator[[1,1]]*determinants[[x,1]] - t.g.numerator[[x,1]]*determinants[[1,1]]
                           },
                           simplify = FALSE)

intersections <- lapply(final.polynomials,
                        function(x){
                          solve(x)
                        }) |>
                        simplify2array()

real.intersections <- apply(intersections,
                            MARGIN = 2,
                            function(x){
                              sort(Re(x)[which(abs(Im(x)) < 1e-7)])
                            })

intersection.data.list <- sapply(1:length(real.intersections),
                                 function(x){
                                   intersections <- real.intersections[[x]]

                                   if(length(intersections) == 0){
                                    return(NULL)
                                   }

                                   repeated.roots <- length(unique(intersections)) != length(intersections)

                                   if(repeated.roots){
                                    stop("Repeated root found! Please contact maintainer via email at ianxu@uchicago.edu.")
                                   } else if(length(intersections) > 2 || length(intersections) == 1){
                                     stop("There are not 2 intersections! Please contact maintainer via email at ianxu@uchicago.edu.")
                                   }
                                   
                                   data.frame(
                                      permNum = x + 1,
                                      intersect.left = intersections[1],
                                      intersect.right = intersections[2],
                                      growth.condition = predict(final.polynomials[[x]], mean(intersections)) > 0
                                   ) 
                                 })

intersection.data <- do.call(rbind, intersection.data.list)

pvals.df <- pval.calculator.iv(intersection.data)


pval.calculator.iv <- function(intersection.data){

  intersect.left <- intersection.data$intersect.left
  intersect.right <- intersection.data$intersect.right
  
  beta0 <- c(intersect.left, intersect.right) |>
    sort()
  
  count_matrix <- outer(beta0, intersect.left, `>=`) & outer(beta0, intersect.right, `<=`)
  
  growth_condition <- intersection.data$growth_condition
  
  count_matrix[, growth_condition] <- !count_matrix[, growth_condition]
  
  pvals <- (apply(count_matrix, MARGIN = 1, sum) + 1)/nPerms
  
  pvals.df <- data.frame(beta0.start = c(-Inf, beta0),
                         beta0.end = c(beta0, Inf),
                         pvals = c(sum(growth_condition, 1)/nPerms, 
                                   pvals[-length(pvals)/2], 
                                   sum(growth_condition, 1)/nPerms))
  
  return(pvals.df)
}

intersection.data$discriminant