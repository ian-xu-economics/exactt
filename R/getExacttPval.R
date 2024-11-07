#' Calculate p-values for two sided case.
#'
#' @param Y.temp The response vector for which the test is being performed.
#' @param X1.temp A numeric column vector of the primary variable.
#' @param X2.temp A numeric matrix of the secondary variables.
#' @param permIndices A matrix of permutation indices used in the test.
#' @param GX.indices A matrix of permutation indices to create GX matrices.
#' @param Q.X1.temp A numeric column vector of the primary variable annihilated by GX2.
#' @param studentize A boolean indicating whether to studentize the randomization statistics
#' @param side A character to indicate the side of the test.
#' @param denominator Character argument indicating how to calculate epsilon hat.
#'
#' @importFrom polynom polynomial
#' @importFrom stats predict coefficients lm
#' @importFrom cli cli_abort
exactt.pval.new.reg <- function(Y.temp, X1.temp, X2.temp, permIndices, GX.indices, Q.X1.temp, studentize, side, denominator){
  
  n <- nrow(Y.temp)
  
  if(denominator == "GX1" || studentize == FALSE){
    # We no longer store X1.temp.permuted, Y.temp.permuted, or eps_hat.permuted; it is a RAM nightmare. 
    if(studentize == TRUE){
      # 1 x nPerms matrix
      sigma.hat <- sqrt(t(Q.X1.temp^2) %*% 
                          matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                             build_GX(X1.temp, GX.indices) + build_GX(X2.temp, GX.indices),
                                           model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                 ncol = ncol(permIndices))^2)
    } else{
      sigma.hat <- 1
    }
  
    Q.X1.temp.dot.X1.temp <- t(Q.X1.temp) %*% matrix(X1.temp[permIndices,], ncol = ncol(permIndices))
    Q.X1.temp.dot.Y.temp <- t(Q.X1.temp) %*% matrix(Y.temp[permIndices], ncol = ncol(permIndices))
    
    lineParams <- data.frame(m = -c(Q.X1.temp.dot.X1.temp/sigma.hat), 
                             b = c(Q.X1.temp.dot.Y.temp/sigma.hat))
    
    m.identity <- lineParams$m[1]
    b.identity <- lineParams$b[1]
    
    if(side == "both"){
      intersect.data <- cbind((lineParams$b[-1] - b.identity)/(m.identity - lineParams$m[-1]), 
                              (-b.identity - lineParams$b[-1])/(m.identity + lineParams$m[-1])) |> 
        apply(MARGIN = 1, function(x){ sort(x, na.last = TRUE) }) |> 
        t() |>
        data.frame()
      
      names(intersect.data) <- c("intersectLeft", "intersectRight")
      
      intersect.data <- cbind(slope = lineParams$m[-1], 
                              intersect = lineParams$b[-1], 
                              intersect.data)
    } else{
      intersect.data <- data.frame(slope = lineParams$m[-1],
                                   yintercept = lineParams$b[-1],
                                   intersections = (lineParams$b[-1] - b.identity)/(m.identity - lineParams$m[-1]))
    }
    
    pvals.df <- pvalCalculator(intersect.data, check.identity = m.identity, intercept = b.identity, iv = FALSE, side = side)
  } else{ # Denominator = X1
    
    Q.X1.GX2.dot.Y.temp.permuted <- matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                                       X1.temp + build_GX(X2.temp, GX.indices),
                                                     model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                           ncol = ncol(permIndices))
    
    Q.X1.GX2.dot.X1.temp.permuted <- matrix(stats::lm(matrix(X1.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                                        X1.temp + build_GX(X2.temp, GX.indices),
                                                      model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                            ncol = ncol(permIndices))
    
    sigma.hat.sq.polynomials <- lapply(1:ncol(permIndices),
                                       function(x){
                                         c(t(Q.X1.GX2.dot.Y.temp.permuted[,x]) %*%
                                             diag(c(Q.X1.temp^2)) %*%
                                             Q.X1.GX2.dot.Y.temp.permuted[,x],
                                           -2*Q.X1.GX2.dot.X1.temp.permuted[,x] %*%
                                             diag(c(Q.X1.temp^2)) %*%
                                             Q.X1.GX2.dot.Y.temp.permuted[,x],
                                           t(Q.X1.GX2.dot.X1.temp.permuted[,x]) %*%
                                             diag(c(Q.X1.temp^2)) %*%
                                             Q.X1.GX2.dot.X1.temp.permuted[,x]) |>
                                           polynom::polynomial()
                                       })

    t.sq.polynomials <- apply(permIndices,
                              MARGIN = 2,
                              function(x){
                                c(t(Y.temp[x,]) %*%
                                    Q.X1.temp %*%
                                    t(Q.X1.temp) %*%
                                    Y.temp[x,],
                                  -2*t(X1.temp[x,]) %*%
                                    Q.X1.temp %*%
                                    t(Q.X1.temp) %*%
                                    Y.temp[x,],
                                  t(X1.temp[x,]) %*%
                                    Q.X1.temp %*%
                                    t(Q.X1.temp) %*%
                                    X1.temp[x,]
                                ) |>
                                  polynom::polynomial()
                              },
                              simplify = FALSE)

    roots <- sapply(2:ncol(permIndices),
                    function(x){
                      polyroot(stats::coefficients(t.sq.polynomials[[1]] *
                                                   sigma.hat.sq.polynomials[[x]] -
                                                   t.sq.polynomials[[x]] *
                                                   sigma.hat.sq.polynomials[[1]]))
                    })

    real.roots <- apply(roots,
                        MARGIN = 2,
                        function(x){
                          sort(Re(x)[!abs(Im(x)) > 1e-5])
                        },
                        simplify = FALSE)
    
    # Need to check if the roots are valid in the original problem. Especially if we are dealing with multiple sides.
    # Issue is that we don't multiply by both sides by least common multiple of the denominators.
    # Doing this would be tricky because we don't know what the new polynomial after dividing LCM by sigma.hat.sq.polynomials[[x]]
  
    # This code focuses on the two sided case.
    
    intersect.data.list <- sapply(1:length(real.roots),
                                  function(x){
                                    real.roots.temp <- real.roots[[x]]
                                   
                                    values.at.roots.test <- stats::predict(t.sq.polynomials[[1]], real.roots.temp) / stats::predict(sigma.hat.sq.polynomials[[1]], real.roots.temp)
                                    values.at.roots.rand <- stats::predict(t.sq.polynomials[[1+x]], real.roots.temp) / stats::predict(sigma.hat.sq.polynomials[[1+x]], real.roots.temp)
              
                                    valid.real.roots <- real.roots.temp[abs(values.at.roots.test - values.at.roots.rand) < 1e-5]
                                    
                                    test.values <- c(valid.real.roots[1] - 1,
                                                     (valid.real.roots[-1] + valid.real.roots[-length(valid.real.roots)])/2,
                                                     valid.real.roots[length(valid.real.roots)] + 1)
                                    
                                    values.at.test.vals.test <- stats::predict(t.sq.polynomials[[1]], test.values) / stats::predict(sigma.hat.sq.polynomials[[1]], test.values)
                                    values.at.test.vals.rand <- stats::predict(t.sq.polynomials[[1+x]], test.values) / stats::predict(sigma.hat.sq.polynomials[[1+x]], test.values)
                                    
                                    return(data.frame(beta0.start = c(-Inf, valid.real.roots),
                                                      beta0.end = c(valid.real.roots, Inf),
                                                      test.stat.smaller = values.at.test.vals.test < values.at.test.vals.rand))
                                  },
                                  simplify = FALSE)
    
    intersect.data <- do.call('rbind', intersect.data.list)
    
    beta0 <- c(intersect.data$beta0.start, intersect.data$beta0.end) |>
      unique() |>
      sort()
    
    intersect.data.final <- data.frame(beta0.start = beta0[-length(beta0)],
                                       beta0.end = beta0[-1])
    
    pvals.df <- pvalCalculator.V2(intersect.data.final, 
                                  intersect.data, 
                                  nPerms = ncol(permIndices))

  }
  
  return(pvals.df)
}

exactt.pval.new.iv <- function(Y.temp, X1.temp, X2.temp, permIndices, GX.indices, Q.Z.temp, studentize){
  
  n <- nrow(Y.temp)
  
  Q.Z.temp.dot.X1.temp <- t(Q.Z.temp) %*% matrix(X1.temp[permIndices,], ncol = ncol(permIndices)) # X1.temp.permuted
  Q.Z.temp.dot.Y.temp <- t(Q.Z.temp) %*% matrix(Y.temp[permIndices], ncol = ncol(permIndices)) # Y.temp.permuted
  
  if(studentize == TRUE){
    eps_hat.permuted <- matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                           build_GX(X1.temp, GX.indices) + build_GX(X2.temp, GX.indices),
                                         model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                               ncol = ncol(permIndices))
    
    # nBlocks! x 1 matrix
    Sigma.hat.inverse <- apply(eps_hat.permuted,
                               MARGIN = 2,
                               function(x){
                                 # We don't divide by n here for numerical precision reasons
                                 solve(crossprod(Q.Z.temp * x)) 
                               },
                               simplify = FALSE) |>
      simplify2array()
    
    dim(Sigma.hat.inverse) <- c(ncol(Q.Z.temp), ncol(Q.Z.temp), ncol(permIndices))
  } else{
    identity.matrix <- diag(ncol(Q.Z.temp))
    Sigma.hat.inverse <- array(identity.matrix, dim = c(nrow(identity.matrix),
                                                        ncol(identity.matrix),
                                                        ncol(permIndices)))
  }
  
  a <- sapply(seq_len(ncol(permIndices)), function(g) {
    t(Q.Z.temp.dot.X1.temp[, g]) %*%
      Sigma.hat.inverse[,,g] %*%
      Q.Z.temp.dot.X1.temp[, g]
  })
  
  b <- -2 * sapply(seq_len(ncol(permIndices)), function(g) {
    t(Q.Z.temp.dot.Y.temp[, g]) %*%
      Sigma.hat.inverse[,,g] %*%
      Q.Z.temp.dot.X1.temp[, g]
  })
  
  c <- sapply(seq_len(ncol(permIndices)), function(g) {
    t(Q.Z.temp.dot.Y.temp[, g]) %*%
      Sigma.hat.inverse[,,g] %*%
      Q.Z.temp.dot.Y.temp[, g]
  })
  
  a.identity <- a[1]
  b.identity <- b[1]
  c.identity <- c[1]
  
  a.component <- a.identity - a[-1]
  b.component <- b.identity - b[-1]
  c.component <- c.identity - c[-1]
  
  discriminant <- b.component^2 - 4*a.component*c.component
  
  intersect.data <- data.frame(a.component,
                               b.component,
                               c.component,
                               discriminant,
                               intersectLeft = NA,
                               intersectRight = NA)
  
  intersect.plus <- suppressWarnings((-b.component + sqrt(discriminant))/(2*a.component))
  
  intersect.minus <- suppressWarnings((-b.component - sqrt(discriminant))/(2*a.component))
  
  intersections <- cbind(intersect.plus, 
                         intersect.minus) |> 
    apply(MARGIN = 1, sort) |>
    unlist() |>
    matrix(byrow = TRUE, ncol = 2)
  
  intersect.data[discriminant > 0, c("intersectLeft", "intersectRight")] <- intersections
  
  pvals.df <- pvalCalculator(intersect.data, check.identity = a.identity, intercept = NULL, iv = TRUE, side = "both")
  
  return(pvals.df)
}

pvalCalculator <- function(intersect.data, check.identity, intercept, iv, side){
  
  nPerms <- nrow(intersect.data) + 1
  
  if(iv){
    intersectLeft <- intersect.data$intersectLeft[intersect.data$discriminant >= 0]
    intersectRight <- intersect.data$intersectRight[intersect.data$discriminant >= 0]
    
    beta0 <- c(intersectLeft, intersectRight) |>
      sort()
    
    count_matrix <- outer(beta0, intersectLeft, `>=`) & outer(beta0, intersectRight, `<=`)
    
    growth_condition <- !check.identity > intersect.data$a.component[intersect.data$discriminant >= 0]
    
    count_matrix[, growth_condition] <- !count_matrix[, growth_condition]
    
    pvals <- (apply(count_matrix, MARGIN = 1, sum) + 1)/nPerms
    
    pvals.df <- data.frame(beta0.start = c(-Inf, beta0),
                           beta0.end = c(beta0, Inf),
                           pvals = c(sum(growth_condition, 1)/nPerms, 
                                     pvals[-length(pvals)/2], 
                                     sum(growth_condition, 1)/nPerms))
  } else{
    if(side == "both"){
      NaN.index <- which(is.nan(intersect.data$intersectLeft) | 
                           is.nan(intersect.data$intersectRight))
      
      if(length(NaN.index) > 0){
        filtered.intersect.data <- intersect.data[-NaN.index,]
      } else{
        filtered.intersect.data <- intersect.data
      }
      
      # Identify rows where either intersectLeft or intersectRight is infinite
      left_inf <- is.infinite(filtered.intersect.data$intersectLeft)
      right_inf <- is.infinite(filtered.intersect.data$intersectRight)
      
      # Swap and multiply Inf by -1 for rows where intersectLeft is Inf
      temp <- filtered.intersect.data$intersectLeft[left_inf]
      filtered.intersect.data$intersectLeft[left_inf] <- filtered.intersect.data$intersectRight[left_inf]
      filtered.intersect.data$intersectRight[left_inf] <- -temp
      
      # Swap and multiply Inf by -1 for rows where intersectRight is Inf
      temp <- filtered.intersect.data$intersectRight[right_inf]
      filtered.intersect.data$intersectRight[right_inf] <- filtered.intersect.data$intersectLeft[right_inf]
      filtered.intersect.data$intersectLeft[right_inf] <- -temp
      
      beta0 <- c(filtered.intersect.data$intersectLeft, filtered.intersect.data$intersectRight)
      beta0 <- beta0[is.finite(beta0)] |>
        sort() |>
        unique()
      
      pvals.df <- data.frame(beta0.start = c(-Inf, beta0),
                             beta0.end = c(beta0, Inf))
      
      pvals.df$pvals <- (apply(pvals.df,
                               MARGIN = 1,
                               function(x){
                                 if(all(is.finite(x))){
                                   test.point = mean(x)
                                 } else if(x[1] == -Inf && x[2] == Inf){
                                   test.point = 0
                                 } else if(x[1] == -Inf){
                                   test.point = x[2] - 1
                                 } else if(x[2] == Inf){
                                   test.point = x[1] + 1
                                 }
                                 sum(test.point >= filtered.intersect.data$intersectLeft & 
                                       test.point <= filtered.intersect.data$intersectRight)
                               }) + 1 + length(NaN.index))/nPerms
    } else{
      not.real.intersects.index <- which(!is.finite(intersect.data$intersections) | 
                                           is.nan(intersect.data$intersections))
      
      # Create the count_matrix without explicit conditional check for index length
      if(length(not.real.intersects.index) > 0) {
        filtered.intersections <- intersect.data$intersections[-not.real.intersects.index]
        filtered.slopes <- intersect.data$slope[-not.real.intersects.index]
        
        ## Check non-real intersects
        if(side == "right"){
          # check always greater
          extra <- sum(intercept <= intersect.data$yintercept[not.real.intersects.index])
        } else{
          # check always smaller
          extra <- sum(intercept >= intersect.data$yintercept[not.real.intersects.index])
        }
      } else {
        filtered.intersections <- intersect.data$intersections
        filtered.slopes <- intersect.data$slope
      }
      
      beta0 <- sort(filtered.intersections)
      
      count_matrix <- outer(beta0, 
                            filtered.intersections, 
                            ifelse(side == "right", 
                                   yes = `<=`, 
                                   no = `>=`)
                            )
      
      slope_condition <- !check.identity < filtered.slopes
      
      count_matrix[, slope_condition] <- !count_matrix[, slope_condition]
  
      pvals <- (apply(count_matrix, MARGIN = 1, sum) + 1 + extra)/nPerms
      
      if(side == "right"){
        pval.left <- sum(check.identity < filtered.slopes, 1, extra)/nPerms
        
        pvals.complete <- c(pval.left, pvals)
      } else{
        pval.right <- sum(check.identity > filtered.slopes, 1, extra)/nPerms
        
        pvals.complete <- c(pvals, pval.right)
      }
        
      pvals.df <- data.frame(beta0.start = c(-Inf, beta0),
                             beta0.end = c(beta0, Inf),
                             pvals = pvals.complete)
    }
  }
  
  return(pvals.df)
}

pvalCalculator.V2 <- function(intersect.data.final, intersect.data, nPerms){
  
  pvals <- apply(intersect.data.final,
                 MARGIN = 1,
                 function(x){
                   x <- unlist(x)
                   if(x[1] == -Inf){
                     x[1] <- x[2] - 1
                   } else if(x[2] == Inf){
                     x[2] <- x[1] + 1 
                   }
                   
                   count <- sum(intersect.data$test.stat.smaller[
                     intersect.data$beta0.start < mean(x) & 
                       intersect.data$beta0.end > mean(x)
                     ])
                   
                   pvalues <- (count+1)/nPerms
                   return(pvalues)
                 })
  
  intersect.data.final$pvals <- pvals
 
  return(intersect.data.final)
}

