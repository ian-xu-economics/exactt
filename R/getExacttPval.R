#' Calculate p-values for two sided case.
#'
#' @param Y.temp The response vector for which the test is being performed.
#' @param X1.temp A numeric column vector of the primary variable.
#' @param X2.temp A numeric matrix of the secondary variables.
#' @param indep.X2.index The indices of variables in X2 that are independent to X1.
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
exactt.pval.new.reg <- function(Y.temp, X1.temp, X2.temp, indep.X2.index, permIndices, GX.indices, Q.X1.temp, studentize, side, denominator){
  
  n <- nrow(Y.temp)
  
  if(denominator == "GX1" || studentize == FALSE){
    # We no longer store X1.temp.permuted, Y.temp.permuted, or eps_hat.permuted; it is a RAM nightmare. 
    if(studentize == TRUE){
      # 1 x nPerms matrix
      sigma.hat <- sqrt(t(Q.X1.temp^2) %*% 
                          matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                             build_GX(X1.temp, GX.indices) + build_GX(X2.temp, GX.indices, indep.X2.index),
                                           model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                 ncol = ncol(permIndices))^2)
    } else{
      sigma.hat <- 1
    }
  
    m <- c(t(Q.X1.temp) %*% matrix(X1.temp[permIndices,], ncol = ncol(permIndices))/sigma.hat)
    b <- c(t(Q.X1.temp) %*% matrix(Y.temp[permIndices], ncol = ncol(permIndices))/sigma.hat)
    
    if(side == "both"){
      line.data <- data.frame(a = abs(m[-1]), h = b[-1]/m[-1])
      
      a.identity <- abs(m[1])
      h.identity <- b[1]/m[1]
    } else{
      line.data <- data.frame(m = m[-1], b = b[-1])
      
      m.identity <- m[1]
      b.identity <- b[1]
    }
    
    if(side == "both"){
      line.data <- line.data |> 
        cbind(cbind((line.data$a*line.data$h - a.identity*h.identity)/(line.data$a - a.identity), 
                    (a.identity*h.identity + line.data$a*line.data$h)/(a.identity + line.data$a)) |> 
                apply(MARGIN = 1, \(x){ sort(x, na.last = TRUE) }) |> 
                t() |> 
                data.frame() |> 
                stats::setNames(c("intersectLeft", "intersectRight")))
      
      pvals.df <- pvalCalculator(line.data, 
                                 check.identity = a.identity, 
                                 intercept = h.identity, 
                                 iv = FALSE, 
                                 side = side)
    } else{
      line.data$intersections <- (line.data$b - b.identity)/(m.identity - line.data$m)
      
      pvals.df <- pvalCalculator(line.data, 
                                 check.identity = m.identity, 
                                 intercept = b.identity, 
                                 iv = FALSE, 
                                 side = side)
    }
    
  } else{ # Denominator = X1
    
    Q.X1.GX2.dot.Y.temp.permuted <- matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                                       X1.temp + build_GX(X2.temp, GX.indices, indep.X2.index),
                                                     model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                           ncol = ncol(permIndices))
    
    Q.X1.GX2.dot.X1.temp.permuted <- matrix(stats::lm(matrix(X1.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                                        X1.temp + build_GX(X2.temp, GX.indices, indep.X2.index),
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

exactt.pval.new.iv <- function(Y.temp, X1.temp, X2.temp, indep.X2.index, permIndices, GX.indices, Q.Z.temp, studentize){
  
  n <- nrow(Y.temp)
  
  Q.Z.temp.dot.X1.temp <- t(Q.Z.temp) %*% matrix(X1.temp[permIndices,], ncol = ncol(permIndices)) # X1.temp.permuted
  Q.Z.temp.dot.Y.temp <- t(Q.Z.temp) %*% matrix(Y.temp[permIndices], ncol = ncol(permIndices)) # Y.temp.permuted
  
  if(studentize == TRUE){
    eps_hat.permuted <- matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                           build_GX(X1.temp, GX.indices) + build_GX(X2.temp, GX.indices, indep.X2.index),
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


#' Calculate p-values
#'
#' @param line.data The data.frame containing line data.
#' @param check.identity The slope (or other parameter) of the identity line.
#' @param intercept The x or y intercept of the identity line.
#' @param iv Is this for instrumental variable regression analysis?
#' @param side Which side is the alternative?
#'
#' @importFrom stats reshape
pvalCalculator <- function(line.data, check.identity, intercept, iv, side){
  
  nPerms <- nrow(line.data) + 1
  
  if(iv){
    intersectLeft <- line.data$intersectLeft[line.data$discriminant >= 0]
    intersectRight <- line.data$intersectRight[line.data$discriminant >= 0]
    
    beta0 <- c(intersectLeft, intersectRight) |>
      sort()
    
    count_matrix <- outer(beta0, intersectLeft, `>=`) & outer(beta0, intersectRight, `<=`)
    
    growth_condition <- !check.identity > line.data$a.component[line.data$discriminant >= 0]
    
    count_matrix[, growth_condition] <- !count_matrix[, growth_condition]
    
    pvals <- (apply(count_matrix, MARGIN = 1, sum) + 1)/nPerms
    
    pvals.df <- data.frame(beta0.start = c(-Inf, beta0),
                           beta0.end = c(beta0, Inf),
                           pvals = c(sum(growth_condition, 1)/nPerms, 
                                     pvals[-length(pvals)/2], 
                                     sum(growth_condition, 1)/nPerms))
  } else{
    if(side == "both"){
      NaN.index <- which(is.nan(line.data$intersectLeft) | 
                           is.nan(line.data$intersectRight))
      
      if(length(NaN.index) > 0){
        filtered.line.data <- line.data[-NaN.index,]
      } else{
        filtered.line.data <- line.data
      }
      
      # Identify rows where either intersectLeft or intersectRight is infinite
      left_inf <- is.infinite(filtered.line.data$intersectLeft)
      right_inf <- is.infinite(filtered.line.data$intersectRight)
      
      # Swap and multiply Inf by -1 for rows where intersectLeft is Inf
      temp <- filtered.line.data$intersectLeft[left_inf]
      filtered.line.data$intersectLeft[left_inf] <- filtered.line.data$intersectRight[left_inf]
      filtered.line.data$intersectRight[left_inf] <- -temp
      
      # Swap and multiply Inf by -1 for rows where intersectRight is Inf
      temp <- filtered.line.data$intersectRight[right_inf]
      filtered.line.data$intersectRight[right_inf] <- filtered.line.data$intersectLeft[right_inf]
      filtered.line.data$intersectLeft[right_inf] <- -temp
      
      filtered.line.data.long <- stats::reshape(filtered.line.data,
                                                varying = c("intersectLeft", "intersectRight"),
                                                v.names = "beta0",
                                                timevar = "type",
                                                times = c("left", "right"),
                                                direction = "long")
      
      rownames(filtered.line.data.long) <- NULL
      filtered.line.data.long$id <- NULL
      filtered.line.data.long <- filtered.line.data.long[order(filtered.line.data.long$beta0),]
      filtered.line.data.long <- filtered.line.data.long[!is.infinite(filtered.line.data.long$beta0),]
      
      slope.intersect.booleans <- ifelse(test = (filtered.line.data.long$a < check.identity & 
                                                   filtered.line.data.long$type == "left") |
                                           (filtered.line.data.long$a > check.identity & 
                                              filtered.line.data.long$type == "right") | 
                                           (filtered.line.data.long$a == check.identity &
                                              filtered.line.data.long$h < intercept),
                                         yes = 1,
                                         no = -1)
        
      starter.count <- 1 + length(NaN.index) + sum((filtered.line.data$a == check.identity & filtered.line.data$h > intercept) | 
                                                     filtered.line.data$a > check.identity)
        
      pvals.df <- data.frame(beta0.start = c(-Inf, filtered.line.data.long$beta0),
                             beta0.end = c(filtered.line.data.long$beta0, Inf),
                             pvals = cumsum(c(starter.count, slope.intersect.booleans))/nPerms)
    } else{
      not.real.intersects.index <- which(!is.finite(line.data$intersections) | 
                                           is.nan(line.data$intersections))
      
      # Create the count_matrix without explicit conditional check for index length
      if(length(not.real.intersects.index) > 0){
        filtered.line.data <- line.data[-not.real.intersects.index,]
        filtered.line.data <- filtered.line.data[order(filtered.line.data$intersections),]
        
        ## Check non-real intersects
        same.slope.unusual <- ifelse(side == "right",
                                     yes = sum(intercept <= line.data$b[not.real.intersects.index]), # check always greater
                                     no = sum(intercept >= line.data$b[not.real.intersects.index])) # check always smaller
        
      } else {
        filtered.line.data <- line.data[order(line.data$intersections),]
        
        same.slope.unusual <- 0
      }
      
      slope.booleans <- check.identity < filtered.line.data$m
      count.unusual <- nPerms - sum(slope.booleans) + same.slope.unusual
      slope.booleans.converted <- ifelse(slope.booleans, yes = 1L, no = -1L)
      
      if(side == "left"){
        counts <- c(count.unusual, 
                    slope.booleans.converted) |> 
          cumsum()
      } else{
        counts <- c(count.unusual,
                    rev(slope.booleans.converted)) |> 
          cumsum() |>
          rev()
      }
      
      pvals.df <- data.frame(beta0.start = c(-Inf, filtered.line.data$intersections),
                             beta0.end = c(filtered.line.data$intersections, Inf),
                             pvals = counts/nPerms) 
    }
  }
  
  
  pvals.df <- pvals.df[pvals.df$beta0.start != pvals.df$beta0.end, ]
  rownames(pvals.df) = NULL
  return(pvals.df)
}

pvalCalculator.V2 <- function(line.data.final, line.data, nPerms){
  
  pvals <- apply(line.data.final,
                 MARGIN = 1,
                 function(x){
                   x <- unlist(x)
                   if(x[1] == -Inf){
                     x[1] <- x[2] - 1
                   } else if(x[2] == Inf){
                     x[2] <- x[1] + 1 
                   }
                   
                   count <- sum(line.data$test.stat.smaller[
                     line.data$beta0.start < mean(x) & 
                       line.data$beta0.end > mean(x)
                   ])
                   
                   pvalues <- (count+1)/nPerms
                   return(pvalues)
                 })
  
  line.data.final$pvals <- pvals
  
  return(line.data.final)
}
