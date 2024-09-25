# Two-sided case

exactt.pval.new.reg <- function(Y.temp, X1.temp, permIndices, Q.X1.temp, QGX1GX2.temp, side){
  
  #if(GX1){
    Y.temp.permuted <- matrix(Y.temp[permIndices], ncol = ncol(permIndices))
    X1.temp.permuted <- matrix(X1.temp[permIndices,], ncol = ncol(permIndices))
    
    if(!is.null(QGX1GX2.temp)){
      eps_hat.permuted <- QGX1GX2.temp %*% Y.temp.permuted
      # nBlocks! x 1 matrix
      sigma.hat <- sqrt(t(Q.X1.temp^2) %*% eps_hat.permuted^2/nrow(Y.temp))
    } else{
      sigma.hat <- 1
    }
  
    Q.X1.temp.dot.X1.temp <- t(Q.X1.temp) %*% X1.temp.permuted
    Q.X1.temp.dot.Y.temp <- t(Q.X1.temp) %*% Y.temp.permuted
    
    lineParams <- data.frame(m = -c(Q.X1.temp.dot.X1.temp/sigma.hat), 
                             b = c(Q.X1.temp.dot.Y.temp/sigma.hat))
    
    m.identity <- lineParams$m[1]
    b.identity <- lineParams$b[1]
    
    if(side == "both"){
      intersect.data <- cbind((lineParams$b[-1] - b.identity)/(m.identity - lineParams$m[-1]), 
                              (-b.identity - lineParams$b[-1])/(m.identity + lineParams$m[-1])) |> 
        apply(MARGIN = 1, sort) |> 
        t() |>
        data.frame()
      
      names(intersect.data) <- c("intersectLeft", "intersectRight")
      
      intersect.data <- cbind(slope = lineParams$m[-1], intersect.data)
      
      pvals.df <- pvalCalculator(intersect.data, m.identity, iv = FALSE, side = side)
    } else{
      intersect.data <- data.frame(slope = lineParams$m[-1],
                                   intersections = (lineParams$b[-1] - b.identity)/(m.identity - lineParams$m[-1]))
      
      pvals.df <- pvalCalculator(intersect.data, m.identity, iv = FALSE, side = side)
    }
  # } else{
  #   Q.X1GX2.temp <- build_QGX1GX2(X1.temp, GX2.temp, GX.indices, GX1 = FALSE)
  #   
  #   sigma.hat.sq.polynomials <- apply(permIndices,
  #                                     MARGIN = 2,
  #                                     function(x){
  #                                       c(t(Y.temp[x,]) %*% 
  #                                           t(Q.X1GX2.temp) %*%
  #                                           diag(X1^2) %*%
  #                                           Q.X1GX2.temp %*%
  #                                           Y.temp[x,],
  #                                         -2*t(X1.temp[x,]) %*% 
  #                                           t(Q.X1GX2.temp) %*%
  #                                           diag(X1^2) %*%
  #                                           Q.X1GX2.temp %*%
  #                                           Y.temp[x,],
  #                                         t(X1.temp[x,]) %*% 
  #                                           t(Q.X1GX2.temp) %*%
  #                                           diag(X1^2) %*%
  #                                           Q.X1GX2.temp %*%
  #                                           X1.temp[x,]) |>
  #                                         polynom::polynomial()
  #                                     },
  #                                     simplify = FALSE)
  #   
  #   t.sq.polynomials <- apply(permIndices,
  #                             MARGIN = 2,
  #                             function(x){
  #                               c(t(Y.temp[x,]) %*% 
  #                                   Q.X1.temp %*% 
  #                                   t(Q.X1.temp) %*% 
  #                                   Y.temp[x,],
  #                                 -2*t(X1.temp[x,]) %*%
  #                                   Q.X1.temp %*% 
  #                                   t(Q.X1.temp) %*% 
  #                                   Y.temp[x,],
  #                                 t(X1.temp[x,]) %*% 
  #                                   Q.X1.temp %*% 
  #                                   t(Q.X1.temp) %*% 
  #                                   X1.temp[x,]
  #                               ) |>
  #                                 polynom::polynomial()
  #                             },
  #                             simplify = FALSE)
  #   
  #   roots <- sapply(2:ncol(permIndices),
  #                   function(x){
  #                     polyroot(coefficients(t.sq.polynomials[[1]] *
  #                                             sigma.hat.sq.polynomials[[x]] - 
  #                                             t.sq.polynomials[[x]] *
  #                                             sigma.hat.sq.polynomials[[1]]))
  #                   })
  #   
  #   real_roots <- apply(roots,
  #                       MARGIN = 2,
  #                       function(x){
  #                         Re(x)[!abs(Im(x)) > 1e-5]
  #                       })
  # }
  
  return(pvals.df)
}

exactt.pval.new.iv <- function(Y.temp, X1.temp, permIndices, Q.Z.temp, QGX1GX2.temp){
  
  Y.temp.permuted <- matrix(Y.temp[permIndices], ncol = ncol(permIndices))
  X1.temp.permuted <- matrix(X1.temp[permIndices,], ncol = ncol(permIndices))
  Q.Z.temp.dot.X1.temp <- t(Q.Z.temp) %*% X1.temp.permuted
  Q.Z.temp.dot.Y.temp <- t(Q.Z.temp) %*% Y.temp.permuted
  
  if(!is.null(QGX1GX2.temp)){
    eps_hat.permuted <- QGX1GX2.temp %*% Y.temp.permuted
    # nBlocks! x 1 matrix
    Sigma.hat.inverse <- apply(eps_hat.permuted,
                               MARGIN = 2,
                               function(x){
                                 t(Q.Z.temp) %*% (Q.Z.temp*x^2/nrow(Y.temp)) |> 
                                   solve()
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
  
  pvals.df <- pvalCalculator(intersect.data, a.identity, iv = TRUE)
  
  return(pvals.df)
}

pvalCalculator <- function(intersect.data, check.identity, iv, side){
  
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
      beta0 <- sort(c(intersect.data$intersectLeft, intersect.data$intersectRight))
      
      count_matrix <- outer(beta0, intersect.data$intersectLeft, `>=`) &
        outer(beta0, intersect.data$intersectRight, `<=`)
      
      slope_condition <- !check.identity < intersect.data$slope
      
      count_matrix[, slope_condition] <- !count_matrix[, slope_condition]
    } else{
      beta0 <- sort(intersect.data$intersections)
      
      if(side == "right"){
        count_matrix <- outer(beta0, intersect.data$intersections, `>=`)
      } else{
        count_matrix <- outer(beta0, intersect.data$intersections, `<=`)
      }
      
      slope_condition <- !check.identity < intersect.data$slope
      
      count_matrix[, slope_condition] <- !count_matrix[, slope_condition]
    }
  
    pvals <- (apply(count_matrix, MARGIN = 1, sum) + 1)/nPerms
    
    if(side == "both"){
      pval.left <- sum(check.identity > intersect.data$slope, 1)/nPerms
      pval.right <- pval.left
      
      pvals.complete <- c(pval.left, pvals[-length(pvals)/2], pval.right)
    } else if(side == "right"){
      pval.left <- sum(check.identity > intersect.data$slope, 1)/nPerms
      
      pvals.complete <- c(pval.left, pvals)
    } else{
      pval.right <- sum(check.identity > intersect.data$slope, 1)/nPerms
      
      pvals.complete <- c(pvals, pval.right)
    }
      
    pvals.df <- data.frame(beta0.start = c(-Inf, beta0),
                           beta0.end = c(beta0, Inf),
                           pvals = pvals.complete)
    
  }
  
  return(pvals.df)
}
