# Two-sided case

exactt.pval.new.reg <- function(Y.temp, X1.temp, permIndices, Q.X1.temp, QGX1GX2.temp){

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
  
  lineParams <- data.frame(m = c(Q.X1.temp.dot.X1.temp/sigma.hat), 
                           l = c(Q.X1.temp.dot.Y.temp/Q.X1.temp.dot.X1.temp))
  
  m.identity <- lineParams$m[1]
  l.identity <- lineParams$l[1]
  
  intersect.data <- cbind((m.identity*l.identity - lineParams$m[-1] * lineParams$l[-1])/(m.identity - lineParams$m[-1]), 
                          (m.identity*l.identity + lineParams$m[-1] * lineParams$l[-1])/(m.identity + lineParams$m[-1])) |> 
    apply(MARGIN = 1, sort) |> 
    t() |>
    data.frame()
  
  names(intersect.data) <- c("intersectLeft", "intersectRight")
  
  intersect.data <- cbind(slope = lineParams$m[-1], intersect.data)
  
  pvals.df <- pvalCalculator(intersect.data, m.identity, iv = FALSE, side = "both")
  
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

pvalCalculator <- function(intersect.data, check.identity, iv, side = "both"){
  
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
    beta0 <- c(intersect.data$intersectLeft, intersect.data$intersectRight) |>
      sort()
    
    if(side == "both"){
      count_matrix <- outer(beta0, intersect.data$intersectLeft, `>=`) &
        outer(beta0, intersect.data$intersectRight, `<=`)
      
      slope_condition <- !check.identity > intersect.data$slope
      
      count_matrix[, slope_condition] <- !count_matrix[, slope_condition]
    } else if(side == "left"){
      
    }
    
    pvals <- (apply(count_matrix, MARGIN = 1, sum) + 1)/nPerms
    
    pvals.df <- data.frame(beta0.start = c(-Inf, beta0),
                           beta0.end = c(beta0, Inf),
                           pvals = c(sum(check.identity <= intersect.data$slope, 1)/nPerms, 
                                     pvals[-length(pvals)/2], 
                                     sum(check.identity <= intersect.data$slope, 1)/nPerms))
  }
  
  return(pvals.df)
}
