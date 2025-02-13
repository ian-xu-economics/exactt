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
#' @param root.tolerance Tolerance for determining real and extraneous roots (when denominator = "X1" or "noX1").
#'
#' @importFrom polynom polynomial
#' @importFrom stats predict coefficients lm
#' @importFrom cli cli_abort
exactt.pval.new.reg <- function(Y.temp, X1.temp, X2.temp, indep.X2.index, permIndices, GX.indices, Q.X1.temp, studentize, side, denominator, root.tolerance){
  
  n <- nrow(Y.temp)
  
  if(denominator == "GX1" || studentize == FALSE){
    # We no longer store X1.temp.permuted, Y.temp.permuted, or eps_hat.permuted; it is a RAM nightmare. 
    if(studentize == TRUE){
      if(ncol(X2.temp) == 0){
        sigma.hat <- sqrt(t(Q.X1.temp^2) %*% 
                            matrix(stats::lm(matrix(Y.temp[c(permIndices)], ncol = ncol(permIndices)) ~ 
                                               build_GX(X1.temp, GX.indices),
                                             model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                   ncol = ncol(permIndices))^2)
      } else{
        # 1 x nPerms matrix
        sigma.hat <- sqrt(t(Q.X1.temp^2) %*% 
                            matrix(stats::lm(matrix(Y.temp[c(permIndices)], ncol = ncol(permIndices)) ~ 
                                               build_GX(X1.temp, GX.indices) + build_GX(X2.temp, GX.indices, indep.X2.index),
                                             model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                   ncol = ncol(permIndices))^2)
      }
    } else{
      sigma.hat <- 1
    }
    
    line.data.num <- data.frame(b = as.numeric(t(Q.X1.temp) %*% matrix(Y.temp[c(permIndices)], ncol = ncol(permIndices))/sigma.hat),
                                m = as.numeric(-t(Q.X1.temp) %*% matrix(X1.temp[permIndices,], ncol = ncol(permIndices))/sigma.hat))
    
    sigma.hat.sq.polynomials <- NULL
    
    if(side == "both"){
      line.data <- data.frame(a = abs(line.data.num$m[-1]), 
                              h = line.data.num$b[-1]/line.data.num$m[-1])
      
      a.identity <- abs(line.data.num$m[1])
      h.identity <- line.data.num$b[1]/line.data.num$m[1]
    } else{
      line.data <- data.frame(m = line.data.num$m[-1], b = line.data.num$b[-1])
      
      m.identity <- line.data.num$m[1]
      b.identity <- line.data.num$b[1]
    }
    
    if(side == "both"){
      line.data <- line.data |> 
        cbind(cbind(-(line.data$a*line.data$h - a.identity*h.identity)/(line.data$a - a.identity), 
                    -(a.identity*h.identity + line.data$a*line.data$h)/(a.identity + line.data$a)) |> 
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
    
  } else if(denominator == "X1"){ # Denominator = X1
    
    if(ncol(X2.temp) == 0){
      Q.X1.GX2.dot.Y.temp.permuted <- matrix(stats::lm(matrix(Y.temp[c(permIndices)], ncol = ncol(permIndices)) ~ 
                                                         X1.temp,
                                                       model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                             ncol = ncol(permIndices))
      
      Q.X1.GX2.dot.X1.temp.permuted <- matrix(stats::lm(matrix(X1.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                                          X1.temp,
                                                        model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                              ncol = ncol(permIndices))
    } else{
      Q.X1.GX2.dot.Y.temp.permuted <- matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                                         X1.temp + build_GX(X2.temp, GX.indices, indep.X2.index),
                                                       model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                             ncol = ncol(permIndices))
      
      Q.X1.GX2.dot.X1.temp.permuted <- matrix(stats::lm(matrix(X1.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                                          X1.temp + build_GX(X2.temp, GX.indices, indep.X2.index),
                                                        model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                              ncol = ncol(permIndices))
    }
    
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
    
    line.data.num <- apply(permIndices,
                           MARGIN = 2,
                           function(x){
                             c(b = t(Q.X1.temp) %*% Y.temp[x,],
                               m = -t(Q.X1.temp) %*% X1.temp[x,]) 
                           }) |>
      t() |>
      data.frame()
    
    t.num.polynomials <- apply(line.data.num,
                               MARGIN = 1,
                               function(x){
                                 polynom::polynomial(x)
                               },
                               simplify = FALSE)
    
    check.polynomials <- lapply(2:ncol(permIndices),
                                function(x){
                                  t.num.polynomials[[1]]^2 *
                                    sigma.hat.sq.polynomials[[x]] -
                                    t.num.polynomials[[x]]^2 *
                                    sigma.hat.sq.polynomials[[1]]
                                })
    
    real.root.count <- sapply(check.polynomials,
                              function(p){
                                
                                sturm.sequence <- sturm_sequence(p)
                                
                                cauchy.bound <- cauchy_bound(p)
                                
                                real_roots_in_interval(sturm.sequence, 
                                                       -cauchy.bound, 
                                                       cauchy.bound)
                              })
    
    possible.real.roots <- lapply(1:length(check.polynomials),
                                  function(x){
                                 
                                    isolate.real.roots(check.polynomials[[x]], 
                                                       real.root.count[x]) |>
                                      Re()
                                 
                                  })
      
    real.roots <- lapply(1:length(possible.real.roots),
                         function(x){
                           
                           roots <- possible.real.roots[[x]]
                           check.roots.1 <- roots[stats::predict(sigma.hat.sq.polynomials[[1]], roots) > 0 & 
                                                    stats::predict(sigma.hat.sq.polynomials[[1+x]], roots) > 0]
                           
                           if(side == "both"){
                             return(check.roots.1)
                           } else{
                             check.roots.2 <- check.roots.1[abs(stats::predict(t.num.polynomials[[1]], check.roots.1) / sqrt(stats::predict(sigma.hat.sq.polynomials[[1]], check.roots.1)) -
                                                                   stats::predict(t.num.polynomials[[1+x]], check.roots.1) / sqrt(stats::predict(sigma.hat.sq.polynomials[[1+x]], check.roots.1))) < root.tolerance]
                           
                             return(check.roots.2)
                           }
                           
                         })
    # Need to check if the roots are valid in the original problem. Especially if we are dealing with multiple sides.
    # Issue is that we don't multiply by both sides by least common multiple of the denominators.
    # Doing this would be tricky because we don't know what the new polynomial after dividing LCM by sigma.hat.sq.polynomials[[x]]
    
    # This code focuses on the two sided case.
    # Check if the denominator is non-zero.
    intersect.data.list <- lapply(1:length(real.roots),
                                  function(x){
                                    
                                    valid.roots <- real.roots[[x]]
                                    
                                    first <- valid.roots[1] - 999
                                    middle <- (valid.roots[-1] + valid.roots[-length(valid.roots)])/2
                                    last <- valid.roots[length(valid.roots)] + 999
                                    
                                    test.values <- c(first, middle, last)
                                    
                                    values.at.test.vals.test <- stats::predict(t.num.polynomials[[1]], test.values) / sqrt(stats::predict(sigma.hat.sq.polynomials[[1]], test.values))
                                    values.at.test.vals.rand <- stats::predict(t.num.polynomials[[1+x]], test.values) / sqrt(stats::predict(sigma.hat.sq.polynomials[[1+x]], test.values))
                                    
                                    if(side == "both"){
                                      return(data.frame(beta0.start = c(-Inf, valid.roots),
                                                        beta0.end = c(valid.roots, Inf),
                                                        test.stat.smaller = abs(values.at.test.vals.test) < abs(values.at.test.vals.rand)))
                                    } else if(side == "left"){
                                      return(data.frame(beta0.start = c(-Inf, valid.roots),
                                                        beta0.end = c(valid.roots, Inf),
                                                        test.stat.smaller = values.at.test.vals.test < values.at.test.vals.rand))
                                    } else{
                                      return(data.frame(beta0.start = c(-Inf, valid.roots),
                                                        beta0.end = c(valid.roots, Inf),
                                                        test.stat.smaller = values.at.test.vals.test > values.at.test.vals.rand))
                                    }
                                    
                                  })
      
    intersect.data <- do.call('rbind', intersect.data.list)
    
    beta0 <- c(intersect.data$beta0.start, 
               intersect.data$beta0.end) |>
      unique() |>
      sort()
    
    intersect.data.final <- data.frame(beta0.start = beta0[-length(beta0)],
                                       beta0.end = beta0[-1])
    
    pvals.df <- pvalCalculator.V2(intersect.data.final, 
                                  intersect.data, 
                                  nPerms = ncol(permIndices))
    } else if(denominator == "noX1"){
      
      if(ncol(X2.temp) == 0){
        Q.GX2.dot.Y.temp.permuted <- matrix(Y.temp[permIndices], ncol = ncol(permIndices)) - mean(Y.temp)
        
        Q.GX2.dot.X1.temp.permuted <- matrix(X1.temp[permIndices], ncol = ncol(permIndices)) - mean(X1.temp)
      } else{
        Q.GX2.dot.Y.temp.permuted <- matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                                        build_GX(X2.temp, GX.indices, indep.X2.index),
                                                         model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                               ncol = ncol(permIndices))
        
        Q.GX2.dot.X1.temp.permuted <- matrix(stats::lm(matrix(X1.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                                         build_GX(X2.temp, GX.indices, indep.X2.index),
                                                          model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                                ncol = ncol(permIndices))
      }
        
      sigma.hat.sq.polynomials <- lapply(1:ncol(permIndices),
                                         function(x){
                                           c(t(Q.GX2.dot.Y.temp.permuted[,x]) %*%
                                               diag(c(Q.X1.temp^2)) %*%
                                               Q.GX2.dot.Y.temp.permuted[,x],
                                             -2*Q.GX2.dot.X1.temp.permuted[,x] %*%
                                               diag(c(Q.X1.temp^2)) %*%
                                               Q.GX2.dot.Y.temp.permuted[,x],
                                             t(Q.GX2.dot.X1.temp.permuted[,x]) %*%
                                               diag(c(Q.X1.temp^2)) %*%
                                               Q.GX2.dot.X1.temp.permuted[,x]) |>
                                             polynom::polynomial()
                                         })
      
      line.data.num <- apply(permIndices,
                             MARGIN = 2,
                             function(x){
                               c(b = t(Q.X1.temp) %*% Y.temp[x,],
                                 m = -t(Q.X1.temp) %*% X1.temp[x,]) 
                             }) |>
        t() |>
        data.frame()
      
      t.num.polynomials <- apply(line.data.num,
                                 MARGIN = 1,
                                 function(x){
                                   polynom::polynomial(x)
                                 },
                                 simplify = FALSE)
      
      check.polynomials <- lapply(2:ncol(permIndices),
                                  function(x){
                                    t.num.polynomials[[1]]^2 *
                                      sigma.hat.sq.polynomials[[x]] -
                                      t.num.polynomials[[x]]^2 *
                                      sigma.hat.sq.polynomials[[1]]
                                  })
      
      real.root.count <- sapply(check.polynomials,
                                function(p){
                                  
                                  sturm.sequence <- sturm_sequence(p)
                                  
                                  cauchy.bound <- cauchy_bound(p)
                                  
                                  real_roots_in_interval(sturm.sequence, 
                                                         -cauchy.bound, 
                                                         cauchy.bound)
                                })
      
      possible.real.roots <- lapply(1:length(check.polynomials),
                                    function(x){
                                      
                                      isolate.real.roots(check.polynomials[[x]], 
                                                         real.root.count[x]) |>
                                        Re()
                                      
                                    })
      
      real.roots <- lapply(1:length(possible.real.roots),
                           function(x){
                             
                             roots <- possible.real.roots[[x]]
                             check.roots.1 <- roots[stats::predict(sigma.hat.sq.polynomials[[1]], roots) > 0 & 
                                                      stats::predict(sigma.hat.sq.polynomials[[1+x]], roots) > 0]
                             
                             if(side == "both"){
                               return(check.roots.1)
                             } else{
                               check.roots.2 <- check.roots.1[abs(stats::predict(t.num.polynomials[[1]], check.roots.1) / sqrt(stats::predict(sigma.hat.sq.polynomials[[1]], check.roots.1)) -
                                                                    stats::predict(t.num.polynomials[[1+x]], check.roots.1) / sqrt(stats::predict(sigma.hat.sq.polynomials[[1+x]], check.roots.1))) < root.tolerance]
                               
                               return(check.roots.2)
                             }
                             
                           })
      # Need to check if the roots are valid in the original problem. Especially if we are dealing with multiple sides.
      # Issue is that we don't multiply by both sides by least common multiple of the denominators.
      # Doing this would be tricky because we don't know what the new polynomial after dividing LCM by sigma.hat.sq.polynomials[[x]]
      
      # This code focuses on the two sided case.
      # Check if the denominator is non-zero.
      intersect.data.list <- lapply(1:length(real.roots),
                                    function(x){
                                      
                                      valid.roots <- real.roots[[x]]
                                      
                                      first <- valid.roots[1] - 999
                                      middle <- (valid.roots[-1] + valid.roots[-length(valid.roots)])/2
                                      last <- valid.roots[length(valid.roots)] + 999
                                      
                                      test.values <- c(first, middle, last)
                                      
                                      values.at.test.vals.test <- stats::predict(t.num.polynomials[[1]], test.values) / sqrt(stats::predict(sigma.hat.sq.polynomials[[1]], test.values))
                                      values.at.test.vals.rand <- stats::predict(t.num.polynomials[[1+x]], test.values) / sqrt(stats::predict(sigma.hat.sq.polynomials[[1+x]], test.values))
                                      
                                      if(side == "both"){
                                        return(data.frame(beta0.start = c(-Inf, valid.roots),
                                                          beta0.end = c(valid.roots, Inf),
                                                          test.stat.smaller = abs(values.at.test.vals.test) < abs(values.at.test.vals.rand)))
                                      } else if(side == "left"){
                                        return(data.frame(beta0.start = c(-Inf, valid.roots),
                                                          beta0.end = c(valid.roots, Inf),
                                                          test.stat.smaller = values.at.test.vals.test < values.at.test.vals.rand))
                                      } else{
                                        return(data.frame(beta0.start = c(-Inf, valid.roots),
                                                          beta0.end = c(valid.roots, Inf),
                                                          test.stat.smaller = values.at.test.vals.test > values.at.test.vals.rand))
                                      }
                                      
                                    })
      
      intersect.data <- do.call('rbind', intersect.data.list)
      
      beta0 <- c(intersect.data$beta0.start, 
                 intersect.data$beta0.end) |>
        unique() |>
        sort()
      
      intersect.data.final <- data.frame(beta0.start = beta0[-length(beta0)],
                                         beta0.end = beta0[-1])
      
      pvals.df <- pvalCalculator.V2(intersect.data.final, 
                                    intersect.data, 
                                    nPerms = ncol(permIndices))
  } else{
    cli::cli_abort("Inputted value into `denominator` parameter is not recognized.")
  }
  
  return(list(pvals.df = pvals.df,
              line.data.num = line.data.num,
              line.data.denom.sq = sigma.hat.sq.polynomials))
}

exactt.pval.new.iv <- function(Y.temp, X1.temp, X2.temp, indep.X2.index, permIndices, GX.indices, Q.Z.temp, studentize){
  
  n <- nrow(Y.temp)
  
  Q.Z.temp.dot.X1.temp <- t(Q.Z.temp) %*% matrix(X1.temp[permIndices,], ncol = ncol(permIndices)) # X1.temp.permuted
  Q.Z.temp.dot.Y.temp <- t(Q.Z.temp) %*% matrix(Y.temp[permIndices], ncol = ncol(permIndices)) # Y.temp.permuted
  
  if(studentize == TRUE){
    if(ncol(X2.temp) == 0){
      eps_hat.permuted <- matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                             build_GX(X1.temp, GX.indices),
                                           model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                 ncol = ncol(permIndices))
    } else{
      eps_hat.permuted <- matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                             build_GX(X1.temp, GX.indices) + build_GX(X2.temp, GX.indices, indep.X2.index),
                                           model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                 ncol = ncol(permIndices))
    }
    
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

exactt.pval.wald <- function(Y.temp, X1.temp, X2.temp, permIndices, GX.indices, Q.Z.temp, studentize, beta.null.matrix){
  
  n <- nrow(Y.temp)
  
  if(studentize == TRUE){
    if(ncol(X2.temp) == 0){
      eps_hat.permuted <- matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                             build_GX(X1.temp, GX.indices),
                                           model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                 ncol = ncol(permIndices))
    } else{
      eps_hat.permuted <- matrix(stats::lm(matrix(Y.temp[permIndices], ncol = ncol(permIndices)) ~ 
                                             build_GX(X1.temp, GX.indices) + build_GX(X2.temp, GX.indices, indep.X2.index),
                                           model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals,
                                 ncol = ncol(permIndices))
    }  
      
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
    Sigma.hat.inverse <- array(identity.matrix, 
                               dim = c(nrow(identity.matrix),
                                       ncol(identity.matrix),
                                       ncol(permIndices)))
  }
  
  if(!is.atomic(beta.null.matrix) && 
     length(beta.null.matrix) > 1 && 
     !is.null(dim(beta.null.matrix)) &&
     nrow(beta.null.matrix) == 1 &&
     all(c(as.matrix(beta.null.matrix)) == 0)){
    
    Q.Z.temp.dot.Y.temp <- apply(permIndices,
                                 MARGIN = 2,
                                 FUN = function(x){
                                   t(Q.Z.temp) %*% Y.temp[x, drop = FALSE]
                                 },
                                 simplify = FALSE) |>
      simplify2array()
    
    randomization.stats <- sapply(1:ncol(permIndices),
                                  FUN = function(i) {
                                    t(matrix(Q.Z.temp.dot.Y.temp[,,i, drop = FALSE],
                                             nrow = ncol(Q.Z.temp),
                                             ncol = ncol(Y.temp))) %*% 
                                      Sigma.hat.inverse[,,i] %*%
                                      matrix(Q.Z.temp.dot.Y.temp[,,i, drop = FALSE],
                                             nrow = ncol(Q.Z.temp),
                                             ncol = ncol(Y.temp))
                                  })
    
    p.value <- mean(randomization.stats[1] <= randomization.stats)
    
    result <- cbind(beta.null.matrix,
                    p.value) |>
      data.frame() 
    
  } else{
    
    Q.Z.temp.dot.X1.Y.temp <- apply(permIndices,
                                    MARGIN = 2,
                                    FUN = function(x){
                                      t(Q.Z.temp) %*% cbind(X1.temp, Y.temp)[x,, drop = FALSE]
                                    },
                                    simplify = FALSE) |>
      simplify2array() |>
      array(dim = c(ncol(Q.Z.temp), ncol(X1.temp) + ncol(Y.temp), ncol(permIndices)))
    
    beta.null.matrix <- cbind(as.matrix(beta.null.matrix), 1)
    
    randomization.stats <- get.wald.randomization.stats(Q.Z.temp.dot.X1.Y.temp, 
                                                        Sigma.hat.inverse, 
                                                        beta.null.matrix)
    
    p.value <- apply(randomization.stats[,1] <= randomization.stats,
                     MARGIN = 1,
                     mean)
    
    result <- cbind(beta.null.matrix[,-ncol(beta.null.matrix), drop = FALSE],
                    p.value) |>
      data.frame() 
    
  }

  return(result)
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
                   
                   count <- sum(line.data$test.stat.smaller[line.data$beta0.start < mean(x) & 
                                                              line.data$beta0.end > mean(x)])
                   
                   pvalues <- (count+1)/nPerms
                   return(pvalues)
                 })
  
  line.data.final$pvals <- pvals
  
  return(line.data.final)
}


# find_CI <- function(df, var, alpha = 0.05) {
#   # 1) Filter out points where p.value <= alpha
#   df_accepted <- df[df$p.value > alpha, ]
#   
#   # 2) If nothing remains, the entire interval is "accepted" => (-Inf, Inf).
#   if(nrow(df_accepted) == 0) {
#     return(c(-Inf, Inf))
#   }
#   
#   # 3) Identify the range of the *entire* grid for the variable
#   var_all_min <- min(df[[var]])
#   var_all_max <- max(df[[var]])
#   
#   # 4) Among the accepted points, get the min and max for that variable
#   var_acc_min <- min(df_accepted[[var]])
#   var_acc_max <- max(df_accepted[[var]])
#   
#   # 5) Decide if the lower bound is -Inf or the actual min
#   #    - If the accepted min is the same as the absolute min of the entire grid,
#   #      we interpret that as going to -Inf on that side.
#   ci_lower <- if(var_acc_min == var_all_min) -Inf else var_acc_min
#   
#   # 6) Similarly for the upper bound
#   ci_upper <- if(var_acc_max == var_all_max) Inf else var_acc_max
#   
#   # 7) Return the confidence interval
#   c(ci_lower, ci_upper)
# }


##### EXTRACT THE REAL ROOTS FROM COMPLEX ROOTS
# Function to perform polynomial division returning quotient and remainder
poly_div <- function(f, g) {
  f_coef <- coef(f)
  g_coef <- coef(g)
  
  # If degree of f < degree of g, quotient = 0, remainder = f
  if (length(f_coef) < length(g_coef)) {
    return(list(quotient = polynomial(0), remainder = f))
  }
  
  q_coef <- rep(0, length(f_coef) - length(g_coef) + 1)
  r_coef <- f_coef
  deg_g <- length(g_coef) - 1
  
  for (i in seq_along(q_coef)) {
    q_coef[i] <- r_coef[1]/g_coef[1]
    for (j in seq_along(g_coef)) {
      r_coef[j] <- r_coef[j] - q_coef[i] * g_coef[j]
    }
    r_coef <- r_coef[-1]
  }
  
  q <- polynomial(q_coef)
  r <- if (length(r_coef) == 0) polynomial(0) else polynomial(r_coef)
  list(quotient = q, remainder = r)
}

# Construct the Sturm sequence
sturm_sequence <- function(p) {
  seq_list <- list(p, stats::deriv(p))
  i <- 2
  repeat {
    division <- poly_div(seq_list[[i-1]], seq_list[[i]])
    next_poly <- -division$remainder
    if (all(coef(next_poly) == 0)) break
    seq_list[[i+1]] <- next_poly
    i <- i + 1
  }
  seq_list
}

# Sign changes count function
sign_changes <- function(seq, x) {
  vals <- sapply(seq, function(poly) predict(poly, x))
  # Remove zero evaluations
  vals <- vals[vals != 0]
  if (length(vals) < 2) return(0)
  sum(diff(sign(vals)) != 0)
}

# Counting number of real roots in an interval [a,b]
real_roots_in_interval <- function(seq, a, b) {
  sign_changes_a <- sign_changes(seq, a)
  sign_changes_b <- sign_changes(seq, b)
  sign_changes_a - sign_changes_b
}

# Compute Cauchy's bound
cauchy_bound <- function(p) {
  coefs <- coef(p)
  # Leading coefficient (highest degree term)
  a_n <- coefs[length(coefs)]
  # Exclude the leading coefficient for computing the max ratio
  other_coefs <- coefs[-length(coefs)]
  M <- 1 + max(abs(other_coefs / a_n))
  M
}

# Isolate and approximate the roots
# We will scan the interval [-M, M] in smaller steps and find subintervals containing a single root
isolate.real.roots <- function(polynomial, real.root_count) {
  
  all.roots <- solve(polynomial)
  all.roots.im <- Im(all.roots)
  
  real.roots.im <- all.roots.im |>
    abs() |>
    sort() |>
    utils::head(real.root_count)
  
  real.roots <- all.roots[which(abs(all.roots.im) %in% real.roots.im)] |>
    Re()
  
  return(real.roots)
  
}

# Alternative method to find real roots
# Idea is to divide cauchy bounds up into little chunks, then use Sturm's Theorem
# To find how many roots there are in that little interval. If it is non-zero, use 
# uniroot() to locate the root.
# isolate_and_find_roots <- function(seq, p, lower, upper, n_intervals = 1000) {
#   xs <- seq(lower, upper, length.out = n_intervals + 1)
#   found_roots <- numeric(0)
#   for (i in seq_len(n_intervals)) {
#     count <- real_roots_in_interval(seq, xs[i], xs[i+1])
#     if (count == 1) {
#       # Exactly one root in (xs[i], xs[i+1])
#       # Use uniroot to find it
#       root <- uniroot(function(x) predict(p, x), interval = c(xs[i], xs[i+1]))$root
#       found_roots <- c(found_roots, root)
#     }
#   }
#   unique(found_roots)
# }
