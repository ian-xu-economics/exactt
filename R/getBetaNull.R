#' Calculate Null Hypothesis Beta Values (Internal Function)
#'
#' This function computes a vector of beta values for null hypothesis testing
#' centered around an estimated value obtained via iterative root finding.
#' The function attempts to find the center where the p-value of the test statistic
#' is closest to 0.8, and then calculates lower and upper bounds around this center
#' where the p-value is equal to the specified alpha level. The function is robust
#' to failure in finding roots by using error handling and fallbacks.
#'
#' @param Y.temp The response vector for which the test is being performed.
#' @param X1.temp A numeric column vector of the primary variable.
#' @param X2.temp A matrix of secondary variables.
#' @param alpha The significance level used to define the confidence interval for
#'        the test statistic.
#' @param nBlocks The number of blocks in the permutation test.
#' @param permIndices A matrix of permutation indices used in the test.
#' @param beta_hat The estimated coefficient on the primary variable for the hypothesis test.
#' @param se The standard error associated with beta_hat.
#' @param precisionToUse The precision level for rounding off the beta values in the
#'        resulting sequence.
#'
#' @return Returns a sorted unique vector of beta values for null hypothesis testing.
#'         These values represent the sequence around the estimated beta center, adjusted
#'         according to the specified alpha and precision.
#'
#' @importFrom rootSolve uniroot.all
#' @importFrom dplyr %>%
#' @importFrom stats uniroot
#'
#' @noRd
getBetaNull = function(Y.temp, X1.temp, X2.temp, Z.temp = NULL, alpha, nBlocks, permIndices, beta_hat, se, precisionToUse){
  
  pvalCenters <- sort(rep(seq(0.85, 0, -0.05), 2),
                      decreasing = TRUE) 
  
  for(i in seq_along(pvalCenters)){
    
    if(is.null(Z.temp)){
      estimatedCenter <- rootSolve::uniroot.all(f = function(x) exactt_pval(x, Y.temp, X1.temp, X2.temp, nBlocks, permIndices, GX1 = TRUE)$pval - pvalCenters[i], 
                                                lower = beta_hat - 3*se, 
                                                upper = beta_hat + 3*se, 
                                                maxiter = 35,
                                                tol = .Machine$double.eps^0.25,
                                                trace = 0,
                                                n = 500) %>%
        try(silent = TRUE)
    } else{
      estimatedCenter <- rootSolve::uniroot.all(f = function(x) exactt_pval_IV(x, Y.temp, X1.temp, X2.temp, Z.temp, nBlocks, permIndices)$pval - pvalCenters[i], 
                                                lower = beta_hat - 3*se, 
                                                upper = beta_hat + 3*se, 
                                                maxiter = 35,
                                                tol = .Machine$double.eps^0.25,
                                                trace = 0,
                                                n = 500) %>%
        try(silent = TRUE)
    }
    
    if(class(estimatedCenter)[1] != "try-error" 
       && length(estimatedCenter) >= 1
       && !all(is.nan(estimatedCenter))){
      estimatedCenter = mean(estimatedCenter)
      break
    } else if(i == length(pvalCenters)){
      estimatedCenter = beta_hat
      #stop("Unable to find estimated center.")
    }
  }
  
  if(is.null(Z.temp)){
    LB <- stats::uniroot(f = function(x) exactt_pval(x, Y.temp, X1.temp, X2.temp, nBlocks, permIndices, GX1 = TRUE)$pval - alpha, 
                       lower = estimatedCenter - 4*se, 
                       upper = estimatedCenter, 
                       extendInt = "upX",
                       maxiter = 35,
                       tol = .Machine$double.eps^0.25,
                       trace = 0) %>%
      try(silent = TRUE)
  } else{
    LB <- stats::uniroot(f = function(x) exactt_pval_IV(x, Y.temp, X1.temp, X2.temp, Z.temp, nBlocks, permIndices)$pval - alpha, 
                   lower = estimatedCenter - 4*se, 
                   upper = estimatedCenter, 
                   extendInt = "upX",
                   maxiter = 35,
                   tol = .Machine$double.eps^0.25,
                   trace = 0) %>%
      try(silent = TRUE)
  }
  
  LBroot <- ifelse(class(LB)[1] == "try-error", 
                   yes = estimatedCenter - 10*se, 
                   no = LB$root)
  
  if(is.null(Z.temp)){
    UB <- stats::uniroot(f = function(x) exactt_pval(x, Y.temp, X1.temp, X2.temp, nBlocks, permIndices, GX1 = TRUE)$pval - alpha, 
                       lower = estimatedCenter, 
                       upper = estimatedCenter + 3*se, 
                       extendInt = "downX",
                       maxiter = 35,
                       tol = .Machine$double.eps^0.25,
                       trace = 0) %>%
      try(silent = TRUE)
  } else{
    UB <- stats::uniroot(f = function(x) exactt_pval_IV(x, Y.temp, X1.temp, X2.temp, Z.temp, nBlocks, permIndices)$pval - alpha, 
                   lower = estimatedCenter, 
                   upper = estimatedCenter + 3*se, 
                   extendInt = "downX",
                   maxiter = 35,
                   tol = .Machine$double.eps^0.25,
                   trace = 0) %>%
      try(silent = TRUE)
  }
  
  UBroot <- ifelse(class(UB)[1] == "try-error", 
                   yes = estimatedCenter + 10*se, 
                   no = UB$root)
  
  betaNull <- seq(round(LBroot - se, digits = -precisionToUse),
                     round(UBroot + se, digits = -precisionToUse),
                     10^(precisionToUse-1)) %>%
    try(silent = TRUE) 
  
  if(class(betaNull)[1] == "try-error"){
    betaNull <- seq(round(estimatedCenter - 12*se, digits = -precisionToUse),
                       round(estimatedCenter + 12*se, digits = -precisionToUse),
                       10^(precisionToUse-1)) 
  }
  
  betaNull <- betaNull %>%
    c(0) %>%
    unique() %>%
    sort()
  
  return(betaNull)
}
