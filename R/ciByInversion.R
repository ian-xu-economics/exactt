#' Confidence Interval by Inversion
#'
#' This function calculates the confidence interval for a parameter by inverting the p-values.
#' It provides an option to compute a weighted interval based on taking a convex combination
#' of the two \beta^0 values that generate p-values above and below alpha.
#'
#' @param pvals A numeric vector of p-values.
#' @param betaNullVec A numeric vector of null values for the parameter of interest, corresponding to the p-values.
#' @param alpha The significance level used for the hypothesis tests.
#' @param weighted A logical value indicating whether to compute a weighted interval. If TRUE, the interval is weighted based on the proximity of p-values to alpha.
#'
#' @return A numeric vector of length 2, representing the lower and upper bounds of the confidence interval.
#'
#' @noRd

ciByInversion <- function(betaNullVec, pvals, alpha, weighted){
  
  lowerTopIndex <- suppressWarnings(min(which(pvals >= alpha)))
  upperTopIndex <- suppressWarnings(max(which(pvals >= alpha)))
  
  betaLowerTop <- ifelse(lowerTopIndex != -Inf, yes = betaNullVec[lowerTopIndex], no = -Inf)
  betaUpperTop <- ifelse(upperTopIndex != Inf, yes = betaNullVec[upperTopIndex], no = Inf)
  
  # If the top element for the lower bound is equal to alpha
  lowerBound <- ifelse(pvals[lowerTopIndex] == alpha, betaLowerTop, NA)
  upperBound <- ifelse(pvals[upperTopIndex] == alpha, betaUpperTop, NA)
  
  if(is.na(lowerBound)){
    
    lowerBotIndex <- ifelse(is.infinite(lowerTopIndex), 
                            yes = lowerTopIndex,
                            no = lowerTopIndex - 1)
    
    betaLowerBot <- ifelse(is.infinite(lowerBotIndex) || lowerBotIndex == 0,
                           yes = -Inf,
                           no = betaNullVec[lowerBotIndex])
    
    if(weighted && is.finite(betaLowerBot)){
      
      pvalLowerTop = pvals[lowerTopIndex]
      pvalLowerBot = pvals[lowerBotIndex]
      
      slope <- (pvalLowerTop - pvalLowerBot)/(betaLowerTop - betaLowerBot)
      
      lowerBound <- betaLowerTop - (pvalLowerTop - alpha)/slope
      
    } else{
      lowerBound <- betaLowerBot
    }
  }
  
  if(is.na(upperBound)){
    
    upperBotIndex <- ifelse(is.infinite(upperTopIndex), 
                            yes = upperTopIndex, 
                            no = upperTopIndex + 1)
    
    betaUpperBot <- ifelse(is.infinite(upperBotIndex) || upperBotIndex == (length(betaNullVec) + 1),
                           yes = Inf,
                           no = betaNullVec[upperBotIndex])
    
    if(weighted && is.finite(betaUpperBot)){
      
      pvalUpperTop = pvals[upperTopIndex]
      pvalUpperBot = pvals[upperBotIndex]
      
      slope <- (pvalUpperTop - pvalUpperBot)/(betaUpperTop - betaUpperBot)
      
      upperBound <- betaUpperTop - (pvalUpperTop - alpha)/slope
      
    } else{
      upperBound <- betaUpperBot
    }
  }
  
  return(c(lowerBound, upperBound))
}
