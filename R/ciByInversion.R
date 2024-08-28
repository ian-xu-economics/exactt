#' Confidence Interval by Inversion
#'
#' This function calculates the confidence interval for a parameter by inverting the p-values.
#' It provides an option to compute a weighted interval based on taking a convex combination
#' of the two \beta^0 values that generate p-values above and below alpha.
#'
#' @param beta0.df$beta0.pval A numeric vector of p-values.
#' @param beta0.df$beta0 A numeric vector of null values for the parameter of interest, corresponding to the p-values.
#' @param alpha The significance level used for the hypothesis tests.
#' @param weighted A logical value indicating whether to compute a weighted interval. If TRUE, the interval is weighted based on the proximity of p-values to alpha.
#'
#' @return A numeric vector of length 2, representing the lower and upper bounds of the confidence interval.
#'
#' @noRd

ciByInversion <- function(beta0.df, alpha, weighted){
  
  lowerTopIndex <- suppressWarnings(min(which(beta0.df$beta0.pval >= alpha)))
  upperTopIndex <- suppressWarnings(max(which(beta0.df$beta0.pval >= alpha)))
  
  betaLowerTop <- ifelse(lowerTopIndex != -Inf, yes = beta0.df$beta0[lowerTopIndex], no = -Inf)
  betaUpperTop <- ifelse(upperTopIndex != Inf, yes = beta0.df$beta0[upperTopIndex], no = Inf)
  
  # If the top element for the lower bound is equal to alpha
  lowerBound <- ifelse(beta0.df$beta0.pval[lowerTopIndex] == alpha, betaLowerTop, NA)
  upperBound <- ifelse(beta0.df$beta0.pval[upperTopIndex] == alpha, betaUpperTop, NA)
  
  if(is.na(lowerBound)){
    
    lowerBotIndex <- ifelse(is.infinite(lowerTopIndex), 
                            yes = lowerTopIndex,
                            no = lowerTopIndex - 1)
    
    betaLowerBot <- ifelse(is.infinite(lowerBotIndex) || lowerBotIndex == 0,
                           yes = -Inf,
                           no = beta0.df$beta0[lowerBotIndex])
    
    if(weighted && is.finite(betaLowerBot)){
      
      pvalLowerTop = beta0.df$beta0.pval[lowerTopIndex]
      pvalLowerBot = beta0.df$beta0.pval[lowerBotIndex]
      
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
    
    betaUpperBot <- ifelse(is.infinite(upperBotIndex) || upperBotIndex == (length(beta0.df$beta0) + 1),
                           yes = Inf,
                           no = beta0.df$beta0[upperBotIndex])
    
    if(weighted && is.finite(betaUpperBot)){
      
      pvalUpperTop = beta0.df$beta0.pval[upperTopIndex]
      pvalUpperBot = beta0.df$beta0.pval[upperBotIndex]
      
      slope <- (pvalUpperTop - pvalUpperBot)/(betaUpperTop - betaUpperBot)
      
      upperBound <- betaUpperTop - (pvalUpperTop - alpha)/slope
      
    } else{
      upperBound <- betaUpperBot
    }
  }
  
  return(c(lowerBound, upperBound))
}
