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
#' @param nBlocks The number of blocks to use for block permutations.
#' @param permIndices A matrix of permutation indices used in the test.
#' @param beta_hat The estimated coefficient on the primary variable for the hypothesis test.
#' @param se The standard error associated with beta_hat.
#' @param precisionToUse The precision level for rounding off the beta values in the
#'        resulting sequence.
#' @param GX1 GX1 Logical indicating whether to use GX1 or X1 when constructing eps_hat. 
#'        Using X1 has slightly more power at slightly more computational expense.
#'
#' @return Returns a sorted unique vector of beta values for null hypothesis testing.
#'         These values represent the sequence around the estimated beta center, adjusted
#'         according to the specified alpha and precision.
#'
#' @importFrom dplyr %>%
#'
#' @noRd
getBetaNull <- function(Y.temp, X1.temp, X2.temp, Z.temp = NULL, alpha, nBlocks, permIndices, GX.indices, beta_hat, se, studentize, precisionToUse, GX1){
  
  beta0 <- seq(from = beta_hat - 25*se,
               to = beta_hat + 25*se,
               by = se/4)
  
  if(is.null(Z.temp)){
    beta0.pvals <- exactt_pval(beta0, Y.temp, X1.temp, X2.temp, nBlocks, permIndices, GX.indices, studentize, GX1)$pval
  } else{
    beta0.pvals <- exactt_pval_IV(beta0, Y.temp, X1.temp, X2.temp, Z.temp, nBlocks, permIndices, GX.indices, studentize, GX1)$pval
  }
  
  # for(i in 1:20){
  #   if(beta0.pvals[1] < alpha &&
  #      beta0.pvals[length(beta0.pvals)] < alpha &&
  #      any(beta0.pvals >= alpha)){
  #     break
  #   }
  # 
  #   # Expand left
  #   beta0.left <- seq(from = beta0[1] - 10*se,
  #                            to = beta0[1] - se/4, # so there are no repeats
  #                            se/2)
  # 
  #   # Expand right
  #   beta0.right <- seq(from = beta0[length(beta0)] + se/4, # so there are no repeats
  #                      to = beta0[length(beta0)] + 10*se,
  #                      se/2)
  # 
  #   if(is.null(Z.temp)){
  #     beta0.pvals.left <- exactt_pval(beta0.left, Y.temp, X1.temp, X2.temp, nBlocks, permIndices, studentize, GX1)$pval
  #     beta0.pvals.right <- exactt_pval(beta0.right, Y.temp, X1.temp, X2.temp, nBlocks, permIndices, studentize, GX1)$pval
  #   } else{
  #     beta0.pvals.left <- exactt_pval_IV(beta0.left, Y.temp, X1.temp, X2.temp, Z.temp, nBlocks, permIndices)$pval
  #     beta0.pvals.right <- exactt_pval_IV(beta0.right, Y.temp, X1.temp, X2.temp, Z.temp, nBlocks, permIndices)$pval
  #   }
  # 
  #   beta0 <- c(beta0.left, beta0, beta0.right)
  #   beta0.pvals <- c(beta0.pvals.left, beta0.pvals, beta0.pvals.right)
  # }
  
  # Check if left and right most beta0 are below alpha
  if(!(beta0.pvals[1] <= alpha && 
       beta0.pvals[length(beta0.pvals)] <= alpha)){
    
     beta0.final <- seq(round(beta_hat - 20*se, digits = -precisionToUse),
                        round(beta_hat + 20*se, digits = -precisionToUse),
                        10^(precisionToUse + 1)) %>%
       c(0) %>%
       unique() %>%
       sort()
  } else{
  
    peakIndex <- which.max(beta0.pvals)
    beta0.left <- beta0[1:peakIndex]
    beta0.right <- beta0[peakIndex:length(beta0.pvals)]
    beta0.pvals.left <- beta0.pvals[1:peakIndex]
    beta0.pvals.right <- beta0.pvals[peakIndex:length(beta0.pvals)]
   
    # Find left point alpha root
    # Check if there is any zero in the vector
    zero_indices.left <- which(beta0.pvals.left == alpha)
    beta0.left.bound <- ifelse(length(zero_indices.left) > 0,
                               yes = beta0.left[min(zero_indices.left)],
                               no = beta0.left[beta0.pvals.left < alpha][which.max(beta0.pvals.left[beta0.pvals.left < alpha])]) + se/4
    
    # Find right point alpha root
    # Check if there is any zero in the vector
    zero_indices.right <- which(beta0.pvals.right == alpha)
    beta0.right.bound <- ifelse(length(zero_indices.right) > 0,
                                yes = beta0.right[max(zero_indices.right)],
                                no = beta0.right[beta0.pvals.right < alpha][which.max(beta0.pvals.right[beta0.pvals.right < alpha])]) - se/4
    
    beta0.final <- c(seq(round(beta0.left.bound - se/2, digits = -precisionToUse),
                         round(beta0.left.bound, digits = -precisionToUse),
                         10^(precisionToUse-1)),
                     seq(round(beta0.left.bound, digits = -precisionToUse),
                         round(beta0.right.bound, digits = -precisionToUse),
                         se/4),
                     seq(round(beta0.right.bound, digits = -precisionToUse),
                         round(beta0.right.bound + se/2, digits = -precisionToUse),
                         10^(precisionToUse-1)),
                     0) %>%
      unique() %>%
      sort()
  }
  
  return(beta0.final)
}
