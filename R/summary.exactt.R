#' Summarize 'exactt' Objects
#'
#' @param object An object of class 'exactt' from the output from `exactt` function.
#' @param ... Additional arguments passed to the summary function.
#' - `alpha`: A value between 0 and 1 indicating the level of significance.
#'
#' @return An object of class 'exactt', which includes:
#'   - `summary`: A matrix summarizing the test results for each variable.
#'   - `detailed`: A list containing detailed test results for each variable.
#'   - `gaResults`: Optional; a list of results from the `GA::ga()` function, included only when
#'     power optimization is performed via genetic algorithm parameters. Each element of the list
#'     corresponds to results for one of the tested variables, containing details like the best
#'     permutations found, fitness scores, and other GA diagnostics.
#'   - `call`: The matched call.
#'
#' @method summary exactt
#' @export
summary.exactt <- function(object, ...){
  
  dots <- list(...)
  alpha <- dots$alpha
  
  if(is.null(alpha)){
    return(object)
  } else{
    object$call$alpha <- alpha
    
    summaryTableList <- vector("list", length = length(object$detailed))
    
    for(i in 1:length(object$detailed)){
      pvals.df <- object$detailed[[i]]
      
      ci.lower.index <- min(which(pvals.df$pvals > alpha))
      ci.upper.index <- max(which(pvals.df$pvals > alpha))
      
      summaryTableList[[i]] <- matrix(data = c(pvals.df$beta0.start[ci.lower.index], 
                                               pvals.df$beta0.end[ci.upper.index]),
                                      nrow = 1, 
                                      ncol = 2)
    }
    
    object$summary[,3:4] <- do.call('rbind', summaryTableList)
    
    return(object)
  }
}
