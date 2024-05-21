#' Exact Testing of Linear Model Coefficients
#'
#' Performs exact tests on specified coefficients of a linear model object
#' using permutation tests to generate p-values and confidence intervals around
#' the estimated coefficients. This function can handle both studentized and
#' non-studentized test statistics, and allows the user to specify various
#' parameters for the test.
#'
#' @param lmObject A linear model object typically created using \code{\link[stats]{lm}}.
#' @param alpha The significance level used for the hypothesis tests; defaults to 0.05.
#' @param variables Optional; a character vector of predictor names to test.
#'        If NULL, all predictors in the model are tested.
#' @param betaNull Optional; a numeric vector of null hypothesis values for the coefficients.
#'        Must be the same length as `variables` if not NULL.
#' @param nBlocks The number of blocks to use for permutations.
#' @param nPerms Optional; the number of permutations to perform.
#'        If NULL or greater than the number of possible permutations, all permutations are used.
#' @param studentize Logical indicating whether to use studentized residuals for the test.
#' @param permutation Optional; a specific permutation vector to use.
#' @param ... Additional arguments passed to `GA::ga()` for optimizing power. 
#' This can include parameters like `popSize`, `maxiter`, `parallel`, etc., 
#' that are used to configure the genetic algorithm. Note that when sample size is large
#' optimizing is computationally expensive and has little effect.
#'
#' @return An object of class 'et', which includes:
#'   - `summary`: A matrix summarizing the test results for each variable.
#'   - `detailed`: A list containing detailed test results for each variable.
#'   - `gaResults`: Optional; a list of results from the `GA::ga()` function, included only when
#'     power optimization is performed via genetic algorithm parameters. Each element of the list
#'     corresponds to results for one of the tested variables, containing details like the best
#'     permutations found, fitness scores, and other GA diagnostics.
#'   - `call`: The matched call.
#'
#' @details
#' The function divides the data into blocks specified by `nBlocks` and performs permutations
#' within these blocks to generate the null distribution of the test statistic. The user can
#' specify a set number of permutations with `nPerms`, or allow the function to calculate all
#' possible permutations if `nPerms` is unspecified or too large.
#'
#' If `studentize` is TRUE, studentized residuals are used to adjust the test statistics,
#' potentially leading to more robust inference under model misspecification.
#'
#' The function allows for a high degree of customization through its parameters and can
#' handle large datasets and complex model structures efficiently.
#'
#' @importFrom stats median
#'
#' @export
exactt <- function(lmObject, 
                   alpha = 0.05,
                   variables = NULL,
                   betaNull = NULL,
                   nBlocks = 5, 
                   nPerms = NULL,
                   studentize = TRUE,
                   permutation = NULL,
                   ...) {
  
  # If betaNull entered, check if equal length as variables
  if(!is.null(betaNull) 
     && is.null(variables)){
    warning("'betaNull' will be ignored since 'variables' is NULL.")
  } else if(!is.null(betaNull) 
            && length(betaNull) != length(variables)){
    stop("Length of 'betaNull' must match length of 'variables'.")
  }
  
  data.in <- lmObject$model
  data.in_n <- nrow(data.in)
  
  n <- floor(data.in_n/nBlocks)*nBlocks
  data <- data.in[1:n, ]
  
  # Demean each column of data
  data <- apply(data, MARGIN = 2, function(x) scale(x, scale = FALSE))
  
  Y <- data[,1, drop = FALSE]
  X <- data[,-1, drop = FALSE]
  
  # Check if desired variables in X
  if(is.null(variables)){
    variables <- colnames(X)
  } else{
    variablesInData <- variables %in% colnames(X)
  
    if(all(!variablesInData)){
      stop("None of the names in 'variables' exist in the set of regressor names.")
    } else if(!all(variablesInData)){
      # Construct message for warning
      nonexistentVariables <- variables[which(!variablesInData)]
      if (length(nonexistentVariables) > 1) {
        formatted_vars <- sprintf('"%s"', nonexistentVariables)
        var_string <- paste(formatted_vars[-length(formatted_vars)], collapse = ", ", sep = "")
        last_var <- paste("and", formatted_vars[length(formatted_vars)])
        message <- paste(var_string, last_var, "do not exist in the set of regressor names. They will be skipped.")
      } else {
        message <- sprintf('"%s" does not exist in the set of regressor names. It will be skipped.', nonexistentVariables)
      }
      
      warning(message)
      
      variables <- variables[which(variablesInData)]
    }
  }
  
  # Construct matrix of block indices
  blockSize <- n/nBlocks
  blockIndexMatrix <- matrix(1:n, 
                             nrow = blockSize, 
                             ncol = nBlocks, 
                             byrow = FALSE)
  
  # If `nPerms` is unspecified or greater than the number of possible permutations, then use all possible permutations. 
  # When number of possible permutations is bigger than MG, then we need to randomly sample.
  if(is.null(nPerms) || nPerms >= factorial(nBlocks)){
    blockPermutations = do.call(rbind, combinat::permn(1:nBlocks))
    permIndices <- apply(blockPermutations, MARGIN = 1, function (x) c(blockIndexMatrix[, x]))
  } else{
    permIndices <- cbind(1:n, replicate(nPerms, c(blockIndexMatrix[, sample(1:nBlocks)])))
  }
  
  # Convert '...' to a list to inspect and possibly modify
  gaArgs <- list(...)
  
  if(length(gaArgs) == 1 && is.null(gaArgs[[1]])){ # Case 1: don't optimize
    if(!is.null(permutation)){
      data <- data[permutation,]
    }
  } else{ # Case 2: optimize
    if(!is.null(permutation)){ 
      warning("Optimization will not be implemented becasue 'permutation' specified. If optimization is desired, either remove 'permutation' or include 'permutation' in matrix of suggestions and pass to function.")
      data <- data[permutation,]
      gaArgs = list()
    } else if("type" %in% names(gaArgs)){
      warning("Custom 'type' value is ignored in this function.")
      gaArgs$type <- NULL
    } else if("fitness" %in% names(gaArgs)){
      warning("Custom 'fitness' value is ignored in this function.")
      gaArgs$fitness <- NULL
    } else if ("lower" %in% names(gaArgs) || "upper" %in% names(gaArgs)) {
      warning("Custom 'lower' and 'upper' values are ignored in this function.")
      gaArgs$lower <- NULL
      gaArgs$upper <- NULL
    }
  }
  
  detailed.out <- vector("list", length = length(variables))
  
  summaryTable.out <- matrix(data = NA_real_,
                             nrow = length(variables), 
                             ncol = 6, 
                             dimnames = list(variables, 
                                             c("Estimate", 
                                               "Pr(>|t|)", 
                                               paste0(alpha*100/2, "%", " W"), 
                                               paste0(100-alpha*100/2, "%", " W"),
                                               paste0(alpha*100/2, "%"), 
                                               paste0(100-alpha*100/2, "%"))))
  
  final_results <- vector("list")
  
  if(length(gaArgs) > 0 && ncol(X) > 1){
    final_results[["gaResults"]] <- vector("list")
  }
  
  if(attr(lmObject$terms, "intercept")){
    summaryTablelm <- summary(lmObject)$coefficients[-1,, drop = FALSE]
  } else{
    summaryTablelm <- summary(lmObject)$coefficients
  }
  
  for(i in 1:nrow(summaryTablelm)){
    if(!row.names(summaryTablelm)[i] %in% variables){
      next
    }
    
    beta_hat <- summaryTablelm[i, 1]
    se <- summaryTablelm[i, 2]
    precisionToUse <- ifelse(se > 0, yes = floor(log(se, base = 10)) - 1, no = -5)
    
    X1.temp <- X[,i, drop = FALSE]
    X2.temp <- X[,-i, drop = FALSE]
    Y.temp <- Y
    
    if(length(gaArgs) > 0){
      
      gaArgs$type <- "permutation"
      gaArgs$fitness <- function(permutation){ -fitness_function(permutation, X1.temp, X2.temp, blockIndexMatrix, blockPermutations) }
      gaArgs$lower <- rep(1, n)
      gaArgs$upper <- rep(n, n)
      
      if(!is.null(gaArgs$parallel) && gaArgs$parallel != FALSE){
        
        if(gaArgs$parallel == TRUE){
          numCores <- parallel::detectCores()
        } else if(is.numeric(gaArgs$parallel) && gaArgs$parallel >= 2){
          numCores <- gaArgs$parallel
        }
        
        cl <- parallel::makeCluster(numCores)
        doParallel::registerDoParallel(cl)
        
        parallel::clusterExport(cl, varlist = c("X1.temp", "X2.temp", "blockIndexMatrix", "blockPermutations", "fitness_function", "build_GX2", "build_GX", "block_permute", "build_QGX2"), envir = environment())
        parallel::clusterCall(cl, library, package = "Matrix", character.only = TRUE)
        parallel::clusterCall(cl, library, package = "MASS", character.only = TRUE)
        parallel::clusterCall(cl, library, package = "combinat", character.only = TRUE)
        parallel::clusterCall(cl, library, package = "dplyr", character.only = TRUE)
        
        gaArgs$parallel = cl
      } 
      
      gaResults <- do.call(GA::ga, gaArgs)
      
      # Close cluster and reset gaArgs
      parallel::stopCluster(cl)
      gaArgs <- list(...)
      
      X1.temp <- X1.temp[gaResults@solution,, drop = FALSE]
      X2.temp <- X2.temp[gaResults@solution,, drop = FALSE]
      Y.temp <- Y.temp[gaResults@solution,, drop = FALSE]
      
      final_results[["gaResults"]][[colnames(X)[i]]] <- gaResults
    }
    
    if(!is.null(betaNull) 
       && !is.null(betaNull[[which(variables == row.names(summaryTablelm)[i])]])){
        betaNullVec <- betaNull[[which(variables == row.names(summaryTablelm)[i])]]
    } else{
      # Find LB and UB roots
      betaNullVec <- get_betaNullVec(Y.temp, X1.temp, X2.temp, alpha, nBlocks, permIndices, beta_hat, se, precisionToUse)
    }
    
    exacttResults <- exactt_pval(betaNullVec, 
                                 Y.temp, 
                                 X1.temp,
                                 X2.temp, 
                                 nBlocks, 
                                 permIndices)
    
    detailed.out[[i]] <- tibble::tibble("betaNull" = betaNullVec, 
                                        "pval" = exacttResults$pval)
                                        #"randomizationDistribution" = list(c(exacttResults$randomizationDist)))

    final_results[["detailed"]][[colnames(X)[i]]] <- detailed.out[[i]]
    
    # Adding weights to find midpoint
    
    lower_bot_idx <- max(which(detailed.out[[i]]$pval <= alpha & betaNullVec < stats::median(betaNullVec)))
    diff_pval_lower_bot <- abs(detailed.out[[i]]$pval[lower_bot_idx] - alpha)
    beta_lower_bot <- betaNullVec[lower_bot_idx]
    
    lower_top_idx <- lower_bot_idx + 1
    diff_pval_lower_top <- abs(detailed.out[[i]]$pval[lower_top_idx] - alpha)
    beta_lower_top <- betaNullVec[lower_top_idx]
    
    weight_lower_top <- 1-diff_pval_lower_top/(diff_pval_lower_bot + diff_pval_lower_top)
    beta_lower <- (1-weight_lower_top)*beta_lower_bot + weight_lower_top*beta_lower_top
    
    upper_bot_idx <- min(which(detailed.out[[i]]$pval <= alpha & betaNullVec > stats::median(betaNullVec)))
    diff_pval_upper_bot <- abs(detailed.out[[i]]$pval[upper_bot_idx] - alpha)
    beta_upper_bot <- betaNullVec[upper_bot_idx]
    
    upper_top_idx <- upper_bot_idx - 1
    diff_pval_upper_top <- abs(detailed.out[[i]]$pval[upper_top_idx] - alpha)
    beta_upper_top <- betaNullVec[upper_top_idx]
    
    weight_upper_top <- 1-diff_pval_upper_top/(diff_pval_upper_bot + diff_pval_upper_top)
    beta_upper <- (1-weight_upper_top)*beta_upper_bot + weight_upper_top*beta_upper_top
    
    idx0 <- which(betaNullVec == 0)
    
    summaryTable.out[i,] <- c(round(beta_hat, digits = -precisionToUse),
                              round(detailed.out[[i]]$pval[idx0], digits = -precisionToUse),
                              round(beta_lower, digits = -precisionToUse),
                              round(beta_upper, digits = -precisionToUse),
                              round(max(betaNullVec[detailed.out[[i]]$pval <= alpha & (betaNullVec <= stats::median(betaNullVec))]), digits = -precisionToUse),
                              round(min(betaNullVec[detailed.out[[i]]$pval <= alpha & (betaNullVec >= stats::median(betaNullVec))]), digits = -precisionToUse))
  }
  
  final_results[["summary"]] <- summaryTable.out
  final_results[["call"]] <- match.call()
  
  class(final_results) <- "et"
  
  return(final_results) 
}


