#' @title Exact t-Test Regression
#' @description Performs an exact t-test regression with optional optimization using a genetic algorithm.
#' @param formula An object of class \code{\link[stats]{formula}} specifying the model.
#' @param data A data frame containing the variables in the model.
#' @param n An integer specifying the number of observations to use.
#' @param Y.var The name of the dependent variable.
#' @param X.var The name of the independent variable(s).
#' @param blockIndexMatrix A matrix specifying the block indices for permutations.
#' @param variables A vector of variable indices to be tested.
#' @param nBlocks An integer specifying the number of blocks.
#' @param alpha A numeric value for the significance level.
#' @param betaNull A list of null hypothesis values for the regression coefficients.
#' @param nPerms An integer specifying the number of permutations. If \code{NULL}, all possible permutations are used.
#' @param studentize A logical value indicating whether to use studentized test statistics.
#' @param permutation A vector specifying a permutation of the data. If \code{NULL}, no permutation is applied.
#' @param optimize A logical value indicating whether to use a genetic algorithm for optimization.
#' @param gaArgs A list of arguments to be passed to the genetic algorithm function.
#' @return A list containing the results of the exact t-test regression, including detailed results and summary tables.
#' @importFrom stats as.formula model.matrix lm summary.lm
#' @importFrom combinat permn
#' @importFrom cli cli_warn
#' @importFrom parallel makeCluster detectCores stopCluster clusterExport clusterCall
#' @importFrom doParallel registerDoParallel
#' @importFrom GA ga
#' @importFrom tibble tibble
#' @importFrom dplyr filter
#' @keywords internal
exacttReg <- function(formula, 
                      data, 
                      n,
                      Y.var, 
                      X.var, 
                      blockIndexMatrix,
                      variables,
                      nBlocks,
                      alpha,
                      betaNull,
                      nPerms,
                      studentize,
                      permutation,
                      optimize,
                      gaArgs){
  
  # Convert the data to a model matrix, which will handle one-hot encoding
  YFormulaString <- paste("~", paste(c(Y.var, 0), collapse = " + "))
  Y <- model.matrix(stats::as.formula(YFormulaString), data)
  
  X <- stats::model.matrix(formula, data)
  assign <- attr(X, "assign")
  
  Y.use <- Y[1:n, , drop = FALSE]
  X.use <- X[1:n, , drop = FALSE]
  
  # If `nPerms` is unspecified or greater than the number of possible permutations, then use all possible permutations. 
  # When number of possible permutations is bigger than MG, then we need to randomly sample.
  if(is.null(nPerms) || nPerms >= factorial(nBlocks)){
    blockPermutations <- do.call(rbind, combinat::permn(1:nBlocks))
    permIndices <- apply(blockPermutations, MARGIN = 1, function (x) c(blockIndexMatrix[, x]))
  } else{
    permIndices <- cbind(1:n, replicate(nPerms, c(blockIndexMatrix[, sample(1:nBlocks)])))
  }
  
  if(!optimize){ # Case 1: don't optimize
    if(!is.null(permutation)){
      X.use <- X.use[permutation,, drop = FALSE]
      Y.use <- Y.use[permutation,, drop = FALSE]
    }
  } else{ # Case 2: optimize
    
    if(!is.null(permutation)){ 
      cli::cli_warn("Optimization will not be implemented because 'permutation' specified. If optimization is desired, either remove 'permutation' or include 'permutation' in matrix of suggestions and pass to function.")
      X.use <- X.use[permutation,, drop = FALSE]
      Y.use <- Y.use[permutation,, drop = FALSE]
      optimize <- FALSE
    } else{
      if("type" %in% names(gaArgs)){
        warning("Custom 'type' value is ignored in this function.")
        gaArgs$type <- NULL
      } 
      if("fitness" %in% names(gaArgs)){
        warning("Custom 'fitness' value is ignored in this function.")
        gaArgs$fitness <- NULL
      } 
      if ("lower" %in% names(gaArgs) || "upper" %in% names(gaArgs)) {
        warning("Custom 'lower' and 'upper' values are ignored in this function.")
        gaArgs$lower <- NULL
        gaArgs$upper <- NULL
      }
      
      gaArgs$type <- "permutation"
      gaArgs$fitness <- function(permutation){ fitness_function(permutation, X1.temp, X2.temp, Z.temp = NULL, blockIndexMatrix, blockPermutations) }
      gaArgs$lower <- rep(1, n)
      gaArgs$upper <- rep(n, n)
    }
  }
  
  lmObject <- stats::lm(formula, data = data)
  
  detailed.out <- vector("list", length = sum(assign %in% variables))
  
  summaryTable.out <- matrix(data = NA_real_,
                             nrow = sum(assign %in% variables), 
                             ncol = 6, 
                             dimnames = list(colnames(X)[assign %in% variables], 
                                             c("Estimate", 
                                               "Pr(>|t|)", 
                                               paste0(alpha*100/2, "%", " W"), 
                                               paste0(100-alpha*100/2, "%", " W"),
                                               paste0(alpha*100/2, "%"), 
                                               paste0(100-alpha*100/2, "%"))))
  
  final_results <- vector("list")
  
  if(optimize && ncol(X) > 1){
    final_results[["gaResults"]] <- vector("list")
  }
  
  summaryTablelm <- summary(lmObject)$coefficients
  
  rowCounter <- 1
  
  for(i in seq_along(assign)){
    if(assign[i] == 0 | !assign[i] %in% variables){
      next
    }
    
    beta_hat <- summaryTablelm[i, 1]
    se <- summaryTablelm[i, 2]
    precisionToUse <- ifelse(se > 0, yes = floor(log(se, base = 10)) - 1, no = -5)
    
    X1.temp <- X.use[,i, drop = FALSE]
    X2.temp <- X.use[,-i, drop = FALSE]
    Y.temp <- as.matrix(Y.use)
    
    if(optimize){
      
      if(!is.null(gaArgs$parallel) && gaArgs$parallel != FALSE){
        
        ogParArg <- gaArgs$parallel
        
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
        
        gaArgs$parallel <- cl
      } 
      
      gaResults <- do.call(GA::ga, gaArgs)
      
      # Close cluster if parallel is true
      if(!is.null(gaArgs$parallel) && ogParArg != FALSE){
        parallel::stopCluster(cl)
        gaArgs$parallel <- ogParArg
      }
      
      X1.temp <- X1.temp[gaResults@solution,, drop = FALSE]
      X2.temp <- X2.temp[gaResults@solution,, drop = FALSE]
      Y.temp <- Y.temp[gaResults@solution,, drop = FALSE]
      
      final_results[["gaResults"]][[colnames(X)[i]]] <- gaResults
    }
    
    if(!is.null(betaNull) && !is.null(betaNull[[which(variables == row.names(summaryTablelm)[i])]])){
      betaNullVec <- betaNull[[which(variables == row.names(summaryTablelm)[i])]]
    } else{
      # Find LB and UB roots
      betaNullVec <- getBetaNull(Y.temp, X1.temp, X2.temp, Z.temp = NULL, alpha, nBlocks, permIndices, beta_hat, se, precisionToUse)
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
    
    summaryTable.out[rowCounter,] <- c(beta_hat,
                                       exacttResults$pval[which(betaNullVec == 0)],
                                       ciByInversion(betaNullVec, exacttResults$pval, alpha, weighted = TRUE),
                                       ciByInversion(betaNullVec, exacttResults$pval, alpha, weighted = FALSE))
    rowCounter <- rowCounter + 1
  }
  
  final_results[["summary"]] <- summaryTable.out
  
  return(final_results)
}
