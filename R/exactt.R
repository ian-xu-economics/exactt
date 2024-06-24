#' Exact Testing of Linear Model Coefficients
#'
#' Performs exact tests on specified coefficients of a linear model object
#' using permutation tests to generate p-values and confidence intervals around
#' the estimated coefficients. This function can handle both studentized and
#' non-studentized test statistics, and allows the user to specify various
#' parameters for the test.
#'
#' @param model A formula specifying the model..
#' @param data A data frame or matrix containing the variables used in the model.
#' @param alpha The significance level used for the hypothesis tests; defaults to 0.05.
#' @param variables Optional; a character vector of predictor names to test.
#'        If NULL, all predictors in the model are tested.
#' @param betaNull Optional; a numeric vector of null hypothesis values for the coefficients.
#'        Must be the same length as `variables` if not NULL.
#' @param nBlocks The number of blocks to use for permutations.
#' @param nPerms Optional; the number of permutations to perform.
#'        If NULL or greater than the number of possible permutations, all permutations are used.
#' @param studentize Logical indicating whether to use studentized residuals for the test.
#' @param permutation Optional; a specific permutation vector to rearrange order of data.
#' @param optimize Logical indicating whether to optimize the ordering of the data (default is FALSE).
#' @param GX1 Logical indicating whether to use GX1 or X1 when constructing eps_hat.
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
#' @importFrom stats median formula model.matrix
#' @importFrom Formula Formula
#' @importFrom cli cli_abort
#' 
#' @export
exactt <- function(model,
                   data,
                   alpha = 0.05,
                   variables = NULL,
                   betaNull = NULL,
                   nBlocks = 5,
                   nPerms = NULL,
                   studentize = TRUE,
                   permutation = NULL,
                   optimize = FALSE,
                   GX1 = TRUE,
                   ...) {
  
  call <- match.call(expand.dots = TRUE)
  
  ####### Do checks #######
  
  # Check if `betaNull` provided. If so, check if length equal to number of variables
  if(!is.null(betaNull) && is.null(variables)){
    warning("'betaNull' will be ignored since 'variables' is NULL.")
  } else if(!is.null(betaNull) 
            && length(betaNull) != length(variables)){
    cli::cli_abort("Length of 'betaNull' must match length of 'variables'.")
  }
  
  # Check if `model` provided and is formula with LHS
  if(is.null(model)){
    stop("The 'model' parameter must be provided.")
  } else if(!rlang::is_formula(model, lhs = TRUE)){
    stop("The 'model' parameter must be a formula with a LHS.")
  }
  
  # Check if length of permutation equals number of observations
  if(!is.null(permutation) && length(permutation) != nrow(data)){
    cli::cli_abort("The length of 'permutation' does not equal the number of observations in 'data'.")
  }
  
  data.n <- nrow(data)
  n <- floor(data.n/nBlocks)*nBlocks
  
  # Construct matrix of block indices
  blockSize <- n/nBlocks
  blockIndexMatrix <- matrix(1:n, 
                             nrow = blockSize, 
                             ncol = nBlocks, 
                             byrow = FALSE)
  
  # Check if length of permutation equals number of observations
  if(length(permutation) > n){
    cli::cli_warn("After dropping remainder observations, 'permutation' is longer than the number of observations in the data. Remainder observations will be dropped from 'permutation' as well.")
    permutation <- permutation[permutation %in% 1:n]
  }
  
  ivregObject <- ivreg::ivreg(model,
                              data = data,
                              model = TRUE,
                              x = TRUE)
  
  regressors <- as.character(unlist(attr(ivregObject$terms$regressors, "variables")))[-1]
  Y.var <- regressors[1]
  X.var <- regressors[-1]
  
  endogenous.var <- names(ivregObject$endogenous)
  exogenous.var <- names(ivregObject$exogenous)

  if(attr(ivregObject$terms$regressors, "intercept")){
    exogenous.var <- exogenous.var[-1]
  }
  
  Y <- matrix(ivregObject$y)
  X <- ivregObject$x$regressors
  
  Y.use <- Y[1:n,, drop = FALSE]
  X.use <- X[1:n,, drop = FALSE]
  
  Z.var <- names(ivregObject$instruments)
  
  if(!is.null(Z.var)){
    IV <- TRUE
    Z <- ivregObject$x$instruments[,Z.var]
    Z.use <- Z[1:n,, drop = FALSE]
  } else{
    IV <- FALSE
  }
  
  assign <- attr(X, "assign")
  
  summaryTableIvreg <- summary(ivregObject)$coefficients
  
  gaArgs = list(...)
  
  if(is.null(variables)){
    variables <- 1:length(X.var)
  }
  
  # If `nPerms` is unspecified or greater than the number of possible permutations, then use all possible permutations. 
  # When number of possible permutations is bigger than MG, then we need to randomly sample.
  if(is.null(nPerms) || nPerms >= factorial(nBlocks)){
    blockPermutations <- do.call(rbind, combinat::permn(1:nBlocks))
    permIndices <- apply(blockPermutations, MARGIN = 1, function (x) c(blockIndexMatrix[, x]))
  } else{
    permIndices <- cbind(1:n, replicate(nPerms, c(blockIndexMatrix[, sample(1:nBlocks)])))
  }
  
  final_results <- vector("list")
  
  if(!optimize){ # Case 1: don't optimize
    if(!is.null(permutation)){
      X.use <- X.use[permutation,, drop = FALSE]
      Y.use <- Y.use[permutation,, drop = FALSE]
      if(IV){
        Z.use <- Z.use[permutation,, drop = FALSE]
      }
    }
  } else{ # Case 2: optimize
    
    if(!is.null(permutation)){ 
      cli::cli_warn("Optimization will not be implemented because 'permutation' specified. If optimization is desired, either remove 'permutation' or include 'permutation' in matrix of suggestions and pass to function.")
      X.use <- X.use[permutation,, drop = FALSE]
      Y.use <- Y.use[permutation,, drop = FALSE]
      if(IV){
        Z.use <- Z.use[permutation,, drop = FALSE]
      }
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
      gaArgs$fitness <- function(permutation){ fitness_function(permutation, X1.temp, X2.temp, Z.temp, blockIndexMatrix, blockPermutations) }
      gaArgs$lower <- rep(1, n)
      gaArgs$upper <- rep(n, n)
    }
  }
  
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
  
  if(optimize && ncol(X) > 1){
    final_results[["gaResults"]] <- vector("list")
  }
  
  rowCounter <- 1
  
  for(i in seq_along(attr(X, "assign"))){
    
    if(assign[i] == 0 | !assign[i] %in% variables){
      next
    } 
    
    exacttIV <- !colnames(X)[i] %in% exogenous.var
    
    beta_hat <- summaryTableIvreg[i, 1]
    se <- summaryTableIvreg[i, 2]
    precisionToUse <- ifelse(se > 0, yes = floor(log(se, base = 10)) - 1, no = -5)
    
    Y.temp <- as.matrix(Y.use)
    X1.temp <- X.use[,i, drop = FALSE]
    X2.temp <- X.use[,-i, drop = FALSE]
    
    if(exacttIV){
      Z.temp <- Z.use
    }
    
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
        
        if(!exacttIV){
          Z.temp <- NULL
        } 
        
        parallel::clusterExport(cl, varlist = c("X1.temp", "X2.temp", "Z.temp", "blockIndexMatrix", "blockPermutations", "fitness_function", "build_GX2", "build_GX", "block_permute", "build_QGX2"), envir = environment())
        
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
      
      Y.temp <- Y.temp[gaResults@solution,, drop = FALSE]
      X1.temp <- X1.temp[gaResults@solution,, drop = FALSE]
      X2.temp <- X2.temp[gaResults@solution,, drop = FALSE]
      
      if(exacttIV){
        Z.temp <- Z.temp[gaResults@solution,, drop = FALSE]
      }
      
      final_results[["gaResults"]][[colnames(X)[i]]] <- gaResults
    }
    
    if(!is.null(betaNull) && !is.null(betaNull[[i-1]])){
      betaNullVec <- betaNull[[i-1]]
    } else{
      # Find LB and UB roots
      if(exacttIV){
        betaNullVec <- getBetaNull(Y.temp, X1.temp, X2.temp, Z.temp, alpha, nBlocks, permIndices, beta_hat, se, precisionToUse)
      } else{
        betaNullVec <- getBetaNull(Y.temp, X1.temp, X2.temp, Z.temp = NULL, alpha, nBlocks, permIndices, beta_hat, se, precisionToUse)
      }
    }
    
    if(exacttIV){
      exacttResults <- exactt_pval_IV(betaNullVec, 
                                      Y.temp, 
                                      X1.temp,
                                      X2.temp, 
                                      Z.temp,
                                      nBlocks, 
                                      permIndices,
                                      studentize)
      
    } else{
      exacttResults <- exactt_pval(betaNullVec, 
                                   Y.temp, 
                                   X1.temp,
                                   X2.temp, 
                                   nBlocks, 
                                   permIndices,
                                   studentize,
                                   GX1)
    }

    final_results[["detailed"]][[colnames(X)[i]]] <- data.frame("betaNull" = betaNullVec, 
                                                                "pval" = exacttResults$pval)
    
    #"randomizationDistribution" = list(c(exacttResults$randomizationDist)))
    
    
    pvalBetaNull0 <- ifelse(length(exacttResults$pval[which(betaNullVec == 0)]) == 0, 
                            yes = NA,
                            no = exacttResults$pval[which(betaNullVec == 0)])
    
    summaryTable.out[rowCounter,] <- c(beta_hat,
                                       pvalBetaNull0,
                                       ciByInversion(betaNullVec, exacttResults$pval, alpha, weighted = TRUE),
                                       ciByInversion(betaNullVec, exacttResults$pval, alpha, weighted = FALSE))
    
    rowCounter <- rowCounter + 1

  }
  
  final_results[["summary"]] <- summaryTable.out
  
  final_results[["call"]] <- call
  
  class(final_results) <- "et"
  
  return(final_results) 
}
