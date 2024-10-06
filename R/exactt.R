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
#' @param side A character to indicate the side of the test.
#' @param alpha The significance level used for the hypothesis tests; defaults to 0.05.
#' @param variables Optional; a character vector of predictor names to test.
#'        If NULL, all predictors in the model are tested.
#' @param beta0 Optional; a numeric vector of null hypothesis values for the coefficients.
#'        Must be the same length as `variables` if not NULL.
#' @param nBlocks The number of blocks to use for block permutations.
#' @param nPerms Optional; the number of permutations to perform.
#'        If NULL or greater than the number of possible permutations, all permutations are used.
#' @param studentize Logical indicating whether to use studentized residuals for the test.
# @param precisionToUse The precision level for rounding off the beta values in the
#        resulting sequence.
# @param randomizationDist Logical indicating whether to return randomization distribution for each null hypothesis value.
#' @param optimize Logical indicating whether to optimize the ordering of the data.
# @param GX1 Logical indicating whether to use GX1 or X1 when constructing eps_hat. 
#        Using X1 has slightly more power at slightly more computational expense.
#' @param seed Seed used when optimizing using `GA::ga()`. Default is 31740.
#' @param Q.X1 Optional argument used for power testing purposes.
#' @param denominator Character argument indicating how to calculate epsilon hat.
#' @param ... Additional arguments passed to `GA::ga()` for optimizing power. 
#' This can include parameters like `popSize`, `maxiter`, `parallel`, etc., 
#' that are used to configure the genetic algorithm. Note that when sample size is large
#' optimizing is computationally expensive and has little effect.
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
#' @importFrom cli cli_abort cli_alert_info
#' @importFrom doRNG registerDoRNG
#' @importFrom GA gaControl
#' @importFrom combinat permn
#' 
#' @export
exactt <- function(model,
                   data,
                   side = "both",
                   alpha = 0.05,
                   variables = NULL,
                   beta0 = NULL,
                   nBlocks = 5,
                   nPerms = NULL,
                   studentize = TRUE,
                   optimize = FALSE,
                   seed = 31740,
                   Q.X1 = NULL,
                   denominator = "GX1",
                   ...) {
  
  call <- match.call(expand.dots = TRUE)
  
  # Evaluate the arguments
  call$alpha <- eval(call$alpha, envir = parent.frame())
  
  ####### Do checks #######
  
  # Check if `model` provided and is formula with LHS
  if(is.null(model)){
    stop("The 'model' parameter must be provided.")
  } else if(!rlang::is_formula(model, lhs = TRUE)){
    stop("The 'model' parameter must be a formula with a LHS.")
  }
  
  ivregObject <- ivreg::ivreg(model,
                              data = data,
                              model = TRUE,
                              x = TRUE)
  
  # Check if `beta0` provided. If so, check if length equal to number of variables
  if(!is.null(beta0) && is.null(variables)){
    warning("'beta0' will be ignored since 'variables' is NULL.")
  } else if(!is.null(beta0) 
            && length(beta0) != length(variables)){
    cli::cli_abort("Length of 'beta0' must match length of 'variables'.")
  }
  
  data <- ivregObject$model
  
  data.n <- nrow(data)
  n <- floor(data.n/nBlocks)*nBlocks
  
  # Construct matrix of block indices
  blockSize <- n/nBlocks
  blockIndexMatrix <- matrix(1:n, 
                             nrow = blockSize, 
                             ncol = nBlocks, 
                             byrow = FALSE)
  
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
  
  gaArgs <- list(seed = seed, ...)
  
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
  
  GX.indices <- build_GX(blockIndexMatrix)
  
  if(optimize){ # Case 1: don't optimize
    if("type" %in% names(gaArgs)){
      cli::cli_warn("Custom 'type' value is ignored in this function.")
      gaArgs$type <- NULL
    } 
    if("fitness" %in% names(gaArgs)){
      cli::cli_warn("Custom 'fitness' value is ignored in this function.")
      gaArgs$fitness <- NULL
    } 
    if ("lower" %in% names(gaArgs) || "upper" %in% names(gaArgs)) {
      cli::cli_warn("Custom 'lower' and 'upper' values are ignored in this function.")
      gaArgs$lower <- NULL
      gaArgs$upper <- NULL
    }
    if("crossover" %in% names(gaArgs) && gaArgs$crossover != "gaperm_oxCrossover_R"){
      cli::cli_warn("'crossover' is restricted to 'gaperm_oxCrossover_R' due to Rcpp issues.")
    }
    
    gaArgs$type <- "permutation"
    gaArgs$fitness <- function(permutation){ fitness_function(permutation, X1.temp, X2.temp, Z.temp, blockIndexMatrix, permIndices, GX.indices, blockPermutations) }
    gaArgs$lower <- rep(1, n)
    gaArgs$upper <- rep(n, n)
    gaArgs$crossover = "gaperm_oxCrossover_R"
  }
  
  summaryTableList <- vector("list")
  detailedList <- vector("list")
  gaResultsList <- vector("list")
  
  for(i in seq_along(attr(X, "assign"))){
    
    if(assign[i] == 0 | !assign[i] %in% variables){
      next
    } 
    
    exacttIV <- !colnames(X)[i] %in% exogenous.var
    
    beta_hat <- summaryTableIvreg[i, 1]
    
    Y.temp <- as.matrix(Y.use)
    X1.temp <- X.use[,i, drop = FALSE]
    X2.temp <- X.use[,-i, drop = FALSE]
    
    if(exacttIV){
      Z.temp <- Z.use
    }
    
    if(optimize){
      
      if(!exacttIV){
        Z.temp <- NULL
      } 
      
      if(!is.null(gaArgs$parallel) && gaArgs$parallel == TRUE){
        
        ogParArg <- gaArgs$parallel
        
        if(gaArgs$parallel == TRUE){
          numCores <- parallel::detectCores()
        } else if(is.numeric(gaArgs$parallel) && gaArgs$parallel >= 2){
          numCores <- gaArgs$parallel
        }
        
        # Create the appropriate cluster
        if (.Platform$OS.type == "windows") {
          # Use socket cluster on Windows or if forking is not desired
          cl <- parallel::makeCluster(numCores, type = "PSOCK")
          
          # Export variables and functions only if using a socket cluster
          parallel::clusterExport(cl, varlist = c("X1.temp", 
                                                  "X2.temp", 
                                                  "Z.temp", 
                                                  "blockIndexMatrix", 
                                                  "permIndices", 
                                                  "GX.indices", 
                                                  "blockPermutations",
                                                  "fitness_function", 
                                                  "build_GX2", 
                                                  "build_GX", 
                                                  "block_permute", 
                                                  "build_QGX2"), 
                                  envir = environment())
          parallel::clusterCall(cl, library, package = "Matrix", character.only = TRUE)
          parallel::clusterCall(cl, library, package = "MASS", character.only = TRUE)
          parallel::clusterCall(cl, library, package = "combinat", character.only = TRUE)
          parallel::clusterCall(cl, library, package = "dplyr", character.only = TRUE)
        } else {
          # Unix-based system and forking is enabled
          cl <- parallel::makeCluster(numCores, type = "FORK")
        }
        
        # Register the parallel backend
        doParallel::registerDoParallel(cl, cores = numCores)
       
        gaArgs$parallel <- cl
      } else{
        ogParArg <- FALSE
      }
      
      cli::cli_alert_success("Optimizing ordering for `{colnames(X)[i]}`.")
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
      
      gaResultsList[[colnames(X)[i]]] <- gaResults
    }
    
    GX2.temp <- build_GX2(X2.temp, GX.indices)
    #GX2.reduced.temp <- remove_dependent_columns(GX2.temp)
    QGX2.temp <- build_QGX2(GX2.temp)
    if(exacttIV){
      Q.Z.temp <- QGX2.temp %*% Z.temp
    } else{
      if(is.null(Q.X1)){
        Q.X1.temp <- QGX2.temp %*% X1.temp
      } else{
        Q.X1.temp <- Q.X1
      }
    }
    
    if(studentize){
      QGX1GX2.temp <- build_QGX1GX2(X1.temp, GX2.temp, GX.indices, denominator)
    } else{
      QGX1GX2.temp <- NULL
    }
    
    if(TRUE){
      if(exacttIV){
        pvals.df <- exactt.pval.new.iv(Y.temp, X1.temp, permIndices, Q.Z.temp, QGX1GX2.temp)
      } else{
        pvals.df <- exactt.pval.new.reg(Y.temp, X1.temp, GX2.temp, permIndices, GX.indices, Q.X1.temp, QGX2.temp, QGX1GX2.temp, side = side, denominator)
      }
      
      attr(pvals.df, "assign") = assign[i]
      detailedList[[colnames(X)[i]]] <- pvals.df
      
      pvalBeta0.index <- which(0 >= pvals.df$beta0.start & 0 <= pvals.df$beta0.end)
      
      if(side == "both"){
        ci.lower.index <- min(which(pvals.df$pvals > alpha))
        ci.upper.index <- max(which(pvals.df$pvals > alpha))
      } else if (side == "left"){
        ci.lower.index <- 1
        ci.upper.index <- max(which(pvals.df$pvals > alpha))
      } else{
        ci.lower.index <- min(which(pvals.df$pvals > alpha))
        ci.upper.index <- ncol(permIndices)
      }
      
      summaryTableList[[i]] <- matrix(data = c(beta_hat,
                                               pvals.df$pvals[pvalBeta0.index],
                                               pvals.df$beta0.start[ci.lower.index], 
                                               pvals.df$beta0.end[ci.upper.index]),
                                      nrow = 1, 
                                      ncol = 4, 
                                      dimnames = list(colnames(X)[i], 
                                                      c("Estimate", 
                                                        "P-value",
                                                        "Lower Bound",
                                                        "Upper Bound")
                                                      ))
    } 
  }

  result <- structure(list(call = call,
                           summary = do.call('rbind', summaryTableList),
                           detailed = detailedList,
                           gaResults = gaResultsList,
                           ivregResults = ivregObject),
                      class = "exactt")
  
  if(length(gaResultsList) > 0){
    result$gaResults <- gaResultsList
  } 
  
  if(exacttIV){
    result$Q.Z <- Q.Z.temp
  } else{
    result$Q.X1 <- Q.X1.temp
  }
  
  return(result) 
}
