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
                   ...) {
  
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
  
  # Parse the formula using the Formula package
  formula <- Formula::Formula(model)
  
  # Extract Y (response), X (explanatory variables), and Z (instrumental variables)
  Y.var <- as.character(attr(stats::terms(formula), "variables")[[2]])
  X.var <- all.vars(stats::formula(formula, lhs = 0, rhs = 1))
  suppressWarnings(Z.var <- all.vars(formula(formula, lhs = 0, rhs = 2)))
  
  IV <- ifelse(length(Z.var) > 0, TRUE, FALSE)
  
  # Assume !IV for now
  
  # Convert the data to a model matrix, which will handle one-hot encoding
  X <- stats::model.matrix(formula, data)
  assign <- attr(X, "assign")
  
  Y <- data[Y.var]
  
  data.n <- nrow(X)
  n <- floor(data.n/nBlocks)*nBlocks
  
  X.use <- X[1:n, , drop = FALSE]
  Y.use <- Y[1:n, , drop = FALSE]
  data.use <- data[1:n, , drop = FALSE]
  
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
  
  if(!optimize){ # Case 1: don't optimize
    if(!is.null(permutation)){
      X.use <- X.use[permutation,, drop = FALSE]
      Y.use <- Y.use[permutation,, drop = FALSE]
    }
  } else{ # Case 2: optimize
    # Convert '...' to a list to inspect and possibly modify
    gaArgs <- list(...)
    
    if(!is.null(permutation)){ 
      warning("Optimization will not be implemented because 'permutation' specified. If optimization is desired, either remove 'permutation' or include 'permutation' in matrix of suggestions and pass to function.")
      X.use <- X.use[permutation,, drop = FALSE]
      Y.use <- Y.use[permutation,, drop = FALSE]
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
    
    gaArgs$type <- "permutation"
    gaArgs$fitness <- function(permutation){ fitness_function(permutation, X1.temp, X2.temp, blockIndexMatrix, blockPermutations) }
    gaArgs$lower <- rep(1, n)
    gaArgs$upper <- rep(n, n)
  }

  if(IV){
    Z <- data.use[,Z.var, drop = FALSE]
  } else{
    lmObject <- stats::lm(model, data = data)
  }
  
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
    
    summaryTable.out[rowCounter,] <- c(beta_hat,
                                       detailed.out[[i]]$pval[idx0],
                                       beta_lower,
                                       beta_upper,
                                       max(betaNullVec[detailed.out[[i]]$pval <= alpha & (betaNullVec <= stats::median(betaNullVec))]),
                                       min(betaNullVec[detailed.out[[i]]$pval <= alpha & (betaNullVec >= stats::median(betaNullVec))]))
    
    rowCounter <- rowCounter + 1
  }
  
  final_results[["summary"]] <- summaryTable.out
  final_results[["call"]] <- match.call()
  
  class(final_results) <- "et"
  
  return(final_results) 
}
