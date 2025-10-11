#' Exact Testing of Linear Model Coefficients
#'
#' Performs exact tests on specified coefficients of a linear model object
#' using permutation tests to generate p-values and confidence intervals around
#' the estimated coefficients. This function can handle both studentized and
#' non-studentized test statistics, and allows the user to specify various
#' parameters for the test.
#'
#' @param model A formula specifying the model.
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
#' @param optimize Logical indicating whether to optimize the ordering of the data.
#' @param seed Seed used when optimizing using `GA::ga()`. Default is 31740.
#' @param denominator Character argument indicating how to calculate epsilon hat.
#' @param Q.X1 Use custom QX1 value.
#' @param GX.indices Indices for max rank GX. Used for warm start.
#' @param root.tolerance Tolerance for determining real and extraneous roots (when denominator = "X1" or "noX1").
#' @param ... Additional arguments passed to `GA::ga()` for optimizing power. 
#' This can include parameters like `popSize`, `maxiter`, `parallel`, etc., 
#' that are used to configure the genetic algorithm. Note that when sample size is large
#' optimizing is computationally expensive and has little effect.
#'
#' @return The p-value of the test that the null values of beta are 0.
#'
#' @details
#' The function divides the data into blocks specified by `nBlocks` and performs permutations
#' within across blocks to generate the null distribution of the test statistic. The user can
#' specify a set number of permutations with `nPerms`, or allow the function to calculate all
#' possible permutations if `nPerms` is unspecified or too large.
#'
#' If `studentize` is TRUE, studentized residuals are used to adjust the test statistics,
#' potentially leading to more robust inference under model misspecification.
#'
#' The function allows for a high degree of customization through its parameters and can
#' handle large datasets and complex model structures efficiently.
#'
#' @importFrom stats median formula model.matrix lm coef
#' @importFrom Formula Formula
#' @importFrom cli cli_abort cli_alert_info
#' @importFrom doRNG registerDoRNG
#' @importFrom GA gaControl
#' @importFrom combinat permn
#' @importFrom utils tail
#' 
#' @export
exactt <- function(model,
                   data,
                   side = "both",
                   alpha = 0.05,
                   variables = NULL,
                   beta0 = NULL,
                   nBlocks = NULL,
                   nPerms = NULL,
                   studentize = TRUE,
                   optimize = FALSE,
                   seed = 31740,
                   denominator = "GX1",
                   Q.X1 = NULL,
                   GX.indices = NULL,
                   root.tolerance = 1e-9,
                   ...) {
  
  call <- match.call(expand.dots = TRUE)
  
  # Evaluate the arguments
  call$alpha <- eval(call$alpha, envir = parent.frame())
  call$model <- eval(call$model, envir = parent.frame())
  call$side <- eval(call$side, envir = parent.frame())
  
  ####### Do checks #######
  
  # Check if `model` provided and is formula with LHS
  if(is.null(model)){
    stop("The 'model' parameter must be provided.")
  } else if(!rlang::is_formula(model, lhs = TRUE)){
    stop("The 'model' parameter must be a formula with a LHS.")
  }
  
  if(!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1){
    stop("The 'alpha' parameter must be numeric, have length one, and between 0 and 1.")
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
  
  endogenous.var <- names(ivregObject$endogenous)
  exogenous.var <- names(ivregObject$exogenous)

  if(attr(ivregObject$terms$regressors, "intercept")){
    exogenous.var <- exogenous.var[-1]
  }
  
  Y.use <- matrix(ivregObject$y)
  X.use <- ivregObject$x$regressors
  
  Z.var <- names(ivregObject$instruments)
  
  if(!is.null(Z.var)){
    IV <- TRUE
    Z.use <- ivregObject$x$instruments[,Z.var]
  } else{
    IV <- FALSE
  }
  
  X.assign <- attr(X.use, "assign")
  
  beta.hats <- coef(ivregObject)
  
  if(is.null(nBlocks)){
    p2 <- ncol(X.use) - 1  
    
    if(data.n < 400*sqrt(p2)){
      nBlocks <- 5
    } else if(data.n < 1000*sqrt(p2)){
      nBlocks <- 6
    } else if(data.n < 1850*sqrt(p2)){
      nBlocks <- 7
    } else if(data.n < 3000*sqrt(p2)){
      nBlocks <- 8
    } else if(data.n < 4500*sqrt(p2)){
      nBlocks <- 9
    } else if(data.n < 6000*sqrt(p2)){
      nBlocks <- 10
    } else{
      nBlocks <- 11
    }
  }
  
  n <- floor(data.n/nBlocks)*nBlocks
  n.remainder.indices <- utils::tail(1:data.n, data.n - n)
  
  # Construct matrix of block indices
  blockSize <- n/nBlocks
  blockIndexMatrix <- matrix(1:n,
                             nrow = blockSize,
                             ncol = nBlocks,
                             byrow = FALSE)
  
  X.var <- as.character(unlist(attr(ivregObject$terms$regressors, "variables")))[-(1:2)]
  
  gaArgs <- list(seed = seed, ...)
  
  if(is.null(variables)){
    variables.construct <- 1:length(X.var)
  } else if(is.list(variables)){
    variables.construct <- sapply(variables, 
                                  FUN = function(x) x[[2]][[2]])
  } else{
    variables.construct <- variables
  }
  
  # If `nPerms` is unspecified or greater than the number of possible permutations, then use all possible permutations. 
  # When number of possible permutations is bigger than MG, then we need to randomly sample.
  if(is.null(nPerms) || nPerms >= factorial(nBlocks)){
    blockPermutations <- do.call(rbind, combinat::permn(1:nBlocks))
    permIndices <- apply(blockPermutations, 
                         MARGIN = 1, 
                         function(x){
                           c(blockIndexMatrix[, x], n.remainder.indices)
                         })
    
    nPerms <- factorial(nBlocks)
  } else{
    permIndices <- cbind(1:n, 
                         unique(replicate(nPerms, c(blockIndexMatrix[, sample(1:nBlocks)])),
                                MARGIN = 2))
  }
  
  
  # If studentize is true, construct GX.indices
  # IF there are 3 or more columns in X.use, construct X.use
  # If there are 2 columns in X.use & neither of the column names is equal to "(Intercept)", construct X.use
  if(studentize == TRUE ||
     ncol(X.use) >= 3 || 
     (ncol(X.use) == 2 && !"(Intercept)" %in% colnames(X.use))){
    
    if(!is.null(GX.indices)){
      if(ncol(GX.indices) != nBlocks * (nBlocks - 2) + 2){
        cli::cli_warn("Number of columns in `GX.indices` is not nBlocks * (nBlocks - 2) + 2. Recalculating GX.indices now.")
        GX.indices <- build_GX.indices(blockIndexMatrix, n.remainder.indices)
      } else if(nrow(GX.indices) != data.n){
        cli::cli_warn("Number of rows in `GX.indices` does not match number of rows in data. Recalculating GX.indices now.")
        GX.indices <- build_GX.indices(blockIndexMatrix, n.remainder.indices)
      }
    } else{
      GX.indices <- build_GX.indices(blockIndexMatrix, n.remainder.indices)
    }
    
  }

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
    gaArgs$fitness <- function(permutation){ fitness_function(permutation = permutation, 
                                                              X1.temp = X1.temp, 
                                                              X2.temp = X2.temp, 
                                                              Z.temp = Z.temp, 
                                                              indep.X2.index = indep.X2.index,
                                                              blockIndexMatrix = blockIndexMatrix, 
                                                              GX.indices = GX.indices, 
                                                              permIndices = permIndices) }
    gaArgs$lower <- rep(1, data.n)
    gaArgs$upper <- rep(data.n, data.n)
    gaArgs$crossover = "gaperm_oxCrossover_R"
  }
  
  summaryTableList <- vector("list")
  detailedList <- vector("list")
  geometryList <- vector("list")
  gaResultsList <- vector("list")
  Q.X1.Z.List <- vector("list")
  
  for(i in seq_along(X.assign)){
    
    if(X.assign[i] == 0 | !X.assign[i] %in% variables.construct){
      next
    } 
    
    exacttIV <- !colnames(X.use)[i] %in% exogenous.var
    
    beta.hats.i <- beta.hats[i]
    
    Y.temp <- as.matrix(Y.use)
    X1.temp <- X.use[,i, drop = FALSE]
    X2.temp <- X.use[,-i, drop = FALSE]
    attr(X2.temp, "assign") <- attr(X.use, "assign")[-i]
    
    if(is.list(variables)){
      variables.index <- which(variables.construct == X.assign[i])
      
      split_formula <- strsplit(as.character(variables[[variables.index]])[[2]], "\\|")[[1]]
      
      indep.X2.index <- strsplit(split_formula[2], "\\+")[[1]] |>
        trimws() |>
        as.numeric()
    } else{
      indep.X2.index <- NULL
    }
    
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
                                                  "indep.X2.index",
                                                  "blockIndexMatrix", 
                                                  "permIndices", 
                                                  "GX.indices", 
                                                  "blockPermutations",
                                                  "n",
                                                  "fitness_function", 
                                                  "build_GX", 
                                                  "build_GX.indices"), 
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
      
      cli::cli_alert_success("Optimizing ordering for `{colnames(X.use)[i]}`.")
      gaResults <- do.call(GA::ga, gaArgs)
      
      # Close cluster if parallel is true
      if(!is.null(gaArgs$parallel) && ogParArg != FALSE){
        parallel::stopCluster(cl)
        gaArgs$parallel <- ogParArg
      }
      
      Y.temp <- Y.temp[gaResults@solution[1,],, drop = FALSE]
      X1.temp <- X1.temp[gaResults@solution[1,],, drop = FALSE]
      X2.temp <- X2.temp[gaResults@solution[1,],, drop = FALSE]
      
      if(exacttIV){
        Z.temp <- Z.temp[gaResults@solution[1,],, drop = FALSE]
      }
      
      gaResultsList[[colnames(X.use)[i]]] <- gaResults
    }
    
    if(!is.null(Q.X1)){
      if(exacttIV){
        Q.Z.temp <- Q.X1
      } else{
        Q.X1.temp <- Q.X1
      }
    } else{
      if(exacttIV){
        if(ncol(X2.temp) == 0){
          Q.Z.temp <- Z.temp
        } else if(ncol(X2.temp) == 1 && 
                  all(X2.temp == 1)){
          Q.Z.temp <- Z.temp - mean(Z.temp)
        } else{
          Q.Z.temp <- stats::lm(Z.temp ~ 0 + build_GX(X2.temp, GX.indices, indep.X2.index),
                                model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals |>
            matrix(nrow = nrow(Z.temp))
        }
        Q.X1.Z.List[[colnames(X.use)[i]]] <- Q.Z.temp
      } else{
        if(ncol(X2.temp) == 0){
          Q.X1.temp <- X1.temp
        } else if(ncol(X2.temp) == 1 && 
                  all(X2.temp == 1)){
          Q.X1.temp <- X1.temp - mean(X1.temp)
        } else{
          Q.X1.temp <- stats::lm(X1.temp ~ 0 + build_GX(X2.temp, GX.indices, indep.X2.index),
                                 model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals |>
            matrix(ncol = ncol(X1.temp))
        }
        Q.X1.Z.List[[colnames(X.use)[i]]] <- Q.X1.temp
      }
    }
    
    if(exacttIV){
      pvals.df <- exactt.pval.new.iv(Y.temp, X1.temp, X2.temp, indep.X2.index, permIndices, GX.indices, Q.Z.temp, studentize)
    } else{
      exactt.pval.reg <- exactt.pval.new.reg(Y.temp, X1.temp, X2.temp, indep.X2.index, permIndices, GX.indices, Q.X1.temp, studentize, side = side, denominator, root.tolerance)
      
      pvals.df <- exactt.pval.reg$pvals.df
      geometryList[[colnames(X.use)[i]]] <- list(line.data.num = exactt.pval.reg$line.data.num,
                                                 line.data.denom.sq = exactt.pval.reg$line.data.denom.sq)
    }
    
    attr(pvals.df, "assign") = X.assign[i]
    detailedList[[colnames(X.use)[i]]] <- pvals.df
    
    pvalBeta0.index <- which(0 >= pvals.df$beta0.start & 0 <= pvals.df$beta0.end)
    
    ci.lower.index <- min(which(pvals.df$pvals > alpha))
    ci.upper.index <- max(which(pvals.df$pvals > alpha))
    
    summaryTableList[[i]] <- matrix(data = c(beta.hats.i,
                                             max(pvals.df$pvals[pvalBeta0.index]),
                                             pvals.df$beta0.start[ci.lower.index], 
                                             pvals.df$beta0.end[ci.upper.index]),
                                    nrow = 1, 
                                    ncol = 4, 
                                    dimnames = list(colnames(X.use)[i], 
                                                    c("Estimate", 
                                                      "P-value",
                                                      "Lower Bound",
                                                      "Upper Bound"))
                                    )
  }

  result <- structure(list(call = call,
                           summary = do.call('rbind', summaryTableList),
                           nBlocks = nBlocks,
                           detailed = detailedList,
                           gaResults = gaResultsList,
                           ivregResults = ivregObject,
                           geometry = geometryList),
                      class = "exactt")
  
  if(length(gaResultsList) > 0){
    result$gaResults <- gaResultsList
  } 
  
  if(!is.null(GX.indices)){
    result$GX.indices <- GX.indices
  }
  
  if(exacttIV){
    result$Q.Z <- Q.X1.Z.List
  } else{
    result$Q.X1 <- Q.X1.Z.List
  }
  
  return(result) 
}

#' Exact Wald-test
#'
#' @param model A formula specifying the model.
#' @param data A data frame or matrix containing the variables used in the model.
#' @param alpha The significance level used for the hypothesis tests; defaults to 0.05.
#' @param variables Indices of variables of interest.
#' @param beta0 If 0, test whether the variables of interest are equal to the zero vector.
#'        If NULL, creates a grid of beta0 values (for confidence intervals) and tests the variables of interest.
#' @param nBlocks The number of blocks to use for block permutations.
#' @param nPerms Optional; the number of permutations to perform.
#'        If NULL or greater than the number of possible permutations, all permutations are used.
#' @param studentize Logical indicating whether to use studentized residuals for the test.
#' @param optimize Logical indicating whether to optimize the ordering of the data.
#' @param seed Seed used when optimizing using `GA::ga()`. Default is 31740.
#' @param GX.indices Indices for max rank GX. Used for warm start.
#' @param confidence.intervals Do you want exact confidence intervals?
#' @param gurobi.params Parameters to pass through to `gurobi::gurobi()`.
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
#' within across blocks to generate the null distribution of the test statistic. The user can
#' specify a set number of permutations with `nPerms`, or allow the function to calculate all
#' possible permutations if `nPerms` is unspecified or too large.
#'
#' If `studentize` is TRUE, studentized residuals are used to adjust the test statistics,
#' potentially leading to more robust inference under model misspecification.
#'
#' The function allows for a high degree of customization through its parameters and can
#' handle large datasets and complex model structures efficiently.
#'
#' @importFrom stats median formula model.matrix lm
#' @importFrom Formula Formula
#' @importFrom cli cli_abort cli_alert_info
#' @importFrom combinat permn
#' @importFrom utils tail
#' 
#' @export
exactt.wald <- function(model,
                        data,
                        alpha = 0.05,
                        variables = NULL,
                        beta0 = c(0, NULL),
                        nBlocks = NULL,
                        nPerms = NULL,
                        studentize = TRUE,
                        optimize = FALSE,
                        seed = 31740,
                        confidence.intervals = TRUE,
                        gurobi.params = list(FeasibilityTol = 1e-9,
                                             OptimalityTol = 1e-9,
                                             IntFeasTol = 1e-9,
                                             OutputFlag = 0),
                        GX.indices = NULL,
                        ...) {

  call <- match.call(expand.dots = TRUE)

  # Evaluate the arguments
  call$alpha <- eval(call$alpha, envir = parent.frame())
  call$model <- eval(call$model, envir = parent.frame())
  
  ####### Do checks #######

  # Check if `model` provided and is formula with LHS
  if(is.null(model)){
    stop("The 'model' parameter must be provided.")
  } else if(!rlang::is_formula(model, lhs = TRUE)){
    stop("The 'model' parameter must be a formula with a LHS.")
  }
  
  if(!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1){
    stop("The 'alpha' parameter must be numeric, have length one, and between 0 and 1.")
  }

  ivregObject <- ivreg::ivreg(model,
                              data = data,
                              model = TRUE,
                              x = TRUE)

  # Check if `beta0` provided. If so, check if length equal to number of variables
  # if(!is.null(beta0) && is.null(variables)){
  #   warning("'beta0' will be ignored since 'variables' is NULL.")
  # } else if(!is.null(beta0)
  #           && length(beta0) != length(variables)){
  #   cli::cli_abort("Length of 'beta0' must match length of 'variables'.")
  # }

  data <- ivregObject$model

  data.n <- nrow(data)

  endogenous.var <- names(ivregObject$endogenous)
  exogenous.var <- names(ivregObject$exogenous)

  if(attr(ivregObject$terms$regressors, "intercept")){
    exogenous.var <- exogenous.var[-1]
  }

  Y.use <- matrix(ivregObject$y)
  X.use <- ivregObject$x$regressors

  Z.var <- names(ivregObject$instruments)

  if(!is.null(Z.var)){
    IV <- TRUE
    Z.use <- ivregObject$x$instruments[,Z.var]
  } else{
    IV <- FALSE
  }

  X.assign <- attr(X.use, "assign")

  beta.hats <- coef(ivregObject)
  se <- sqrt(diag(ivregObject$cov.unscaled * ivregObject$sigma^2))
  
  if(is.null(variables)){
    variables <- X.assign[which(X.assign != 0)]
  }
  
  if(is.null(nBlocks)){
    p2 <- ncol(X.use) - length(variables)  
    
    if(data.n < 400*sqrt(p2)){
      nBlocks <- 5
    } else if(data.n < 1000*sqrt(p2)){
      nBlocks <- 6
    } else if(data.n < 1850*sqrt(p2)){
      nBlocks <- 7
    } else if(data.n < 3000*sqrt(p2)){
      nBlocks <- 8
    } else if(data.n < 4500*sqrt(p2)){
      nBlocks <- 9
    } else if(data.n < 6000*sqrt(p2)){
      nBlocks <- 10
    } else{
      nBlocks <- 11
    }
  }
  
  n <- floor(data.n/nBlocks)*nBlocks
  n.remainder.indices <- utils::tail(1:data.n, data.n - n)
  
  # Construct matrix of block indices
  blockSize <- n/nBlocks
  blockIndexMatrix <- matrix(1:n,
                             nrow = blockSize,
                             ncol = nBlocks,
                             byrow = FALSE)
  
  gaArgs <- list(seed = seed, ...)

  # If `nPerms` is unspecified or greater than the number of possible permutations, then use all possible permutations.
  # When number of possible permutations is bigger than MG, then we need to randomly sample.
  if(is.null(nPerms) || nPerms >= factorial(nBlocks)){
    blockPermutations <- do.call(rbind, combinat::permn(1:nBlocks))
    permIndices <- apply(blockPermutations,
                         MARGIN = 1,
                         function(x){
                           c(blockIndexMatrix[, x], n.remainder.indices)
                         })
    
    nPerms <- factorial(nBlocks)
  } else{
    permIndices <- cbind(1:n, 
                         unique(replicate(nPerms, c(blockIndexMatrix[, sample(1:nBlocks)])),
                                MARGIN = 2))
  }

  if(!is.null(GX.indices)){
    if(ncol(GX.indices) != nBlocks * (nBlocks - 2) + 2){
      cli::cli_warn("Number of columns in `GX.indices` is not nBlocks * (nBlocks - 2) + 2. Recalculating GX.indices now.")
      GX.indices <- build_GX.indices(blockIndexMatrix, n.remainder.indices)
    } else if(nrow(GX.indices) != data.n){
      cli::cli_warn("Number of rows in `GX.indices` does not match number of rows in data. Recalculating GX.indices now.")
      GX.indices <- build_GX.indices(blockIndexMatrix, n.remainder.indices)
    }
  } else{
    GX.indices <- build_GX.indices(blockIndexMatrix, n.remainder.indices)
  }

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
    gaArgs$fitness <- function(permutation){ fitness_function(permutation = permutation,
                                                              X1.temp = X1.temp,
                                                              X2.temp = X2.temp,
                                                              Z.temp = Z.temp,
                                                              blockIndexMatrix = blockIndexMatrix,
                                                              GX.indices = GX.indices,
                                                              permIndices = permIndices) }
    gaArgs$lower <- rep(1, data.n)
    gaArgs$upper <- rep(data.n, data.n)
    gaArgs$crossover = "gaperm_oxCrossover_R"
  }

  summaryTableList <- vector("list")
  detailedList <- vector("list")
  gaResultsList <- vector("list")
  Q.X1.Z.List <- vector("list")
  
  i <- which(X.assign %in% variables)
  
  # Change this to be based on the formula, if they include another |
  exacttIV <- any(!colnames(X.use)[i] %in% exogenous.var)

  beta.hats.i <- beta.hats[i]
  se.i <- se[i]
  
  if(is.null(beta0)){
    beta0.matrix <- lapply(1:length(beta.hats.i),
                           function(x){
                             precisionToUse <- ifelse(se.i[x] > 0,
                                                      yes = floor(log(se.i[x], base = 10)) - 1,
                                                      no = -1)
  
                             return(round(beta.hats.i[x],
                                          -precisionToUse) +
                                      seq(-30*(10^(precisionToUse+1)),
                                           30*(10^(precisionToUse+1)),
                                          10^precisionToUse)
                                    )
                           })

    names(beta0.matrix) <- names(beta.hats.i)

    beta0.matrix <- do.call('expand.grid', beta0.matrix) |>
      rbind(0)
  } else
  if(beta0 == 0){
    beta0.matrix <- matrix(0, 
                           nrow = 1, 
                           ncol = length(beta.hats.i),
                           dimnames = list(NULL, names(beta.hats.i))) |>
      data.frame()
  } else{
    cli::cli_abort("The `beta0` parameter currently only supports '0' or 'NULL'.")
  }
  
  Y.temp <- as.matrix(Y.use)
  X1.temp <- X.use[,i, drop = FALSE]
  X2.temp <- X.use[,-i, drop = FALSE]

  if(exacttIV){
    Z.temp <- Z.use
  }

  if(optimize){
    
    X1.temp <- scale(X1.temp, center = FALSE, scale = TRUE)

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
                                                "n",
                                                "fitness_function",
                                                "build_GX",
                                                "build_GX.indices"),
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

    cli::cli_alert_success("Optimizing ordering for Wald test.")
    gaResults <- do.call(GA::ga, gaArgs)

    # Close cluster if parallel is true
    if(!is.null(gaArgs$parallel) && ogParArg != FALSE){
      parallel::stopCluster(cl)
      gaArgs$parallel <- ogParArg
    }

    Y.temp <- Y.temp[gaResults@solution[1,],, drop = FALSE]
    X1.temp <- X.use[,i, drop = FALSE][gaResults@solution[1,],, drop = FALSE]
    X2.temp <- X2.temp[gaResults@solution[1,],, drop = FALSE]

    if(exacttIV){
      Z.temp <- Z.temp[gaResults@solution[1,],, drop = FALSE]
    }
  }

  if(exacttIV){
    if(ncol(X2.temp) == 0){
      Q.Z.temp <- Z.temp
    } else if(ncol(X2.temp) == 1 && "(Intercept)" %in% colnames(X2.temp)){
      Q.Z.temp <- apply(Z.temp,
                        MARGIN = 2,
                        function(x) scale(x, center = TRUE))
    } else{ 
      Q.Z.temp <- stats::lm(Z.temp ~ 0 + build_GX(X2.temp, GX.indices),
                            model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals |>
        matrix(nrow = nrow(Z.temp))
    }
    
    Q.X1.Z.List <- Q.Z.temp
  } else{
    if(ncol(X2.temp) == 0){
      Q.X1.temp <- X1.temp
    } else if(ncol(X2.temp) == 1 && "(Intercept)" %in% colnames(X2.temp)){
      Q.X1.temp <- apply(X1.temp,
                         MARGIN = 2,
                         function(x) scale(x, center = TRUE))
    } else{
      Q.X1.temp <- stats::lm(X1.temp ~ 0 + build_GX(X2.temp, GX.indices),
                             model = FALSE, x = FALSE, y = FALSE, qr = FALSE)$residuals |>
        matrix(ncol = ncol(X1.temp))
    }
    
    Q.X1.Z.List <- Q.X1.temp
  }

  p.values.and.extra <- exactt.pval.wald(Y.temp, X1.temp, X2.temp, permIndices, GX.indices, Q.X1.temp, studentize, beta.null.matrix = beta0.matrix)
  p.values <- p.values.and.extra$p.values
  
  result <- structure(list(call = call,
                           nBlocks = nBlocks,
                           detailed = p.values,
                           gaResults = gaResultsList,
                           ivregResults = ivregObject,
                           omega.g = p.values.and.extra$omega.g,
                           GX.indices = GX.indices),
                      class = "exactt.wald")
  
  if(confidence.intervals){
    
    if(alpha >= 1/nPerms & alpha < 1){
      # Create exact projected confidence intervals using Gurobi
      conf.ints <- exact.wald.projected.confidence.intervals(omega.g = p.values.and.extra$omega.g,
                                                             gurobi.params = gurobi.params,
                                                             X1.names = names(beta.hats.i),
                                                             alpha = alpha)
      
      summaryTable <- matrix(data = cbind(beta.hats.i, conf.ints$confidence.intervals),
                             nrow = length(beta.hats.i), 
                             ncol = 3, 
                             dimnames = list(names(beta.hats.i), 
                                             c("Estimate",
                                               "Lower Bound",
                                               "Upper Bound")))
      
      result$gurobi.results <- conf.ints$gurobi.results
    } else if(alpha < 1/nPerms){
      cli::cli_warn("'alpha' is less than 1/nPerms. Confidence interval bounds will be [-Inf, Inf].")
      
      summaryTable <- matrix(data = cbind(beta.hats.i, -Inf, Inf),
                             nrow = length(beta.hats.i), 
                             ncol = 3, 
                             dimnames = list(names(beta.hats.i), 
                                             c("Estimate",
                                               "Lower Bound",
                                               "Upper Bound")))
    } else{
      cli::cli_warn("'alpha' is equal to 1. Confidence interval is empty set because we reject everywhere.")
      
      summaryTable <- matrix(data = beta.hats.i,
                             nrow = length(beta.hats.i), 
                             ncol = 1, 
                             dimnames = list(names(beta.hats.i), 
                                             "Estimate"))
    }
  } else{
    summaryTable <- matrix(data = beta.hats.i,
                           nrow = length(beta.hats.i), 
                           ncol = 1, 
                           dimnames = list(names(beta.hats.i), 
                                           "Estimate"))
  }
  
  result$summary.table <- summaryTable
  
  if(length(gaResultsList) > 0){
    result$gaResults <- gaResultsList
  } 
  
  if(exacttIV){
    result$Q.Z <- Q.X1.Z.List
  } else{
    result$Q.X1 <- Q.X1.Z.List
  }
  
  return(result)
}

