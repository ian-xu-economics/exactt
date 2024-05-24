#' Check Model and Data Compatibility
#'
#' This function verifies that the `model` is a valid formula with a left-hand side
#' and that the `data` is a data frame or matrix containing all the variables specified
#' in the `model`.
#'
#' @param model A formula with a left-hand side specifying the model.
#' @param data A data frame or matrix with column names containing the variables used in the `model`.
#' 
#' @return Returns `TRUE` if all variables specified in the `model` are present in the `data`.
#' 
#' @importFrom rlang is_formula
#' 
#' @noRd
check_model_data <- function(model, data) {
  # Check if `model` provided and is formula with LHS
  if(is.null(model)){
    stop("The 'model' parameter must be provided.")
  } else if(!rlang::is_formula(model, lhs = TRUE)){
    stop("The 'model' parameter must be a formula with a LHS.")
  }
  
  # Check if `data` provided and is data.frame or matrix with column names.
  if(is.null(data)){
    stop("The 'data' parameter must be provided.'")
  } else if(!is.matrix(data) && !is.data.frame(data)){
    stop("The 'data' parameter must be a data.frame or matrix with column names.")
  } else if(is.matrix(data)){
    data <- data.frame(data)
  }
  
  # Extract variable names from `model`
  terms <- all.vars(model)
  
  # Check if all the variables in the `model` are columns in the data
  missing_vars <- setdiff(terms, names(data))
  
  if (length(missing_vars) > 0) {
    stop(paste("The following variables in the 'model' parameter are not present in the data:", paste(missing_vars, collapse = ", ")))
  }
  
  if (is.matrix(data)) {
    return(data) # Return modified data if it was a matrix
  } else {
    return(NULL) # Return NULL if no modification is needed
  }
}
