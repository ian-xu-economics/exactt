#' Check if Desired Variables are in Regressor Names
#'
#' This internal function checks if the specified variables are present in the set of regressor names (`X.var`).
#' It stops execution if none of the variables are present and issues a warning if some of the variables are not found,
#' skipping those variables in the returned result.
#'
#' @param variables A character vector of variable names to check.
#' @param X.var A character vector of regressor names.
#' 
#' @return A character vector of variable names that exist in `X.var`. If none of the variables are present,
#' the function stops execution with an error message.
#' 
#' @noRd
check_variables_X <- function(variables, X.var){
  
  variablesInData <- variables %in% X.var
  
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
  return(variables)
}
