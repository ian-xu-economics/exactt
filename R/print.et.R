#' Print Method for Objects of Class 'et'
#'
#' This function provides a customized print method for objects of class 'et'.
#' It formats and displays the call and summary stored within such objects.
#' The output is intended to provide a clear, readable presentation of the essential
#' information contained in the object.
#'
#' @param x An object of class 'et', typically containing at least a 'call' and a 'summary' component.
#'            The 'call' is expected to be an expression representing how the object was created.
#'            The 'summary' should contain a summary of the object, which can be any form that can be
#'            formatted and printed.
#' @param ... Additional arguments passed to the print function, potentially modifying 
#'            the output or print behavior (e.g., controlling the number of digits printed).
#'
#' @return Invisible NULL. The function is used for its side effect of printing to the console.
#'
#' @export
print.et = function(x, ...){
  
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("\nSummary:\n")
  print.default(format(x$summary),
                print.gap = 2L, quote = FALSE)
  
}
