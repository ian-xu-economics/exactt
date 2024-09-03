#' Plot P-values for 'exactt' Objects
#'
#' This function generates plots of P-values against beta-null (\eqn{\beta_0}) values for each variable
#' specified in an 'exactt' object. The function uses ggplot2 for plotting, showing the significance level
#' with a horizontal line and the estimate from the summary with a vertical line.
#'
#' @param x An object of class 'exactt', typically the output from `exactt` function,
#'   containing elements 'detailed' for plotting data, 'summary' for vertical lines at estimates,
#'   and 'call' from which the significance level 'alpha' is extracted.
#' @param variables A character vector specifying which variables to plot.
#'   If NULL, plots are generated for all variables contained in the 'exactt' object.
#' @param ... Additional arguments passed to the plot function (not used currently, but included for consistency).
#'
#' @return A list of ggplot objects, one for each variable specified. Each plot represents
#'   the relationship between P-values and beta-null values for that variable,
#'   highlighted with reference lines for the estimated coefficient and significance level.
#'
#' @importFrom graphics abline polygon lines points axis legend
#' @importFrom grDevices rgb
#'
#' @method plot exactt
#' @export
plot.exactt = function(x, variables = NULL, ...){
  
  dots = list(...)
  
  if("alpha" %in% names(dots)) alpha <- dots$alpha
  else alpha <- ifelse(is.null(x$call$alpha), yes = 0.05, no = x$call$alpha)
  
  if(is.null(variables)){
    variables <- 1:length(x$detailed)
  }
  
  for(i in 1:length(x$detailed)){
    if(!i %in% variables){
      next
    }
    
    point_estimate <- x$summary[i, 1]
    variable_name <- names(x$detailed)[i]
    
    data <- x$detailed[[i]]
    data$point_estimate = point_estimate
    data$alpha = alpha
    
    if(data$beta0[1] == 0 
       && data$beta0[2] - data$beta0[1] != data$beta0[3] - data$beta0[2]){
      data <- data[-1, ]
    }
    
    beta0 <- data$beta0
    beta0.pval <- data$beta0.pval
    
    # Find the indices where p-values are above alpha
    above_alpha <- which(beta0.pval > alpha)
    
    # Plot the basic setup without any data first
    plot(beta0, beta0.pval, type = "n",
         xlab = bquote(beta[.(variable_name)]^0),
         ylab = "P-value",
         ylim = c(0, 1),
         xlim = range(beta0),
         xaxt = "n", yaxt = "n", # Remove default axes for custom control
         cex.axis = 1, cex.lab = 1.4, family = "Helvetica"
    )
    
    # Lighter gridlines
    graphics::abline(h = seq(0, 1, 0.1), col = "gray90", lty = "dashed")
    graphics::abline(v = pretty(beta0, n = 8), col = "gray90", lty = "dashed")
    
    # Add the green shaded region where p-values are above alpha
    if(length(above_alpha) > 0) {
      graphics::polygon(c(beta0[above_alpha], 
                          rev(beta0[above_alpha])), 
                        c(rep(alpha, length(above_alpha)), 
                          rev(beta0.pval[above_alpha])),
              col = grDevices::rgb(0, 1, 0, alpha = 0.3), border = NA)
    }
    
    # Add the points and line with improved styling
    graphics::lines(beta0, 
                    beta0.pval, 
                    col = "black", 
                    lwd = 1.5)
    graphics::points(beta0, 
                     beta0.pval, 
                     col = grDevices::rgb(0, 0, 0, alpha = 0.70),
                     pch = 1, 
                     cex = 0.8)
    
    # Add horizontal and vertical reference lines
    graphics::abline(h = alpha, 
                     col = "blue", 
                     lty = "dashed", 
                     lwd = 1.5)
    graphics::abline(v = point_estimate, 
                     col = "red", 
                     lty = "dotted", 
                     lwd = 1.5)
    
    # Customize the axes
    graphics::axis(1, 
                   at = pretty(beta0, n = 8), 
                   labels = TRUE, 
                   las = 1, 
                   cex.axis = 1, 
                   family = "Helvetica")
    
    # Add ticks at 0.1 intervals, but labels at 0.2 intervals
    graphics::axis(2, 
                   at = seq(0, 1, 0.1), 
                   labels = NA, 
                   tck = -0.02)
    graphics::axis(2, 
                   at = seq(0, 1, 0.2), 
                   labels = seq(0, 1, 0.2), 
                   cex.axis = 1, 
                   family = "Helvetica")
    
    # Position the legend inside the plot area with a transparent background
    graphics::legend("topright", inset = c(0.04, 0.04),
                     legend = c(expression(alpha), expression(hat(beta))),
                     col = c("blue", "red"), 
                     lty = c("dashed", "dotted"),
                     bty = "o", cex = 1.2, box.col = "black")
  }
}
