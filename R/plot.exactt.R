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
  
  for(var.num in 1:length(x$detailed)){
    if(!var.num %in% variables){
      next
    }
    
    variable_name <- names(x$detailed)[var.num]
    
    data <- x$detailed[[var.num]]
    
    ci.lower <- data$beta0.start[min(which(data$pvals >= alpha))]
    ci.upper <- data$beta0.end[max(which(data$pvals >= alpha))]
    
    # Extend the bottom margin to make space for the legend
    #old_par <- par(no.readonly = TRUE)  # Save old par settings
    graphics::par(mar = c(7, 4.5, 1, 1))  # Increase bottom margin
    
    # Determine finite data for calculating plot limits
    finite_beta0 <- c(data$beta0.start, data$beta0.end)
    finite_beta0 <- finite_beta0[is.finite(finite_beta0)]
    
    # Calculate the finite range
    finite_range <- range(finite_beta0)
    
    # Extend the finite range by 5% on both sides for plotting limits
    x_extension <- 0.05 * diff(finite_range)
    x_limits <- c(finite_range[1] - x_extension, finite_range[2] + x_extension)
    
    # Replace -Inf and Inf with finite values beyond the plotting limits
    data$beta0.start[!is.finite(data$beta0.start)] <- x_limits[1] - x_extension
    data$beta0.end[!is.finite(data$beta0.end)] <- x_limits[2] + x_extension
    
    # Initialize the plot
    graphics::plot(1, 
                   type = "n",
                   xlim = x_limits,
                   ylim = c(0, 1),
                   xlab = bquote(beta[.(variable_name)]^0),
                   ylab = "P-value",
                   xaxt = "n", yaxt = "n",
                   cex.axis = 1, cex.lab = 1.4, family = "Helvetica")
    
    # Add gridlines
    graphics::abline(h = seq(0, 1, 0.1), col = "gray90", lty = "dotted", lwd = 0.75)
    graphics::abline(v = pretty(x_limits, n = 8), col = "gray90", lty = "dotted", lwd = 0.75)
    
    # Shade regions where p-values are above alpha
    above_alpha <- which(data$pvals > alpha)
    if(length(above_alpha) > 0) {
      for(i in above_alpha) {
        graphics::rect(xleft = max(data$beta0.start[i], x_limits[1]),
                       xright = min(data$beta0.end[i], x_limits[2]),
                       ybottom = alpha,
                       ytop = data$pvals[i],
                       col = rgb(0, 1, 0, alpha = 0.3),
                       border = NA)
      }
    }
    
    num.rows <- nrow(data)
    
    # Plot horizontal lines for each interval
    for(i in 1:num.rows) {
      graphics::segments(x0 = max(data$beta0.start[i], x_limits[1]),
                         x1 = min(data$beta0.end[i], x_limits[2]),
                         y0 = data$pvals[i],
                         y1 = data$pvals[i],
                         col = "black",
                         lwd = 1.25)
    }
    
    # Arrowhead at the start
    graphics::arrows(x0 = x_limits[1],
                     y0 = data$pvals[1],
                     x1 = data$beta0.end[1],
                     y1 = data$pvals[1],
                     col = "black",
                     lwd = 1.1,
                     length = 0.05,
                     angle = 30,
                     code = 1)
    
    # Arrowhead at the end
    graphics::arrows(x0 = data$beta0.start[num.rows],
                     y0 = data$pvals[num.rows],
                     x1 = x_limits[2],
                     y1 = data$pvals[num.rows],
                     col = "black",
                     lwd = 1.1,
                     length = 0.05,
                     angle = 30,
                     code = 2)
    
    # Add horizontal and vertical reference lines
    graphics::abline(h = alpha, col = "orange", lty = "dashed", lwd = 0.75)
    graphics::abline(v = x$summary[var.num, 1], col = "red", lty = "dashed", lwd = 0.75)
    graphics::abline(v = c(ci.lower, ci.upper), col = "cornflowerblue", lty = "dashed", lwd = 0.75)
    
    # Calculate offset for text positioning (2% of the x-axis range)
    x_offset <- 0.02 * diff(x_limits)
    
    # Add text label for ci.lower (Lower Bound)
    graphics::text(x = ci.lower - x_offset,
                   y = 0.5,  # Vertical position (middle of the plot)
                   labels = paste0("Lower bound of ", 
                                   (1-alpha)*100, 
                                   "% CI (", 
                                   round(ci.lower, 3), 
                                   ")"),
                   srt = 90,  # Rotate text 90 degrees
                   adj = c(0.5, 0.25),
                   col = "cornflowerblue",
                   cex = 0.9,
                   family = "Helvetica")
    
    # Add text label for ci.upper (Upper Bound)
    graphics::text(x = ci.upper + x_offset,
                   y = 0.5,  # Vertical position (middle of the plot)
                   labels = paste0("Upper bound of ", 
                                   (1-alpha)*100, 
                                   "% CI (", 
                                   round(ci.upper, 3), 
                                   ")"),
                   srt = -90,  # Rotate text 90 degrees
                   adj = c(0.5, 0.25),
                   col = "cornflowerblue",
                   cex = 0.9,
                   family = "Helvetica")
    
    # Customize the axes
    graphics::axis(1, at = pretty(x_limits, n = 8), labels = TRUE, las = 1, cex.axis = 1, family = "Helvetica")
    
    # Add ticks at 0.1 intervals, but labels at 0.2 intervals
    graphics::axis(2, at = seq(0, 1, 0.1), labels = NA, tck = -0.02)
    graphics::axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), cex.axis = 1, family = "Helvetica")
    
    # Position the legend inside the plot area with a transparent background
    graphics::legend("bottom", 
                     inset=c(0, -0.28),
                     legend = c(expression(alpha), expression(hat(beta)), "CI bounds"),
                     col = c("orange", "red", "cornflowerblue"), lty = c("dashed", "dashed", "dashed"),
                     cex = 0.9, 
                     bty = "n", 
                     bg = "transparent",
                     horiz = TRUE,
                     xpd = TRUE)
    #par(old_par)
    
  }
}

# Previous code
# point_estimate <- x$summary[i, 1]
# variable_name <- names(x$detailed)[i]
# 
# pvalue.segments <- x$detailed[[i]]
# data$point_estimate = point_estimate
# data$alpha = alpha
# 
# if(data$beta0[1] == 0 
#    && data$beta0[2] - data$beta0[1] != data$beta0[3] - data$beta0[2]){
#   data <- data[-1, ]
# }
# 
# beta0 <- data$beta0
# beta0.pval <- data$beta0.pval
# 
# # Find the indices where p-values are above alpha
# above_alpha <- which(beta0.pval > alpha)
# 
# # Plot the basic setup without any data first
# plot(beta0, beta0.pval, type = "n",
#      xlab = bquote(beta[.(variable_name)]^0),
#      ylab = "P-value",
#      ylim = c(0, 1),
#      xlim = range(beta0),
#      xaxt = "n", yaxt = "n", # Remove default axes for custom control
#      cex.axis = 1, cex.lab = 1.4, family = "Helvetica"
# )
# 
# # Lighter gridlines
# graphics::abline(h = seq(0, 1, 0.1), col = "gray90", lty = "dashed")
# graphics::abline(v = pretty(beta0, n = 8), col = "gray90", lty = "dashed")
# 
# # Add the green shaded region where p-values are above alpha
# if(length(above_alpha) > 0) {
#   graphics::polygon(c(beta0[above_alpha], 
#                       rev(beta0[above_alpha])), 
#                     c(rep(alpha, length(above_alpha)), 
#                       rev(beta0.pval[above_alpha])),
#                     col = grDevices::rgb(0, 1, 0, alpha = 0.3), border = NA)
# }
# 
# # Add the points and line with improved styling
# graphics::lines(beta0, 
#                 beta0.pval, 
#                 col = "black", 
#                 lwd = 1.5)
# graphics::points(beta0, 
#                  beta0.pval, 
#                  col = grDevices::rgb(0, 0, 0, alpha = 0.75),
#                  pch = 1, 
#                  cex = 0.8)
# 
# # Add horizontal and vertical reference lines
# graphics::abline(h = alpha, 
#                  col = "blue", 
#                  lty = "dashed", 
#                  lwd = 1.5)
# graphics::abline(v = point_estimate, 
#                  col = "red", 
#                  lty = "dotted", 
#                  lwd = 1.5)
# 
# # Customize the axes
# graphics::axis(1, 
#                at = pretty(beta0, n = 8), 
#                labels = TRUE, 
#                las = 1, 
#                cex.axis = 1, 
#                family = "Helvetica")
# 
# # Add ticks at 0.1 intervals, but labels at 0.2 intervals
# graphics::axis(2, 
#                at = seq(0, 1, 0.1), 
#                labels = NA, 
#                tck = -0.02)
# graphics::axis(2, 
#                at = seq(0, 1, 0.2), 
#                labels = seq(0, 1, 0.2), 
#                cex.axis = 1, 
#                family = "Helvetica")
# 
# # Position the legend inside the plot area with a transparent background
# graphics::legend("topright", inset = c(0.04, 0.04),
#                  legend = c(expression(alpha), expression(hat(beta))),
#                  col = c("blue", "red"), 
#                  lty = c("dashed", "dotted"),
#                  bty = "o", cex = 1.2, box.col = "black")
