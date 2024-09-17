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
#' @param pointEstimate A Boolean indicating whether to include the point estimate in the plot.
#' @param ciBounds a Boolean indicating whether to include confidence interval bounds in the plot.
#' @param ... Additional arguments passed to the legend function.
#' - `legend_position`: A character indicating where to place the legend ("topright" by default).
#' - `legend_inset`: A vector of length 2 indicating the inset distance(s) from the margins (c(0.025, 0.025) by default).
#' - `legend_cex`: A double indicating the character expansion factor (0.9 by default).
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
plot.exactt = function(x, variables = NULL, pointEstimate = TRUE, ciBounds = TRUE, xlimits = NULL, ylimits = NULL, ...){
  
  dots = list(...)
  
  if("alpha" %in% names(dots)){
    alpha <- dots$alpha
  } else{
    alpha <- ifelse(is.null(x$call$alpha), yes = 0.05, no = x$call$alpha)
  }
  
  if(is.null(variables)){
    variables <- 1:length(x$detailed)
  }
  
  for(var.num in 1:length(x$detailed)){
    if(!attr(x$detailed[[var.num]], "assign") %in% variables){
      next
    }
    
    variable_name <- names(x$detailed)[var.num]
    
    data <- x$detailed[[var.num]]
    
    # Extend the bottom margin to make space for the legend
    graphics::par(mar = c(4.5, 4.5, 1, 1))  # Increase bottom margin
    
    # Determine finite data for calculating plot limits
    finite_beta0 <- c(data$beta0.start, data$beta0.end)
    finite_beta0 <- finite_beta0[is.finite(finite_beta0)]
    
    # Calculate the finite range
    finite_range <- range(finite_beta0)
    
    # Extend the finite range by 5% on both sides for plotting limits
    x_extension <- 0.05 * diff(finite_range)
    x_limits_auto <- c(finite_range[1] - x_extension, finite_range[2] + x_extension)
    
    # Adjust x_limits if user provides xlimits
    if(!is.null(xlimits)){
      x_limits <- xlimits
    } else {
      x_limits <- x_limits_auto
    }
    
    # Ensure x_extension is defined
    if(diff(x_limits) == 0){
      x_extension <- 0.05  # A small fixed extension
    } else {
      x_extension <- 0.05 * diff(x_limits)
    }
    
    # Replace -Inf and Inf with x_limits[1] and x_limits[2]
    data$beta0.start[!is.finite(data$beta0.start)] <- x_limits[1]
    data$beta0.end[!is.finite(data$beta0.end)] <- x_limits[2]
    
    # Set y_limits
    if(is.null(ylimits)){
      y_limits <- c(0, 1)
    } else {
      y_limits <- ylimits
    }
    
    # Initialize the plot
    graphics::plot(1, 
                   type = "n",
                   xlim = x_limits,
                   ylim = y_limits,
                   xlab = bquote(beta[.(variable_name)]^0),
                   ylab = "P-value",
                   xaxt = "n", yaxt = "n",
                   cex.axis = 1, cex.lab = 1.4, family = "Helvetica")
    
    # Add gridlines
    graphics::abline(h = seq(y_limits[1], y_limits[2], 0.1), col = "gray90", lty = "dotted", lwd = 0.75)
    graphics::abline(v = pretty(x_limits, n = 8), col = "gray90", lty = "dotted", lwd = 0.75)
    
    # Shade regions where p-values are above alpha
    above_alpha <- which(data$pvals > alpha)
    if(length(above_alpha) > 0) {
      for(i in above_alpha) {
        shade_start <- max(data$beta0.start[i], x_limits[1])
        shade_end <- min(data$beta0.end[i], x_limits[2])
        if(shade_end > shade_start) {
          graphics::rect(xleft = shade_start,
                         xright = shade_end,
                         ybottom = alpha,
                         ytop = data$pvals[i],
                         col = rgb(0, 1, 0, alpha = 0.3),
                         border = NA)
        }
      }
    }
    
    num.rows <- nrow(data)
    
    # Plot horizontal lines for each interval within x_limits
    for(i in 1:num.rows) {
      segment_start <- max(data$beta0.start[i], x_limits[1])
      segment_end <- min(data$beta0.end[i], x_limits[2])
      if(segment_end > segment_start) {
        graphics::segments(x0 = segment_start,
                           x1 = segment_end,
                           y0 = data$pvals[i],
                           y1 = data$pvals[i],
                           col = "black",
                           lwd = 1.25)
      }
    }
    
    # Adjust arrow plotting at the start
    if(data$beta0.end[1] >= x_limits[1]) {
      arrow_start <- max(data$beta0.start[1], x_limits[1])
      arrow_end <- min(data$beta0.end[1], x_limits[2])
      if(arrow_end > arrow_start) {
        graphics::arrows(x0 = arrow_start,
                         y0 = data$pvals[1],
                         x1 = arrow_end,
                         y1 = data$pvals[1],
                         col = "black",
                         lwd = 1.1,
                         length = 0.05,
                         angle = 30,
                         code = ifelse(data$beta0.start[1] <= x_limits[1], 1, 0))
      }
    }
    
    # Adjust arrow plotting at the end
    if(data$beta0.start[num.rows] < x_limits[2]) {
      arrow_start <- max(data$beta0.start[num.rows], x_limits[1])
      arrow_end <- min(data$beta0.end[num.rows], x_limits[2])
      if(arrow_end > arrow_start) {
        graphics::arrows(x0 = arrow_start,
                         y0 = data$pvals[num.rows],
                         x1 = arrow_end,
                         y1 = data$pvals[num.rows],
                         col = "black",
                         lwd = 1.1,
                         length = 0.05,
                         angle = 30,
                         code = ifelse(data$beta0.end[num.rows] >= x_limits[2], 2, 0))
      }
    }
    
    legendText <- NULL
    legendCol <- NULL
    legendLty <- NULL
    # Add horizontal and vertical reference lines
    graphics::abline(h = alpha, col = "orange", lty = "dashed", lwd = 0.75)
    legendText <- c(legendText, expression(alpha))
    legendCol <- c(legendCol, "orange")
    legendLty <- c(legendLty, "dashed")
    
    if(pointEstimate){
      # Only plot if point estimate is within x_limits
      if(x$summary[var.num, 1] >= x_limits[1] && x$summary[var.num, 1] <= x_limits[2]){
        graphics::abline(v = x$summary[var.num, 1], col = "red", lty = "dashed", lwd = 0.75)
        
        legendText <- c(legendText, expression(hat(beta)))
        legendCol <- c(legendCol, "red")
        legendLty <- c(legendLty, "dashed")
      }
    }
    
    if(ciBounds){
      ci.lower.index <- min(which(data$pvals >= alpha))
      ci.upper.index <- max(which(data$pvals >= alpha))
      ci.lower <- data$beta0.start[ci.lower.index]
      ci.upper <- data$beta0.end[ci.upper.index]
      
      # Only plot if CI bounds are within x_limits
      ci_bounds_within_limits <- (ci.lower >= x_limits[1] && ci.lower <= x_limits[2]) ||
        (ci.upper >= x_limits[1] && ci.upper <= x_limits[2])
      if(ci_bounds_within_limits){
        graphics::abline(v = c(ci.lower, ci.upper), col = "cornflowerblue", lty = "dashed", lwd = 0.75)
        
        # Calculate offset for text positioning (2% of the x-axis range)
        x_offset <- 0.02 * diff(x_limits)
        
        # Add text label for ci.lower (Lower Bound)
        if(ci.lower >= x_limits[1] && ci.lower <= x_limits[2]){
          graphics::text(x = ci.lower - x_offset,
                         y = 0.5 * (y_limits[1] + y_limits[2]),
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
        }
        
        # Add text label for ci.upper (Upper Bound)
        if(ci.upper >= x_limits[1] && ci.upper <= x_limits[2]){
          graphics::text(x = ci.upper + x_offset,
                         y = 0.5 * (y_limits[1] + y_limits[2]),
                         labels = paste0("Upper bound of ", 
                                         (1-alpha)*100, 
                                         "% CI (", 
                                         round(ci.upper, 3), 
                                         ")"),
                         srt = -90,  # Rotate text -90 degrees
                         adj = c(0.5, 0.25),
                         col = "cornflowerblue",
                         cex = 0.9,
                         family = "Helvetica")
        }
        
        legendText <- c(legendText, "CI Bounds")
        legendCol <- c(legendCol, "cornflowerblue")
        legendLty <- c(legendLty, "dashed")
      }
    }
    
    # Customize the axes
    graphics::axis(1, at = pretty(x_limits, n = 8), labels = TRUE, las = 1, cex.axis = 1, family = "Helvetica")
    
    # Add ticks at 0.1 intervals, but labels at 0.2 intervals
    graphics::axis(2, at = seq(y_limits[1], y_limits[2], 0.1), labels = NA, tck = -0.02)
    graphics::axis(2, at = seq(y_limits[1], y_limits[2], 0.2), labels = seq(y_limits[1], y_limits[2], 0.2), cex.axis = 1, family = "Helvetica")
    
    # Build legend
    legend_position <- ifelse(is.null(dots$legend_position), 
                              yes = "topright",
                              no = dots$legend_position)
    
    if(is.null(dots$legend_inset)){ 
      legend_inset <- c(0.025, 0.025)
    } else{
      legend_inset <- dots$legend_inset
    }
    
    legend_cex <- ifelse(is.null(dots$legend_cex), 
                         yes = 0.9,
                         no = dots$legend_cex)
    
    graphics::legend(legend_position, 
                     inset = legend_inset,
                     legend = legendText,
                     col = legendCol, 
                     lty = legendLty,
                     cex = legend_cex, 
                     bty = "o")
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
