#' Plot P-values for 'et' Objects
#'
#' This function generates plots of P-values against beta-null (\eqn{\beta_0}) values for each variable
#' specified in an 'et' object. The function uses ggplot2 for plotting, showing the significance level
#' with a horizontal line and the estimate from the summary with a vertical line.
#'
#' @param et An object of class 'et', typically the output from `exactt` function,
#'   containing elements 'detailed' for plotting data, 'summary' for vertical lines at estimates,
#'   and 'call' from which the significance level 'alpha' is extracted.
#' @param variables A character vector specifying which variables to plot.
#'   If NULL, plots are generated for all variables contained in the 'et' object.
#'
#' @return A list of ggplot objects, one for each variable specified. Each plot represents
#'   the relationship between P-values and beta-null values for that variable,
#'   highlighted with reference lines for the estimated coefficient and significance level.
#'
#' @importFrom rlang sym
#' @export
exacttPlot = function(et, variables = NULL){
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function. Please install it using `install.packages('ggplot2')`.")
  }
  if (!requireNamespace("latex2exp", quietly = TRUE)) {
    stop("Package 'latex2exp' is required for this function. Please install it using `install.packages('latex2exp')`.")
  }
  
  alpha <- ifelse(is.null(et$call$alpha), yes = 0.05, no = et$call$alpha)
  
  if(is.null(variables)){
    variables <- 1:length(et$detailed)
  }
  
  plots <- vector("list", length = length(variables))
  plotsCounter <- 1
  
  for(i in 1:length(et$detailed)){
    if(!i %in% variables){
      next
    }
    
    point_estimate <- et$summary[i, 1]
    variable_name <- names(et$detailed)[i]
    
    data <- et$detailed[[i]]
    data$point_estimate = point_estimate
    data$alpha = alpha
    
    if(data$beta0[1] == 0 
       && data$beta0[2] - data$beta0[1] != data$beta0[3] - data$beta0[2]){
      data <- data[-1, ]
    }

    # Create plot with correct settings
    plots[[plotsCounter]] <- ggplot2::ggplot(data = data) + 
      ggplot2::theme(legend.position = "bottom",
                     panel.grid.major = element_line(color = "cyan"),
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "transparent", color = NA),
                     axis.line = element_line(color = "black"),
                     plot.title = element_text(hjust = 0.5),
                     text = element_text(size = 12)) +
      ggplot2::geom_line(ggplot2::aes(x = !!rlang::sym("beta0"), y = !!rlang::sym("pval"))) + 
      ggplot2::geom_hline(ggplot2::aes(yintercept = alpha, color = "alpha"), linetype = "dashed") + 
      ggplot2::geom_vline(ggplot2::aes(xintercept = point_estimate, color = "estimate"), linetype = "dotted") + 
      ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.1), labels = seq(0, 1, 0.1)) + 
      ggplot2::scale_x_continuous(breaks = ggplot2::waiver(), n.breaks = 10) + 
      ggplot2::labs(x = latex2exp::TeX(paste0("$\\beta^0_{", variable_name, "}$")), y = "P-value") + 
      ggplot2::ggtitle(label = latex2exp::TeX(paste0("P-value vs. $\\beta^0_{", variable_name, "}$"))) + 
      ggplot2::scale_color_manual(name = "", 
                                  values = c("blue", "red"), 
                                  labels = c(latex2exp::TeX("$\\alpha$"), 
                                             latex2exp::TeX("$\\hat{\\beta}$")))
    
    plotsCounter <- plotsCounter + 1
  }
  
  return(plots)
}
