% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exacttPlot.R
\name{exacttPlot}
\alias{exacttPlot}
\title{Plot P-values for 'et' Objects}
\usage{
exacttPlot(et, variables = NULL)
}
\arguments{
\item{et}{An object of class 'et', typically the output from \code{exactt} function,
containing elements 'detailed' for plotting data, 'summary' for vertical lines at estimates,
and 'call' from which the significance level 'alpha' is extracted.}

\item{variables}{A character vector specifying which variables to plot.
If NULL, plots are generated for all variables contained in the 'et' object.}
}
\value{
A list of ggplot objects, one for each variable specified. Each plot represents
the relationship between P-values and beta-null values for that variable,
highlighted with reference lines for the estimated coefficient and significance level.
}
\description{
This function generates plots of P-values against beta-null (\eqn{\beta_0}) values for each variable
specified in an 'et' object. The function uses ggplot2 for plotting, showing the significance level
with a horizontal line and the estimate from the summary with a vertical line.
}
