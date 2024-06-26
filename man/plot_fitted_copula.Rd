% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bivariate-modelling.R
\name{plot_fitted_copula}
\alias{plot_fitted_copula}
\title{Plot a fitted bivariate copula}
\usage{
plot_fitted_copula(joint_model, levels = c(0.5, 1, 1.5), ...)
}
\arguments{
\item{joint_model}{A list containing two nonstationary marginal distriutions and a joint copula, as returned by 'fit_copula'.}

\item{levels}{Vector of levels at which contours of the empirical and fitted contours should be plotted.}

\item{...}{Other graphical parameters to be passed to the plotting function.}
}
\description{
Plot the fitted copula against the observed sample to check the goodness of fit. Formal goodness of fit testing is not yet implemented.
}
