# BIVARIATE MODELLING USING COPULAS
#
#' Fit a copula to two nonstationary marginal distributions
#'
#' @param mdl_x List of attributes & parameters defining first nonstationary model, as returned by 'fit_ns'
#' @param mdl_y List of attributes & parameters defining second nonstationary model, as returned by 'fit_ns'
#'
#' @description
#' Fit a 2d copula to model the dependence between two nonstationary marginal distributions. Currently only the t-copula is implemented.
#'
#' @export
#'
fit_copula <- function(mdl_x, mdl_y) {

  # transform marginals to U
  u_x <- map_to_u(mdl_x)
  u_y <- map_to_u(mdl_y)

  # Fit the copula and extract a copula object - need to generalise this
  # if(is(copulatype, "tCopula")) {
    fitted_copula <- fitCopula(tCopula(), data = cbind(u_x, u_y), hideWarnings = T)
    cfit <- tCopula(coef(fitted_copula)["rho.1"], df = round(coef(fitted_copula)["df"],0), df.fixed = T) # integer DF needed for eg. pCopula & goodness of fit
  # } else {
  #   print("Copula type not yet implemented")
  #   return(NULL)
  # }

  # output
  return(list(mdl_x = mdl_x, mdl_y = mdl_y, copula = cfit))
}


################################################################################################################################
#' Plot a fitted bivariate copula
#'
#' @description
#' Plot the fitted copula against the observed sample to check the goodness of fit. Formal goodness of fit testing is not yet implemented.
#'
#' @param joint_model A list containing two nonstationary marginal distriutions and a joint copula, as returned by 'fit_copula'.
#' @param ... Other graphical parameters to be passed to the plotting function.
#'
#' @export
#'
plot_fitted_copula <- function(joint_model, ...) {

  # transform marginals to U
  u_x <- map_to_u(joint_model$mdl_x)
  u_y <- map_to_u(joint_model$mdl_y)

  # generate sample from the copula
  samp <- rCopula(length(u_x), joint_model$copula)

  plot(u_x, u_y, col = "black", pch = 20, ...)
  points(samp, col = "cornflowerblue", pch = 1)

  contour(kde2d(u_x, u_y), col = "black", add = T, levels = c(0.5,1))
  contour(joint_model$copula, dCopula, add = T, col = "cornflowerblue", lty = 2, levels = c(0.5,1))
}


################################################################################################################################
#' Get joint return period of two marginal events in a specified climate
#'
#' @param joint_model A list containing two nonstationary marginal distributions and a joint copula, as returned by 'fit_copula'.
#' @param fixed_cov Data.frame with rows specifying the factual/counterfactual climates for which the joint return period is to be estimated.
#' @param ev_x (Optional) scalar: the event value to be evaluated for the first marginal distribution. Default is to use the event value specified during model fitting.
#' @param ev_y (Optional) scalar: the event value to be evaluated for the second marginal distribution. Default is to use the event value specified during model fitting.
#'
#' @returns Scalar value per row of 'fixed_cov', givving the joint return period according to the specified bivariate model.
#' @export
#'
joint_returnperiod <- function(joint_model, fixed_cov, ev_x, ev_y) {

  if(missing(ev_x)) ev_x <- joint_model$mdl_x$ev
  if(missing(ev_y)) ev_y <- joint_model$mdl_y$ev

  setNames(sapply(1:nrow(fixed_cov), function(i)  {
    1/pCopula(c(map_to_u(joint_model$mdl_x, ev_x, fixed_cov[i,,drop = F]), map_to_u(joint_model$mdl_y, ev_y, fixed_cov[i,,drop = F])), joint_model$copula)
  }), rownames(fixed_cov))
}
