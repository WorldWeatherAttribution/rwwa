# PLOTTING FUNCTIONS TO SUPPORT BIVARIATE MODELLING USING COPULAS
#
#' Support function to calculate a 2d matrix of joint exceedance probabilities from a bivariate model
#'
#' @param joint_model A list containing two nonstationary marginal distributions and a joint copula, as returned by 'fit_copula'.
#' @param fixed_cov Data.frame with rows specifying the factual/counterfactual climates for which the joint return period is to be estimated.
#' @param xrange (Optional) vector giving max and min values of x for which exceedance probabilities should be computed. If not provided, will be taken from the range of values in the first marginal model.
#' @param yrange (Optional) vector giving max and min values of y for which exceedance probabilities should be computed. If not provided, will be taken from the range of values in the second marginal model.
#' @param n Integer: number of points along each axis at which exceedance probabilities should be evaluated. Default is 32.
#'
#' @keywords internal
#'
#' @export
#'
copula_mesh <- function(joint_model, fixed_cov, xrange, yrange, n = 32) {

  # compute joint exceedances over regular grid for easy plotting
  mdl_x <- joint_model$mdl_x
  mdl_y <- joint_model$mdl_y

  if(missing(xrange)) xrange <- range(pretty(mdl_x$x))
  if(missing(yrange)) yrange <- range(pretty(mdl_y$x))

  # define the regular mesh for plotting
  x_mesh <- seq(xrange[1], xrange[2],length.out = n)
  y_mesh <- seq(yrange[1], yrange[2], length.out = n)

  # convert the regular mesh to U space
  x_umesh <- sapply(x_mesh, map_to_u, mdl = mdl_x, fixed_cov = fixed_cov)
  y_umesh <- sapply(y_mesh, map_to_u, mdl = mdl_y, fixed_cov = fixed_cov)

  return(list("x" = x_mesh, "y" = y_mesh, "z" = sapply(x_umesh, function(x) sapply(y_umesh, function(y) pCopula(cbind(x,y), joint_model$copula)))))
}


################################################################################################################################
#' Plot contours of equal return periods for a bivariate model.
#'
#' @param joint_model A list containing two nonstationary marginal distributions and a joint copula, as returned by 'fit_copula'.
#' @param fixed_cov Data.frame with rows specifying the factual/counterfactual climates for which the joint return period is to be estimated.
#' @param add Boolean: add to an existing plot? Default is F (create new plot).
#' @param rp Vector of scalars: return periods for which contours should be plotted. Default is c(5,10,20,50).
#' @param xlim (Optional) vector giving max and min values of x to be plotted. If not provided, will be taken from the range of values in the first marginal model. Ignored if 'add' is T.
#' @param ylim (Optional) vector giving max and min values of y to be plotted. If not provided, will be taken from the range of values in the second marginal model. Ignored if 'add' is T.
#' @param ... Other graphical parameters to be passed to the plotting function.
#'
#' @export
#'
plot_joint_contour <- function(joint_model, fixed_cov, add = F, rp = c(5,10,20,50), xlim, ylim, ...) {

  if(nrow(fixed_cov) > 1) {
    cat("More than one set of covariates provided - only showing first set, plot others separately")
    fixed_cov <- fixed_cov[1,,drop = F]
  }

  if(missing(xlim)) { if(add) { xlim <- par("usr")[1:2] } else { xlim <- range(pretty(joint_model$mdl_x$x)) }}
  if(missing(ylim)) { if(add) { ylim <- par("usr")[3:4] } else { ylim <- range(pretty(joint_model$mdl_y$x)) }}

  cmesh <- copula_mesh(joint_model, fixed_cov = fixed_cov, xrange = xlim, yrange = ylim)
  contour(cmesh, levels = 1/rp, labels = rp, xaxs = "i", yaxs = "i", add = add, ...)
}


################################################################################################################################
#' Transform the event of interest into a specified factual/counterfactual climate and plot on existing axes.
#'
#' @param joint_model A list containing two nonstationary marginal distributions and a joint copula, as returned by 'fit_copula'.
#' @param fixed_cov Data.frame with rows specifying the factual/counterfactual climates for which the joint return period is to be estimated.
#' @param ... Other graphical parameters to be passed to the plotting function.
#'
#' @export
#'
plot_joint_event <- function(joint_model, fixed_cov, ...) {

  erl_x <- stransform(joint_model$mdl_x, fixed_cov = fixed_cov)[joint_model$mdl_x$ev_idx]
  erl_y <- stransform(joint_model$mdl_y, fixed_cov = fixed_cov)[joint_model$mdl_y$ev_idx]

  points(erl_x, erl_y, ...)
}
