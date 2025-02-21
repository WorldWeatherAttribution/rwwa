# ESTIMATE PARAMETERS OF NONSTATIONARY MODEL
#
# Methods to estimate the parameters of a nonstationary GEV or normal distribution using max-likelihood estimation
#
#' Support function to calculate the log-likelihood of a vector of data 'x' for a nonstationary distribution defined by the parameters 'pars' (used by 'fit_ns' to estimate the model parameters)
#'
#' @param pars Vector of named parameters of the nonstationary model
#' @param cov Data.frame containing values of the covariate(s) at which the location is to be evaluated
#' @param x Vector of data to evaluate log-likelihood
#' @param dist String defining the parametric form to be fitted (currently only implemented for 'gev' and 'norm')
#' @param fittype String defining the type of model fit to use - method is currently implemented for 'shift' and 'fixeddisp'
#'
#' @return Scalar log-likelihood
#'
#' @keywords internal
#' @export
#'
ns_loglik <- function(pars, cov, x, dist, fittype) {

  nspars <- get.ns_pars(fittype = fittype, pars = pars, fixed_cov = cov)

  loc <- nspars$loc
  scale <- nspars$scale
  shape <- nspars$shape

  # constrain variance to be strictly positive
  if(any(scale <= 0)) return(NA)

  # return negative log-likelihood
  if(dist == "norm") {
    return(-sum(dnorm(x, mean = loc, sd = scale, log = T)))
  } else if(dist == "gev") {
    shape = pars["shape"]
    return(-sum(devd(x, loc = loc, scale = scale, shape = shape, log = T)))
  } else {
    print(paste(dist, "not implemented"))
    return()
  }
}

################################################################################################################################
#' Use maximum likelihood estimation to find the parameters of the nonstationary model that best describes a vector of data 'x' using a set of given covariates
#'
#' @param dist String defining the parametric form to be fitted (currently only implemented for GEV and Gaussian)
#' @param type String defining the relationship between the fitted values and the covariate (currently implemented: 'shift' and 'fixeddisp')
#' @param data Data.frame with named columns containing the variable and any covariates of interest
#' @param varnm String identifying the dependent variable (must be a column name in 'data')
#' @param covnm String or vector of strings identifying the predictors (must be column names in 'data'). If not provided, a constant model will be fitted using a dummy covariate.
#' @param lower Boolean indicating whether to evaluate the lower tail of the data or not: default is F (evaluate the upper tail).
#' @param ev_year (optional) Scalar specifying the year of the event of interest; default is to use the value from the last row of 'data'
#' @param ev (optional) Scalar specifying the magnitude of the event of interest; default is to use the value corresponding to 'ev_year'
#' @param method String defining the method to be used by 'optim' to maximise the log-likelihood: default is 'BFGS'
#'
#' @return List of attributes & parameters defining a nonstationary model
#'
#' @details Further details to be added for dist & type
#'
#' @export
#'
fit_ns <- function(dist, type = "fixeddisp", data, varnm, covnm = NA, lower = F, ev_year = NA, ev = NA, method = "BFGS") {

  # if no covariate name provided, add a dummy covariate & fit stationary series instead
  if(is.na(covnm[1])) {
    data <- cbind(data, "const" = 1)
    covnm <- "const"
  }

  # remove extraneous columns from dataframes
  cov <- data[, covnm, drop = F]
  k <- length(covnm)
  x <- data[,varnm]

  # should also add something to handle case with no covariates

  # currently only works for distributions fully specified by mean & sd: only tested for normal, lognormal
  if(! dist %in% c("norm", "gev")) {
    print("Not yet implemented: use norm or gev")
    return()
  }

  # if looking at lower tail with a GEV, necessary to negate data and consider block maxima - add flag to keep track
  minima <- F
  if (lower & (dist %in% c("gev"))) {
    x <- -x
    minima <- T
  }

  # fit model with appropriate number of parameters, pad if necessary
  init <- c("mu0" = mean(x), "sigma0" = sd(x), setNames(rep(0,k), paste0("alpha_", covnm)))

  if(type %in% c("shiftscale")) init <- c(init, setNames(rep(0,k), paste0("beta_", covnm)))

  if(dist %in% c("gev")) init <- c(init, "shape" = 0)
  fitted <- suppressWarnings(optim(par = init, ns_loglik, cov = cov, x = x, dist = dist, fittype = type, method = method))

  # if looking at lower tail with a GEV, necessary to flip data and consider block maxima, so trend & location parameters have been flipped. This may cause some confusion so may have to modify later!
  if(minima) {
    fitted[["NOTE"]] <- "NB: model parameters are estimated for negated values"
    fitted$par["mu0"] <- -fitted$par["mu0"]
    fitted$par[grepl("alpha", names(fitted$par))] <- -fitted$par[grepl("alpha", names(fitted$par))]
    x <- -x
  }

  # attach assorted useful information
  fitted[["dist"]] <- dist
  fitted[["type"]] <- type
  fitted[["varnm"]] <- varnm
  fitted[["covnm"]] <- covnm
  fitted[["data"]] <- data
  fitted[["x"]] <- x
  fitted[["cov"]] <- cov

  fitted[["lower"]] <- lower               # saves having to specify every time later on
  fitted[["minima"]] <- minima             # look at maxima of 0-temps, rather than minima of observed temps

  # event year: assume that event of interest is most recent, unless told otherwise (used in later plotting functions)
  if(is.na(ev_year)) { ev_year <- data$year[length(x)] }
  fitted[["ev_year"]] <- ev_year

  if (ev_year %in% data$year) {
    ev_idx <- which(data$year == ev_year)
  } else {
    ev_idx <- nrow(data)+1
  }

  # event value: assume that event of interest is most recent, unless told otherwise (used in later plotting functions)
  if(is.na(ev)) {
    if(ev_year %in% data$year) {
      ev <- data[data$year == ev_year,varnm]
    } else {
      print("WARNING: Event year not in data, no event value recorded")
      ev <- NA
    }
  }
  fitted[["ev"]] <- ev
  fitted[["ev_idx"]] <- ev_idx

  return(fitted)
}


################################################################################################################################
#' Support function to refit the parameters of a model using new data (used in bootstrapping)
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param new_data Data.frame with named columns containing the variable and any covariates of interest: must contain the same covariates as the data used to fit 'mdl'
#'
#' @return List of attributes & parameters defining a nonstationary model
#'
#' @keywords internal
#'
#' @export
#'
refit <- function(mdl, new_data) {
  fit_ns(dist = mdl$dist, type = mdl$type, data = new_data, varnm = mdl$varnm, covnm = mdl$covnm, lower = mdl$lower, ev = mdl$ev)
}

################################################################################################################################
#' Calculate AIC from model output
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#'
#' @return AIC as scalar
#'
#' @export
#'
aic <- function(mdl) 2 * length(mdl$par) + 2 * mdl$value


################################################################################################################################
#' Calculate R^2 and adjusted R^2 from model output
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#'
#' @return Vector containing R^2 and adjusted R^2
#'
#' @export
#'
rsquared <- function(mdl) {

  tss = sum((mdl$x - mean(mdl$x))^2)
  rss = sum((ns_pars(mdl)$loc - mdl$x)^2)

  n = length(mdl$x)
  p = length(mdl$covnm)

  rsquared = 1 - (rss / tss)
  rsquared_adj = 1 - (1 - rsquared) * (n-1) / (n-p-1)

  return(c(r2 = rsquared, r2.adj = rsquared_adj))
}
