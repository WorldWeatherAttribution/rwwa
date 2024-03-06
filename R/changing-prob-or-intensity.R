# CHANGES IN LIKELIHOOD & INTENSITY
#
#' Support function to map a series of events x to exceedance probabilities under a given model
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param x Vector of event magnitudes to be transformed
#' @param fixed_cov data.frame containing values of the covariate to be used to transform the event magnitudes. Default value is NA, in which case nonstationary parameters are returned for all values of the covariate(s) used to fit the model
#'
#' @return Vector of exceedance probabilities
#'
#' @keywords internal
#'
#' @export
#'
map_to_u <- function(mdl, x, fixed_cov = NA) {

  pars <- ns_pars(mdl, fixed_cov = fixed_cov)
  if(missing(x)) x <- mdl$x

  # retrieve the actual fitted model parameters if they were negated for fitting (ie if looking at lower tails of GEV)
  if(mdl$minima) {
    pars$loc <- -pars$loc
    x = -x
    mdl$lower <- !mdl$lower # also have to look at the opposite tail
  }

  # get exceedance probability
  if(mdl$dist == "norm") {
    pit <- pnorm(x, mean = pars$loc, sd = pars$scale, lower.tail = mdl$lower)
  } else if(mdl$dist == "gev") {
    pit <- sapply(1:length(pars$loc), function(i) pevd(x[i], loc = pars$loc[i], scale = pars$scale[i], shape = pars$shape[i], lower.tail = mdl$lower))
  } else {
    return(NULL)
  }
  return(pit)
}


################################################################################################################################
#' Support function to map a series of exceedance probabilities u to event magnitudes x under a given model
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param u Vector of exceedance probabilities to be transformed
#' @param fixed_cov data.frame containing values of the covariate to be used to transform the probabilities. Default value is NA, in which case nonstationary parameters are returned for all values of the covariate(s) used to fit the model
#'
#' @return Vector of event magnitudes
#'
#' @keywords internal
#'
#' @export
#'
map_from_u <- function(mdl, u, fixed_cov = NA) {

  if(missing(u)) u <- map_to_u(mdl, fixed_cov = fixed_cov)
  pars <- ns_pars(mdl, fixed_cov = fixed_cov)

  # retrieve the actual fitted model parameters if they were flipped for fitting
  if(mdl$minima) {
    pars$loc <- -pars$loc
    mdl$lower <- !mdl$lower       # also have to look at the opposite tail
  }

  # map quantile onto stationary distribution
  if(mdl$dist == "norm") {
    # erl <- sapply(1:length(u), function(j) {
    #     qnorm(u[j], mean = pars$loc, sd = pars$scale, lower.tail = mdl$lower)
    # })
    erl <- qnorm(u, mean = pars$loc, sd = pars$scale, lower.tail = mdl$lower)
  } else if(mdl$dist == "gev") {
    erl <- sapply(1:length(u), function(j) {
      sapply(1:length(pars$loc), function(i) {
        qevd(u[j], loc = pars$loc[i], scale = pars$scale[i], shape = pars$shape[i], lower.tail = mdl$lower)
      })
    })
  } else {
    return(NULL)
  }

  # if parameters flipped for fitting, flip 'em back
  if(mdl$minima) erl <- -erl

  return(erl)
}


################################################################################################################################
#' Transform one nonstationary distribution into another (or into a stationary distribution)
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param fixed_cov data.frame containing values of the covariate defining a stationary distribution
#'
#' @return Vector of event magnitudes corresponding to the stationary covariates
#'
#' @export
#'
stransform <- function(mdl, fixed_cov = NA) {

  # current default is to map NS distribution to itself - could give a stationary default instead
  map_from_u(mdl, u = map_to_u(mdl), fixed_cov = fixed_cov)
}


################################################################################################################################
#' Estimate the return period of an event or events for given covariates
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param x Vector of data to be transformed
#' @param fixed_cov data.frame containing values of the covariate at which the location is to be evaluated. Default value is NA, in which case nonstationary parameters are returned for all values of the covariate(s) used to fit the model
#'
#' @return Vector of return periods
#'
#' @export
#'
return_period <- function(mdl, x, fixed_cov = NA) {
  1 / map_to_u(mdl, x, fixed_cov)
}

################################################################################################################################
#' Estimate the probability ratio of an event under a given model: the factor-change in the probability of exceeding a set threshold when changing from the counterfactual to the factual world
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param ev Scalar: magnitude of event
#' @param cov_f Data.frame of covariate values defining the factual climate
#' @param cov_cf Data.frame of covariate values defining the counterfactual climate
#'
#' @return Scalar representing how much more or less likely it is that we will exceed the value 'ev' in the factual world vs the counterfactual world
#'
#' @export
#'
prob_ratio <- function(mdl, ev, cov_f, cov_cf) {

  if(missing(ev)) ev <- mdl$ev

  ep_f <- map_to_u(mdl, ev, fixed_cov = cov_f)
  ep_cf <- map_to_u(mdl, ev, fixed_cov = cov_cf)

  ep_f / ep_cf
}

################################################################################################################################
#' Estimate the effective return level of a fixed return period under a nonstationary model
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param rp Scalar indicating the return period of interest
#' @param fixed_cov data.frame containing values of the covariate at which the return level is to be evaluated. Default value is NA, in which case effective return levels are returned for all values of the covariate(s) used to fit the model
#'
#' @return Vector of effective return levels corresponding to the covariates provided
#'
#' @export
#'
eff_return_level <- function(mdl, rp, fixed_cov = NA) {

  map_from_u(mdl, 1/rp, fixed_cov = fixed_cov)
}


################################################################################################################################
#' Estimate the change in intensity of events of a specified return period under a given model, between the counterfactual and factual worlds
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param rp Scalar: return period of interest
#' @param cov_f Data.frame of covariate values defining the factual climate
#' @param cov_cf Data.frame of covariate values defining the counterfactual climate
#' @param relative Boolean: return change as relative (F) or absolute (T)? Default is F
#'
#' @return Scalar representing the estimated change in intensity
#'
#' @export
#'
int_change <- function(mdl, rp = NA, cov_f, cov_cf, relative = F) {

  if(is.na(rp)) {
    if(relative) {
      cat("Return period needed to calculate relative change")
      return(NA)
    } else {
      rp <- 10
    }
  }
  # if return period is less than 1, assume it's an exceedance probability and convert to a return period
  if(rp < 1) { rp <- 1/rp } else { rp <- rp }

  # get effective return levels
  rl <- eff_return_level(mdl, rp, fixed_cov = cov_f)
  rl_cf <- eff_return_level(mdl, rp, fixed_cov = cov_cf)

  # if variable is logged, convert to real values first
  if(substr(mdl$varnm, 1, 5) == "log10") {
    rl <- 10^rl
    rl_cf <- 10^rl_cf
  } else if (substr(mdl$varnm, 1, 3) == "log"){
    rl <- exp(rl)
    rl_cf <- exp(rl_cf)
  }

  if(relative) {
    (rl - rl_cf) / rl_cf * 100
  } else {
    rl - rl_cf
  }
}
