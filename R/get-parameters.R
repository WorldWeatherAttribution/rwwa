# EXTRACT PARAMETERS FROM MODEL
#
#' Support function to get location of nonstationary distribution at specified values of the covariates,
#' given the form and parameters of the nonstationary model
#'
#' @param fittype string defining the type of model fit to use - method is currently implemented for 'shift' and 'fixeddisp'
#' @param pars vector of named parameters of the nonstationary model
#' @param fixed_cov data.frame containing values of the covariate(s) at which the location is to be evaluated
#'
#' @return Return a list containing the location, scale and (if applicable) shape parameters corresponding to the given covariates
#'
#' @keywords internal
#' @export
#'
get.ns_pars <- function(fittype, pars, fixed_cov) {

  # support function only - better to use ns_pars with model

  if(missing(fittype)) {
    print("No fit type given")
    return()
  }

  if(missing(pars)) {
    print("No parameters given")
    return()
  }

  if(missing(fixed_cov)) {
    # set all covariates to zero
    fixed_cov <- setNames(data.frame(rep(0, sum(grepl("alpha", names(pars))))),
                          gsub("alpha_", "", names(pars)[(grepl("alpha", names(pars)))]))
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # calculate effect of all alpha parameters
  effect_a <- rowSums(matrix(sapply(names(fixed_cov), function(cnm) pars[paste0("alpha_",cnm)] * fixed_cov[,cnm]), nrow = nrow(fixed_cov)))
  effect_b <- rowSums(matrix(sapply(names(fixed_cov), function(cnm) pars[paste0("beta_",cnm)] * fixed_cov[,cnm]), nrow = nrow(fixed_cov)))

  if(fittype == "fixeddisp") {

    ns_const = exp(effect_a / pars["mu0"])
    loc = pars["mu0"] * ns_const
    scale = pars["sigma0"] * ns_const

  } else if(fittype == "shift") {

    loc = pars["mu0"] + effect_a
    scale = rep(pars["sigma0"], length(loc))

  } else if(fittype == "shiftscale") {

    loc = pars["mu0"] + effect_a
    scale = exp(pars["sigma0"] + effect_b)

  } else {
    cat(fittype, "not implemented")
    return()
  }

  # return the list of named parameters: location, scale, shape (if applicable)
  if("shape" %in% names(pars)) {
    return(lapply(list("loc" = loc, "scale" = scale, "shape" = rep(pars["shape"], length(scale))), unname))
  } else {
    return(lapply(list("loc" = loc, "scale" = scale), unname))
  }
}

####################################################################################################################
#' Get location, scale and (if relevant) shape of fitted nonstationary distribution at specified values of the covariates
#'
#' @param mdl list of parameters and attributes defining a fitted model, as returned by 'fit_ns'
#' @param fixed_cov data.frame containing values of the covariate at which the location is to be evaluated. Default value is NA, in which case nonstationary parameters are returned for all values of the covariate used to fit the model
#'
#' @return Return a list containing the location, scale and (if applicable) shape parameters corresponding to the given covariates
#'
#' @export
#'
ns_pars <- function(mdl, fixed_cov = NA) {

  # if no covariate value given, evaluate at all covariate values
  if(is.na(unlist(fixed_cov)[1])) fixed_cov <- mdl$cov

  # can either provide one covariate value per obs OR one fixed covariate
  if((nrow(fixed_cov) < length(mdl$x)) & (nrow(fixed_cov) > 1)) {
    print("Not implemented for multiple sets of fixed covariates: using first row of fixed_cov")
    fixed_cov <- fixed_cov[1,,drop = F]
  }

  # if fixed_cov contains extra covariates, drop them
  fixed_cov <- fixed_cov[,sapply(colnames(fixed_cov), function(cnm) cnm %in% mdl$covnm), drop = F]

  # if fixed_cov is missing any covariates, they will be assumed to be zero by support function
  if(!all(c(sapply(mdl$covnm, function(cnm) cnm %in% colnames(fixed_cov)), sapply(mdl$covnm, function(cnm) cnm %in% colnames(fixed_cov))))) {
    print("Not all model covariates appear in factual/counterfactual covariates: missing covariates will be assumed to be zero throughout")
  }

  # use support function to get parameters
  return(get.ns_pars(fittype = mdl$type, pars = mdl$par, fixed_cov = fixed_cov))
}
