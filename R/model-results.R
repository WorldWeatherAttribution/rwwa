# SUMMARISE MODEL RESULTS
#
#' Get best estimates of parameters from fitted model
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param cov_f Data.frame with one row containing named covariates defining the factual climate
#' @param cov_cf Data.frame with one or more rows containing named covariates defining the counterfactual climate
#' @param ev (Optional) scalar: magnitude of the event of interest. If not provided, event value is picked up from the fitted model
#' @param rp (Optional) scalar: fixed return period of interest
#'
#' @return Vector containing estimates of all model parameters and quantities of interest
#'
#' @export
#'
mdl_ests <- function(mdl, cov_f, cov_cf, ev, rp = NA) {

  if(nrow(cov_f) > 1) {
    print("cov_f has more than one row: only first row will be used as factual covariates")
    cov_f <- cov_f[1,,drop = F]
  }

  if(!all(c(sapply(mdl$covnm, function(cnm) cnm %in% colnames(cov_f)), sapply(mdl$covnm, function(cnm) cnm %in% colnames(cov_cf))))) {
    print("Not all model covariates appear in factual/counterfactual covariates: missing covariates will be assumed to be zero throughout")
  }

  pars <- mdl$par
  current_pars <- ns_pars(mdl, fixed_cov = cov_f)
  disp <- current_pars$scale / current_pars$loc

  if(is.na(rp)) {
    # if no fixed RP is given, estimate it from the event value
    if(missing(ev)) ev <- mdl$ev
    rp <- return_period(mdl, ev, fixed_cov = cov_f)
  } else {
    # if fixed RP given, use it to bootstrap the expected magnitude of the event value
    ev <- eff_return_level(mdl, rp, fixed_cov = cov_f)
  }

  # loop over counterfactual covariates (if necessary) & get PRs and intensity changes
  if(nrow(cov_cf) == 1) {
      changes <- c("PR" = prob_ratio(mdl, ev, cov_f, cov_cf),
                   "dI_abs" = tryCatch(int_change(mdl, rp, cov_f, cov_cf, relative = F), error = function(cond) {return(NA)}),
                   "dI_rel" = tryCatch(int_change(mdl, rp, cov_f, cov_cf, relative = T), error = function(cond) {return(NA)}))
  } else {
      changes <- unlist(lapply(rownames(cov_cf), function(rnm) {
          pr <- prob_ratio(mdl, ev, cov_f, cov_cf[rnm,,drop = F])
          di_abs <- tryCatch(int_change(mdl, rp, cov_f, cov_cf[rnm,,drop = F], relative = F), error = function(cond) {return(NA)})
          di_rel <- tryCatch(int_change(mdl, rp, cov_f, cov_cf[rnm,,drop = F], relative = T), error = function(cond) {return(NA)})
          setNames(c(pr, di_abs, di_rel), paste0(c("PR", "dI_abs", "dI_rel"), "_",rnm))
      }))
  }

  return(c(mdl$par, "disp" = disp, "event_magnitude" = ev, "return_period" = rp, changes))
}


################################################################################################################################
#' Bootstrapped confidence intervals of parameter estimates for a fitted model
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param cov_f Data.frame with one row containing named covariates defining the factual climate
#' @param cov_cf Data.frame with one or more rows containing named covariates defining the counterfactual climate
#' @param ev (Optional) scalar: magnitude of the event of interest. If not provided, event value is picked up from the fitted model
#' @param rp (Optional) scalar: fixed return period of interest
#' @param seed Scalar: seed to be used to initialise random sample for bootstrapped confidence intervals (if using)
#' @param nsamp Scalar: number of bootstrap samples to be used to estimate confidence intervals for location parameter. Set to NA if no confidence intervals are required. Default is 500.
#' @param ci Scalar from 0 to 1 defining width of confidence interval. Default is 0.95
#' @param return_sample Boolean: return confidence interval (F) or full bootstrap sample (T)? Default is to return the interval (F).
#'
#' @return Data.frame containing estimates of all model parameters and quantities of interest, along with limits of central confidence interval; or matrix of bootstrapped values of each quantity (can be useful if 95% interval is unstable)
#'
#' @export
#'
boot_ci <- function(mdl, cov_f, cov_cf, ev, rp = NA, seed = 42, nsamp = 500, ci = 0.95, return_sample = F) {

  alpha <- 1-ci

  # if event value not provided, get stored value from fitted model
  if(missing(ev)) ev <- mdl$ev

  # remove any extraneous variables from covariate tables
  cov_f <- cov_f[,mdl$covnm, drop = F]
  cov_cf <- cov_cf[,mdl$covnm, drop = F]

  # trim factual covariates if necessary
  if(nrow(cov_f) > 1) {
    print("cov_f has more than one row: only first row will be used as factual covariates")
    cov_f <- cov_f[1,,drop = F]
  }

  # if any covariates not provided, replace with 0 (otherwise will print error message for every bootstrap sample)
  for (cnm in mdl$covnm) {
    if(!cnm %in% colnames(cov_f)) {
      cat(cnm,"missing from factual covariates: assumed to be zero\n\n")
      cols_f <- colnames(cov_f)
      cov_f <- cbind(cov_f, 0)
      colnames(cov_f) <- c(cols_f, cnm)
    }
    if(!cnm %in% colnames(cov_cf)) {
      cat(cnm,"missing from counterfactual covariates: assumed to be zero\n\n")
      cols_cf <- colnames(cov_cf)
      cov_cf <- cbind(cov_cf, 0)
      colnames(cov_cf) <- c(cols_cf, cnm)
    }
  }

  obs_res <- mdl_ests(mdl, cov_f, cov_cf, ev, rp = rp)

  # get bootstrap sample
  set.seed(seed)
  boot_res <- sapply(1:nsamp, function(i) {
    boot_df <- mdl$data[sample(1:nrow(mdl$data), replace = T),]
    tryCatch({
      boot_mdl <- refit(mdl, new_data = boot_df)
      mdl_ests(boot_mdl, cov_f, cov_cf, ev = ev, rp = rp)
    },
    error = function(cond) {return(rep(NA, length(obs_res)))})
  })
  if(return_sample) {
    return(boot_res)
  } else {
    boot_qq <- t(rbind("est" = obs_res, apply(boot_res, 1, quantile, c(alpha/2, 1-(alpha/2)), na.rm = T)))
    return(boot_qq)
  }
}


################################################################################################################################
#' Bootstrapped confidence intervals of parameter estimates for a fitted model, tailored to WWA requirements for evaluation, attribution and projection of results from climate models
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param rp Scalar: fixed return period of interest
#' @param cov_f Data.frame with one row containing named covariates defining the factual climate
#' @param cov_hist Data.frame with one or more rows containing named covariates defining the historic counterfactual climate
#' @param cov_fut (optional) Data.frame with one or more rows containing named covariates defining the future counterfactual climate
#' @param y_start Integer: first year for which observations are available, used to determine the start of the evaluation period. Default is 1979.
#' @param y_now Integer: last year for which observations are available, used to determine the start of the evaluation/attribution period. Default is to use the current year.
#' @param y_fut Integer: last year to be used in the projection model. Default is 2050.
#' @param seed Scalar: seed to be used to initialise random sample for bootstrapped confidence intervals (if using)
#' @param nsamp Scalar: number of bootstrap samples to be used to estimate confidence intervals for location parameter. Set to NA if no confidence intervals are required. Default is 500.
#' @param ci Scalar from 0 to 1 defining width of confidence interval. Default is 0.95
#' @param return_sample Boolean: return confidence interval (F) or full bootstrap sample (T)? Default is to return the interval (F).
#' @param di_relative Boolean: force the return of relative (F) or absolute (T) change in intensity. Default is to choose based on fit type (shift = absolute, fixeddisp = relative)
#'
#' @return Data.frame containing estimates of all model parameters and quantities of interest, along with limits of central confidence interval; or matrix of bootstrapped values of each quantity (can be useful if 95% interval is unstable)
#'
#' @export
#'
cmodel_results <- function(mdl, rp = 10, cov_f, cov_hist, cov_fut,
                           y_start = 1979, y_now = as.integer(substr(Sys.Date(),1,4)), y_fut = 2050,
                           nsamp = 5, seed = 42, ci = 0.95, return_sample = F, di_relative = NA) {

  alpha <- 1-ci
  set.seed(seed)

  # fill in missing parameters
  if(is.na(di_relative)) di_relative <- mdl$type == "fixeddisp"

  # remove any extraneous variables from covariate tables
  cov_f <- cov_f[,mdl$covnm, drop = F]
  cov_hist <- cov_hist[,mdl$covnm, drop = F]

  # trim factual covariates if necessary
  if(nrow(cov_f) > 1) {
    print("cov_f has more than one row: only first row will be used as factual covariates")
    cov_f <- cov_f[1,,drop = F]
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # update model for evaluation & attribution

  df <- mdl$data

  # fit models to subsets of the data
  mdl_eval <- refit(mdl, new_data = df[df$year >= y_start & df$year <= y_now,])
  mdl_attr <- refit(mdl, new_data = df[df$year <= y_now,])

  # get return level to use for analysis
  event_rl <- eff_return_level(mdl_attr, rp = rp, fixed_cov = cov_f)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Bootstrap each set of model results

  if(mdl$type == "fixeddisp") { key_par <- "disp" } else { key_par <- "sigma0" }
  if(mdl$dist == "gev") { key_par <- c(key_par, "shape") }
  if(di_relative) { di_cnm <- "dI_rel" } else { di_cnm <- "dI_abs" }

  # get bootstrapped intervals, select only the elements of interest
  ci_eval <- boot_ci(mdl_eval, cov_f = cov_f, cov_cf = cov_hist, ev = event_rl, rp = rp, nsamp = nsamp)[key_par,,drop = F]
  ci_attr <- boot_ci(mdl_attr, cov_f = cov_f, cov_cf = cov_hist, ev = event_rl, rp = rp, nsamp = nsamp)
  ci_attr <- ci_attr[grepl(paste0("PR|",di_cnm), rownames(ci_attr)),]

  # flatten & rename
  ci_eval <- unlist(lapply(rownames(ci_eval), function(cnm) setNames(ci_eval[cnm,], paste("eval", gsub("_", "-", cnm), c("est", "lower", "upper"), sep = "_"))))
  ci_attr <- unlist(lapply(rownames(ci_attr), function(cnm) setNames(ci_attr[cnm,], paste("attr", gsub("_", "-", cnm), c("est", "lower", "upper"), sep = "_"))))

  # compile results so far
  res <- c(ci_eval, "rp_value" = event_rl, ci_attr)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # bootstrap model results (if future covariates given)
  if(!missing(cov_fut)) {

    # drop any extraneous covariates
    cov_fut <- cov_fut[,mdl$covnm, drop = F]

    # refit the model
    mdl_proj <- refit(mdl, new_data = df[df$year <= y_fut,])

    # bootstrap results
    ci_proj <- boot_ci(mdl_proj, cov_f = cov_f, cov_cf = cov_fut, ev = event_rl, rp = rp, nsamp = nsamp)
    ci_proj <- ci_proj[grepl(paste0("PR|",di_cnm), rownames(ci_proj)),]

    # invert future projections
    ci_proj[grepl("PR", rownames(ci_proj)),] <- 1/ci_proj[grepl("PR", rownames(ci_proj)),c(1,3,2)]
    ci_proj[grepl(di_cnm, rownames(ci_proj)),] <- -ci_proj[grepl(di_cnm, rownames(ci_proj)), c(1,3,2)]

    ci_proj <- unlist(lapply(rownames(ci_proj), function(cnm) setNames(ci_proj[cnm,], paste("proj", gsub("_", "-", cnm), c("est", "lower", "upper"), sep = "_"))))

    res <- c(res, ci_proj)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # reshape & relabel results
  res <- t(data.frame(res))
  rownames(res) <- paste0(mdl$varnm, " ~ ", paste(mdl$covnm, collapse = " + "), " (rp", rp,")")
  return(res)
}
