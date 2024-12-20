# SYNTHESIS OF ATTRIBUTION RESULTS FROM OBSERVATIONS AND CLIMATE MODELS
#
# Methods to synthesise fitted model results from observations and climate models
#
#' Support function to calculate the weighted mean of climate model results
#'
#' @param data Data.frame with one row per climate model, with columns containing best estimate and upper and lower bound of 95pc confidence interval for the quantity of interest
#' @param sig_mod Scalar: model representation error. Default value is 0.
#'
#' @return Vector of best estimate and lower and upper 95pc confidence bounds for weighted model mean, with and without model representation error
#'
#' @keywords internal
#' @export
#'
getsynmean <- function(data, sig_mod = 0) {

  # calculate weight for each model based on inverse variance
  w = 1/(((data$upper - data$lower)/(2*1.96))^2 + sig_mod^2)
  w1 = sum(w)

  # weighted mean
  s1 <- sum(w*data$est) / w1

  # get weighted interval by adding variances
  sig_lower = sqrt(sum(w * (((data$est - data$lower)/1.96)^2 + sig_mod^2)) / w1)
  sig_upper = sqrt(sum(w * (((data$est - data$upper)/1.96)^2 + sig_mod^2)) / w1)

  return(setNames(s1 + c(0, -1.96*sig_lower, +1.96*sig_upper), c("est", "lower", "upper")))
}


################################################################################################################################
#' Support function to calculate the chi^2 function used to estimate sig_mod
#'
#' @param data Data.frame with one row per climate model, with columns containing best estimate and upper and lower bound of 95pc confidence interval for the quantity of interest
#' @param sig_mod Scalar: model representation error. Default value is 0
#'
#' @return Vector of best estimate and 95pcconfidence bounds for weighted model mean, with and without model representation error
#'
#' @keywords internal
#' @export
#'
getsynchi2 <- function(data, sig_mod = 0) {

  # function to be minimized by finding sig_mod such that chi^2/mdof ~= 1

  # get best estimate of weighted mean
  s1 <- getsynmean(data, sig_mod)["est"]

  # compute chi2 by converting model intervals to standard deviations & adding sig_mod adjustment
  chi2 <- sum(apply(data[,c("est", "lower", "upper")], 1, function(x) {
    if(x["est"] > s1) {
      (x["est"] - s1)^2 / (((x["est"] - x["lower"])/1.96)^2 + sig_mod^2)
    } else {
      (x["est"] - s1)^2 / (((x["est"] - x["upper"])/1.96)^2 + sig_mod^2)
    }
  }))
  return(chi2)
}


################################################################################################################################
#' Produce synthesis of attribution results from observations and climate models
#'
#' @param obs_in Data.frame with one row per observational dataset, with columns "est", "lower" and "upper" giving the best estimate and lower and upper bounds for the quantity of interest. Dataset names should be given as rownames. Default is NA (no observations), in which case only the model synthesis is carried out.
#' @param models_in Data.frame with one row per climate model, with columns "est", "lower" and "upper" giving the best estimate and lower and upper bounds for the quantity of interest. Model names should be given as rownames.
#' @param synth_type String defining the type of synthesis to carry out. Options are 'abs' (absolute changes); 'rel' (percentage changes); and 'PR' (probability ratios). Default is 'abs'.
#'
#' @return List containing synth_type; sig_obs (scalar value indicating obs representation error); sig_mod (scalar value indicating model representation error); chi2/dof (initial estimate of ratio of chi^2 to model DOF); and df, a data.frame containing the synthesised results.
#'
#' @export
#'
synthesis <- function(obs_in = NA, models_in, synth_type = "abs") {

  if(is.na(unlist(obs_in))[1]) {
    no_obs <- T
    # create a dummy dataframe to avoid having to rewrite everything twice
    obs_in <- data.frame("est" = 0, "lower" = 0, "upper" = 0)
    rownames(obs_in) <- "dummy"
  } else {
    no_obs <- F
  }

  # relabel the data for easier reference later
  colnames(obs_in) <- colnames(models_in) <- c("est", "lower", "upper")

  if(!("model" %in% colnames(obs_in))) obs_in$model <- rownames(obs_in)
  if(!("model" %in% colnames(models_in))) models_in$model <- rownames(models_in)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(!synth_type %in% c("abs", "rel", "PR")) {
    cat(paste0("Synthesis type '",synth_type,"' not implemented - must be abs, rel or PR"), "\n")
  }

  if(synth_type == "PR") {
    obs_in[,c("est", "lower", "upper")] <- log(obs_in[,c("est", "lower", "upper")])
    models_in[,c("est", "lower", "upper")] <- log(models_in[,c("est", "lower", "upper")])
  } else if(synth_type == "rel") {
    obs_in[,c("est", "lower", "upper")] <- log(1+obs_in[,c("est", "lower", "upper")]/100)
    models_in[,c("est", "lower", "upper")] <- log(1+models_in[,c("est", "lower", "upper")]/100)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # get mean of intervals to estimate natural variability component

  # compute representation error from scatter of mean ($\sigma_{rep}$ in the paper)
  nobs = nrow(obs_in)
  obs <- apply(obs_in[,c("est", "lower", "upper"),drop = F], 2, mean)
  if(nobs == 1) {
    sig_obs = 0
  } else {
    s2 = sum((obs_in$est - obs[1])^2)
    sig_obs = sqrt(s2/(nobs-1))
  }

  # add representation error to individual observations
  obs_in$l_wb <- obs_in$est - sqrt((obs_in$est - obs_in$lower)**2 + (1.96*sig_obs)**2)
  obs_in$u_wb <- obs_in$est + sqrt((obs_in$est - obs_in$upper)**2 + (1.96*sig_obs)**2)

  # apply representation error to obs synthesis
  # we're working with confidence intervals here, so we extend them by adding (1.96sig_obs)^2 in quadrature
  obs[2] <- obs[1] - sqrt( (obs[1] - obs[2])**2 + (1.96*sig_obs)**2 )
  obs[3] <- obs[1] + sqrt( (obs[1] - obs[3])**2 + (1.96*sig_obs)**2 )

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # get initial estimate of model mean & calculate chi^2
  chi2 <- getsynchi2(models_in, sig_mod = 0)
  mdof <- nrow(models_in)-1

  if ( chi2/mdof > 1 ) {
    # find sig_mod such that chi^2/dof = 1
    sig_mod <- optim(0, function(x) {(getsynchi2(models_in, sig_mod = x) - (nrow(models_in)-1))^2},
                     method = "Brent", lower = 0, upper = 5)$par
  } else {
    sig_mod <- 0
  }

  # get weighted model mean
  models <- getsynmean(models_in, sig_mod = sig_mod)

  # add representation error to individual models
  models_in$l_wb <- models_in$est - sqrt((models_in$est - models_in$lower)**2 + (1.96*sig_mod)**2)
  models_in$u_wb <- models_in$est + sqrt((models_in$est - models_in$upper)**2 + (1.96*sig_mod)**2)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # weighted mean of models & obs (coloured bar)
  w_obs <- unname((obs["upper"] - obs["lower"])^{-2})
  w_mod <- unname((models["upper"] - models["lower"])^{-2})

  wmean <- (w_obs * obs["est"] + w_mod * models["est"]) / (w_obs + w_mod)

  # get weighted interval by averaging variances
  sig_lower = sqrt((w_obs * ((obs["est"] - obs["lower"])/1.96)^2 + w_mod * ((models["est"] - models["lower"])/1.96)^2) / (w_obs + w_mod))
  sig_upper = sqrt((w_obs * ((obs["est"] - obs["upper"])/1.96)^2 + w_mod * ((models["est"] - models["upper"])/1.96)^2) / (w_obs + w_mod))
  synth <- setNames(c(wmean, wmean - 1.96*sig_lower, wmean + 1.96*sig_upper), c("est", "lower", "upper"))

  # unweighted mean of obs and models
  umean <- (obs["est"] +  models["est"]) / 2
  synth["l_wb"] <- umean - sqrt((obs["est"]-obs["lower"])^2 + (models["est"]-models["lower"])^2)/2
  synth["u_wb"] <- umean + sqrt((obs["est"]-obs["upper"])^2 + (models["est"]-models["upper"])^2)/2

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # add group labels for easier plotting
  obs_in <- cbind(obs_in, "group" = "obs")
  obs <- data.frame(t(c("model" = "Observations", "group" = "obs_synth", obs)))
  models_in <- cbind(models_in, "group" = "models")
  models <- data.frame(t(c("model" = "Models", "group" = "model_synth", models)))
  synth <- data.frame(t(c("model" = "Synthesis", "group" = "synth", synth)))

  # combine all the data together in one dataframe
  res <- rbind.fill(obs_in, obs, models_in, models, synth)[,c("group", "model", "est", "lower", "upper", "l_wb", "u_wb")]
  for(cnm in c("est", "lower", "upper", "l_wb", "u_wb")) { res[,cnm] <- as.numeric(res[,cnm]) }

  # if only dummy obs, remove
  if(no_obs) {
    # drop all rows that don't relate to models
    res <- res[grepl("model", res$group),]
    sig_obs <- NA
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # reverse any transformations applied
  if(synth_type == "PR") {
    res[,c("est", "lower", "upper", "l_wb", "u_wb")] <- exp(res[,c("est", "lower", "upper", "l_wb", "u_wb")])
    sig_obs <- exp(sig_obs)
    sig_mod <- exp(sig_mod)
    umean <- exp(umean)
  } else if(synth_type == "rel") {
    res[,c("est", "lower", "upper", "l_wb", "u_wb")] <- 100*(exp(res[,c("est", "lower", "upper", "l_wb", "u_wb")])-1)
    sig_obs <- 100*(exp(sig_obs)-1)
    sig_mod <- 100*(exp(sig_mod)-1)
    umean <- 100*(exp(umean)-1)
  }

  return(list(synth_type = synth_type, sig_obs = sig_obs, "chi2/dof" = chi2 / mdof, df = res, uw_mean = umean))
}


################################################################################################################################
#' Plot synthesis of attribution results from observations and climate models
#'
#' @param synth Either a list containing the results of a synthesis carrier out using the 'synthesis' function or a data.frame containing the raw output from a synthesis on the Climate Explorer.
#' @param xlim Vector of length 2 (optional) giving range of x coordinates. If not provided, will be estimated from the range of the synthesised values.
#' @param lwd Scalar: width of line to be used for each bar. Default is 10.
#' @param xlab String with which to label the x axis. Default is an empty string.
#' @param main String: main title of plot. Default is an empty string.
#' @param add_space Boolean: add a blank row between the observations, models and overall synthesis? Default is T.
#' @param log Boolean (optional unless using output from Climate Explorer): should x axis be plotted on a log scale? Default is NA, with value determined by the type of synthesis (PR = T, rel/abs = F).
#' @param hide_labels Boolean: Hide the model names on the y axis? Default is F (show labels)
#'
#' @export
#'
plot_synthesis <- function(synth, xlim, lwd = 10, xlab = "", main = "", add_space = T, log = NA, hide_labels = F) {

  gcols = c("obs" = adjustcolor("blue", 0.5),
            "obs_synth" = "blue",
            "models" = adjustcolor("red", 0.5),
            "model_synth" = "red",
            "synth" = "magenta")

  # determine whether to plot on log axes or not (assume not unless told otherwise)
  if(is.na(log)) {
    if(!is.null(synth$synth)) {
      if (synth$synth_type == "PR") {logaxs <- "x"} else {logaxs <- ""}
    } else {
      logaxs <- ""
    }
  } else {
    if(log) {logaxs <- "x"} else {logaxs <- ""}
  }

  if (is(synth, "list")) synth <- synth$df

  if (missing(xlim)) {
    if(logaxs == "x") {
      xlim <- exp(range(pretty(log(as.numeric(unlist(synth[,c("lower", "upper", "l_wb", "u_wb")]))))))
    } else {
      xlim <- range(pretty(as.numeric(unlist(synth[,c("lower", "upper", "l_wb", "u_wb")]))))
    }
  }

  # relabel groups if needed (eg. is using results from climate explorer)
  if(is.numeric(synth$group)) synth$group <- names(gcols)[synth$group]

  nobs <- sum(synth$group == "obs")
  nmod <- sum(synth$group == "models")

  if(add_space) {
    yy <- c(rev(0:nobs+nmod+4), rev(0:nmod+2), 0)
  } else {
    yy <- nrow(synth):1
  }

  if(logaxs == "x") {vline <- 1} else {vline <- 0}

  plot(0, type = "n", xlim = xlim, ylim = range(yy) + c(-0.5,0.5), log = logaxs,
       yaxt = "n", ylab = "", xlab = xlab, main = main)

  grid(ny = NA, col = adjustcolor("black", 0.1), lty = 1)
  abline(v = vline, lty = 2)

  gcols <- gcols[synth$group]

  segments(y0 = yy, x0 = synth$l_wb, x1 = synth$u_wb, lwd = lwd, col = "black", lend = 2)
  segments(y0 = yy, x0 = synth$l_wb, x1 = synth$u_wb, lwd = lwd-2, col = "white", lend = 2)
  segments(y0 = yy, x0 = synth$lower, x1 = synth$upper, lwd = lwd, col = gcols, lend = 2)
  points(synth$est, yy, pch = 21, bg = gcols, lwd = 2, cex = lwd/10)

  if(!hide_labels) axis(2, at = yy, labels = synth$model, las = 1)
}
