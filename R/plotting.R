# PLOTTING FUNCTIONS
#
#' Plot a fitted trend over time
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param ev (Optional) scalar: magnitude of the event of interest. If not provided, event value is picked up from the fitted model
#' @param ev_year (Optional) scalar: year of the event of interest. If not provided, inferred from the fitted model
#' @param rp (Optional) vector of length two, setting return period for which effective return levels should be plotted. Default is c(6,40)
#' @param add_loess Boolean: add a Loess smoother to the plot? Default is F.
#' @param loess_col String: set colour to be used for Loess smoother (if using). Default is 'forestgreen'.
#' @param xlab (Optional) string: label for x axis: default is 'Year'
#' @param ylab (Optional) string: label for y axis: default is to use variable name
#' @param legend_pos String indicating location of legend: default is 'topleft'. Change to NA to remove legend.
#' @param main String: main title for plot. Default is to leave blank.
#' @param xlim (Optional) vector defining the limits of the x axis
#' @param ylim (Optional) vector defining the lower and upper limits of the y axis
#' @param lwd Scalar: line weight to be used in plotting
#'
#' @export
#'
plot_trend <- function(mdl, ev, ev_year, rp = c(6, 40), add_loess = F, loess_col = "forestgreen",
                       xlab = "Year", ylab = NA, legend_pos = "topleft", main = "", xlim = NA, ylim = NA, lwd = 2) {

  if(is.na(ylab)) {ylab <- mdl$varnm}
  if(is.na(unlist(xlim)[1])) { xlim <- range(mdl$data$year) }
  if(is.na(unlist(ylim)[1])) { ylim <- range(pretty(mdl$x)) }
  if(missing(ev)) { ev <- mdl$ev }
  if(missing(ev_year)) { ev_year <- mdl$data$year[which.min(abs(mdl$x - ev))] }

  # set up legend
  legend_labels = "Fitted value"
  legend_cols = "black"
  legend_lty = "solid"
  legend_lwd = lwd

  # modify legend if adding effective return levels
  rp <- unique(rp[!is.na(rp)])
  if(length(rp) > 0) {
    legend_labels = c(legend_labels, paste0("1-in-",rp,"-year event"))
    legend_cols = c(legend_cols, rep("blue",length(rp)))
    legend_lty = c(legend_lty, rep("solid",length(rp)))
    legend_lwd = c(legend_lwd,c(lwd,max(1,lwd -1))[1:length(rp)])
  }

  plot(mdl$data$year, mdl$x, type = "S", lwd = lwd, col = adjustcolor("black", 0.5), xlab = xlab,
       ylab = ylab, main = main, xlim = xlim, ylim = ylim)

  lines(mdl$data$year-0.5, ns_pars(mdl)$loc, col = adjustcolor("black", 1), lwd = lwd)
  lines(mdl$data$year-0.5, eff_return_level(mdl, rp[1]), type = "l", lty = 1, col = adjustcolor("blue", 1), lwd = lwd)
  lines(mdl$data$year-0.5, eff_return_level(mdl, rp[2]), type = "l", lty = 1, col = adjustcolor("blue", 1), lwd = max(1,lwd -1))

  # add a loess smoother
  if(add_loess) {
    lines(mdl$data$year, fitted(loess(formula(paste0(mdl$varnm," ~ year")), mdl$data)), col = loess_col, lwd = lwd, lty = "22")
    legend_labels <- c(legend_labels, "Loess smoothed")
    legend_cols <- c(legend_cols, loess_col)
    legend_lty <- c(legend_lty, "22")
    legend_lwd <- c(legend_lwd, lwd)
  }

  points(ev_year-0.5, ev, col = "magenta", lwd = 2, pch = 0)

  # add legend
  legend(legend_pos, legend = legend_labels, lty = legend_lty, col = legend_cols, lwd = legend_lwd)
}


################################################################################################################################
#' Plot fitted trend against a single covariate
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param xcov String indicating the covariate to plot on the x-axis. Must appear in mdl$cov.
#' @param trend_cov Data.frame containing values of the covariates to be used to plot the trend. Default value is NA, in which case the trend is estimated at all values of xcov with all other covariates held at their mean value
#' @param ci_cov (Optional) Data.frame containing values of the covariates at which confidence intervals for the location parameter should be estimated. Default is NA, in which case no confidence bounds are plotted.
#' @param ci_col String: set colour to be used for confidence bounds (if using). Default is 'black'.
#' @param ev (Optional) scalar: magnitude of the event of interest. If not provided, event value is picked up from the fitted model
#' @param ev_x (Optional) scalar: x-value against which to plot the event of interest. If not provided, event year is picked up from the fitted model
#' @param rp (Optional) vector of length two, setting return period for which effective return levels should be plotted. Default is c(6,40)
#' @param add_loess Boolean: add a Loess smoother to the plot? Default is F.
#' @param loess_col String: set colour to be used for Loess smoother (if using). Default is 'forestgreen'.
#' @param seed Scalar: seed to be used to initialise random sample for bootstrapped confidence intervals (if using)
#' @param nsamp Scalar: number of bootstrap samples to be used to estimate confidence intervals for location parameter. Set to NA if no confidence intervals are required. Default is 500.
#' @param xlim (Optional) vector defining the left and right limits of the x axes
#' @param ylim (Optional) vector defining the lower and upper limits of the y axes
#' @param xlab (Optional) string: label for x axis. Default is to use 'xcov'
#' @param ylab (Optional) string: label for y axis. Default is to leave blank.
#' @param legend_pos String indicating location of legend: default is 'topleft'. Change to NA to remove legend.
#' @param main String: main title for plot. Default is to leave blank.
#' @param lwd Scalar: line weight to be used in plotting. Default is 3
#'
#' @export
#'
plot_covtrend <- function(mdl, xcov, trend_cov = NA, ci_cov = NA,  ci_col = "black", ev, ev_x, rp = c(6,40), add_loess = F, loess_col = "forestgreen",
                          seed = 42, nsamp = 500, xlim = NA, ylim = NA, xlab = NA, ylab = NA, legend_pos = "topleft", main = "", lwd = 3) {

  if(is.na(xlab)) { xlab <- toupper(xcov)}
  if(is.na(ylab)) { ylab <- mdl$varnm}
  if(is.na(ylim[1])) { ylim <- range(pretty(mdl$x)) }
  if(missing(ev)) { ev <- mdl$ev }

  x <- mdl$data[,xcov]
  if(missing(ev))    { ev <- mdl$ev }
  if(missing(ev_x)) {ev_x <- x[which(mdl$x == ev)]}
  o <- order(x)

  if(is.na(unlist(ci_cov)[1])) {
    # can't plot confidence intervals so set nsamp to NA
    nsamp <- NA
    xlims <- range(x)
  } else {
    # extend the x-axis to accommodate the CI covariates
    xlims <- range(pretty(c(ci_cov[,xcov], x)))

    # if only the x-covariate is provided, set the other covariates to the mean value for plotting
    for(cnm in mdl$covnm) {
      if(!cnm %in% colnames(ci_cov)) ci_cov[,cnm] <- mean(mdl$data[,cnm])
    }
  }
  if(is.na(unlist(xlim)[1])) xlim <- xlims

  if(is.na(unlist(trend_cov)[1])) {
    # if no plotting covariate provided, fix all covariates at mean value except for xcov
    trend_cov <- data.frame(sapply(mdl$covnm, function(cnm) if(cnm == xcov) {mdl$data[,cnm]} else {mean(mdl$data[,cnm])}, simplify = F))
  }

  # set up legend
  legend_labels = "Fitted value"
  legend_cols = "black"
  legend_lty = "solid"
  legend_lwd = lwd

  # modify legend if adding effective return levels
  rp <- unique(rp[!is.na(rp)])
  if(length(rp) > 0) {
    legend_labels = c(legend_labels, paste0("1-in-",rp,"-year event"))
    legend_cols = c(legend_cols, rep("blue",length(rp)))
    legend_lty = c(legend_lty, rep("solid",length(rp)))
    legend_lwd = c(legend_lwd,c(lwd,max(1,lwd -1))[1:length(rp)])
  }

  plot(x, mdl$x, pch = 20, main = main, xlab = "", ylab = "", ylim = ylim, xlim = xlim,
       col = adjustcolor("black", 0.6))
  mtext(xlab, side = 1, line = 2.5, cex = par("cex.lab"))
  mtext(ylab, side = 2, line = 2.5, cex = par("cex.lab"))

  points(ev_x, ev, col = "magenta", lwd = 2, pch = 0)

  # trend lines
  lines(x[o], ns_pars(mdl, fixed_cov = trend_cov)$loc[o], lwd = 3, col = "black", lty = 1)
  lines(x[o], eff_return_level(mdl, rp[1], fixed_cov = trend_cov)[o], col = "blue", lwd = 3, lty = 1)
  lines(x[o], eff_return_level(mdl, rp[2], fixed_cov = trend_cov)[o], col = "blue", lwd = 2, lty = 1)

  # get confidence interval for mu' (if not required, set ci_cov to NA)
  if(!is.na(nsamp)) {
    mdl_df <- mdl$data
    set.seed(seed)
    mu_ci <- apply(sapply(1:nsamp, function(i) {
      boot_df <- mdl_df[sample(1:nrow(mdl_df), nrow(mdl_df), replace = T),]
      tryCatch({
        boot_mdl <- refit(mdl, boot_df)
        sapply(rownames(ci_cov), function(rnm) ns_pars(boot_mdl, ci_cov[rnm,,drop = F])$loc)
      }, error = function(cond) {return(rep(NA, nrow(ci_cov)))})
    }), 1, quantile, c(0.025, 0.975), na.rm = T)

    # confidence interval & markers for mu' at factual & counterfactual covariates
    segments(x0 = ci_cov[,xcov], y0 = mu_ci["2.5%",], y1 = mu_ci["97.5%",], lwd = 3, col = ci_col, lend = 1)
    # matplot(ci_cov[,"gmst"], t(mu_ci), pch = 3, add = T, col = "red3") # line ends: not very elegant, so removed for now
    points(ci_cov[,xcov], sapply(rownames(ci_cov), function(rnm) ns_pars(mdl, ci_cov[rnm,,drop = F])$loc), pch = "_", col = ci_col, lwd = 2)
  }

  # add a loess smoother
  if(add_loess) {
    dfx <- mdl$data[order(mdl$data[,xcov]),]
    lines(dfx[,xcov], fitted(loess(formula(paste0(mdl$varnm," ~ ", xcov)), dfx)), col = loess_col, lwd = lwd, lty = "22")
    legend_labels <- c(legend_labels, "Loess smoothed")
    legend_cols <- c(legend_cols, loess_col)
    legend_lty <- c(legend_lty, "22")
    legend_lwd <- c(legend_lwd, lwd)
  }

  # add legend
  legend(legend_pos, legend = legend_labels, lty = legend_lty, col = legend_cols, lwd = legend_lwd, cex = par()$cex.lab)
}


################################################################################################################################
#' Goodness-of-fit plot comparing return levels in factual and counterfactual climates
#'
#' @param mdl List of attributes & parameters defining a nonstationary model, as returned by 'fit_ns'
#' @param cov_f Data.frame with one row containing named covariates defining the factual climate
#' @param cov_cf Data.frame with one row containing named covariates defining the counterfactual climate
#' @param ev (Optional) scalar: magnitude of the event of interest. If not provided, event value is picked up from the fitted model
#' @param seed Scalar: seed to be used to initialise random sample for bootstrapped confidence intervals (if using)
#' @param nsamp Scalar: number of bootstrap samples to be used to estimate confidence intervals for location parameter. Set to NA if no confidence intervals are required. Default is 500.
#' @param model_desc Boolean: add description of model to plot? Default is T.
#' @param xlim Vector defining the lower and upper limits of the x axes (default is c(1,10000))
#' @param ylim (Optional) vector defining the lower and upper limits of the y axes
#' @param pch Scalar determing the plotting character to be used. Default is 20.
#' @param xlab (Optional) string: label for x axis. Default is "Return period (years)".
#' @param ylab (Optional) string: label for y axis. Default is to use the covariate name.
#' @param main String: main title for plot. Default is to leave blank.
#' @param legend_pos String indicating location of legend: default is 'topleft'. Change to NA to remove legend.
#' @param legend_labels Vector of labels for legend: default is c("Present climate", "Counterfactual climate").
#'
#' @export
#'
plot_returnlevels <- function(mdl, cov_f, cov_cf, ev, seed = 42, nsamp = 500, model_desc = T,
                              xlim = c(1,10000), ylim = NA, pch = 20, xlab = "Return period (years)", ylab = NA, main = "",
                              legend_pos = "topright", legend_labels = c("Present climate", "Counterfactual climate")) {

  x <- mdl$x
  if(missing(ev)) { ev <- mdl$ev }

  rp_x <- unique(c(seq(1.1,2,0.1), seq(2,100,1), seq(100,1000,10), seq(100,1000,100), seq(1000,10000,1000)))     # return periods at which to calculate values for curves
  rp_th <- 1/seq(1,0,length.out = length(x)+2)[2:(length(x)+1)]                                                  # quantiles to map against observations to check fit

  # trim covariates if necessary
  if(nrow(cov_f) > 1) {
    print("cov_f has more than one row: only first row will be used as factual covariates")
    cov_f <- cov_f[1,,drop = F]
  }
  # trim covariates if necessary
  if(nrow(cov_cf) > 1) {
    print("cov_cf has more than one row: only first row will be used as counterfactual covariates")
    cov_cf <- cov_cf[1,,drop = F]
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # calculate return periods & return levels

  rl_curve_pres <- map_from_u(mdl, 1/rp_x, fixed_cov = cov_f)
  rl_curve_cf <- map_from_u(mdl, 1/rp_x, fixed_cov = cov_cf)

  rl_obs_pres <- map_from_u(mdl, map_to_u(mdl), fixed_cov = cov_f)
  rl_obs_cf <- map_from_u(mdl, map_to_u(mdl), fixed_cov = cov_cf)

  rp_event_pres <- 1/map_to_u(mdl, ev, fixed_cov = cov_f)
  rp_event_cf <- 1/map_to_u(mdl, ev, fixed_cov = cov_cf)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # prep axes

  if(is.na(ylim[1])) { ylim <- range(pretty(c(x, rl_curve_pres, rl_curve_cf))) }
  if(is.na(ylab)) {ylab <- mdl$varnm}

  # plot
  plot(0,type = "n", xlim = xlim, ylim = ylim, log = "x", xlab = "", ylab = "", main = main)
  mtext(xlab, side = 1, line = 2.5, cex = par("cex.lab"))
  mtext(ylab, side = 2, line = 2.5, cex = par("cex.lab"))

  # if model description is required, add it as legend title
  if(model_desc) {
    legend_title <- paste0(mdl$varnm, " ~ ", paste0(mdl$covnm, collapse = " + "), " (",mdl$dist, ", ", mdl$type, ")")
  } else {
    legend_title <- ""
  }

  # add legend
  legend(legend_pos, legend = c(legend_labels, "Observed event"), col = c("firebrick", "blue", "magenta"), lty = 1, pch = c(pch,pch,NA),
         bty = "n", cex = par()$cex.lab, title = legend_title)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # return period curves
  lines(rp_x, rl_curve_pres, lwd = 2, col = "firebrick", lty = 1)       # present climate
  lines(rp_x, rl_curve_cf, lwd = 2, col = "blue", lty = 1)              # counterfactual

  # expected return periods vs return levels transformed to stationarity at that covariate value
  points(rp_th, sort(rl_obs_pres, decreasing = mdl$lower), col = "firebrick", pch = pch)      # present
  points(rp_th, sort(rl_obs_cf, decreasing = mdl$lower), col = "blue", pch = pch)             # counterfactual

  # horizontal line showing observed event, plus ticks showing return periods
  abline(h = ev, col = "magenta", lty = 2)
  suppressWarnings(rug(rp_event_pres, lwd = 3, col = "firebrick"))   # present
  suppressWarnings(rug(rp_event_cf, lwd = 3, col = "blue"))          # counterfactual

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Add confidence intervals to return periods

  if(!is.na(nsamp)) {
    x_ci <- c(5,10,20,50,100,200,500,1000,2000,5000,10000)
    set.seed(seed)

    mdl_df <- mdl$data[,c(mdl$varnm, mdl$covnm)]
    boot_res <- sapply(1:nsamp, function(i) {
      boot_df <- mdl_df[sample(1:nrow(mdl_df), nrow(mdl_df), replace = T),]
      tryCatch({
        boot_mdl <- refit(mdl, boot_df)
        # print(boot_mdl$par)
        c(map_from_u(boot_mdl, 1/x_ci, fixed_cov = cov_f), map_from_u(boot_mdl, 1/x_ci, fixed_cov = cov_cf))
      }, error = function(cond) {return(rep(NA, length(x_ci)*2))})
    })
    est_ci <- apply(boot_res, 1, quantile, c(0.025, 0.975), na.rm = T)

    # shaded region for confidence intervals
    polygon(x = c(x_ci, rev(x_ci)), y = c(est_ci[1,1:length(x_ci)], rev(est_ci[2,1:length(x_ci)])), density = NULL, border = NA, col = adjustcolor("firebrick", 0.1))
    polygon(x = c(x_ci, rev(x_ci)), y = c(est_ci[1,-(1:length(x_ci))], rev(est_ci[2,-(1:length(x_ci))])), density = NULL, border = NA, col = adjustcolor("blue", 0.1))
  }
}
