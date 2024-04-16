# ADDITIONAL SUPPORT FUNCTIONS
#
#' Define available plotting region in Jupyter notebook (and set default graphical pars)
#'
#' @param rc Vector setting the number of rows and columns. Default is c(1,1)
#' @param w Scalar setting the width of each of the (r,c) plotting regions. Default is 4
#' @param h Scalar setting the height of each of the (r,c) plotting regions. Default is 4
#' @param res Scalar setting the resolution of the plots. Default is 200.
#' @param pch Set default plotting character. Default is 20 (small filled circles)
#' @param ... Additional graphical parameters to be passed to par()
#'
#' @export
#'
prep_window <- function(rc = c(1,1), w = 4, h = 4, res = 200, pch = 20, ...) {
  options(repr.plot.width = rc[2]*w, repr.plot.height = rc[1]*h, repr.plot.res = res)
  par(mfrow = rc, pch = pch, ...)
}


################################################################################################################################
#' Load a time series in Climate Explorer .dat format
#'
#' @param fnm String defining filename to be loaded
#' @param ... Additional parameters to be passed to read.csv (for example, used to specify col.names)
#'
#' @return Data.frame containing the .csv data
#'
#' @export
#'
load_ts <- function(fnm, ...) {
  read.csv(fnm, comment.char = "#", sep = " ", header = F, ...)
}


################################################################################################################################
#' Support function to show correlations between variables in panels when using pairs()
#'
#' @param x Vector of x values
#' @param y Vector of y values
#' @param digits Integer: number of decimal places to show correlations to. Default is 2.
#' @param ... Additional graphical paramaters to be passed to text() (eg. to set colour, font size)
#'
#' @keywords internal
#'
#' @export
#'
panel.cor <- function(x, y, digits = 2, ...) {
  usr <- par("usr")
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "pairwise.complete.obs")
  txt <- round(r, digits)
  cex.cor <- 1
  text(0.5, 0.5, txt, font = 2, ...)
  par(usr = usr)
}
