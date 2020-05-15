#' Plot the Bivariate Moran's Index
#'
#' @description
#' Create a scatter plot graphic with lagged variables, one variable in each axis,
#' where can be seen the trend of the spatial autocorrelation.
#'
#' @param x the first variable
#' @param y the second variable
#' @param listw a neighbours list with spatial weights. From package spdep:
#' a listw object. Use poly2nb (class nb)
#' and nb2listw (class listw, nb) from package spdep. Can be any type of listw
#' object, for instance, rook contiguity (common edge) or queen contiguity (common
#' edge or common vertex)
#' @param spChk by default = NULL. Similar meaning and values than parameter spChk
#' of \code{\link[spdep]{localmoran}}
#' @param labels by default = NULL, place labels (point ID) of potential influential
#' observations. If NULL or TRUE will place labels, if FALSE will remove.
#' @param xlab by default = NULL, change this string value if a different label for
#' horizontal axis is needed.
#' @param ylab by default = NULL, change this string value if a different label for
#' vertical axis is needed.
#' @param zero.policy by default = NULL, if FALSE stop with error for any empty
#' neighbours sets, if TRUE permit the weights list to be formed with zero-length
#' weights vectors. Similar meaning and values than parameter zero.policy
#' of \code{\link[spdep]{localmoran}}
#' @param quiet by default = NULL. If FALSE a table with potentially
#' influential observations will be displayed. The table will show next columns:
#' dfb.1_, dfb.x, dffit, cov.r, cook.d, and hat. If TRUE, not R command window output
#' will be displayed.
#' @param ...  other parameters similar to original
#' @return a graphic showing Bivariate Moran's Index behavior
#'
#' @details
#' Bivariate spatial autocorrelation can be seen graphically using spatial lags plots
#' of scaled variables. Method used to create spatial lag values in morambi.plot is
#' lag.listw, and is explained in pages 257 and s58 of Edzer J. Pebesma (2008) -- Edzer
#' J. Pebesma, Roger S. Bivand und Virgilio Gómez-Rubio. 2008. “Areal Data and Spatial
#' Autocorrelation.” In Applied Spatial Data Analysis with R, Online-Ausg., 237–72.
#' Springer New York. For graphic analysis, Variables need to have the same relative
#' scale. Beforehand use the 'scale' function of R base to perform the transformation.
#' In the resulting graphic, trend line depicts the spatial autocorrelation.
#' Random pattern represented in horizontal trend line means no spatial
#' autocorrelation. Cluster spatial autocorrelation is represented by a regression line
#' with positive slope. Dispersed spatial autocorrelation (spatial outliers) is shown by
#' a regression line with negative slope. In the graphic it is possible to identify
#' with different symbol and using labels, atypical observations that fall apart the
#' displayed pattern.
#'
#' @seealso
#' \itemize{
#'   \item Bivariate Moran's I: \code{\link{moran.bi}}
#'   \item Plot Bivariate Moran's I: \code{\link{moranbi.plot}}
#'   \item Bivariate Moran's I Test: \code{\link{moranbi.test}}
#'   \item LISA Cluster and Significance Map: \code{\link{moran.cluster}}
#'   \item Create object "nb": \code{\link[spdep]{poly2nb}}
#'   \item Create object "listw"/"nb": \code{\link[spdep]{nb2listw}}
#' }
#'
#' @section Links:
#'
#' \enumerate{
#'    \item \href{https://en.wikipedia.org/wiki/Indicators_of_spatial_association}{Indicators of Spatial Association}
#'    \item Spatial Autocorrelation \href{http://blogs.oregonstate.edu/geo599spatialstatistics/2016/06/08/spatial-autocorrelation-morans/}{(Moran’s I) Test}
#'    \item Moran's I \href{https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/}{Test}
#'    \item \href{http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Confidence_Intervals/BS704_Confidence_Intervals_print.html}{Confidence Intervals}
#' }
#'
#' @import spatialreg
#' @import spdep
#' @import sp
#' @import grDevices
#' @import graphics
#' @import stats
#' @export
#'
#' @examples
#' library(spdep)
#' columbus <- st_read(system.file("extdata/columbus.shp", package="spdepbi")[1], quiet=TRUE)
#' col_nbq <- poly2nb(columbus)
#' a.lw <- nb2listw(col_nbq, style="W")
#' CRIME <- as.vector(scale(columbus$CRIME))
#' INCOME <- as.vector(scale(columbus$INC))
#' moranbi.plot(CRIME,INCOME,quiet=FALSE,zero.policy=FALSE,listw=a.lw)
moranbi.plot <-  function(x, y, listw, spChk = NULL, labels = NULL, xlab = NULL, ylab = NULL,
           zero.policy = NULL, quiet = NULL, ...){ # zero.policy = NULL, quiet = NULL, ...)
  if (!inherits(listw, "listw"))
      stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (is.null(quiet))
    quiet <- !get("verbose", envir = 'spdep' %:::% '.spdepOptions')
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(is.logical(quiet))
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = 'spdep' %:::% '.spdepOptions')
  stopifnot(is.logical(zero.policy))
  xname <- deparse(substitute(x))
  yname <- deparse(substitute(y))
  if (!is.numeric(x))
    stop(paste(xname, "is not a numeric vector"))
  if (any(is.na(x)))
    stop("NA in X")
  if (!is.numeric(y))
    stop(paste(yname, "is not a numeric vector"))
  if (any(is.na(y)))
    stop("NA in Y")
  n <- length(listw$neighbours)
  if (n != length(x))
    stop("objects of different length")
  if (is.null(spChk))
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(array(c(x, y), dim=c(length(x),length(y))),listw))
    stop("Check of data and weights ID integrity failed")
  labs <- TRUE
  if (is.logical(labels) && !labels)
    labs <- FALSE
  if (is.null(labels) || length(labels) != n)
    labels <- as.character(attr(listw, "region.id"))
  wy <- lag.listw(listw, y, zero.policy = zero.policy)
  if (is.null(xlab))
    xlab <- xname
  if (is.null(ylab))
    ylab <- paste("spatially lagged", yname)
  plot(x, wy, xlab = xlab, ylab = ylab, ...)
  if (zero.policy) {
  n0 <- wy == 0
  if (any(n0)) {
  symbols(x[n0], wy[n0], inches = FALSE, circles = rep(diff(range(x))/50,
                                                           length(which(n0))), bg = "grey", add = TRUE)
    }
  }
  xwy.lm <- lm(wy ~ x)
  abline(xwy.lm)
  abline(h = mean(wy), lty = 2)
  abline(v = mean(x), lty = 2)
  infl.xwy <- influence.measures(xwy.lm)
  is.inf <- which(apply(infl.xwy$is.inf, 1, any))
  points(x[is.inf], wy[is.inf], pch = 9, cex = 1.2)
  if (labs)
    text(x[is.inf], wy[is.inf], labels = labels[is.inf],
         pos = 2, cex = 0.7)
  rownames(infl.xwy$infmat) <- labels
  if (!quiet)
  summary(infl.xwy)
  invisible(infl.xwy)
}
