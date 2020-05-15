#' Bivariate Moran's I
#'
#' @description
#' For lattice/polygon objects, sf or SpatialPolygonsDataFrame (spdep), calculate
#' Bivariate Moran's Index. Moran's I is used to calculate global spatial
#' autocorrelation. This statistic evaluate the existence of clusters in the spatial
#' context of two lattice data variables simultaneously.
#'
#' @param X first variable
#' @param Y second variable
#' @param listw a neighbours list with spatial weights. From package spdep:
#' a listw object. Use poly2nb (class nb)
#' and nb2listw (class listw, nb) from package spdep. Can be any type of listw
#' object, for instance, rook contiguity (common edge) or queen contiguity (common
#' edge or common vertex)
#' @param zero.policy by default = NULL, if FALSE stop with error for any empty
#' neighbour sets, if TRUE permit the weights list to be formed with zero-length
#' weights vectors. Similar meaning and values than parameter zero.policy
#' of \code{\link[spdep]{localmoran}}
#' @param adjust.n by default = TRUE. if FALSE the number of observations is not
#' adjusted for no neighbour observations, if TRUE, the number of observations is
#' adjusted. Parameter inherited from the spdep package
#'
#' @return a statistical index ranging from -1 to 1
#'
#' @details
#' Bivariate Moran's Index measure global spatial autocorrelation - how related the
#' values of two variables are based on the locations where they were measured.
#' Negative Spatial Autocorrelation occurs when the index is significant below of
#' -1/(N-1), where N is the number of rows (polygons), and represent dispersed
#' values. Positive Spatial Autocorrelation occurs when the index significant above
#' of -1/(N-1), where N is the number of rows (polygons), and represent values that
#' are clustered. Random (no spatial autocorrelation) is close to -1/(N-1), this is
#' values around zero when N is high (enough rows or polygons)
#'
#' @seealso
#' \itemize{
#'   \item Bivariate Moran's I Test: \code{\link{moranbi.test}}
#'   \item Plot Bivariate Moran's I: \code{\link{moranbi.plot}}
#'   \item Create object "nb": \code{\link[spdep]{poly2nb}}
#'   \item Create object "listw"/"nb": \code{\link[spdep]{nb2listw}}
#' }
#'
#' @section Links:
#'
#' \enumerate{
#'    \item Neighbours files with spatial weights: \href{http://geodacenter.github.io/workbook/4a_contig_weights/lab4a.html}{*.gal}
#'    \item \href{https://www.statisticshowto.datasciencecentral.com/morans-i/}{Univariate Moran's I}
#'    \item Moran's I in \href{https://es.wikipedia.org/wiki/I_de_Moran}{Wikipedia}
#'    \item \href{https://en.wikipedia.org/wiki/Indicators_of_spatial_association}{Autocorrelation Indicators}
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
#' moran.bi(columbus$CRIME,columbus$INC,a.lw,zero.policy =TRUE)
moran.bi <-  function(X,Y,listw,zero.policy = NULL,adjust.n = TRUE){
    if (!inherits(listw, "listw"))
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (!is.numeric(X))
        stop(paste(deparse(substitute(X)), "is not a numeric vector"))
    if (!is.numeric(Y))
        stop(paste(deparse(substitute(Y)), "is not a numeric vector"))
    if (is.null(zero.policy))
        zero.policy <- get("zeroPolicy", envir = 'spdep' %:::% '.spdepOptions')
    stopifnot(is.logical(zero.policy))
   wc <- spweights.constants(listw, zero.policy = zero.policy, adjust.n = adjust.n)
   n <- wc$n
   morans<-(n/(wc$S0))%*%((t(scale(X))%*%as.matrix(as_dgRMatrix_listw(listw))%*%scale(Y))/(t(scale(X))%*%scale(X)))
   return(as.vector(morans))
}
