#' Getis and Ord's Gi* Cluster and Significance Map
#'
#' @description
#' Create the Getis Gi* Cluster Map and the corresponding Significance Map.
#' Maps are done calculating the Local Gi* (localG - spdep) for each
#' spatial unit and testing its significance.
#'
#' @param x variable to create cluster and significance map
#' @param listw a neighbours list with spatial weights. From package spdep:
#' a listw object. Use poly2nb (class nb)
#' and nb2listw (class listw, nb) from package spdep. Can be any type of listw
#' object, for instance, rook contiguity (common edge) or queen contiguity (common
#' edge or common vertex)
#' @param zero.policy by default = NULL, if FALSE stop with error for any empty
#' neighbour sets, if TRUE permit the weights list to be formed with zero-length
#' weights vectors. Parameter inherited from the spdep package.
#' @param shp the spatial dataset: sf or SpatialPolygonsDataFrame (spdep)
#' @param significant by default is TRUE, if FALSE the significant map is not created
#' @param ... other parameters similar to original
#'
#' @return one or two maps
#'
#' @details
#' Using the function localG (spdep) create the Getis Gi* Cluster Map and the
#' corresponding Significance Map.
#' The significance map is done testing the null hypothesis (Ho) of zero spatial
#' autocorrelation for each spatial unit, then plotting a choropleth map with this
#' legend values: (Not Significant, p-value=0.05, p-value= 0.01, p-value=0.001,
#' p-value=0.0001, and Neighborless). Most significant clustered spatial units are
#' those with p-values smaller than 0.0001. Not significant
#' clustered spatial units are those with p-values grather than 0.05. Gi* Cluster Map
#' is done based on the significance map, but the choropleth legend is different (Not
#' - Significant, High-High, Low-Low, Low-High, High-Low, and Neighborless).
#'
#' @seealso
#' \itemize{
#'   \item Bivariate Moran's I: \code{\link{moran.bi}}
#'   \item Plot Bivariate Moran's I: \code{\link{moranbi.plot}}
#'   \item Bivariate Moran's I Test: \code{\link{moranbi.test}}
#'   \item Bivariate Local Moran's I and Test: \code{\link{localmoran.bi}}
#'   \item Create object "nb": \code{\link[spdep]{poly2nb}}
#'   \item Create object "listw"/"nb": \code{\link[spdep]{nb2listw}}
#' }
#'
#' @section Links:
#'
#' \enumerate{
#'    \item \href{https://pysal.readthedocs.io/en/v1.11.0/users/tutorials/autocorrelation.html}{Spatial Autocorrelation}
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
#' getis.cluster(columbus$CRIME, a.lw, zero.policy = FALSE, columbus, significant=TRUE)
getis.cluster <-  function(x, listw, zero.policy = NULL, shp, significant=TRUE, ...){
  get.dat <- localG(x, listw, zero.policy)

  get.dat1<-data.frame(as.vector(get.dat))
  names(get.dat1)<-c("Z.Gi")


  get.dat1$cluster<-"UN"
  get.dat1$cluster[get.dat1[,"Z.Gi"]>=qnorm(0.975,mean=0,sd=1)]<-"HH"     # both z scores are "high"
  get.dat1$cluster[get.dat1[,"Z.Gi"]<=qnorm(0.025,mean=0,sd=1)]<-"LL"     # both z scores are "low"
  get.dat1$cluster[is.na(get.dat1[,"Z.Gi"])]<-"NA"

  cols<-c(brewer.pal(5, "RdBu"),"#BEBEBE")
  get.dat1$col[get.dat1$cluster=="UN"]<-cols[3]
  get.dat1$col[get.dat1$cluster=="HH"]<-cols[1]
  get.dat1$col[get.dat1$cluster=="LL"]<-cols[5]
  get.dat1$col[get.dat1$cluster=="NA"]<-cols[6]
  get.dat1

  par(pty="s",mar=c(0,0,2,0), oma=c(0,0,0,0))
  P1 <- plot(st_geometry(shp), col=get.dat1$col, main= deparse(substitute(x)), ...)
  #
  #legend(locator(1), legend=c(paste("Not Significant  ",  "(",length(get.dat1$cluster[get.dat1$cluster=="UN"]),")",sep="",collapse=""), paste("High-High  ",  "(",length(get.dat1$cluster[get.dat1$cluster=="HH"]),")",sep="",collapse=""), paste("Low-Low  ",  "(",length(get.dat1$cluster[get.dat1$cluster=="LL"]),")",sep="",collapse=""), paste("Neighborless  ","(",length(get.dat1$cluster[get.dat1$cluster=="NA"]),")",sep="",collapse="")), fill=c(cols[3],cols[1],cols[5],cols[6]), title = "Gi* Cluster Map", bty="n", cex=1.2, y.intersp=0.8)
  legend("topleft",
         legend=c(paste("Not Significant  ",  "(",length(get.dat1$cluster[get.dat1$cluster=="UN"]),")",sep="",collapse=""),
                  paste("High-High  ",  "(",length(get.dat1$cluster[get.dat1$cluster=="HH"]),")",sep="",collapse=""),
                  paste("Low-Low  ",  "(",length(get.dat1$cluster[get.dat1$cluster=="LL"]),")",sep="",collapse=""),
                  paste("Neighborless  ","(",length(get.dat1$cluster[get.dat1$cluster=="NA"]),")",sep="",collapse="")),
         fill=c(cols[3],cols[1],cols[5],cols[6]),
         title = "Gi* Cluster Map", bty="n", cex=1.2, y.intersp=0.8)

  if (significant) {
  get.dat1$prob<-"UN"
  get.dat1$prob[get.dat1[,"Z.Gi"]>=qnorm(0.975,mean=0,sd=1)|get.dat1[,"Z.Gi"]<=qnorm(0.025,mean=0,sd=1)]<-"5%"
  get.dat1$prob[get.dat1[,"Z.Gi"]>=qnorm(0.995,mean=0,sd=1)|get.dat1[,"Z.Gi"]<=qnorm(0.005,mean=0,sd=1)]<-"1%"
  get.dat1$prob[get.dat1[,"Z.Gi"]>=qnorm(0.9995,mean=0,sd=1)|get.dat1[,"Z.Gi"]<=qnorm(0.0005,mean=0,sd=1)]<-"0.1%"
  get.dat1$prob[get.dat1[,"Z.Gi"]>=qnorm(0.99995,mean=0,sd=1)|get.dat1[,"Z.Gi"]<=qnorm(0.00005,mean=0,sd=1)]<-"0.01%"
  get.dat1$prob[is.na(get.dat1[,"Z.Gi"])]<-"NA"

  colsp<-c(brewer.pal(5, "Greens"),"#BEBEBE")
  get.dat1$col1[get.dat1$prob=="UN"]<-colsp[1]
  get.dat1$col1[get.dat1$prob=="5%"]<-colsp[2]
  get.dat1$col1[get.dat1$prob=="1%"]<-colsp[3]
  get.dat1$col1[get.dat1$prob=="0.1%"]<-colsp[4]
  get.dat1$col1[get.dat1$prob=="0.01%"]<-colsp[5]
  get.dat1$col1[is.na(get.dat1[,"Z.Gi"])]<-colsp[6]
  get.dat1

  #dev.new()
  par(pty="s",mar=c(0,0,2,0), oma=c(0,0,0,0))
  P2 <- plot(st_geometry(shp), col=get.dat1$col1, main= deparse(substitute(x)), ...)
  #legend(locator(1), legend=c(paste("Not Significant  ", "(",length(get.dat1$prob[get.dat1$prob=="UN"]),")",sep="",collapse=""), paste("p=0.05  ",  "(",length(get.dat1$prob[get.dat1$prob=="5%"]),")",sep="",collapse=""), paste("p=0.01  ",  "(",length(get.dat1$prob[get.dat1$prob=="1%"]),")",sep="",collapse=""),
  #paste("p=0.001  ",  "(",length(get.dat1$prob[get.dat1$prob=="0.1%"]),")",sep="",collapse=""), paste("p=0.0001  ",  "(",length(get.dat1$prob[get.dat1$prob=="0.01%"]),")",sep="",collapse=""), paste("Neighborless  ",   "(",length(get.dat1$cluster[get.dat1$cluster=="NA"]),")",sep="",collapse="")), bty="n", fill=colsp, title = "Gi* Significance Map", cex=1.2, y.intersp=0.8)
  legend("topleft",
         legend=c(paste("Not Significant  ", "(",length(get.dat1$prob[get.dat1$prob=="UN"]),")",sep="",collapse=""),
                  paste("p=0.05  ",  "(",length(get.dat1$prob[get.dat1$prob=="5%"]),")",sep="",collapse=""),
                  paste("p=0.01  ",  "(",length(get.dat1$prob[get.dat1$prob=="1%"]),")",sep="",collapse=""),
                  paste("p=0.001  ",  "(",length(get.dat1$prob[get.dat1$prob=="0.1%"]),")",sep="",collapse=""),
                  paste("p=0.0001  ",  "(",length(get.dat1$prob[get.dat1$prob=="0.01%"]),")",sep="",collapse=""),
                  paste("Neighborless  ",   "(",length(get.dat1$cluster[get.dat1$cluster=="NA"]),")",sep="",collapse="")),
                  bty="n", fill=colsp,
          title = "Gi* Significance Map", cex=1.2, y.intersp=0.8)
  }
}
