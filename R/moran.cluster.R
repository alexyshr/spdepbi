#' LISA Cluster Map and Significance Map
#'
#' @description
#' Using the function localmoran (spdep) create Local Indicators of Spatial
#' Association (LISA) Cluster Map and corresponding Significance Map.
#' Maps are done calculating the Local Moran Index (localmoran- spdep) for each
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
#' @param shp the spatial dataset: sf SpatialPolygonsDataFrame (spdep)
#' @param significant by default is TRUE, if FALSE the significant map is not created
#' @param alternative by default is two.sided. Type of alternative hypothesis test.
#' Other values are 'less' or 'greater'.
#' @param ... other parameters similar to internal function
#'
#' @return one or two maps
#'
#' @details
#' Using the function localmoran (spdep) create the Local Indicators of Spatial
#' Association - LISA Cluster Map and the corresponding Significance Map.
#' The significance map is done testing the null hypothesis (Ho) of zero spatial
#' autocorrelation for each spatial unit and then plotting a choropleth map with this
#' legend values: (Not Significant, p-value=0.05, p-value= 0.01, p-value=0.001,
#' p-value=0.0001, and Neighborless). Maps can represent concentrations of similar (cluster)
#' or dissimilar values (spatial outliers). Most significant clustered spatial units are
#' those with p-values smaller than 0.0001. Not significant clustered spatial units are
#' those with p-values greater than 0.05. LISA Cluster Map is done based on the
#' significance map but the choropleth legend is different (Not - Significant, High-High, Low-Low,
#' Low-High, High-Low, and Neighborless).
#'
#' @seealso
#' \itemize{
#'   \item Bivariate Moran's I: \code{\link{moran.bi}}
#'   \item Plot Bivariate Moran's I: \code{\link{moranbi.plot}}
#'   \item Bivariate Moran's I Test: \code{\link{moranbi.test}}
#'   \item Create object "nb": \code{\link[spdep]{poly2nb}}
#'   \item Create object "listw"/"nb": \code{\link[spdep]{nb2listw}}
#' }
#'
#' @section Links:
#'
#' \enumerate{
#'    \item \href{https://en.wikipedia.org/wiki/Indicators_of_spatial_association}{Indicators of Spatial Association}
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
#'library(spdep)
#'columbus <- st_read(system.file("extdata/columbus.shp", package="spdepbi")[1], quiet=TRUE)
#'col_nbq <- poly2nb(columbus)
#'a.lw <- nb2listw(col_nbq, style="W")
#'moran.cluster(columbus$CRIME, a.lw, zero.policy = FALSE, columbus, significant=TRUE)
moran.cluster <- function(x, listw, zero.policy = NULL, shp, significant=TRUE, alternative="two.sided", ...){
  mor.dat <- localmoran(x, listw, zero.policy=zero.policy, alternative=alternative)
  wx<-lag.listw(listw, x)
  lag.z<-scale(wx, center=T, scale=T)
  dat.z<-scale(x, center=T, scale=T)

  mor.dat1<-data.frame(mor.dat,lag.z,dat.z)
  names(mor.dat1)<-c("Ii","E.Ii","Var.Ii","Z.Ii","Pr.Zi","WY","Y")

  mor.dat1$cluster<-"UN"
  # both z scores are "high"
  mor.dat1$cluster[mor.dat1[,"Z.Ii"]>=qnorm(0.975,mean=0,sd=1)&mor.dat1[,"Y"]>0&mor.dat1[,"WY"]>0]<-"HH"
  # both z scores are "low"
  mor.dat1$cluster[mor.dat1[,"Z.Ii"]>=qnorm(0.975,mean=0,sd=1)&mor.dat1[,"Y"]<0&mor.dat1[,"WY"]<0]<-"LL"
  # one z score "high", the other "low"
  mor.dat1$cluster[mor.dat1[,"Z.Ii"]<=qnorm(0.025,mean=0,sd=1)&mor.dat1[,"Y"]>0&mor.dat1[,"WY"]<0]<-"HL"
  # one z score "low", the other "high"
  mor.dat1$cluster[mor.dat1[,"Z.Ii"]<=qnorm(0.025,mean=0,sd=1)&mor.dat1[,"Y"]<0&mor.dat1[,"WY"]>0]<-"LH"
  mor.dat1$cluster[is.na(mor.dat1[,"Z.Ii"])]<-"NA"

  cols<-c(brewer.pal(5, "RdBu"),"#BEBEBE")
  mor.dat1$col[mor.dat1$cluster=="UN"]<-cols[3]
  mor.dat1$col[mor.dat1$cluster=="HH"]<-cols[1]
  mor.dat1$col[mor.dat1$cluster=="LL"]<-cols[5]
  mor.dat1$col[mor.dat1$cluster=="HL"]<-cols[2]
  mor.dat1$col[mor.dat1$cluster=="LH"]<-cols[4]
  mor.dat1$col[mor.dat1$cluster=="NA"]<-cols[6]
  mor.dat1

  #dev.new()
  par(pty="s",mar=c(0,0,2,0), oma=c(0,0,0,0))
  P1 <- plot(st_geometry(shp), col=mor.dat1$col, main= deparse(substitute(x)), ...)
  #box(which = "plot", lty = "solid", col='red')
  #box(which = "figure", lty = "solid", col='green')
  #box(which = "inner", lty = "solid", col='blue')
  #box(which = "outer", lty = "solid", col='orange')
  #title("Line 0", line = 0)
  #title("Line 1", line = 1)
  #title("Line 2", line = 2)
  #title("Line 0", line = 0)
  #title("Line -1", line = -1)
  #title("Line -2", line = -2)

  #legend(locator(1), legend=c(paste("Not Significant  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="UN"]),")",sep="",collapse=""), paste("High-High  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="HH"]),")",sep="",collapse=""), paste("Low-Low  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="LL"]),")",sep="",collapse=""),
  #paste("Low-High  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="LH"]),")",sep="",collapse=""), paste("High-Low  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="HL"]),")",sep="",collapse=""), paste("Neighborless  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="NA"]),")",sep="",collapse="")), fill=c(cols[3],cols[1],cols[5],cols[4],cols[2],cols[6]), title = "LISA Cluster Map", bty="n", cex=1.2, y.intersp=0.8)

  legend("topleft",
         legend=c(paste("Not Significant  ",
                                    "(",
                                    length(mor.dat1$cluster[mor.dat1$cluster=="UN"]),
                                    ")",
                                    sep="",collapse=""),
                              paste("High-High  ",
                                    "(",
                                    length(mor.dat1$cluster[mor.dat1$cluster=="HH"]),
                                    ")",
                                    sep="",collapse=""),
                              paste("Low-Low  ",
                                    "(",
                                    length(mor.dat1$cluster[mor.dat1$cluster=="LL"]),
                                    ")",
                                    sep="",collapse=""),
                              paste("Low-High  ",
                                    "(",
                                    length(mor.dat1$cluster[mor.dat1$cluster=="LH"]),
                                    ")",
                                    sep="",collapse=""),
                              paste("High-Low  ",
                                    "(",
                                    length(mor.dat1$cluster[mor.dat1$cluster=="HL"]),
                                    ")",
                                    sep="",collapse=""),
                              paste("Neighborless  ",
                                    "(",
                                    length(mor.dat1$cluster[mor.dat1$cluster=="NA"]),
                                    ")",
                                    sep="",collapse="")),
          fill=c(cols[3],cols[1],cols[5],cols[4],cols[2],cols[6]),
          title = "LISA Cluster Map", bty="n", cex=1.2, y.intersp=0.8)

  #legend(0, 1.5, "cambiar")

  if (significant & alternative=="less" | significant & alternative=="greater"){
  mor.dat1$prob<-"UN"
  mor.dat1$prob[2*mor.dat1[,"Pr.Zi"]<=0.05&2*mor.dat1[,"Pr.Zi"]>0.01]<-"5%"
  mor.dat1$prob[2*mor.dat1[,"Pr.Zi"]<=0.01&2*mor.dat1[,"Pr.Zi"]>0.001]<-"1%"
  mor.dat1$prob[2*mor.dat1[,"Pr.Zi"]<=0.001&2*mor.dat1[,"Pr.Zi"]>0.0001]<-"0.1%"
  mor.dat1$prob[2*mor.dat1[,"Pr.Zi"]<=0.0001&2*mor.dat1[,"Pr.Zi"]>=0]<-"0.01%"
  mor.dat1$prob[is.na(mor.dat1[,"Z.Ii"])]<-"NA"


  colsp<-c(brewer.pal(5, "Greens"),"#BEBEBE")
  mor.dat1$col1[mor.dat1$prob=="UN"]<-colsp[1]
  mor.dat1$col1[mor.dat1$prob=="5%"]<-colsp[2]
  mor.dat1$col1[mor.dat1$prob=="1%"]<-colsp[3]
  mor.dat1$col1[mor.dat1$prob=="0.1%"]<-colsp[4]
  mor.dat1$col1[mor.dat1$prob=="0.01%"]<-colsp[5]
  mor.dat1$col1[is.na(mor.dat1[,"Z.Ii"])]<-colsp[6]
  mor.dat1

  #dev.new()
  par(pty="s",mar=c(0,0,2,0), oma=c(0,0,0,0))
  P2 <- plot(st_geometry(shp), col=mor.dat1$col1, main= deparse(substitute(x)), ...)
  if (significant & alternative=="greater") {
  #legend(locator(1), legend=c(paste("Not Significant  ","(",length(mor.dat1$prob[mor.dat1$prob=="UN"]),")",sep="",collapse=""), paste("p=0.05  ","(",length(mor.dat1$prob[mor.dat1$prob=="5%"]),")",sep="",collapse=""), paste("p=0.01 ","(",length(mor.dat1$prob[mor.dat1$prob=="1%"]),")",sep="",collapse=""),
  #paste("p=0.001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.1%"]),")",sep="",collapse=""), paste("p=0.0001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.01%"]),")",sep="",collapse=""), paste("Neighborless  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="NA"]),")",sep="",collapse="")), bty="n", fill=colsp, title = expression(paste("LISA Significance Map, ",'H'['a']:rho>0)), cex=1.2, y.intersp=0.8)
  legend("topleft",
         legend=c(paste("Not Significant  ","(",length(mor.dat1$prob[mor.dat1$prob=="UN"]),")",sep="",collapse=""),
                  paste("p=0.05  ","(",length(mor.dat1$prob[mor.dat1$prob=="5%"]),")",sep="",collapse=""),
                  paste("p=0.01 ","(",length(mor.dat1$prob[mor.dat1$prob=="1%"]),")",sep="",collapse=""),
                  paste("p=0.001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.1%"]),")",sep="",collapse=""),
                  paste("p=0.0001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.01%"]),")",sep="",collapse=""),
                  paste("Neighborless  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="NA"]),")",sep="",collapse="")),
         bty="n", fill=colsp,
         #title = expression(paste("LISA Significance Map, ",'H'['a']:rho>0)), cex=1.2, y.intersp=0.8)
         title = expression(paste("LISA Significance Map")), cex=1.2, y.intersp=0.8)
  }

  if (significant & alternative=="less"){
  #legend(locator(1), legend=c(paste("Not Significant  ","(",length(mor.dat1$prob[mor.dat1$prob=="UN"]),")",sep="",collapse=""), paste("p=0.05  ","(",length(mor.dat1$prob[mor.dat1$prob=="5%"]),")",sep="",collapse=""), paste("p=0.01 ","(",length(mor.dat1$prob[mor.dat1$prob=="1%"]),")",sep="",collapse=""),
  #paste("p=0.001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.1%"]),")",sep="",collapse=""), paste("p=0.0001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.01%"]),")",sep="",collapse=""), paste("Neighborless  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="NA"]),")",sep="",collapse="")), bty="n", fill=colsp, title = expression(paste("LISA Significance Map, ",'H'['a']:rho<0)), cex=1.2, y.intersp=0.8)

  legend("topleft",
         legend=c(paste("Not Significant  ","(",length(mor.dat1$prob[mor.dat1$prob=="UN"]),")",sep="",collapse=""),
                  paste("p=0.05  ","(",length(mor.dat1$prob[mor.dat1$prob=="5%"]),")",sep="",collapse=""),
                  paste("p=0.01 ","(",length(mor.dat1$prob[mor.dat1$prob=="1%"]),")",sep="",collapse=""),
                  paste("p=0.001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.1%"]),")",sep="",collapse=""),
                  paste("p=0.0001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.01%"]),")",sep="",collapse=""),
                  paste("Neighborless  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="NA"]),")",sep="",collapse="")),
        bty="n", fill=colsp,
        #title = expression(paste("LISA Significance Map, ",'H'['a']:rho<0)), cex=1.2, y.intersp=0.8)
        title = expression(paste("LISA Significance Map")), cex=1.2, y.intersp=0.8)
  }
  }

  if (significant & alternative=="two.sided"){
  mor.dat1$prob<-"UN"
  mor.dat1$prob[mor.dat1[,"Pr.Zi"]<=0.05&mor.dat1[,"Pr.Zi"]>0.01]<-"5%"
  mor.dat1$prob[mor.dat1[,"Pr.Zi"]<=0.01&mor.dat1[,"Pr.Zi"]>0.001]<-"1%"
  mor.dat1$prob[mor.dat1[,"Pr.Zi"]<=0.001&mor.dat1[,"Pr.Zi"]>0.0001]<-"0.1%"
  mor.dat1$prob[mor.dat1[,"Pr.Zi"]<=0.0001&mor.dat1[,"Pr.Zi"]>=0]<-"0.01%"
  mor.dat1$prob[is.na(mor.dat1[,"Z.Ii"])]<-"NA"


  colsp<-c(brewer.pal(5, "Greens"),"#BEBEBE")
  mor.dat1$col1[mor.dat1$prob=="UN"]<-colsp[1]
  mor.dat1$col1[mor.dat1$prob=="5%"]<-colsp[2]
  mor.dat1$col1[mor.dat1$prob=="1%"]<-colsp[3]
  mor.dat1$col1[mor.dat1$prob=="0.1%"]<-colsp[4]
  mor.dat1$col1[mor.dat1$prob=="0.01%"]<-colsp[5]
  mor.dat1$col1[is.na(mor.dat1[,"Z.Ii"])]<-colsp[6]
  mor.dat1

  #dev.new()
  par(pty="s",mar=c(0,0,2,0), oma=c(0,0,0,0))
  P2 <- plot(st_geometry(shp), col=mor.dat1$col1, main= deparse(substitute(x)), ...)
  #legend(locator(1), legend=c(paste("Not Significant  ","(",length(mor.dat1$prob[mor.dat1$prob=="UN"]),")",sep="",collapse=""), paste("p=0.05  ","(",length(mor.dat1$prob[mor.dat1$prob=="5%"]),")",sep="",collapse=""), paste("p=0.01 ","(",length(mor.dat1$prob[mor.dat1$prob=="1%"]),")",sep="",collapse=""),
  #paste("p=0.001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.1%"]),")",sep="",collapse=""), paste("p=0.0001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.01%"]),")",sep="",collapse=""), paste("Neighborless  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="NA"]),")",sep="",collapse="")), bty="n", fill=colsp, title = expression(atop("LISA Significance Map",'H'['a']:rho!=0)), cex=1.2, y.intersp=0.8)
  legend("topleft",
         legend=c(paste("Not Significant  ","(",length(mor.dat1$prob[mor.dat1$prob=="UN"]),")",sep="",collapse=""),
                  paste("p=0.05  ","(",length(mor.dat1$prob[mor.dat1$prob=="5%"]),")",sep="",collapse=""),
                  paste("p=0.01 ","(",length(mor.dat1$prob[mor.dat1$prob=="1%"]),")",sep="",collapse=""),
                  paste("p=0.001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.1%"]),")",sep="",collapse=""),
                  paste("p=0.0001  ","(",length(mor.dat1$prob[mor.dat1$prob=="0.01%"]),")",sep="",collapse=""),
                  paste("Neighborless  ","(",length(mor.dat1$cluster[mor.dat1$cluster=="NA"]),")",sep="",collapse="")),
         bty="n", fill=colsp,
         #title = expression(atop("LISA Significance Map",'H'['a']:rho!=0)), cex=1.2, y.intersp=0.8)
         title = expression(paste("LISA Significance Map")), cex=1.2, y.intersp=0.8)
  }
}
