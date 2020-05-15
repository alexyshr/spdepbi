#' Test for Bivariate Moran's I
#'
#' @description
#' Allow to test the null hyphothesis (ho) of zero spatial autocorrelation,
#' comparing p-value with a significance levels of alpha (0.1, or 0.05, or 0.01).
#'
#'
#' @param X first variable
#' @param Y second variable
#' @param listw a neighbours list with spatial weights. From package spdep:
#' a listw object. Use poly2nb (class nb)
#' and nb2listw (class listw, nb) from package spdep. Can be any type of listw
#' object, for instance, rook contiguity (common edge) or queen contiguity (common
#' edge or common vertex)
#' @param zero.policy by default = NULL, if FALSE stop with error for any empty
#' neighbours sets, if TRUE permit the weights list to be formed with zero-length
#' weights vectors. Parameter inherited from the spdep package.
#' @param adjust.n by default = TRUE. if FALSE the number of observations is not
#' adjusted for noneighbour observations, if TRUE, the number of observations is
#' adjusted. Parameter inherited from the spdep package
#' @param N number of random elements for empirical density
#' @param graph by default = FALSE. Use TRUE to create test's graphic.
#' @param print.results by default = TRUE. Use FALSE to hide test results (table).
#' Results are: observed, expected and p-value.
#' @param ... other parameters similar to original
#'
#' @return a list with three values: observed (Bivariate Moran's Index), expected
#' (-1/(N-1), and p-value.
#'
#' @details
#'
#' Compare the observed Bivariate Moran's I (moran.bi function) with the
#' expected value empirical density. The expected value is -1/(N-1), where N is
#' the number of rows/samples (number of polygons), and represents the null hyphothesis
#' (ho) of No Spatial Autocorrelation (random, observed Bivariate Moran's I around zero).
#' This expected value density is constructed with Monte Carlo simulations.
#' Values significant below of -1/(N-1) represents negative spatial autocorrelation
#' (generally negative values of observed bivariate Moran's I), and values
#' significant above of -1/(N-1) represents positive spatial autocorrelation
#' (generally positive values of observed bivariate Moran's I).
#' For hypothesis testing the sample values are compared with empirical density,
#' and p-value is calculated. Significant values of p-value (reject Ho) are smaller
#' than 10\% (if alpha = 10\%), smaller than 5\% (if alpha = 5\%) or smaller than 1\% (
#' if alpha = 1\%). For significant values of p-value (reject Ho), the conclusion
#' of the test could be: "given the value of p-value, there is less than alpha (1\%,
#' or 5\%, or 10\%) likelihood that the pattern (clustered or dispersed) could be the
#' result of random change".
#'
#' @seealso
#' \itemize{
#'   \item Bivariate Moran's I: \code{\link{moran.bi}}
#'   \item Plot Bivariate Moran's I: \code{\link{moranbi.plot}}
#'   \item Create object "nb": \code{\link[spdep]{poly2nb}}
#'   \item Create object "listw"/"nb": \code{\link[spdep]{nb2listw}}
#' }
#'
#' @section Links:
#'
#' \enumerate{
#'    \item Spatial Autocorrelation \href{http://blogs.oregonstate.edu/geo599spatialstatistics/2016/06/08/spatial-autocorrelation-morans/}{(Moran’s I) Test}
#'    \item Moran's I \href{https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/}{Test}
#' }
#'
#' @import spatialreg
#' @import spdep
#' @import sp
#' @import grDevices
#' @import graphics
#' @import stats
#' @import ggplot2
#' @import RColorBrewer
#' @import combinat
#' @import sf
#' @export
#'
#' @examples
#' library(spdep)
#' columbus <- st_read(system.file("extdata/columbus.shp", package="spdepbi")[1], quiet=TRUE)
#' col_nbq <- poly2nb(columbus)
#' a.lw <- nb2listw(col_nbq, style="W")
#' set.seed(123)
#' MBCrime <- moranbi.test(columbus$CRIME,columbus$INC,a.lw,N=999,graph=TRUE,zero.policy =TRUE)
#' moranbi.test(columbus$INC,columbus$HOVAL,a.lw,999,graph=TRUE,zero.policy=TRUE,N=1000)
moranbi.test <-  function(X,Y,listw,zero.policy=NULL,adjust.n=TRUE,N,graph=FALSE,print.results=TRUE,...){
  observed<-moran.bi(X,Y,listw=listw,zero.policy=zero.policy,adjust.n = adjust.n,...)
  DF <- data.frame(1:length(X),X,Y)
  names(DF) <- c("Obs","X","Y")
  if(length(X)<8){
    X1<-unique(permn(DF$Obs))
  }
  else{
    X1<-randomize_vector(DF$Obs,N)
  }
  if(length(Y)<8){
    Y1<-unique(permn(DF$Obs))
  }
  else{
    Y1<-randomize_vector(DF$Obs,N)
  }
  store<-rep(NA,length(X1))
  for(i in 1:length(store)){
    store[i]<-moran.bi(X[X1[[i]]],Y[Y1[[i]]],listw,zero.policy = zero.policy, adjust.n = adjust.n)
  }
  if(observed>=0){
    p.value<-(sum(ifelse(store>observed,1,0))+1)/(length(store)+1)
    expected=(-1/(length(X)-1))
  }
  else if(observed<0){
    p.value<-(sum(ifelse(store<observed,1,0))+1)/(length(store)+1)
    expected=(-1/(length(X)-1))
  }
  if(graph==T){
    tmp.dat<-data.frame(store=store,observed=observed)
    export.graph<-ggplot(tmp.dat,aes(x=store))+
      scale_y_sqrt()+geom_density()+
      geom_vline(aes(xintercept=observed),color="red",size=1)+
      xlab("Bivariate Moran's I Coefficient") +
      ylab("Empirical Density")+theme_bw()
    print(export.graph)
  }
  if(print.results==T){
    print(list(Observed=observed,Expected=expected,p.value=p.value))
  }
  z<-list(Observed=observed,Expected=expected,p.value=p.value,Values=store)
}
