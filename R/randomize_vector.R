#' Generate a random sample from other vector
#'
#' @description
#' Use function 'sample' from base R library
#' to generate a random sample.
#' Function taken from Edzer Pebesma package.
#'
#' @param X vector to choose from
#' @param N number of random elements to select from X
#'
#' @return a list, a vector
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
#' example(columbus)
#' #col_nbq <- poly2nb(columbus)
#' #a.lw <- nb2listw(col_nbq, style="W")
#' #set.seed(123)
#' DF <- data.frame(1:length(columbus$CRIME),columbus$CRIME,columbus$INC)
#' X1<-randomize_vector(DF$Obs,999)
randomize_vector <-  function(X,N){
   lst<-list()
   for(i in 1:N){
     lst[[i]]<-sample(X,length(X))
   }
   return(lst)
}
