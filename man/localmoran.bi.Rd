% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/localmoran.bi.R
\name{localmoran.bi}
\alias{localmoran.bi}
\title{Bivariate Local Moran's Index and Test. LISA - Local Indicators of Spatial Association}
\usage{
localmoran.bi(x, y, listw, zero.policy = NULL, na.action = na.fail,
  alternative = "greater", p.adjust.method = "none", mlvar = TRUE,
  spChk = NULL)
}
\arguments{
\item{x}{the first variable}

\item{y}{the second variable}

\item{listw}{a neighbours list with spatial weights. From package spdep:
a listw object. Use poly2nb (class nb)
and nb2listw (class listw, nb) from package spdep. Can be any type of listw
object, for instance, rook contiguity (common edge) or queen contiguity (common
edge or common vertex)}

\item{zero.policy}{by default = NULL, if FALSE stop with error for any empty
neighbours sets, if TRUE permit the weights list to be formed with zero-length
weights vectors. Parameter inherited from the spdep package. Similar meaning and
values than parameter zero.policy of \code{\link[spdep]{localmoran}}}

\item{na.action}{by default na.fail. In case of NA values in the variables
the function can have two options: "na.fail" and "na.pass". Similar meaning and
values than parameter zero.policy of \code{\link[spdep]{localmoran}}}

\item{alternative}{by default "greater". Type of alternative hypothesis test.
Other values are 'less' or 'two.sided'. Similar meaning and
values than parameter zero.policy of \code{\link[spdep]{localmoran}}}

\item{p.adjust.method}{correction method as defined in
\code{\link[spdep]{p.adjustSP}}. Include "bonferroni", "holm", "hochberg",
"hommel", "fdr". The default value "none" is a pass-through option.}

\item{mlvar}{by default TRUE. Similar meaning and
values than parameter zero.policy of \code{\link[spdep]{localmoran}}}

\item{spChk}{by default NULL. Similar meaning and
values than parameter zero.policy of \code{\link[spdep]{localmoran}}}
}
\value{
a matrix with five columns: 'Ii', 'E.Ii', 'var.Ii, 'Z.Ii', 'Pr (Z !=
0)'. This last column could be: 'Pr (Z != 0)' for alternative = 'two.sided',
'Pr (Z < 0)' for alternative = 'less', and 'Pr (Z > 0)' for alternative =
'greater'
}
\description{
Calculate the Local Moran's Index (LISA), and perform a significance test
for each spatial unit.
}
\details{
For hypothesis testing Moran's values are transformed to z-scores and its
probability is calculated based in parameter "alternative" (less, greater,
two.sided). Critical values for z-scores are 1.65 (confidence interval of 90\%,
alpha significance level of 10\%), 1.96 (confidence interval of 95\%, alpha
significance level of 5\%), and 2.58 (confidence interval of 99\%, alpha
significance level of 1\%)
}
\section{Links}{


\enumerate{
   \item \href{https://en.wikipedia.org/wiki/Indicators_of_spatial_association}{Indicators of Spatial Association}
   \item Spatial Autocorrelation \href{http://blogs.oregonstate.edu/geo599spatialstatistics/2016/06/08/spatial-autocorrelation-morans/}{(Moran’s I) Test}
   \item Moran's I \href{https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/}{Test}
   \item \href{http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Confidence_Intervals/BS704_Confidence_Intervals_print.html}{Confidence Intervals}
}
}

\examples{
library(spdep)
columbus <- st_read(system.file("extdata/columbus.shp", package="spdepbi")[1], quiet=TRUE)
col_nbq <- poly2nb(columbus)
a.lw <- nb2listw(col_nbq, style="W")
mor.dat <- localmoran.bi(columbus$CRIME, columbus$INC,
  a.lw, zero.policy=FALSE, alternative="two.sided")
}
\seealso{
\itemize{
  \item Bivariate Moran's I: \code{\link{moran.bi}}
  \item Plot Bivariate Moran's I: \code{\link{moranbi.plot}}
  \item Bivariate Moran's I Test: \code{\link{moranbi.test}}
  \item LISA Cluster and Significance Map: \code{\link{moran.cluster}}
  \item Create object "nb": \code{\link[spdep]{poly2nb}}
  \item Create object "listw"/"nb": \code{\link[spdep]{nb2listw}}
}
}
