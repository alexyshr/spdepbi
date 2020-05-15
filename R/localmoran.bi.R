#' Bivariate Local Moran's Index and Test. LISA - Local Indicators of Spatial Association
#'
#' @description
#' Calculate the Local Moran's Index (LISA), and perform a significance test
#' for each spatial unit.
#'
#' @param x the first variable
#' @param y the second variable
#' @param listw a neighbours list with spatial weights. From package spdep:
#' a listw object. Use poly2nb (class nb)
#' and nb2listw (class listw, nb) from package spdep. Can be any type of listw
#' object, for instance, rook contiguity (common edge) or queen contiguity (common
#' edge or common vertex)
#' @param zero.policy by default = NULL, if FALSE stop with error for any empty
#' neighbours sets, if TRUE permit the weights list to be formed with zero-length
#' weights vectors. Parameter inherited from the spdep package. Similar meaning and
#' values than parameter zero.policy of \code{\link[spdep]{localmoran}}
#' @param na.action by default na.fail. In case of NA values in the variables
#' the function can have two options: "na.fail" and "na.pass". Similar meaning and
#' values than parameter zero.policy of \code{\link[spdep]{localmoran}}
#' @param alternative by default "greater". Type of alternative hypothesis test.
#' Other values are 'less' or 'two.sided'. Similar meaning and
#' values than parameter zero.policy of \code{\link[spdep]{localmoran}}
#' @param p.adjust.method correction method as defined in
#' \code{\link[spdep]{p.adjustSP}}. Include "bonferroni", "holm", "hochberg",
#' "hommel", "fdr". The default value "none" is a pass-through option.
#' @param mlvar by default TRUE. Similar meaning and
#' values than parameter zero.policy of \code{\link[spdep]{localmoran}}
#' @param spChk by default NULL. Similar meaning and
#' values than parameter zero.policy of \code{\link[spdep]{localmoran}}
#'
#' @return a matrix with five columns: 'Ii', 'E.Ii', 'var.Ii, 'Z.Ii', 'Pr (Z !=
#' 0)'. This last column could be: 'Pr (Z != 0)' for alternative = 'two.sided',
#' 'Pr (Z < 0)' for alternative = 'less', and 'Pr (Z > 0)' for alternative =
#' 'greater'
#'
#' @details
#' For hypothesis testing Moran's values are transformed to z-scores and its
#' probability is calculated based in parameter "alternative" (less, greater,
#' two.sided). Critical values for z-scores are 1.65 (confidence interval of 90\%,
#' alpha significance level of 10\%), 1.96 (confidence interval of 95\%, alpha
#' significance level of 5\%), and 2.58 (confidence interval of 99\%, alpha
#' significance level of 1\%)
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
#' mor.dat <- localmoran.bi(columbus$CRIME, columbus$INC,
#'   a.lw, zero.policy=FALSE, alternative="two.sided")
localmoran.bi <- function (x, y, listw, zero.policy = NULL, na.action = na.fail,
    alternative = "greater", p.adjust.method = "none", mlvar = TRUE,
    spChk = NULL)
{
    stopifnot(is.vector(x))
    if (!inherits(listw, "listw"))
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (is.null(zero.policy))
        zero.policy <- get("zeroPolicy", envir = 'spdep' %:::% '.spdepOptions')
    stopifnot(is.logical(zero.policy))
    if (!is.null(attr(listw$neighbours, "self.included")) &&
        attr(listw$neighbours, "self.included"))
        stop("Self included among neighbours")
    if (is.null(spChk))
        spChk <- get.spChkOption()
    if (spChk && !chkIDs(x, listw))
        stop("Check of data and weights ID integrity failed")
    if (!is.numeric(x))
        stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (!is.numeric(y))
        stop(paste(deparse(substitute(y)), "is not a numeric vector"))
    NAOK <- deparse(substitute(na.action)) == "na.pass"
    x <- na.action(x)
    na.act <- attr(x, "na.action")
    y <- na.action(y)
    na.act1 <- attr(y, "na.action")
    rn <- attr(listw, "region.id")
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
        excl <- class(na.act) == "exclude"
    }
    n <- length(listw$neighbours)
    if (n != length(x))
        stop("Different numbers of observations")
    res <- matrix(nrow = n, ncol = 5)
    if (alternative == "two.sided")
        Prname <- "Pr(z != 0)"
    else if (alternative == "greater")
        Prname <- "Pr(z > 0)"
    else Prname <- "Pr(z < 0)"
    colnames(res) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", Prname)
    z <- as.vector(scale(x))
    yz <- as.vector(scale(y))
    lyz <- lag.listw(listw, yz, zero.policy = zero.policy, NAOK = NAOK)
    if (mlvar)
        s2 <- sum(z^2, na.rm = NAOK)/n
    else s2 <- sum(z^2, na.rm = NAOK)/(n - 1)
    res[, 1] <- (z/s2) * lyz
    Wi <- sapply(listw$weights, sum)
    res[, 2] <- -Wi/(n - 1)
    if (!mlvar)
        s2 <- sum(z^2, na.rm = NAOK)/n
    b2 <- (sum(z^4, na.rm = NAOK)/n)/(s2^2)
    A <- (n - b2)/(n - 1)
    B <- (2 * b2 - n)/((n - 1) * (n - 2))
    C <- Wi^2/((n - 1)^2)
    Wi2 <- sapply(listw$weights, function(x) sum(x^2))
    Wikh2 <- sapply(listw$weights, function(x) {
        ifelse(is.null(x), 0, 1 - crossprod(x, x))
    })
    res[, 3] <- A * Wi2 + B * Wikh2 - C
    res[, 4] <- (res[, 1] - res[, 2])/sqrt(res[, 3])
    if (alternative == "two.sided")
        pv <- 2 * pnorm(abs(res[, 4]), lower.tail = FALSE)
    else if (alternative == "greater")
        pv <- pnorm(res[, 4], lower.tail = FALSE)
    else pv <- pnorm(res[, 4])
    res[, 5] <- p.adjustSP(pv, listw$neighbours, method = p.adjust.method)
    if (!is.null(na.act) && excl) {
        res <- naresid(na.act, res)
    }
    if (!is.null(rn))
        rownames(res) <- rn
    attr(res, "call") <- match.call()
    if (!is.null(na.act))
        attr(res, "na.action") <- na.act
    class(res) <- c("localmoran.bi", class(res))
    res
}
