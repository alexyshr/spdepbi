############################################################
## Funcionamiento moran.bi, moranbi.test y moranbi.plot:  ##
############################################################

#library(spdep)
#library(sp)
#library(spatialreg)
library(sf)
library(spdep)
#sf::st_geometry
#sf::st_read
#example(columbus)
columbus <- st_read(system.file("extdata/columbus.shp", package="spdepbi")[1], quiet=TRUE)
#col.gal.nb <- read.gal(system.file("extdata/columbus.gal", package="spdepbi")[1])
#class(columbus)
#columbus = as(columbus, 'Spatial')
#class(columbus)
#coords <- coordinates(columbus)
coords = st_geometry(st_centroid(columbus$geometry))
#columbus = as(columbus, 'sf')
col_nbq <- poly2nb(columbus)
class(columbus)
par.lags1 <- nblag(col_nbq, 6)                  # Orden 6
plot(st_geometry(columbus), border="grey")
#myIds <-c("1","2","3","4","5")
#par.lags1.sub <- subset(par.lags1, par.lags1$region.id %in% myIds)
#coords.sub <- subset(coords, columbus$POLYID < 5)
#plot(par.lags1[[3]][1:3], coords[1:3,1:2], add=TRUE, col="red", lty=2)
#plot(par.lags1[[3]], coords, add=TRUE, col="red", lty=2)
plot(par.lags1[[2]], coords, add=TRUE, col="red", lty=2)
#zero.policy=T
#including regions with no neighbours.
#weights vectors of zero length are inserted for regions without neighbour in the neighbours list
#e.lw2 <- nb2listw(par.lags1[[6]], style="W",zero.policy=T)
#summary(e.lw2)
#style="W"
#row-standardized form.
#Row-standardization takes the given weights wij (e.g, the binary zero-one weights)
#and divides them by the row sum
#As a result, each row sum of the row-standardized weights equals one.
#Also, the sum of all weights equals n, the total number of observations
a.lw <- nb2listw(col_nbq, style="W")
summary(a.lw)
CRIME <- columbus$CRIME

####################################
########     moran.bi      #########
####################################

#source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moran.bi.R")
W <- as.matrix(as_dgRMatrix_listw(a.lw))
library(spdepbi)
moran.bi(columbus$CRIME,columbus$INC,a.lw,zero.policy =T)

########################################
########     moranbi.test      #########
########################################

#source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moranbi.test.R")
#source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/randomize_vector.R")
set.seed(123)
MBCrime <- moranbi.test(columbus$CRIME,columbus$INC,a.lw,999,graph=T,zero.policy =T,N=1000)
moranbi.test(columbus$INC,columbus$HOVAL,a.lw,999,graph=T,zero.policy =T,N=1000)

########################################
########     moranbi.plot      #########
########################################

#source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moranbi.plot.R")
# Editando las etiquetas de los ejes
CRIME <- as.vector(scale(columbus$CRIME))
INCOME <- as.vector(scale(columbus$INC))
moranbi.plot(CRIME,INCOME,quiet =F,zero.policy =F,listw=a.lw)
# Sin editar la etiqueta de los ejes
moranbi.plot(as.vector(scale(columbus$CRIME)),as.vector(scale(columbus$INC)),quiet =F,zero.policy =F,listw=a.lw)

#########################################
########     moran.cluster      #########
#########################################

#source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moran.cluster.R")
# LISA Cluster Map: COLUMBUS
# Si desean el mapa de significancia deben escribir T, de lo contrario F
#x11()
moran.cluster(CRIME, a.lw, zero.policy = FALSE, columbus, significant=T)
#moran.cluster(CRIME, e.lw2, zero.policy = FALSE, columbus, significant=T)

#########################################
########     getis.cluster      #########
#########################################

#source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/getis.cluster.R")
# Getis Cluster Map: COLUMBUS
# Si desean el mapa de significancia deben escribir T, de lo contrario F
#x11()
getis.cluster(CRIME, a.lw, zero.policy = FALSE, columbus, significant=T)
#getis.cluster(CRIME, e.lw2, zero.policy = FALSE, columbus, significant=T)


# ### Using package sp
# columbus <- st_read(system.file("extdata/columbus.shp", package="spdepbi")[1], quiet=TRUE)
# columbus = as(columbus, 'sf')
# plot(st_geometry(columbus))
# columbus = as(columbus, 'Spatial')
# #coords = st_geometry(st_centroid(columbus$geometry))
# coords <- coordinates(columbus)
# rn <- sapply(slot(columbus, "polygons"), function(x) slot(x, "ID"))
# #rn <- sapply(columbus, function(x) x$POLYID)
#
# columbus@data$ID <- rn
# #x <- locator(1)
# x = list(8.940741, 12.35734)
# names(x) = c("x","y")
# over(SpatialPoints(x),columbus)
# plot(columbus[columbus$ID==over(SpatialPoints(x),columbus)$ID,], col = "green4", add = TRUE)
# #plot(columbus[columbus$ID==over(SpatialPoints(x),columbus)$ID,], col = "white", add = TRUE)
#
# #Rook ploting
# col_nbr <- poly2nb(columbus,queen=F)            # rook
# col.lags10 <- nblag(col_nbr, 10)             # Orden 10
#
# #Plot first and second order polygons for rook
# plot(columbus[col.lags10[[1]][[as.numeric(over(SpatialPoints(x),columbus)$ID)]],], col = "springgreen4", density=20, add = TRUE)
# plot(columbus[col.lags10[[2]][[as.numeric(over(SpatialPoints(x),columbus)$ID)]],], col = "yellow4", density=20, add = TRUE)
#
# # Get centroid of selected polygon and plot its ID
# centroidsp <- coordinates(columbus[columbus$ID==over(SpatialPoints(x),columbus)$ID,])
# colnames(centroidsp) <- c("x","y")
# text(centroidsp[,1],centroidsp[,2],columbus[columbus$ID==over(SpatialPoints(x),columbus)$ID,]$POLYID)
#
# #Get centroid of first order - rook
# #centroidsp <- getSpPPolygonsLabptSlots(columbus[columbus$ID==over(SpatialPoints(x),columbus)$ID,])
# centroids <- coordinates(columbus[col.lags10[[1]][[as.numeric(over(SpatialPoints(x),columbus)$ID)]],])
# colnames(centroids) <- c("x","y")
# #
# #plot text of first order - rook
# text(centroids[,1],centroids[,2],columbus[col.lags10[[1]][[as.numeric(over(SpatialPoints(x),columbus)$ID)]],]$POLYID,cex=0.8)

# #Queen plotting
# columbus <- st_read(system.file("extdata/columbus.shp", package="spdepbi")[1], quiet=TRUE)
# columbus = as(columbus, 'sf')
# plot(st_geometry(columbus))
# columbus = as(columbus, 'Spatial')
# rn <- sapply(slot(columbus, "polygons"), function(x) slot(x, "ID"))
# columbus@data$ID <- rn
# x = list(8.940741, 12.35734)
# names(x) = c("x","y")
# over(SpatialPoints(x),columbus)
# plot(columbus[columbus$ID==over(SpatialPoints(x),columbus)$ID,], col = "green4", add = TRUE)
#
# col_nbq1 <- poly2nb(columbus)                   # queen
# col.lags <- nblag(col_nbq1, 10)
#
# #Plot first and second order polygons for queen
# plot(columbus[col.lags[[1]][[as.numeric(over(SpatialPoints(x),columbus)$ID)]],], col = "red", density=20, add = TRUE)
# plot(columbus[col.lags[[2]][[as.numeric(over(SpatialPoints(x),columbus)$ID)]],], col = "yellow4", density=20, add = TRUE)
#
# # Get centroid of selected polygon and plot its ID
# centroidsp <- coordinates(columbus[columbus$ID==over(SpatialPoints(x),columbus)$ID,])
# colnames(centroidsp) <- c("x","y")
# text(centroidsp[,1],centroidsp[,2],columbus[columbus$ID==over(SpatialPoints(x),columbus)$ID,]$POLYID)
#
# #Get centroid of first order - queen
# #centroidsp <- getSpPPolygonsLabptSlots(columbus[columbus$ID==over(SpatialPoints(x),columbus)$ID,])
# centroids <- coordinates(columbus[col.lags[[1]][[as.numeric(over(SpatialPoints(x),columbus)$ID)]],])
# colnames(centroids) <- c("x","y")
# #
# #plot text of first order - queen
# text(centroids[,1],centroids[,2],columbus[col.lags[[1]][[as.numeric(over(SpatialPoints(x),columbus)$ID)]],]$POLYID,cex=0.8)

# ### Using package sf
#Rook
columbus <- st_read(system.file("extdata/columbus.shp", package="spdepbi")[1], quiet=TRUE)
plot(st_geometry(columbus))
coords = st_geometry(st_centroid(columbus$geometry))
p1 = st_point(c(8.940741,12.35734))
(p1 = st_sfc(p1))
(polcontain = st_within(p1,columbus))
(myPol = which(as.matrix(polcontain)))
result = columbus[myPol, ]
#plot(result, col = "blue", add = TRUE)
plot(p1, add=T, pch= 16, cex= 0.5, col="red")
#Get centroid of selected polygon and plot its ID
centroidsp <- st_geometry(st_centroid(result$geometry))
text(st_coordinates(centroidsp)[,1],st_coordinates(centroidsp)[,2],columbus[myPol,]$POLYID)
# #Rook ploting
col_nbr <- poly2nb(columbus,queen=F)            # rook
col.lags10 <- nblag(col_nbr, 10)             # Orden 10
# #Plot first and second order polygons for rook
plot(st_geometry(columbus[col.lags10[[1]][[myPol]],]), col = "springgreen4", add = TRUE)
plot(st_geometry(columbus[col.lags10[[2]][[myPol]],]), col = "yellow4", add = TRUE)
# #Get centroid of first order - rook
centroids <- st_geometry(st_centroid(columbus[col.lags10[[1]][[myPol]],]$geometry))
#plot text of first order - rook
text(st_coordinates(centroids)[,1],st_coordinates(centroids)[,2],columbus[col.lags10[[1]][[myPol]],]$POLYID ,cex=0.8)
#
#Queen
#
plot(st_geometry(columbus))
coords = st_geometry(st_centroid(columbus$geometry))
p1 = st_point(c(8.940741,12.35734))
(p1 = st_sfc(p1))
(polcontain = st_within(p1,columbus))
(myPol = which(as.matrix(polcontain)))
result = columbus[myPol, ]
#plot(result, col = "blue", add = TRUE)
plot(p1, add=T, pch= 16, cex= 0.5, col="red")
#Get centroid of selected polygon and plot its ID
centroidsp <- st_geometry(st_centroid(result$geometry))
text(st_coordinates(centroidsp)[,1],st_coordinates(centroidsp)[,2],columbus[myPol,]$POLYID)
#Queen neighbours
col_nbq1 <- poly2nb(columbus)                   # queen
col.lags10 <- nblag(col_nbq1, 10)
# #Plot first and second order polygons for rook
plot(st_geometry(columbus[col.lags10[[1]][[myPol]],]), col = "springgreen4", add = TRUE)
plot(st_geometry(columbus[col.lags10[[2]][[myPol]],]), col = "yellow4", add = TRUE)
# #Get centroid of first order - rook
centroids <- st_geometry(st_centroid(columbus[col.lags10[[1]][[myPol]],]$geometry))
#plot text of first order - rook
text(st_coordinates(centroids)[,1],st_coordinates(centroids)[,2],columbus[col.lags10[[1]][[myPol]],]$POLYID ,cex=0.8)



