
###########################################################################
library("maptools")
library("spatstat")
library("raster")
library("rgdal")
library("tidyr")
library ("tidyverse")
library("maps")
#setwd("C:/Users/dlburnett/Desktop/shape files from Raph")
setwd("~/Desktop/project files /shape files from Raph")

#get points into spatial points
#from desktop
oneday <-read_csv("~/Desktop/files on the data set/oneday.csv")
farm1 <- subset(oneday, farm == 1)
write.csv(farm1, file = "farm1_oneday.csv")
locs_f1 <- subset(farm1, select = c("lat2", "lon2"))


#define projection, using the custom projection that Jeff wrote
crs.geo <- CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs ") 

#and get the lat/lon into a spatialPointsDataFrame
coordinates(locs_f1) <- c("lon2", "lat2") 
proj4string(locs_f1) <- crs.geo
summary(locs_f1)

#turn into UTM 9N projection
locs_latlon <- spTransform(locs_f1, CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs"))
proj4string(locs_latlon)
plot(locs_latlon)

#Load sites and water contour
#points <-readShapePoints("")
#just renaming so I don't have to change the code
points <- readOGR("correct projection.shp")
poly <-readOGR("Waterway-BC_UTM9N.shp")
plot(poly)
map.axes()
points(points)
plot(points)
map.axes()
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points=extent(points)[2]-extent(points)[1]
yrange_points=extent(points)[4]-extent(points)[3]
bandwidth <-min(xrange_points,yrange_points)/8

gridsize<-500
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx<- xrange_poly/gridsize
cellsy<- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf = as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+init=epsg:2953")
#Convert from m2 to km2
raster_out<-raster_out*1000000

writeRaster(raster_out,"Kdensity_ALL_sites_NB.tif","GTiff", overwrite=T)

###########################################################################
raster_in<-raster("Kdensity_ALL_sites_NB.tif")

plot(raster_in)
plot(points, add=T, pch=16, cex=0.5)

farm_densities<-as.data.frame(extract(raster_in,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="NB_site_densities.csv")