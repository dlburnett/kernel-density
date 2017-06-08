
###########################################################################
library("maptools")
library("spatstat")
library("raster")
library("rgdal")
library("tidyr")
library ("tidyverse")
library("maps")
#setwd("C:/Users/dlburnett/Desktop/shape files from Raph")
setwd("~/Desktop/project files /")
#####################################################################################
#use this to filter the larger csv's down to which farms you want
#can filter out based on maturity/mortality as well 
#then load into QGIS to save in correct projection
#get points into spatial points
#from desktop
oneday <- read_csv("~/Desktop/project files /march_all_final.csv")

summary(oneday)
plot(oneday$mat)
plot(oneday$mor)
mature <- subset(oneday, mat >= 0.8)
mature2 <- subset(mature, mor <= 0.5)
write.csv(mature2, file = "allfarms_march_mature.csv") 

farm1 <- subset(oneday, farm == 1)
write.csv(farm1, file = "farm1_march.csv")

farm2 <- subset(oneday, farm == 2)
write.csv(farm2, file = "farm2_march.csv")

farm3 <- subset(oneday, farm == 3)
write.csv(farm3, file = "farm3_march.csv")

farm4 <- subset(oneday, farm == 4)
write.csv(farm4, file = "farm4_march.csv")

farm5 <- subset(oneday, farm == 5)
write.csv(farm5, file = "farm5_march.csv")

farm6 <- subset(oneday, farm == 6)
write.csv(farm6, file = "farm6_march.csv")

farm7 <- subset(oneday, farm == 7)
write.csv(farm7, file = "farm7_march.csv")

farm8 <- subset(oneday, farm == 8)
write.csv(farm8, file = "farm8_march.csv")

farm9 <- subset(oneday, farm == 9)
write.csv(farm9, file = "farm9_march.csv")

farm10 <- subset(oneday, farm == 10)
write.csv(farm10, file = "farm10_march.csv")

farm11 <- subset(oneday,farm == 11)
write.csv(farm11, file = "farm11_march.csv")

farm12 <- subset(oneday,farm == 12)
write.csv(farm12, file = "farm12_march.csv")

farm13 <- subset(oneday, farm == 13)
write.csv(farm13, file = "farm13_march.csv")

farm14 <- subset(oneday, farm == 14)
write.csv(farm14, file = "farm14_march.csv")

farm15 <- subset(oneday, farm == 15)
write.csv(farm15, file = "farm15_march.csv")

farm16 <- subset(oneday, farm == 16)
write.csv(farm16, file = "farm16_march.csv")

farm17 <- subset(oneday, farm == 17)
write.csv(farm17, file = "farm17_march.csv")

farm18 <- subset(oneday, farm == 18)
write.csv(farm18, file = "farm18_march.csv")

farm19 <- subset(oneday, farm == 19)
write.csv(farm19, file = "farm19_march.csv")

farm20 <- subset(oneday, farm == 20)
write.csv(farm20, file = "farm20_march.csv")

#############################################################################################
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
########################################################################################################

#START HERE
#Load sites and water contour

###########################all farms, sorted by maturity >= 0.8, and mor <= 0.5############################
setwd("~/Desktop/project files /sorted points in correct projection")
points <- readOGR("march_all_mature.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
plot(points)
map.axes()
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_allfarms_march_mature.tif","GTiff", overwrite=T)
plot(raster_out)
plot(points, add=T, pch=16, cex=0.5)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/allfarms_march_mature_densities.csv")

###############all farms for march###########################
setwd("~/Desktop/project files /sorted points in correct projection")
list.files()
points <- readOGR("march_all.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
plot(points)
map.axes()
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_allfarms_march_allparticles.tif","GTiff", overwrite=T)
plot(raster_out)
plot(points, add=T, pch=16, cex=0.5)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/allfarms_march_mature_densities.csv")

################farm1################################################################### 
setwd("~/Desktop/project files /sorted points in correct projection")
points <- readOGR("farm1_march_ALL.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
plot(points)
map.axes()
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm1_march.tif","GTiff", overwrite=T)


#bw selection using the diggle method (mean square cross validation)
out.diggle <- density(dat.ppp, sigma = bw.diggle, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf.diggle<- as(as(out.diggle, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf.diggle))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm1_march_digglebw.tif","GTiff", overwrite=T)
raster_in <- raster("Kdensity_farm1_march_digglebw.tif")

plot(raster_in)
plot(points, add=T, pch=16, cex=0.5)

b <- bw.diggle(dat.ppp)
b2 <- bw.ppl(dat.ppp)
farm_densities<-as.data.frame(raster::extract(raster_in,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv /farm1_densities.csv")


#bw selection using the ppl method
b2 <- bw.ppl(dat.ppp)

b2 <- 853.3

out.den.ppl <- density(dat.ppp, sigma = b2, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf.ppl <- as(as(out.den.ppl, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out_ppl <- brick(raster(out.spdf.ppl))
raster_out_ppl <- raster_out_ppl*1000000
writeRaster(raster_out_ppl, file = "~/Desktop/tif_files /Kdensity_farm1_march_pplbw.tif","GTiff", overwrite=T)
raster_in <- raster("~/Desktop/tif_files /Kdensity_farm1_march_pplbw.tif")

plot(raster_in)
plot(points, add=T, pch=16, cex=0.5)

farm_densities<-as.data.frame(raster::extract(raster_in,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv /farm2b_densities.csv")


#########################farm2################################################
setwd("~/Desktop/project files /sorted points in correct projection")
points <- readOGR("farm2_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
plot(points)
map.axes()
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm2_march.tif","GTiff", overwrite=T)


#bw selection using the diggle method (mean square cross validation)
out.diggle <- density(dat.ppp, sigma = bw.diggle, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf.diggle<- as(as(out.diggle, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf.diggle))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm2_march_digglebw.tif","GTiff", overwrite=T)
raster_in <- raster("Kdensity_farm2_march_digglebw.tif")

plot(raster_in)
plot(points, add=T, pch=16, cex=0.5)






#########################farm3################################################
setwd("~/Desktop/project files /sorted points in correct projection")
points <- readOGR("farm3_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
plot(points)
map.axes()
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm3_march.tif","GTiff", overwrite=T)


#bw selection using the diggle method (mean square cross validation)
out.diggle <- density(dat.ppp, sigma = bw.diggle, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf.diggle<- as(as(out.diggle, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf.diggle))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm3_march_digglebw.tif","GTiff", overwrite=T)
raster_in <- raster("Kdensity_farm3_march_digglebw.tif")

plot(raster_in)
plot(points, add=T, pch=16, cex=0.5)









#########################farm4################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm4_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
plot(points)
map.axes()
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm4_march.tif","GTiff", overwrite=T)
raster_out <- raster_out*1000000
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm4_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm4_march_densities.csv")

#bw selection using the diggle method (mean square cross validation)
out.diggle <- density(dat.ppp, sigma = bw.diggle, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf.diggle<- as(as(out.diggle, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf.diggle))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm4_march_digglebw.tif","GTiff", overwrite=T)
raster_in <- raster("~/Desktop/tif_files /Kdensity_farm4_march_digglebw.tif")

plot(raster_in)
plot(points, add=T, pch=16, cex=0.5)
farm_densities<-as.data.frame(raster::extract(raster_in,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm4_march_densities.csv")







#########################farm5################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm5_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm5_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm5_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm5_march_densities.csv")






#########################farm7################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm7_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm7_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm7_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm7_march_densities.csv")








#########################farm8################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm8_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm8_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm8_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm8_march_densities.csv")












#########################farm10################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm10_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm10_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm10_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm10_march_densities.csv")













#########################farm11################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm11_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm11_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm11_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm11_march_densities.csv")








#########################farm12################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm12_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm12_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm12_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm12_march_densities.csv")






#########################farm13################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm13_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm13_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm13_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm13_march_densities.csv")










#########################farm14################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm14_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm14_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm14_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm14_march_densities.csv")









#########################farm15################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm15_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm15_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm15_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm15_march_densities.csv")




#########################farm16################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm16_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm16_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm16_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm16_march_densities.csv")








#########################farm17################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm17_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm17_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm17_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm17_march_densities.csv")










#########################farm18################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm18_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm18_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm18_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm18_march_densities.csv")







#########################farm19################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm19_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm19_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm19_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm19_march_densities.csv")



#########################farm20################################################
setwd ("/Volumes/My Passport for Mac/sorted points in correct projection")
points <- readOGR("farm20_march.shp")
setwd("~/Desktop/project files /shape files from Raph")
poly <-readOGR("Broughton area_waterway_UTM9N.shp")
plot(poly)
map.axes()
points(points)
proj4string(points)
proj4string(poly)

#Get bandwidth information and determine cell size (based on desired dimensions)
#Note both are in EPSG:2953 - projected in metres
xrange_points = extent(points)[2]-extent(points)[1]
yrange_points = extent(points)[4]-extent(points)[3]
bandwidth <- min(xrange_points,yrange_points)/8

gridsize <- 100
xrange_poly <- extent(poly)[2]-extent(poly)[1]
yrange_poly <- extent(poly)[4]-extent(poly)[3]
cellsx <- xrange_poly/gridsize
cellsy <- yrange_poly/gridsize

#Generate kernel density surface (using Diggle's edge correction)
dimuse <- c(cellsy, cellsx)
dat.ow <- as(as(poly, "SpatialPolygons"), "owin")

dat.pp <- as(points, "SpatialPoints")
dat.ppp <- ppp(x = coordinates(dat.pp)[,1], y = coordinates(dat.pp)[,2], window = dat.ow)
out.den <- density(dat.ppp, sigma = bandwidth, dimyx = dimuse, edge=TRUE, diggle=TRUE)
out.spdf <- as(as(out.den, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
raster_out <- brick(raster(out.spdf))
#Set projection (optional - default is WGS84) (2953 = NB projection)
proj4string(raster_out) = CRS("+proj=utm +zone=9 +datum=WGS84 +units=m +no_defs")
#Convert from m2 to km2
raster_out <- raster_out*1000000
writeRaster(raster_out, file = "~/Desktop/tif_files /Kdensity_farm20_march.tif","GTiff", overwrite=T)
raster_out <- raster("~/Desktop/tif_files /Kdensity_farm20_march.tif")
plot(raster_out)

farm_densities<-as.data.frame(raster::extract(raster_out,points))
farm_densities$siteid<-points$SiteID

write.csv(farm_densities,file="~/Desktop/kernel density csv/farm20_march_densities.csv")
