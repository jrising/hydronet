setwd("~/projects/hydronet/snowtest/survey")

library(maps)
library(mapdata)

glaciers <- data.frame(id=c(1048, 3928, 3929, 3991), lat=c(31.28, 31.4, 31.37, 30.45), lon=c(78.33, 78.5, 78.49, 81.333), year0=c(1981, 1974, 1976, 2005), year1=c(1991, 1983, 1985, 2010))

basin <- as.matrix(read.csv("basinmask.csv"))
lats <- seq(29.74583, 33.9125, by=0.008333334)
lons <- seq(74.77084, 85.17084, by=0.008333334)
lons <- c(lons, 85.17084)
lats <- lats[-length(lats)]

stations <- read.csv("globalsod.csv")
stations2 <- read.csv("ghcnv2.csv")

image(lons, lats, t(basin), xlab="", ylab="", main="Glaciers and Stations in the Bhakra Basin")
map("worldHires", xlim=c(78, 82), ylim=c(30, 32), add=T)
points(glaciers$lon, glaciers$lat)
points(stations$longitude, stations$latitude, pch=19)
points(stations2$lon, stations2$lat, pch=18, col=3)
