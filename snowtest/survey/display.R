setwd("~/projects/hydronet/snowtest/survey")

library(maps)
library(mapdata)

glaciers <- data.frame(id=c(1048, 3928, 3929, 3991), lat=c(31.28, 31.4, 31.37, 30.45), lon=c(78.33, 78.5, 78.49, 81.333), year0=c(1981, 1974, 1976, 2005), year1=c(1991, 1983, 1985, 2010), col=c("#f35e5a", "#6ba204", "#00bfc4", "#b95eff"))

basin <- as.matrix(read.csv("basinmask.csv"))
lats <- seq(29.74583, 33.9125, by=0.008333334)
lons <- seq(74.77084, 85.17084, by=0.008333334)
lons <- c(lons, 85.17084)
lats <- lats[-length(lats)]

stations <- read.csv("globalsod.csv")
stations2 <- read.csv("ghcnv2.csv")

image(lons, lats, t(basin), xlab="", ylab="", main="Glaciers and Stations in the Bhakra Basin", col=c("#808080", "#ffffff"))
map("worldHires", xlim=c(78, 82), ylim=c(30, 32), add=T)
points(glaciers$lon, glaciers$lat, col=as.character(glaciers$col), pch=19)
points(stations$longitude, stations$latitude, pch=1)
points(stations2$lon, stations2$lat, pch=18)

## Display the changes in mass balance
fog1048 <- read.table("FoG_MB_1048.csv", skip=13, sep=';', header=T)
fog3928 <- read.table("FoG_MB_3928.csv", skip=13, sep=';', header=T)
fog3929 <- read.table("FoG_MB_3929.csv", skip=13, sep=';', header=T)
fog3991 <- read.table("FoG_MB_3991.csv", skip=13, sep=';', header=T)
allfog <- list("1048"=fog1048, "3928"=fog3928, "3929"=fog3929, "3991"=fog3991)

data <- data.frame(year=c(), diff=c(), group=c())
for (name in c("1048", "3928", "3929", "3991")) {
    data <- rbind(data, data.frame(year=c(allfog[[name]]$REFERENCE_YEAR[1], allfog[[name]]$SURVEY_YEAR), diff=c(0, cumsum(allfog[[name]]$ANNUAL_BALANCE)), group=name))
}

library(ggplot2)

ggplot(data, aes(x=year, y=diff, colour=group)) +
    geom_line() + xlab("") + ylab("Relative mass balance (mm w.e.)") + scale_colour_discrete(name="FoG ID") + ggtitle("Relative mass balances of Bhakra basin glaciers")
