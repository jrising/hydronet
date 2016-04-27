setwd("~/projects/hydronet/snowtest/survey")

do.plots <- F

glaciers <- data.frame(id=c(1048, 3928, 3929, 3991, 1047, 1049, 1516, 2921, 3454, 3666, 3694, 732, 906, 3641),
                       lat=c(31.28, 31.4, 31.37, 30.45, 30.55, 30.73, 28.82, 32.2349, 30.885, 35.47, 28.45, 39.62, 27.7, 31.398),
                       lon=c(78.33, 78.5, 78.49, 81.333, 79.9, 79.68, 83.49, 77.5157, 78.818, 77.04, 85.75, 71.56, 86.57, 78.407),
                       year0=c(1981, 1974, 1976, 2005, 1984, 1981, 1998, 1986, 1992, 1986, 1991, 1967, 1978, 2000),
                       year1=c(1991, 1983, 1985, 2010, 1990, 1988, 2013, 2014, 2000, 1991, 2010, 2014, 1999, 2003),
                       col=c("#f35e5a", "#6ba204", "#00bfc4", "#b95eff", rep("#FFFF80", 10)))

basin <- as.matrix(read.csv("basinmask.csv"))
lats <- seq(29.74583, 33.9125, by=0.008333334)
lons <- seq(74.77084, 85.17084, by=0.008333334)
lons <- c(lons, 85.17084)
lats <- lats[-length(lats)]

if (do.plots) {
    library(maps)
    library(mapdata)

    stations <- read.csv("globalsod.csv")
    stations2 <- read.csv("ghcnv2.csv")

    image(lons, lats, t(basin), xlab="", ylab="", main="Glaciers and Stations in the Bhakra Basin", col=c("#808080", "#ffffff"))
    map("worldHires", xlim=c(78, 82), ylim=c(30, 32), add=T)
    points(glaciers$lon, glaciers$lat, col=as.character(glaciers$col), pch=19)
    points(stations$longitude, stations$latitude, pch=1)
    points(stations2$lon, stations2$lat, pch=18)
}

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

if (do.plots) {
    library(ggplot2)

    ggplot(data, aes(x=year, y=diff, colour=group)) +
        geom_line() + xlab("") + ylab("Relative mass balance (mm w.e.)") + scale_colour_discrete(name="FoG ID") + ggtitle("Relative mass balances of Bhakra basin glaciers")

    ggsave("balances.pdf", width=6, height=6)
}
