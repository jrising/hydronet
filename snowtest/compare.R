source("~/projects/hydronet/snowtest/survey/display.R")
setwd("~/projects/hydronet/snowtest")

modeled <- read.csv("glacierrecord.csv")
lats <- modeled$value[is.na(modeled$time) & modeled$variable == "latitude"]
lons <- modeled$value[is.na(modeled$time) & modeled$variable == "longitude"]

## Check that these are in the same order
sum(lats - glaciers$lat) == 0
sum(lons - glaciers$lon) == 0

modeled <- modeled[!is.na(modeled$time),]
modeled$date <- as.POSIXct(modeled$time, origin="1970-01-01")

radius <- 6371e3 # m
latspan <- .25 * 2*pi*radius / 180
lonspans <- .25 * 2*pi*radius*cos(lats * pi / 180) / 180
areas <- latspan * lonspans # m^2

for (ii in 1:nrow(glaciers)) {
    ## Translate into mm w.e.
    modeled$value[modeled$location == ii] <- 1000 * modeled$value[modeled$location == ii] / areas[ii]

    ## Recenter to year0
    if (ii < 4)
        volume1 <- modeled$value[modeled$location == ii & modeled$variable == 'volume' & substr(modeled$date, 1, 10) == paste0(glaciers$year0[ii], "-01-01")]
    else
        volume1 <- modeled$value[modeled$location == ii & modeled$variable == 'volume' & substr(modeled$date, 1, 10) == "2004-12-30"]
    modeled$value[modeled$location == ii & modeled$variable == 'volume'] <- modeled$value[modeled$location == ii & modeled$variable == 'volume'] - volume1
    
    ## Based on lats and lons
    modeled$location[modeled$location == ii] <- glaciers$id[ii]
}

modeled$location <- as.factor(modeled$location)

alldata <- data.frame(date=c(modeled$date[modeled$variable == "volume"], as.POSIXct(paste0(data$year, "-01-01"))),
                      diff=c(modeled$value[modeled$variable == "volume"], data$diff),
                      group=c(as.character(modeled$location[modeled$variable == "volume"]), as.character(data$group)),
                      source=c(rep("modeled", sum(modeled$variable == "volume")),
                          rep("observed", nrow(data))))

library(ggplot2)

ggplot(alldata, aes(x=date, y=diff, colour=group, linetype=source)) +
    geom_line() + xlab("") + ylab("Relative mass balance (mm w.e.)") + scale_colour_discrete(name="FoG ID") + ggtitle("Relative mass balances of Bhakra basin glaciers") + scale_y_continuous(limits=c(-1e4, 1e3), expand=c(0, 0)) + scale_linetype_discrete(name="Data source", breaks=c("observed", "modeled"), labels=c("Observed", "Modeled"))

ggsave("comparison.pdf", width=8, height=6)
