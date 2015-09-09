setwd("~/projects/hydronet/snowtest")

snowcover <- read.csv("snowcover.txt", header=F)
precip <- read.csv("precip.txt", header=F)
temps <- read.csv("temps.txt", header=F)

data <- snowcover
names(data) <- c("date0", "date1", "covered")
data$precip <- NA
data$meantemp <- NA
data$snows <- NA
data$melts <- NA
for (ii in 1:nrow(data)) {
    precipvalues <- precip[precip[,1] >= data$date0[ii] & precip[,1] <= data$date1[ii], 2]
    tempsvalues <- temps[temps[,1] >= data$date0[ii] & temps[,1] <= data$date1[ii], 2]
    data$precip[ii] <- sum(precipvalues)
    data$meantemp[ii] <- mean(tempsvalues)
    if (length(precipvalues) == length(tempsvalues))
        data$snows[ii] <- sum(precipvalues * (tempsvalues < 273.15))
    data$melts[ii] <- sum((tempsvalues - 273.15) * (tempsvalues > 273.15))
}

data <- data[!is.nan(data$meantemp),]
data$coverdelay <- c(NA, data$covered[-nrow(data)])

summary(lm(covered ~ meantemp + snows + coverdelay, data))

summary(lm(covered ~ precip + meantemp + precip*meantemp + coverdelay, data))

data$year0 <- floor(data$date0 / 10000)
data$month0 <- floor((data$date0 %% 10000) / 100)
data$day0 <- data$date0 %% 100
data$mydate0 <- as.Date(paste(data$year0, data$month0, data$day0, sep='-'))

library(ggplot2)

ggplot(data, aes(x=mydate0)) +
    geom_line(aes(y=snows), col=3) + geom_line(aes(y=melts), col=4) +
        geom_rect(data=data, aes(xmin=mydate0, xmax=mydate0+7, ymin=0, ymax=250, fill=covered), alpha=.1) + scale_y_continuous(limits=c(0, 250), expand = c(0,0)) +
            scale_x_date(limits=c(as.Date("1966-01-01"), as.Date("2012-03-01")), expand=c(0, 0)) + theme(legend.position="top")

ggplot(data, aes(x=mydate0)) +
    geom_line(aes(y=precip), col=1) + geom_line(aes(y=meantemp - 273.15), col=2) +
        geom_rect(data=data, aes(xmin=mydate0, xmax=mydate0+7, ymin=0, ymax=250, fill=covered), alpha=.1) + scale_y_continuous(limits=c(0, 250), expand = c(0,0)) +
                scale_x_date(limits=c(as.Date("1966-01-01"), as.Date("2012-03-01")), expand=c(0, 0)) + theme(legend.position="top")
