setwd("~/projects/hydronet/bhakra/")

snows <- read.delim("snowcompare.tsv")

snows2 <- subset(snows, Row == 10 & Col == 12)
plot(density(snows2$True))

plot(snows2$True, snows2$Pred)

