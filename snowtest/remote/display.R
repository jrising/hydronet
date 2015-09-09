setwd("~/projects/hydronet/snowtest")

map <- read.csv("output.csv", header=F)
map <- as.matrix(map)
image(map)

image(map[361:720, 361:720])

image(map[400:420, 580:600])
