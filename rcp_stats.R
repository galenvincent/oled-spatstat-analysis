# file for analyzing RCP pattern to make model of K function for RCP clusters

rcp <- read.rcp('~/Research/point_patterns/Final/FinalConfig23','~/Research/point_patterns/Final/system23',scaleUp = TRUE,newRadius = 0.5)
rcp2 <- stitch.size(rcp, boxSize = c(60,60,60))

b <- nndist.pp3(rcp2,k=1:15)
a <- bdist.points3(rcp2)
c <- b[a<2,]

m <- apply(c,2,mean)