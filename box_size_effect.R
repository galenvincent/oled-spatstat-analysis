# Effect of Box Size on Function Results
# Written to be tested locally, run on HPC for full analysis
library(rapt)
library(parallel)

#### Some tests to start ####

# upload 2 rcps
under <- read.rcp('~/scratch/Rcode/FinalConfigs/FinalConfig1','~/scratch/Rcode/systems/system1',scaleUp = TRUE,newRadius = 0.5)
#over <- read.rcp('~/Research/point_patterns/Final/FinalConfig2','~/Research/point_patterns/Final/system2',scaleUp = TRUE,newRadius = 0.5)

#under.rrl <- pK3est(0.06,under,2000,nrval=200,correction="trans",anom=TRUE)
#envPlot(under.rrl[[1]],ylim=c(-9,9))

#cluster.1 <- makecluster(under,over,0.5,0.5,type = "cr",cr = 2,pic=.75,toPlot=FALSE)
#cluster.1.k <- anomK3est(cluster.1[[1]],under.rrl,correction="trans")
#lines(cluster.1.k,col="black",lwd=2)

#### 50,000 RRL on all of the different box conditions ####
# 1. Create all of the different box sizes
cube.15x15 <- stitch.size(under,boxSize=c(15,15,15))
a <- pK3est(0.06,cube.15x15,50000,nrval=200,correction="trans",anom=TRUE)
write.csv(a[[1]],file="~/scratch/Rcode/cube15x15.csv")
write.csv(a[[2]],file="~/scratch/Rcode/cube15x15toSub.csv")
rm(cube.15x15, a)

cube.20x20 <- stitch.size(under,boxSize=c(20,20,20))
a <- pK3est(0.06,cube.20x20,50000,nrval=200,correction="trans",anom=TRUE)
write.csv(a[[1]],file="~/scratch/Rcode/cube20x20.csv")
write.csv(a[[2]],file="~/scratch/Rcode/cube20x20toSub.csv")
rm(cube.20x20, a)

cube.30x30 <- stitch.size(under,boxSize=c(30,30,30))
a <- pK3est(0.06,cube.30x30,50000,nrval=200,correction="trans",anom=TRUE)
write.csv(a[[1]],file="~/scratch/Rcode/cube30x30.csv")
write.csv(a[[2]],file="~/scratch/Rcode/cube30x30toSub.csv")
rm(cube.30x30, a)

cube.40x40 <- stitch.size(under,boxSize=c(40,40,40))
a <- pK3est(0.06,cube.40x40,50000,nrval=200,correction="trans",anom=TRUE)
write.csv(a[[1]],file="~/scratch/Rcode/cube40x40.csv")
write.csv(a[[2]],file="~/scratch/Rcode/cube40x40toSub.csv")
rm(cube.40x40, a)

cube.60x60 <- stitch.size(under,boxSize=c(60,60,60))
a <- pK3est(0.06,cube.60x60,50000,nrval=200,correction="trans",anom=TRUE)
write.csv(a[[1]],file="~/scratch/Rcode/cube60x60.csv")
write.csv(a[[2]],file="~/scratch/Rcode/cube60x60toSub.csv")
rm(cube.60x60, a)

box.60x60x10 <- stitch.size(under,boxSize=c(60,60,10))
a <- pK3est(0.06,box.60x60x10,50000,nrval=200,correction="trans",anom=TRUE)
write.csv(a[[1]],file="~/scratch/Rcode/box.60x60x10.csv")
write.csv(a[[2]],file="~/scratch/Rcode/box.60x60x10toSub.csv")
rm(box.60x60x10, a)

box.60x60x15 <- stitch.size(under,boxSize=c(60,60,15))
a <- pK3est(0.06,box.60x60x15,50000,nrval=200,correction="trans",anom=TRUE)
write.csv(a[[1]],file="~/scratch/Rcode/box.60x60x15.csv")
write.csv(a[[2]],file="~/scratch/Rcode/box.60x60x15toSub.csv")
rm(box.60x60x15, a)

box.60x60x20 <- stitch.size(under,boxSize=c(60,60,20))
a <- pK3est(0.06,box.60x60x20,50000,nrval=200,correction="trans",anom=TRUE)
write.csv(a[[1]],file="~/scratch/Rcode/box.60x60x20.csv")
write.csv(a[[2]],file="~/scratch/Rcode/box.60x60x20toSub.csv")
rm(box.60x60x20, a)

box.60x60x30 <- stitch.size(under,boxSize=c(60,60,30))
a <- pK3est(0.06,box.60x60x30,50000,nrval=200,correction="trans",anom=TRUE)
write.csv(a[[1]],file="~/scratch/Rcode/box.60x60x30.csv")
write.csv(a[[2]],file="~/scratch/Rcode/box.60x60x30toSub.csv")
rm(box.60x60x30, a)

box.60x60x40 <- stitch.size(under,boxSize=c(60,60,40))
a <- pK3est(0.06,box.60x60x40,50000,nrval=200,correction="trans",anom=TRUE)
write.csv(a[[1]],file="~/scratch/Rcode/box.60x60x40.csv")
write.csv(a[[2]],file="~/scratch/Rcode/box.60x60x40toSub.csv")
rm(box.60x60x40, a)
