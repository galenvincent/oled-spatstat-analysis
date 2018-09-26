# Effect of Box Size on Function Results
# Written to be tested locally, run on HPC for full analysis
library(rapt)
library(parallel)
library(data.table)

#### Some tests to start ####

# upload 2 rcps
under <- read.rcp('~/Research/point_patterns/Final/FinalConfig1','~/Research/point_patterns/Final/system1',scaleUp = TRUE,newRadius = 0.5)
#under <- read.rcp('~/scratch/Rcode/FinalConfigs/FinalConfig1','~/scratch/Rcode/systems/system1',scaleUp = TRUE,newRadius = 0.5)
#under.rrl <- pK3est(0.06,under,2000,nrval=200,correction="trans",anom=TRUE)
#envPlot(under.rrl[[1]],ylim=c(-9,9))

#cluster.1 <- makecluster(under,over,0.5,0.5,type = "cr",cr = 2,pic=.75,toPlot=FALSE)
#cluster.1.k <- anomK3est(cluster.1[[1]],under.rrl,correction="trans")
#lines(cluster.1.k,col="black",lwd=2)

#### create a 60x60x60 pattern of clusters, subsample it into different sized
# patterns, and run the functions on them (with correct subtractions) ####

# Doing 40 realizations of the same type of clustering, split the rcp patterns
# into 2 different sets. The under patterns and the over patterns
under.nums <- seq(1,40,1)
over.nums <- seq(41,80,1)

toSub <- fread('~/Research/box_size_effect/RRLtoSub1.csv',drop=1)
sample <- fread('~/Research/box_size_effect/RRL1.csv',select=2)
env <- fread('~/Research/box_size_effect/RRL1.csv',drop=1)


#loop
for(i in(1:40)){
i <- 1
#upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)

  under.big <- stitch(under,c(3,3,3))
  over.big <- stitch(over,c(3,3,3))
#make cluster
  cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",superfast = TRUE,cr=2)
  
#test on each of the different cluster sizes below here
  
#cut it down to size
  cluster.cut <- subSquare(cluster[[1]],c(60,60,60))
  
#test
  result <- anomK3est(cluster.cut,toSub,max(sample),200,correction="trans")

#output
  fwrite(result,'~/Research/box_size_effect/result1_1.csv',sep=',')

#clear memory
  rm(under,over,under.big,over.big,cluster,cluster,cluster.cut,result)
  gc()

#repeat
}


