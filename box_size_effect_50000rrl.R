# Effect of Box Size on Function Results
# Written to be tested locally, run on HPC for full analysis
library(rapt)
library(parallel)



# K Function --------------------------------------------------------------

# upload rcp
under <- read.rcp('~/Research/point_patterns/Final/FinalConfig1','~/Research/point_patterns/Final/system1',scaleUp = TRUE,newRadius = 0.5)
#under <- read.rcp('~/scratch/Rcode/FinalConfigs/FinalConfig1','~/scratch/Rcode/systems/system1',scaleUp = TRUE,newRadius = 0.5)

#under.rrl <- pK3est(0.06,under,2000,nrval=200,correction="trans",anom=TRUE)
#envPlot(under.rrl[[1]],ylim=c(-9,9))

#cluster.1 <- makecluster(under,over,0.5,0.5,type = "cr",cr = 2,pic=.75,toPlot=FALSE)
#cluster.1.k <- anomK3est(cluster.1[[1]],under.rrl,correction="trans")
#lines(cluster.1.k,col="black",lwd=2)

# 1. Create all of the different box sizes
t1 <- Sys.time()
cube.15x15 <- stitch.size(under,boxSize=c(15,15,15))
a <- pK3est(0.06,cube.15x15,50000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print('15x15x15 cube: ')
print(t2-t1)
write.csv(a[[1]],file="~/scratch/Rcode/cube15x15.csv")
write.csv(a[[2]],file="~/scratch/Rcode/cube15x15toSub.csv")
rm(cube.15x15, a)
gc()

t1 <- Sys.time()
cube.20x20 <- stitch.size(under,boxSize=c(20,20,20))
a <- pK3est(0.06,cube.20x20,50000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print('20x20x20 cube: ')
print(t2-t1)
write.csv(a[[1]],file="~/scratch/Rcode/cube20x20.csv")
write.csv(a[[2]],file="~/scratch/Rcode/cube20x20toSub.csv")
rm(cube.20x20, a)
gc()

t1 <- Sys.time()
cube.30x30 <- stitch.size(under,boxSize=c(30,30,30))
a <- pK3est(0.06,cube.30x30,50000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print('30x30x30 cube: ')
print(t2-t1)
write.csv(a[[1]],file="~/scratch/Rcode/cube30x30.csv")
write.csv(a[[2]],file="~/scratch/Rcode/cube30x30toSub.csv")
rm(cube.30x30, a)
gc()

t1 <- Sys.time()
cube.40x40 <- stitch.size(under,boxSize=c(40,40,40))
a <- pK3est(0.06,cube.40x40,50000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print('40x40x40 cube: ')
print(t2-t1)
write.csv(a[[1]],file="~/scratch/Rcode/cube40x40.csv")
write.csv(a[[2]],file="~/scratch/Rcode/cube40x40toSub.csv")
rm(cube.40x40, a)
gc()

t1 <- Sys.time()
cube.60x60 <- stitch.size(under,boxSize=c(60,60,60))
a <- pK3est(0.06,cube.60x60,50000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print('60x60x60 cube: ')
print(t2-t1)
write.csv(a[[1]],file="~/scratch/Rcode/cube60x60.csv")
write.csv(a[[2]],file="~/scratch/Rcode/cube60x60toSub.csv")
rm(cube.60x60, a)
gc()

t1 <- Sys.time()
box.60x60x10 <- stitch.size(under,boxSize=c(60,60,10))
a <- pK3est(0.06,box.60x60x10,50000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print('60x60x10 box: ')
print(t2-t1)
write.csv(a[[1]],file="~/scratch/Rcode/box.60x60x10.csv")
write.csv(a[[2]],file="~/scratch/Rcode/box.60x60x10toSub.csv")
rm(box.60x60x10, a)
gc()

t1 <- Sys.time()
box.60x60x15 <- stitch.size(under,boxSize=c(60,60,15))
a <- pK3est(0.06,box.60x60x15,50000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print('60x60x15 box: ')
print(t2-t1)
write.csv(a[[1]],file="~/scratch/Rcode/box.60x60x15.csv")
write.csv(a[[2]],file="~/scratch/Rcode/box.60x60x15toSub.csv")
rm(box.60x60x15, a)
gc()

t1 <- Sys.time()
box.60x60x20 <- stitch.size(under,boxSize=c(60,60,20))
a <- pK3est(0.06,box.60x60x20,50000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print('60x60x20 box: ')
print(t2-t1)
write.csv(a[[1]],file="~/scratch/Rcode/box.60x60x20.csv")
write.csv(a[[2]],file="~/scratch/Rcode/box.60x60x20toSub.csv")
rm(box.60x60x20, a)
gc()

t1 <- Sys.time()
box.60x60x30 <- stitch.size(under,boxSize=c(60,60,30))
a <- pK3est(0.06,box.60x60x30,50000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print('60x60x30 box: ')
print(t2-t1)
write.csv(a[[1]],file="~/scratch/Rcode/box.60x60x30.csv")
write.csv(a[[2]],file="~/scratch/Rcode/box.60x60x30toSub.csv")
rm(box.60x60x30, a)
gc()

t1 <- Sys.time()
box.60x60x40 <- stitch.size(under,boxSize=c(60,60,40))
a <- pK3est(0.06,box.60x60x40,50000,nrval=200,correction="trans",anom=TRUE)
t2 <- Sys.time()
print('60x60x40 box: ')
print(t2-t1)
write.csv(a[[1]],file="~/scratch/Rcode/box.60x60x40.csv")
write.csv(a[[2]],file="~/scratch/Rcode/box.60x60x40toSub.csv")
rm(box.60x60x40, a)
gc()



# G Function --------------------------------------------------------------

# upload background points  
rcp <- read.rcp('~/Research/point_patterns/Final/FinalConfig1','~/Research/point_patterns/Final/system1',scaleUp = TRUE,newRadius = 0.5)

#set test parameters
perc <- 0.1
k <- 10
nEvals <- 200
rmax <- 8
stepsize <- 0.05
anom <- TRUE

sizes <- matrix(c(20,20,20,
                  40,40,40,
                  60,60,60,
                  20,60,60,
                  40,60,60), nrow = 5, ncol = 3, byrow = TRUE)

names <- c("20x20x20",
           "40x40x40",
           "60x60x60",
           "20x60x60",
           "40x60x60")

# get RRL envelopes
for(i in 1:1){
  t1 <- Sys.time()
  samp <- stitch.size(rcp,boxSize=sizes[i,])
  a <- gbsa(perc = perc, pattern = samp, k = k, nEvals = nEvals, rmax = rmax, stepsize = stepsize, anom = anom)
  t2 <- Sys.time()
  print(paste(names[i], "vol: ", sep = " "))
  print(t2-t1)
  save(a, file = paste(names[i],".RData",sep = ""))
  rm(samp, a)
  gc()
}




gbsa <- function(perc, pattern, k, nEvals, rmax, stepsize, anom = FALSE, sorted = FALSE){
  # set up cluster export
  cl <- makePSOCKcluster(detectCores())
  clusterExport(cl,c("percentSelect"))
  clusterExport(cl, c("pattern","k","stepsize","rmax","perc"), envir = environment())
  clusterEvalQ(cl,library(spatstat))
  
  plist <- seq(1, nEvals)
  
  # manual G calculation - get out G for k nns
  hapl <- parLapply(cl, plist, function(x){
    X <- percentSelect(perc, pattern, x)
    nn <- nndist(X, k = 1:k)
    h <- apply(nn, 2, function(x){hist(x, breaks = seq(0,rmax, stepsize), plot = FALSE)})
    g <- sapply(h, function(x){cumsum((x$breaks[2]-x$breaks[1])*x$density)})
    return(g)
  })
  
  stopCluster(cl)
  
  nns <- list()
  for(i in 1:k){
    nns[[i]] <- matrix(0, nrow = nrow(hapl[[1]]), ncol = nEvals)
    for(j in 1:nEvals){
      nns[[i]][,j] <- hapl[[j]][,i]
    }
  }
  
  if(sorted == TRUE){
    nnss <- lapply(nns, function(x){t(apply(x,1,sort))})
    nns <- nnss
  }
  if(anom == TRUE){
    if(!exists("nnss")){
      nnss <- lapply(nns, function(x){t(apply(x,1,sort))})
    }
    mids <- lapply(nnss, function(x){x[,nEvals/2]})
    nnssanom <- lapply(nnss, function(x){mid <- x[,nEvals/2]
      apply(x, 2, function(y, mid){y - mid}, mid)})
    nns <- nnssanom
  }
  
  nns[[k+1]] = seq(0,rmax-stepsize,stepsize)
  
  names <- vector("character", k+1)
  for(i in 1:k){
    names[i] <- paste("nn",toString(i),sep = "")
  }
  names[k+1] <- "r"
  
  names(nns) <- names
  
  if(anom == TRUE){
    names(mids) <- names[1:k]
    return(list(nns,mids))
  }
  
  return(nns)
}

