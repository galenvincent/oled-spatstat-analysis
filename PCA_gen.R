# Script to output random data points for PCA
# This will only work on a linux machine

library(data.table)
library(zoo)
library(rapt)
library(parallel)

rm(list=ls())
gc()

source("~/Documents/rapt/R/analysis_functions.R")

under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)
p <- length(under.nums)

cluster.size <- c(60,60,60)

# Number of random simulations desired
nrand <- 10

# Set up a 5x5x5x5 grid of simulation conditions
grdsize <-2
r <- seq(2, 7, len = grdsize)
den <- seq(0.15, 1, len = grdsize)
rbp <- seq(0, 0.6, len = grdsize)
gbp <- seq(0, 0.6, len = grdsize)
tot <- list()
cnt <- 1
#make the grid
for(i in 1:grdsize){ 
  for(j in 1:grdsize){
    for(k in 1:grdsize){
      for(l in 1:grdsize){
        tot[[cnt]] <- c(r[i], den[j], rbp[k], gbp[l])
        cnt <- cnt + 1
      }
    }
  }
}

#add in the random values
if(nrand != 0){
  set.seed(10)
  rr <- runif(nrand, min = 2, max = 7)
  denr <- runif(nrand, min = 0.15, max = 1) 
  rbpr <- runif(nrand, min = 0, max = 0.6)
  gbpr <- runif(nrand, min = 0, max = 0.6)
  for(i in 1:nrand){
    tot[[cnt]] <- c(rr[i], denr[i], rbpr[i], gbpr[i])
    cnt <- cnt + 1
  }
}

#Set up the final matrix
outfin <- matrix(NA, nrow = length(tot), ncol = 14)
for(i in 1:length(tot)){
  outfin[i,1:4] <- tot[[i]]  
}

outtemp <- list()
outtemp2 <- list()

# Upload cube RRL files
toSub <- fread('~/Documents/Research/cubetoSub_big.csv',drop=1)
env.r <- fread('~/Documents/Research/cube_big.csv',select=2)
#toSub <- fread('~/scratch/Rcode/cubetoSub_big.csv',drop=1)
#env.r <- fread('~/scratch/Rcode/cube_big.csv',select=2)

# Max r value and number of r values to go to in the k tests
maxr <- max(env.r)
nr <- nrow(env.r)

cores2use <- detectCores() - 1 
# cl <- makePSOCKcluster(cores2use)
# clusterEvalQ(cl, library(rapt))
# cluserExport(cl, "kseries", "p", "maxr", "nr", "toSub")

p <- 2

print("starting")
t1 <- Sys.time()
for(i in 1:p){
  #upload
  under <- read.rcp(paste('~/Documents/Research/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Documents/Research/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Documents/Research/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Documents/Research/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  #clusterExport(cl, "under.big", "over.big")
  
  RNGkind("L'Ecuyer-CMRG")
  #out <- matrix(unlist(parLapply(cl,tot, kseries, maxr, nr, under.big, over.big, toSub)),nrow = length(tot), ncol = 5, byrow = TRUE)
  outtemp[[i]] <- matrix(unlist(mclapply(tot, kseries, maxr = maxr, nr = nr, under = under.big, over = over.big, toSub = toSub,
                  mc.set.seed = TRUE,
                  mc.cores = cores2use)),
                  nrow = length(tot), ncol = 6, byrow = TRUE)
  
  rm(over, under, over.big, under.big)
  print(toString(i))
}
t2 <- Sys.time()
print(t2-t1)

outfin[,5:9] <- apply(simplify2array(outtemp), 1:2, mean, na.rm = TRUE)[,1:5]
outfin[,10:14] <- apply(simplify2array(outtemp), 1:2, sd, na.rm = TRUE)[,1:5]

outfin <- as.data.table(outfin)

fwrite(outfin, "~/Documents/Research/pca.csv")
##fwrite(outfin, "~/scratch/Rcode/pca.csv")