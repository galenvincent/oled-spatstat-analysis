# Script to output random data points for PCA
# This will only work on a linux machine

library(data.table)
library(zoo)
library(rapt)
library(parallel)

rm(list=ls())
gc()

#source("~/Research/rapt/R/analysis_functions.R")

p <- 101

cluster.size <- c(60,60,60)

# Number of random simulations desired
nrand <- 0

# Set up a 5x5x5x5 grid of simulation conditions
grdsize <-2
r <- seq(2, 6.5, len = grdsize)
den <- seq(0.15, 1, len = grdsize)
rbp <- seq(0, 0.5, len = grdsize)
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

# Upload cube RRL files
toSub <- fread('~/Research/K_cluster_series/cubetoSub_r35.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube_r35.csv',select=2)

#toSub <- fread('~/Documents/Research/cubetoSub_big.csv',drop=1)
#env.r <- fread('~/Documents/Research/cube_big.csv',select=2)

#toSub <- fread('~/scratch/Rcode/cubetoSub_big.csv',drop=1)
#env.r <- fread('~/scratch/Rcode/cube_big.csv',select=2)

# Max r value and number of r values to go to in the k tests
maxr <- max(env.r)
nr <- nrow(env.r)

cores2use <- detectCores() 

print("starting")

cl <- makePSOCKcluster(cores2use) 
clusterEvalQ(cl, library(rapt))

fin <- 101

t1 <- Sys.time()
outtemp <- parLapply(cl,1:fin, kseries2, p ,tot, maxr, nr, toSub)
t2 <- Sys.time()
print(t2-t1)

stopCluster(cl)

for(i in 1:fin){
  fwrite(as.data.frame(outtemp[[i]]), paste("~/full_",toString(i),".csv",sep = ""))
  #fwrite(outtemp[[i]], paste("~/scratch/Rcode/pca_full/full_",toString(i),".csv"), sep = "")
}

outfin[,5:9] <- apply(simplify2array(outtemp), 1:2, mean, na.rm = TRUE)[,1:5]
outfin[,10:14] <- apply(simplify2array(outtemp), 1:2, sd, na.rm = TRUE)[,1:5]

outfin <- as.data.table(outfin)

fwrite(outfin, "~/pca.csv")
##fwrite(outfin, "~/scratch/Rcode/pca.csv")
