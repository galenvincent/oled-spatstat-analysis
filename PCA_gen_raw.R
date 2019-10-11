# Script to output random data points for PCA
# This will only work on a linux machine

library(data.table)
library(zoo)
library(rapt)
library(parallel)

rm(list=ls())
gc()

cluster.size <- c(60,60,60)

# Number of random simulations desired
nrand <- 2 # Number of random parameters to generate
nclust <- 10 # Number of iterations to test each parameter set on
s <- 104 # random seed (change between runs)

#add in the random values
set.seed(s)
params <- list()
r <- runif(nrand, min = 2, max = 6.5)
den <- runif(nrand, min = 0.15, max = 1) 
rbp <- runif(nrand, min = 0, max = 0.5)
gbp <- runif(nrand, min = 0, max = 0.6)
for(i in 1:nrand){
  params[[i]] <- c(r[i], den[i], rbp[i], gbp[i])
}

#Set up the final data structure
sim.res <- list()
for(i in 1:nrand){
  sim.res[[i]] <- matrix(NA, nrow = nclust, ncol = 5)
}

save(params, file = 'params.RData')
print("Parameter set initialized...")
print(paste("Total parameters to test: ", toString(nrand), sep = ''))

# Upload cube RRL files
toSub <- fread('~/Research/K_cluster_series/cubetoSub_r35.csv', drop=1)
env.r <- fread('~/Research/K_cluster_series/cube_r35.csv', select=2)

# HPC
#toSub <- fread('cubetoSub_r35.csv', drop=1)
#env.r <- fread('cube_r35.csv', select=2)

# Max r value and number of r values to go to in the k tests
maxr <- max(env.r)
nr <- nrow(env.r)

cores2use <- detectCores() 

print(paste("Using ", toString(cores2use), " cores.", sep = ''))
print("starting simulations...")

cl <- makePSOCKcluster(cores2use)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))

t1 <- Sys.time()
outtemp <- parLapply(cl, 1:nclust, kseries2, nclust, params, maxr, nr, toSub, pcp = 0.051, verbose = FALSE)
t2 <- Sys.time()
print(t2-t1)

stopCluster(cl)

for(i in 1:nclust){
  for(j in 1:nrand){
    sim.res[[j]][i,] <- outtemp[[i]][j,] 
  }
}

save(sim.res, file = "sim.res.RData")
