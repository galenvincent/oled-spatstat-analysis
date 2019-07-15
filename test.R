library(data.table)
library(zoo)
library(rapt)
library(parallel)

rm(list = ls())
gc()

cluster.size <- c(60,60,60)

nrand <- 5
s <- 4767

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

nclust <- sample(1:523, nrand, replace = TRUE)
seeds <- sample((1+3452):(nrand+3452), nrand, replace = FALSE)

#save(params, file = 'ml.test.params.RData')
print('Parameter set initialized...')

# Upload cube RRL files
toSub <- fread('~/Research/K_cluster_series/cubetoSub_r35.csv', drop=1)
env.r <- fread('~/Research/K_cluster_series/cube_r35.csv', select=2)

# HPC

#toSub <- fread('~/scratch/Rcode/cubetoSub_r35.csv', drop = 1)
#env.r <- fread('~/scratch/Rcode/cube_r35.csv', select = 2)

maxr <- max(env.r)
nr <- nrow(env.r)

cores2use <- detectCores()

print(paste('Using', toString(cores2use), " cores.", sep = ''))
print('Starting...')

cl <- makePSOCKcluster(cores2use)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))
clusterExport(cl, c('seeds', 'nclust','params'))

t1 <- Sys.time()
out <- parLapply(cl, 1:nrand, function(x, maxr, nr, toSub){
  kseries(nclust[x], params[[x]],
          maxr = maxr, nr = nr,
          toSub = toSub,
          #rcp_path = '~/scratch/Rcode/RCP',
          s = seeds[x])},
  maxr, nr, toSub)
t2 <- Sys.time()
print(t2-t1)

test.set <- matrix(NA, nrow = nrand, ncol = 9)
for(i in 1:nrand){
  test.set[i,1:4] <- params[[i]]
  test.set[i,5:9] <- out[[i]]
}

save(test.set, file = 'ml.test.results.RData')

