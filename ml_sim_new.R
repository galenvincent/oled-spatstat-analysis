# File for getting data for ML fitting:
library(data.table)
library(zoo)
library(rapt)
library(parallel)

rm(list=ls())
gc()

cluster.size <- c(60,60,60)

# Number of random simulations desired
nrand <- 8 # Number of random parameters to generate
nclust <- 4 # Number of iterations to test each parameter set on
s <- 122 # random seed (change between runs)

#add in the random values
set.seed(s)
params <- list()
cr <- runif(nrand, min = 2, max = 6.5)
rho1 <- runif(nrand, min = 0.2, max = 1)
rho2 <- runif(nrand, min = 0, max = 0.03)
rb <- runif(nrand, min = 0, max = 0.5)
pb <- runif(nrand, min = 0, max = 0.2)
for(i in 1:nrand){
  params[[i]] <- c(cr[i], rho1[i], rho2[i], rb[i], pb[i])
}

save(params, file = 'params.RData')
print("Parameter set initialized...")
print(paste("Total parameters to test: ", toString(nrand), sep = ''))

# Upload cube RRL files
#toSub <- fread('~/Research/K_cluster_series/cubetoSub_r35_6percent.csv', drop=1)
#env.r <- fread('~/Research/K_cluster_series/cube_r35.csv', select=2)

toSub <- fread('RCP_RRL_toSub_5.1percent.csv')

# HPC
#toSub <- fread('cubetoSub_r35.csv', drop=1)
#env.r <- fread('cube_r35.csv', select=2)

# Max r value and number of r values to go to in the k tests
maxr <- max(env.r)
nr <- nrow(env.r)
rm(env.r)

# Upload RCP files to use
under.nums <- seq(2,(nclust+1),1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,nclust,1)
rcp_path <- '~/Research/point_patterns/Final'
#rcp_path <- '~/scratch/Rcode/RCP' 

under.list <- list()
over.list <- list()
print('Uploading RCP files:')
for(i in 1:nclust){
  print(i)
  under <- read.rcp(paste(rcp_path, '/FinalConfig', toString(under.nums[i]), sep=''),
                    paste(rcp_path, '/system', toString(under.nums[i]), sep=''),
                    scaleUp = TRUE, newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, '/FinalConfig', toString(over.nums[i]), sep=''),
                   paste(rcp_path, '/system', toString(over.nums[i]), sep=''),
                   scaleUp = TRUE, newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  under.list[[i]] <- under.big
  over.list[[i]] <- over.big
}
rm(over, under, over.big, under.big)
gc()

#cores2use <- 28
cores2use <- 4 

print(paste("Using ", toString(cores2use), " cores.", sep = ''))

cl <- makePSOCKcluster(cores2use)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))
clusterExport(cl, c('over.list', 'under.list', 'params'))
print("starting simulations...")

t1 <- Sys.time()
sim.res <- parLapply(cl, 1:nrand, kmet_extract, nclust, maxr, nr, toSub, pcp = 0.05114235, tol = 0.005, 
                     verbose = TRUE, junk_path = '~/Research/junk/')
t2 <- Sys.time()
print(t2-t1)

stopCluster(cl)

save(sim.res, file = 'sim.res.RData')


kmet_extract <- function(i, nper, maxr, nr, toSub, pcp = 0.05114235, tol = 0.005,
                         verbose = FALSE, junk_path = '~/Research/junk/'){
  res <- matrix(NA, nrow = nper, ncol = 5)
  colnames(res) <- c('Km', 'Rm', 'Rdm', 'Rddm', 'Kdm')
  set.seed(i)
  seeds <- round(runif(nper, 1, 1e8))
  for(j in 1:nper){
    cluster <- clustersim(under.list[[j]], over.list[[j]], 0.5,
                          pcp = pcp,
                          cr = params[[i]][1],
                          rho1 = params[[i]][2],
                          rho2 = params[[i]][3],
                          rb = params[[i]][4],
                          pb = params[[i]][5],
                          tol = tol,
                          s = seeds[j])
    if(is.numeric(cluster)){
      res[j,] <- c(NA, NA, NA, NA, NA)
      next
    }
    result <- anomK3est(cluster[[1]], toSub, maxr, nr)
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    rvals.new <- rvals[15:length(rvals)]
    tvals.new <- tvals[15:length(rvals)]
    
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new, F)
    
    res[j,] <- c(metrics[[1]], metrics[[2]], metrics[[3]], metrics[[4]], metrics[[5]])
    
    rm(cluster, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  a <- as.data.frame(seeds)
  if(verbose == TRUE){
    write.csv(a, file = paste(junk_path, '/', toString(i), '.csv', sep = ''))
  }
  return(res)
}
