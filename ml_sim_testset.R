library(data.table)
library(zoo)
library(rapt)
library(parallel)

rm(list = ls())
gc()

#source('~/scratch/Rcode/kmet_extract_testset.R')

cluster.size <- c(60,60,60)

#nrand <- 25000
nrand <- 100
s <- 5235

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

nclust <- sample(1:509, nrand, replace = TRUE)
seeds <- sample((1+6543):(nrand+6543), nrand, replace = FALSE)

#save(params, file = 'test.params.RData')
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
clusterEvalQ(cl, library(data.table))
clusterExport(cl, c('seeds', 'nclust','params','kmet_extract_testset'))

t1 <- Sys.time()
out <- parLapply(cl, 1:nrand, function(x, maxr, nr, toSub){
#out <- lapply(1:nrand, function(x, maxr, nr, toSub){
  #print(x)
  res <- kmet_extract_testset(nclust[x], params[[x]], 
                              maxr = maxr, nr = nr, 
                              toSub = toSub, 
                              #rcp_path = '~/scratch/Rcode/RCP', 
                              s = seeds[x])
  if(x %% 1000 == 0){
    a <- data.frame(1)
    #fwrite(a, file = paste('~/scratch/Rcode/junk/',toString(x),'.csv',sep=''))
  }
  return(res)}, 
  maxr, nr, toSub)
t2 <- Sys.time()
print(t2-t1)

test.set <- matrix(NA, nrow = nrand, ncol = 10)
for(i in 1:nrand){
  test.set[i,1:5] <- params[[i]]
  test.set[i,6:10] <- out[[i]]
}

save(test.set, file = 'test.results.RData')




kmet_extract_testset <- function(j, params, maxr, nr, toSub, pcp = 0.05114235, tol = 0.005,
                                 rcp_path = '~/Research/point_patterns/Final',
                                 s = NULL){
  #upload
  under <- read.rcp(paste(rcp_path, '/FinalConfig', toString(j), sep=''),
                    paste(rcp_path, '/system', toString(j), sep=''),
                    scaleUp = TRUE, newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, '/FinalConfig', toString(j), sep=''),
                   paste(rcp_path, '/system', toString(j), sep=''),
                   scaleUp = TRUE, newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  cluster <- clustersim(under.big, over.big, 0.5,
                        pcp = pcp,
                        cr = params[1],
                        rho1 = params[2],
                        rho2 = params[3],
                        rb = params[4],
                        pb = params[5],
                        tol = tol,
                        s = s)
  if(is.numeric(cluster)){
    return(c(NA, NA, NA, NA, NA))
  }
  
  result <- anomK3est(cluster[[1]],toSub,maxr,nr)
  rvals <- result$r
  tvals <- result$trans
  
  # get out that peak info son
  rvals.new <- rvals[15:length(rvals)]
  tvals.new <- tvals[15:length(rvals)]
  
  #get those metrics out
  metrics <- k3metrics(rvals.new, tvals.new, FALSE)
  
  out <- c(metrics[[1]], metrics[[2]], metrics[[3]], metrics[[4]], metrics[[5]])
  
  rm(cluster, result, rvals, tvals, rvals.new, tvals.new, over, under, over.big, under.big)
  gc()
  
  return(out)
}
