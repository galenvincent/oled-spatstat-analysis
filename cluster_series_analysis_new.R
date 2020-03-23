# Moving on from the box size effect analysis, this script begins to analyze
# different returns of K-functions when cluster parameters are changed

library(data.table)
library(zoo)
library(rapt)
library(parallel)

# Detector Efficiency -----------------------------------------------------

# Detector efficiency series - select out [10, 20, 30, 40, 50, 60]% of points
# from clusters at random, see how the K-function changes
rm(list = ls())
gc()

nrcp <- 8

# detector efficiency series values
de <- list(0.9, 0.8, 0.7, 0.6, 0.5, 0.4)
n <- length(de)

set.seed(10)

toSub <- fread('RCP_RRL_toSub_5.1percent.csv')$`rcp_rrl[[2]]`
maxr <- 35
nr <- 400

cl <- makePSOCKcluster(8)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))
clusterExport(cl, c('de', 'n', 'toSub', 'maxr', 'nr'))

par.out <- parLapply(cl, 1:nrcp, function(i){
  rcp_path <- 'C:/Users/galen/Documents/Research/point_patterns/Final/'
  #rcp_path <- '~/scratch/Rcode/RCP/'
  under <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i), sep = ''),
                    paste(rcp_path, 'system', toString(i), sep = ''),
                    scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i+1), sep = ''),
                   paste(rcp_path, 'system', toString(i+1), sep = ''),
                   scaleUp = TRUE,newRadius = 0.5)
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  clust <- clustersim(under.big, over.big, 0.5,
                      pcp = 0.05114235,
                      cr = 3,
                      rho1 = 1,
                      rho2 = 0,
                      rb = 0,
                      pb = 0,
                      tol = 0.005,
                      s = i)
  if(is.numeric(clust)){
    while(is.numeric(clust)){
      clust <- clustersim(under.big, over.big, 0.5,
                          pcp = 0.05114235,
                          cr = 3,
                          rho1 = 1,
                          rho2 = 0,
                          rb = 0,
                          pb = 0,
                          tol = 0.005,
                          s = sample(1:10000, 1))
    }
  }
  
  met.out <- matrix(NA, nrow = n, ncol = 5)

  #test on different detector efficiencies below here
  for(j in 1:n){
    #browser()
    #remove the proper number of points
    cluster.de <- percentSelect(de[[j]], clust[[1]])
    #test
    result <- anomK3est(cluster.de, toSub, maxr, nr, correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    #plot(rvals,tvals,type = "n",xlab = "", ylab = "")
    rvals.new <- rvals[15:length(rvals)]
    tvals.new <- tvals[15:length(rvals)]
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new, FALSE)
    
    met.out[j,] <- unlist(metrics)
    
    rm(cluster.de, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  return(met.out)
  #clear memory
  rm(under,over,under.big,over.big, clust)
  gc()
  #repeat
})

stopCluster(cl)

Rm <- matrix(NA,nrcp,n)
Km <- matrix(NA,nrcp,n)
Rdm <- matrix(NA,nrcp,n)
Kdm <- matrix(NA,nrcp,n)
Rddm <- matrix(NA,nrcp,n)

for(i in 1:n){
  for(j in 1:nrcp){
    Km[j,i] <- par.out[[j]][i,1]
    Rm[j,i] <- par.out[[j]][i,2]
    Rdm[j,i] <- par.out[[j]][i,3]
    Rddm[j,i] <- par.out[[j]][i,4]
    Kdm[j,i] <- par.out[[j]][i,5]
  }
}

de.res <- list(Rm, Km, Rdm, Kdm, Rddm)
names(de.res) <- c('Rm', 'Km', 'Rdm', 'Kdm', 'Rddm')

save(de.res, file = 'de.res.RData')

# Cluster Radius ----------------------------------------------------------
rm(list=ls())
gc()
# see how the K-function changes

set.seed(11)

nrcp <- 8

# series values
r <- list(1, 2, 3, 4, 5, 8)
n <- length(r)

toSub <- fread('RCP_RRL_toSub_5.1percent.csv')$`rcp_rrl[[2]]`
maxr <- 35
nr <- 400

cl <- makePSOCKcluster(8)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))
clusterExport(cl, c('r', 'n', 'toSub', 'maxr', 'nr'))

par.out <- parLapply(cl, 1:nrcp, function(i){
  rcp_path <- 'C:/Users/galen/Documents/Research/point_patterns/Final/'
  #rcp_path <- '~/scratch/Rcode/RCP/'
  under <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i), sep = ''),
                    paste(rcp_path, 'system', toString(i), sep = ''),
                    scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i+1), sep = ''),
                   paste(rcp_path, 'system', toString(i+1), sep = ''),
                   scaleUp = TRUE,newRadius = 0.5)
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  met.out <- matrix(NA, nrow = n, ncol = 5)
  
  #test on different radius values below here
  for(j in 1:n){
    clust <- clustersim(under.big, over.big, 0.5,
                        pcp = 0.05114235,
                        cr = r[[j]],
                        rho1 = 1,
                        rho2 = 0,
                        rb = 0,
                        pb = 0,
                        tol = 0.005,
                        s = i)
    
    if(is.numeric(clust)){
      met.out[j,] <- c(NA, NA, NA, NA, NA)
      rm(clust, result, rvals, tvals, rvals.new, tvals.new)
      next
    }

    #test
    result <- anomK3est(clust[[1]], toSub, maxr, nr, correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    rvals.new <- rvals[15:length(rvals)]
    tvals.new <- tvals[15:length(rvals)]
    
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new, FALSE)
    
    met.out[j,] <- unlist(metrics)

    rm(clust, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  return(met.out)
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
})

stopCluster(cl)

Km <- matrix(NA,nrcp,n)
Rm <- matrix(NA,nrcp,n)
Rdm <- matrix(NA,nrcp,n)
Rddm <- matrix(NA,nrcp,n)
Kdm <- matrix(NA,nrcp,n)

for(i in 1:n){
  for(j in 1:nrcp){
    Km[j,i] <- par.out[[j]][i,1]
    Rm[j,i] <- par.out[[j]][i,2]
    Rdm[j,i] <- par.out[[j]][i,3]
    Rddm[j,i] <- par.out[[j]][i,4]
    Kdm[j,i] <- par.out[[j]][i,5]
  }
}

cr.res <- list(Rm, Km, Rdm, Kdm, Rddm)
names(cr.res) <- c('Rm', 'Km', 'Rdm', 'Kdm', 'Rddm')

save(cr.res, file = 'cr.res.RData')


# Rho1  ---------------------------------------------------------
rm(list=ls())
gc()

set.seed(12)

nrcp <- 8

rho1 <- list(1, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)
n <- length(rho1)

toSub <- fread('RCP_RRL_toSub_5.1percent.csv')$`rcp_rrl[[2]]`
maxr <- 35
nr <- 400

cl <- makePSOCKcluster(8)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))
clusterExport(cl, c('rho1', 'n', 'toSub', 'maxr', 'nr'))

par.out <- parLapply(cl, 1:nrcp, function(i){
  rcp_path <- 'C:/Users/galen/Documents/Research/point_patterns/Final/'
  #rcp_path <- '~/scratch/Rcode/RCP/'
  under <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i), sep = ''),
                    paste(rcp_path, 'system', toString(i), sep = ''),
                    scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i+1), sep = ''),
                   paste(rcp_path, 'system', toString(i+1), sep = ''),
                   scaleUp = TRUE,newRadius = 0.5)
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  met.out <- matrix(NA, nrow = n, ncol = 5)
  
  #test on different densities
  for(j in 1:n){
    clust <- clustersim(under.big, over.big, 0.5,
                        pcp = 0.05114235,
                        cr = 3,
                        rho1 = rho1[[j]],
                        rho2 = 0,
                        rb = 0,
                        pb = 0,
                        tol = 0.005,
                        s = i)
    if(is.numeric(clust)){
      met.out[j,] <- c(NA, NA, NA, NA, NA)
      rm(clust, result, rvals, tvals, rvals.new, tvals.new)
      next
    }
    
    #test
    result <- anomK3est(clust[[1]], toSub, maxr, nr, correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    rvals.new <- rvals[15:length(rvals)]
    tvals.new <- tvals[15:length(rvals)]
    
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new, FALSE)
    
    met.out[j,] <- unlist(metrics)
    
    rm(clust, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  return(met.out)
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
})

stopCluster(cl)

Rm <- matrix(NA,nrcp,n)
Km <- matrix(NA,nrcp,n)
Rdm <- matrix(NA,nrcp,n)
Kdm <- matrix(NA,nrcp,n)
Rddm <- matrix(NA,nrcp,n)

for(i in 1:n){
  for(j in 1:nrcp){
    Km[j,i] <- par.out[[j]][i,1]
    Rm[j,i] <- par.out[[j]][i,2]
    Rdm[j,i] <- par.out[[j]][i,3]
    Rddm[j,i] <- par.out[[j]][i,4]
    Kdm[j,i] <- par.out[[j]][i,5]
  }
}

rho1.res <- list(Rm, Km, Rdm, Kdm, Rddm)
names(rho1.res) <- c('Rm', 'Km', 'Rdm', 'Kdm', 'Rddm')

save(rho1.res, file = 'rho1.res.RData')

# Rho2  ---------------------------------------------------------
rm(list=ls())
gc()

set.seed(13)

nrcp <- 8

rho2 <- list(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035)
n <- length(rho2)

toSub <- fread('RCP_RRL_toSub_5.1percent.csv')$`rcp_rrl[[2]]`
maxr <- 35
nr <- 400

cl <- makePSOCKcluster(8)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))
clusterExport(cl, c('rho2', 'n', 'toSub', 'maxr', 'nr'))

par.out <- parLapply(cl, 1:nrcp, function(i){
  rcp_path <- 'C:/Users/galen/Documents/Research/point_patterns/Final/'
  #rcp_path <- '~/scratch/Rcode/RCP/'
  under <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i), sep = ''),
                    paste(rcp_path, 'system', toString(i), sep = ''),
                    scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i+1), sep = ''),
                   paste(rcp_path, 'system', toString(i+1), sep = ''),
                   scaleUp = TRUE,newRadius = 0.5)
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  met.out <- matrix(NA, nrow = n, ncol = 5)
  
  #test on different densities
  for(j in 1:n){
    clust <- clustersim(under.big, over.big, 0.5,
                        pcp = 0.05114235,
                        cr = 3,
                        rho1 = 1,
                        rho2 = rho2[[j]],
                        rb = 0,
                        pb = 0,
                        tol = 0.005,
                        s = i)
    if(is.numeric(clust)){
      met.out[j,] <- c(NA, NA, NA, NA, NA)
      rm(clust, result, rvals, tvals, rvals.new, tvals.new)
      next
    }
    
    #test
    result <- anomK3est(clust[[1]], toSub, maxr, nr, correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    rvals.new <- rvals[15:length(rvals)]
    tvals.new <- tvals[15:length(rvals)]
    
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new, FALSE)
    
    met.out[j,] <- unlist(metrics)
    
    rm(clust, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  return(met.out)
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
})

stopCluster(cl)

Rm <- matrix(NA,nrcp,n)
Km <- matrix(NA,nrcp,n)
Rdm <- matrix(NA,nrcp,n)
Kdm <- matrix(NA,nrcp,n)
Rddm <- matrix(NA,nrcp,n)

for(i in 1:n){
  for(j in 1:nrcp){
    Km[j,i] <- par.out[[j]][i,1]
    Rm[j,i] <- par.out[[j]][i,2]
    Rdm[j,i] <- par.out[[j]][i,3]
    Rddm[j,i] <- par.out[[j]][i,4]
    Kdm[j,i] <- par.out[[j]][i,5]
  }
}

rho2.res <- list(Rm, Km, Rdm, Kdm, Rddm)
names(rho2.res) <- c('Rm', 'Km', 'Rdm', 'Kdm', 'Rddm')

save(rho2.res, file = 'rho2.res.RData')

# Radius Blur -------------------------------------------------------------
rm(list=ls())
gc()

set.seed(14)

nrcp <- 7

# percent of cluster radius to set blur sds
rb <- list(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6)
n <- length(rb)

toSub <- fread('RCP_RRL_toSub_5.1percent.csv')$`rcp_rrl[[2]]`
maxr <- 35
nr <- 400

cl <- makePSOCKcluster(7)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))
clusterExport(cl, c('rb', 'n', 'toSub', 'maxr', 'nr'))

par.out <- parLapply(cl, 1:nrcp, function(i){
  rcp_path <- 'C:/Users/galen/Documents/Research/point_patterns/Final/'
  #rcp_path <- '~/scratch/Rcode/RCP/'
  under <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i), sep = ''),
                    paste(rcp_path, 'system', toString(i), sep = ''),
                    scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i+1), sep = ''),
                   paste(rcp_path, 'system', toString(i+1), sep = ''),
                   scaleUp = TRUE,newRadius = 0.5)
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  met.out <- matrix(NA, nrow = n, ncol = 5)
  
  for(j in 1:n){
    clust <- clustersim(under.big, over.big, 0.5,
                        pcp = 0.05114235,
                        cr = 3,
                        rho1 = 1,
                        rho2 = 0,
                        rb = rb[[j]],
                        pb = 0,
                        tol = 0.005,
                        s = i)
    if(is.numeric(clust)){
      met.out[j,] <- c(NA, NA, NA, NA, NA)
      rm(clust, result, rvals, tvals, rvals.new, tvals.new)
      next
    }
    
    #test
    result <- anomK3est(clust[[1]], toSub, maxr, nr, correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    rvals.new <- rvals[15:length(rvals)]
    tvals.new <- tvals[15:length(rvals)]
    
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new, FALSE)
    
    met.out[j,] <- unlist(metrics)
    
    rm(clust, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  return(met.out)
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
})

stopCluster(cl)

Rm <- matrix(NA,nrcp,n)
Km <- matrix(NA,nrcp,n)
Rdm <- matrix(NA,nrcp,n)
Kdm <- matrix(NA,nrcp,n)
Rddm <- matrix(NA,nrcp,n)

for(i in 1:n){
  for(j in 1:nrcp){
    Km[j,i] <- par.out[[j]][i,1]
    Rm[j,i] <- par.out[[j]][i,2]
    Rdm[j,i] <- par.out[[j]][i,3]
    Rddm[j,i] <- par.out[[j]][i,4]
    Kdm[j,i] <- par.out[[j]][i,5]
  }
}

rb.res <- list(Rm, Km, Rdm, Kdm, Rddm)
names(rb.res) <- c('Rm', 'Km', 'Rdm', 'Kdm', 'Rddm')

save(rb.res, file = 'rb.res.RData')


# Position Blur -----------------------------------------------------------
rm(list=ls())
gc()

set.seed(15)

nrcp <- 8

# percent of seperation distance for blur series
pb <- list(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
n <- length(pb)

toSub <- fread('RCP_RRL_toSub_5.1percent.csv')$`rcp_rrl[[2]]`
maxr <- 35
nr <- 400

cl <- makePSOCKcluster(8)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))
clusterExport(cl, c('pb', 'n', 'toSub', 'maxr', 'nr'))

par.out <- parLapply(cl, 1:nrcp, function(i){
  rcp_path <- 'C:/Users/galen/Documents/Research/point_patterns/Final/'
  #rcp_path <- '~/scratch/Rcode/RCP/'
  under <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i), sep = ''),
                    paste(rcp_path, 'system', toString(i), sep = ''),
                    scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i+1), sep = ''),
                   paste(rcp_path, 'system', toString(i+1), sep = ''),
                   scaleUp = TRUE,newRadius = 0.5)
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  met.out <- matrix(NA, nrow = n, ncol = 5)
  
  #test on different gaussian blurs
  for(j in 1:n){
    clust <- clustersim(under.big, over.big, 0.5,
                        pcp = 0.05114235,
                        cr = 3,
                        rho1 = 1,
                        rho2 = 0,
                        rb = 0,
                        pb = pb[[j]],
                        tol = 0.005,
                        s = i)
    if(is.numeric(clust)){
      met.out[j,] <- c(NA, NA, NA, NA, NA)
      rm(clust, result, rvals, tvals, rvals.new, tvals.new)
      next
    }
    
    #test
    result <- anomK3est(clust[[1]], toSub, maxr, nr, correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    rvals.new <- rvals[15:length(rvals)]
    tvals.new <- tvals[15:length(rvals)]
    
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new, FALSE)
    
    met.out[j,] <- unlist(metrics)
    
    rm(cluster, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  return(met.out)
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
})

stopCluster(cl)

Rm <- matrix(NA,nrcp,n)
Km <- matrix(NA,nrcp,n)
Rdm <- matrix(NA,nrcp,n)
Kdm <- matrix(NA,nrcp,n)
Rddm <- matrix(NA,nrcp,n)

for(i in 1:n){
  for(j in 1:nrcp){
    Km[j,i] <- par.out[[j]][i,1]
    Rm[j,i] <- par.out[[j]][i,2]
    Rdm[j,i] <- par.out[[j]][i,3]
    Rddm[j,i] <- par.out[[j]][i,4]
    Kdm[j,i] <- par.out[[j]][i,5]
  }
}

pb.res <- list(Rm, Km, Rdm, Kdm, Rddm)
names(pb.res) <- c('Rm', 'Km', 'Rdm', 'Kdm', 'Rddm')

save(pb.res, file = 'pb.res.RData')




