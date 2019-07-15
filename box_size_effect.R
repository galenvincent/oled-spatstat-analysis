# Effect of Box Size on Function Results
# Written to be tested locally, run on HPC for full analysis
library(rapt)
library(parallel)
library(data.table)
library(RColorBrewer)
library(zoo)

# K function --------------------------------------------------------------

# Doing n realizations of the same type of clustering, split the rcp patterns
# into 2 different sets. The under patterns and the over patterns

#### 100% CLUSTER DENSITY ####
n <- 523 # number of RCP patterns you have

under.nums <- seq(2,(n+1),1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,n,1)

cluster.sizes <- list(c(15,15,15),
                      c(20,20,20),
                      c(30,30,30),
                      c(40,40,40),
                      c(60,60,60),
                      c(60,60,10),
                      c(60,60,15),
                      c(60,60,20),
                      c(60,60,30),
                      c(60,60,40))

toSub <- vector("list",10)
env.r <- vector("list",10)
for(i in(1:10)){
  toSub[[i]] <- fread(paste('~/Research/box_size_effect/RRLtoSub',toString(i),'.csv',sep=""),drop=1)
  env.r[[i]] <- fread(paste('~/Research/box_size_effect/RRL',toString(i),'.csv',sep=""),select=2)
}

p <- length(under.nums)

p <- 4

cl <- makePSOCKcluster(detectCores())
clusterExport(cl, c("under.nums","over.nums","cluster.sizes", "toSub", "env.r"), envir = environment())
clusterEvalQ(cl,library(rapt))

#loop
parLapply(cl, 1:p, function(i){
  if(exists('cluster.sizes')){
    print('we good')
  }else{
    print('no good')
  }
  
#upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)

  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
#make cluster
  cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast", cr=3)
  
#test on each of the different cluster sizes below here
  
  for(j in(1:10)){
    #cut it down to size
    cluster.cut <- subSquare(cluster[[1]],cluster.sizes[[j]])
    
    #test
    result <- anomK3est(cluster.cut, toSub[[j]], max(env.r[[j]]), 200, correction="trans")
    
    #output
    #write(result,paste('~/Research/box_size_effect/density_1/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    write(1, paste('~/Research/oled-spatstat-analysis/',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    print(paste(toString(i),"_",toString(j)))
    rm(cluster.cut)
    gc()
  }
#clear memory
  rm(under,over,under.big,over.big,cluster,result)
  gc()
#repeat
})
stopCluster(cl)

#### 50% CLUSTER DENSITY ####
# Doing 101 realizations of the same type of clustering, split the rcp patterns
# into 2 different sets. The under patterns and the over patterns
under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

cluster.sizes <- list(c(15,15,15),
                      c(20,20,20),
                      c(30,30,30),
                      c(40,40,40),
                      c(60,60,60),
                      c(60,60,10),
                      c(60,60,15),
                      c(60,60,20),
                      c(60,60,30),
                      c(60,60,40))

toSub <- vector("list",10)
env.r <- vector("list",10)
for(i in(1:10)){
  toSub[[i]] <- fread(paste('~/Research/box_size_effect/RRLtoSub',toString(i),'.csv',sep=""),drop=1)
  env.r[[i]] <- fread(paste('~/Research/box_size_effect/RRL',toString(i),'.csv',sep=""),select=2)
}

#loop
for(i in(1:length(under.nums))){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch(under,c(3,3,3))
  over.big <- stitch(over,c(3,3,3))
  #make cluster
  cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=2,den=0.5)
  
  #test on each of the different cluster sizes below here
  
  for(j in(1:10)){
    #cut it down to size
    cluster.cut <- subSquare(cluster[[1]],cluster.sizes[[j]])
    
    #test
    result <- anomK3est(cluster.cut,toSub[[j]],max(env.r[[j]]),200,correction="trans")
    
    #output
    fwrite(result,paste('~/Research/box_size_effect/density_0.5/result',toString(j),'_',toString(i),'.csv',sep=""),sep=',')
    rm(cluster.cut)
    gc()
  }
  #clear memory
  rm(under,over,under.big,over.big,cluster,result)
  gc()
  #repeat
}


#### 50% CLUSTER DENSITY - Metric tracking ####
rm(list = ls())
gc()
n <- 8

under.nums <- seq(2,n+1,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,n,1)

cluster.sizes <- list(c(15,15,15),
                      c(20,20,20),
                      c(30,30,30),
                      c(40,40,40),
                      c(50,50,50),
                      c(60,60,60))

toSub <- vector("list",length(cluster.sizes))
env.r <- vector("list",length(cluster.sizes))
for(i in 1:length(cluster.sizes)){
  if(i == 5){
    toSub[[i]] <- 'placeholder'
    env.r[[i]] <- 'placeholder'
    
    #toSub[[i]] <- fread('Z:/Galen/Box\ Size\ RRLs/RRLtoSub11.csv', drop=1)
    #env.r[[i]] <- fread('Z:/Galen/Box\ Size\ RRLs/RRL11.csv', select=2)
    
    #HPC
    #toSub[[i]] <- fread('~/scratch/Rcode/RRLs/RRLtoSub11.csv', drop=1)
    #env.r[[i]] <- fread('~/scratch/Rcode/RRLs/RRL11.csv', select=2)
    
    next
  }
  toSub[[i]] <- fread(paste('Z:/Galen/Box\ Size\ RRLs/RRLtoSub',toString(i),'.csv',sep=''), drop=1)
  env.r[[i]] <- fread(paste('Z:/Galen/Box\ Size\ RRLs/RRL',toString(i),'.csv',sep=''), select=2)
  
  #HPC
  #toSub[[i]] <- fread(paste('~/scratch/Rcode/RRLs/RRLtoSub',toString(i),'.csv',sep=''), drop=1)
  #env.r[[i]] <- fread(paste('~/scratch/Rcode/RRLs/RRL',toString(i),'.csv',sep=''), select=2)
}

cl <- makePSOCKcluster(detectCores())
clusterEvalQ(cl, c(library(rapt)))
clusterExport(cl, c('cluster.sizes','toSub','env.r','under.nums','over.nums'))


#loop
res <- parLapply(cl, 1:n, function(i){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  #HPC
  #under <- read.rcp(paste('~/scratch/Rcode/RCP/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/RCP/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  #over <- read.rcp(paste('~/scratch/Rcode/RCP/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/RCP/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  #make cluster
  cluster <- makecluster(under.big, over.big, 0.5, 0.5, type="cr", speed="superfast", cr=3, den=0.5)
  
  rvals <- matrix(NA, nrow = 200, ncol = length(cluster.sizes))
  kres <- matrix(NA, nrow = 200, ncol = length(cluster.sizes))
  metricres <- matrix(NA, nrow = 5, ncol = length(cluster.sizes))
  
  #test on each of the different cluster sizes below here
  for(j in(1:length(cluster.sizes))){
    if(j == 5){
      next
    }
    
    #cut it down to size
    cluster.cut <- subSquare(cluster[[1]], cluster.sizes[[j]])
    
    #test
    result <- anomK3est(cluster.cut, toSub[[j]], max(env.r[[j]]), 200, correction="trans")
    rvals[,j] <- result$r
    kres[,j] <- result$trans
    
    rvals.new <- result$r[10:length(result$r)]
    tvals.new <- result$trans[10:length(result$trans)]
    
    metricres[,j] <- unlist(k3metrics(rvals.new, tvals.new, toplot = FALSE))
    
    rm(cluster.cut, result, rvals.new, tvals.new)
    gc()
    #print(paste(toString(i),'_',toString(j),sep = ''))
  }
  
  #clear memory
  rm(under,over,under.big,over.big,cluster)
  gc()
  
  return(list(rvals = rvals, kres = kres, mets = metricres))
})
stopCluster(cl)

kres <- list()
metricres <- list()
for(i in 1:length(cluster.sizes)){
  kres[[i]] <- matrix(NA, nrow = 200, ncol = n+1)
  metricres[[i]] <- matrix(NA, nrow = 5, ncol = n)
  
  kres[[i]][,1] <- res[[1]]$rvals[,i]
  for(j in 1:n){
    kres[[i]][,j+1] <- res[[j]]$kres[,i]
    metricres[[i]][,j] <- res[[j]]$mets[,i]
  }
} 

save('kres', file = 'kres_cluster_RRL.RData')
save('metricres', file= 'metricres_cluster_RRL.RData')


# G Function --------------------------------------------------------------
# Updload RRL data from HPC
cube.data <- list()
load("C:/Users/galen/Documents/Research/box_size_effect_G/20x20x20.RData")
cube.data[[1]] <- a
load("C:/Users/galen/Documents/Research/box_size_effect_G/40x40x40.RData")
cube.data[[2]] <- a
load("C:/Users/galen/Documents/Research/box_size_effect_G/60x60x60.RData")
cube.data[[3]] <- a
load("C:/Users/galen/Documents/Research/box_size_effect_G/20x60x60.RData")
cube.data[[4]] <- a
load("C:/Users/galen/Documents/Research/box_size_effect_G/40x60x60.RData")
cube.data[[5]] <- a
names(cube.data) <- c("20x20x20",
                      "40x40x40",
                      "60x60x60",
                      "20x60x60",
                      "40x60x60")
rm(a)

# Run RRL on 101 clustered data sets
under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

cluster.sizes <- list(c(20,20,20),
                      c(40,40,40),
                      c(60,60,60),
                      c(60,60,20),
                      c(60,60,40))

toSub <- vector("list",5)
env.r <- vector("list",5)

res <- list()
for(i in 1:5){
  res[[i]] <- list()
}

#set G parameters
k <- 10
rmax <- 8
stepsize <- 0.05

for(i in(1:5)){
  toSub[[i]] <- list()
  for(j in 1:k){
    toSub[[i]][[j]] <- cube.data[[i]][[2]][[j]]
  }
  env.r[[i]] <- cube.data[[i]][[1]]$r
}

#loop
for(i in(1:length(under.nums))){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  #make cluster
  cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=2, den = 0.5, pcp = 0.1)
  
  #test on each of the different cluster sizes below here
  
  for(j in(1:5)){
    #cut it down to size
    cluster.cut <- subSquare(cluster[[1]],cluster.sizes[[j]])
    
    #test
    obs <- nndist(cluster.cut, k = 1:k)
    hobs <- apply(obs, 2, function(x){hist(x, breaks = seq(0,rmax, stepsize), plot = FALSE)})
    gobs <- sapply(hobs, function(x){cumsum((x$breaks[2]-x$breaks[1])*x$density)})
    
    obsanom <- matrix(NA, ncol = ncol(gobs), nrow = nrow(gobs))
    for(m in 1:k){
      obsanom[,m] <- gobs[,m]- toSub[[j]][[m]]
    }
    
    res[[j]][[i]] <- obsanom
    
    print(paste(toString(i), "_", toString(j), sep = ""))
  }
  #clear memory
  rm(under,over,under.big,over.big,cluster, cluster.cut, obs, hobs, gobs, obsanom)
  gc()
  #repeat
}

res.sorted <- list()
for(i in 1:5){
  res.sorted[[i]] <- list()
  for(j in 1:k){
    res.sorted[[i]][[j]] <- matrix(NA, ncol = length(over.nums), nrow = rmax/stepsize)
    for(m in 1:length(over.nums)){
      res.sorted[[i]][[j]][,m] <- res[[i]][[m]][,j]
    }
  }
}

names(res.sorted) <- c("20x20x20",
                       "40x40x40",
                       "60x60x60",
                       "20x60x60",
                       "40x60x60")

for(i in 1:5){
  names(res.sorted[[i]]) <- c("nn1",
                              "nn2",
                              "nn3",
                              "nn4",
                              "nn5",
                              "nn6",
                              "nn7",
                              "nn8",
                              "nn9",
                              "nn10")  
}

# plot these things
for(l in 1:10){
  percentile_env <- .999
  percentile_clust <- .95
  nn <- l
  
  env.plot.big <- matrix(NaN,160,5)
  env.plot.small <- matrix(NaN,160,5)
  clust.plot.big <- matrix(NaN,160,5)
  clust.plot.small <- matrix(NaN,160,5)
  rvals <- matrix(NaN,160,5)
  
  for(i in (1:5)){
    rtemp <- cube.data[[i]][[1]]$r # r values
    tvals <- cube.data[[i]][[1]][[nn]] # results 
    nTests <- ncol(tvals) # number of tests done
    prange <- percentile_env*nTests # get the range of indices for which each percentile spans
    
    sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
    
    ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
    ind.small <- round(nTests/2)-floor(prange/2) # do the same for the low end
    
    env.plot.big[,i] <- sortedtVals[,ind.big]
    env.plot.small[,i] <- sortedtVals[,ind.small]
    rvals[,i] <- as.numeric(rtemp)
    rm(rtemp,tvals,nTests,prange,sortedtVals,ind.big,ind.small)
    ##
    tvals <- res.sorted[[i]][[nn]]
    nTests <- ncol(tvals) # number of tests done
    prange <- percentile_clust*nTests # get the range of indices for which each percentile spans
    
    sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
    
    ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
    ind.small <- round(nTests/2)-floor(prange/2) # do the same for the low end
    
    clust.plot.big[,i] <- sortedtVals[,ind.big]
    clust.plot.small[,i] <- sortedtVals[,ind.small]
    
    rm(tvals,nTests,prange,sortedtVals,ind.big,ind.small)
    gc()
  }
  
  
  color = c("red", "blue", "aquamarine4", "lightpink", "skyblue")
  xlim = c(0,max(rvals))
  ylim = c(-0.1,0.375)
  titles <- c(paste("20x20x20"," - nn", toString(l), sep = ""),
              paste("40x40x40"," - nn", toString(l), sep = ""),
              paste("60x60x60"," - nn", toString(l), sep = ""),
              paste("20x60x60"," - nn", toString(l), sep = ""),
              paste("40x60x60"," - nn", toString(l), sep = ""))
  
  
  par(mfcol=c(3,2),mar=c(3.5,4.25,2,1),mgp = c(2,1,0))#,mar = c(3.5,3.5,3.5,2.5))
  for (i in 1:5){
    plot(rvals[,1], env.plot.big[,1], type="n", main=titles[i],
         xlab="r", ylab=expression('G'[3]*'(r)'*'  Anomaly'),
         ylim=ylim, xlim=xlim,
         cex.main = 1.75, cex.lab = 1.75, cex.axis = 1.25)
    
    polygon(c(rvals[,i],rev(rvals[,i])),c(env.plot.big[,i],rev(env.plot.small[,i])),col=color[i])
    lines(rvals[,i],clust.plot.big[,i],col="black",lwd=2)
    lines(rvals[,i],clust.plot.small[,i],col="black",lwd=2)
    abline(h=0,lty=2,lwd=1,col="black")
    if(i == 1){
      text(3.5,0.35,"Random: 99.9% AI", pos=4,cex = 1.25)
      text(3.5,0.3,"Clusters: 95% AI",pos=4,cex = 1.25)
      text(3.75,0.25,"50% Den",pos=4,cex = 1.25)
      text(3.75,0.2,"R = 2",pos=4,cex = 1.25)
      
      #legend(20,14,c("Random - 99.9% AI", "Clusters - 95% AI"), col = c(color[2],"black"), lty = c(1,1,1), lwd = c(10,2), bty = "n")
    }
  }
}



