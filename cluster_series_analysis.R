# Moving on from the box size effect analysis, this script begins to analyze
# different returns of K-functions when cluster parameters are changed

library(data.table)
library(zoo)
library(rapt)
library(parallel)

# First, establish all of the information that we will collect...
  # First peak radius: [K_max] for longhand, Rm for shorthand
    #Info for each cluster realization + Mean and SD
  # First peak height: K_max for longhand, Km for shorthand
    #Info for each cluster realization + Mean and SD
  # First radius where K'(r) is minimized: [K_max*] for longhand, Rdm for shorthand
    #Info for each cluster realization + Mean and SD
  # Value of K at [K_max*]: Kdm
    # Info for each cluster realization + Mean and SD
  # Cluster separation: csep
  # First radius where K''(r) is maximized: Rddm
 

# Detector Efficiency -----------------------------------------------------

# Detector efficiency series - select out [10, 20, 30, 40, 50, 60]% of points
# from clusters at random, see how the K-function changes
rm(list = ls())
gc()

under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

# detector efficiency series values
de <- list(0.9, 0.8, 0.7, 0.6, 0.5, 0.4)
n <- length(de)

Rm <- matrix(0,101,n)
Km <- matrix(0,101,n)
Rdm <- matrix(0,101,n)
Kdm <- matrix(0,101,n)
Rddm <- matrix(0,101,n)

csep <- list()
for(i in 1:n){
  csep[[i]] <- list()
}

set.seed(10)

cluster.size <- c(60,60,60)

toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)
##toSub <- fread('~/scratch/Rcode/cubetoSub.csv',drop=1)
##env.r <- fread('~/scratch/Rcode/cube.csv',select=2)

p <- length(under.nums)

#loops
for(i in(1:p)){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  #make cluster
  cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=2)
  
  #test on different detector efficiencies below here
  for(j in 1:n){
    #browser()
    #remove the proper number of points
    cluster.de <- percentSelect(de[[j]],cluster[[1]])
    #test
    result <- anomK3est(cluster.de,toSub,max(env.r),nrow(env.r),correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    #plot(rvals,tvals,type = "n",xlab = "", ylab = "")
    rvals.new <- rvals[13:length(rvals)]
    tvals.new <- tvals[13:length(rvals)]
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new)
    
    Km[i,j] <- metrics[[1]]
    Rm[i,j] <- metrics[[2]]
    Rdm[i,j] <- metrics[[3]]
    Rddm[i,j] <- metrics[[4]]
    Kdm[i,j] <- metrics[[5]]
    
    csep[[j]][[i]] <- cluster[[4]]
    
    #output
    #fwrite(peak,paste('~/Research/K_cluster_series/de/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    ##fwrite(peak,paste('~/scratch/Rcode/de/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    
    print(paste(toString(i),"_",toString(j)))
    rm(cluster.de, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  #clear memory
  rm(under,over,under.big,over.big,cluster)
  gc()
  #repeat
}

Rm <- as.data.frame(Rm)
Km <- as.data.frame(Km)
Rdm <- as.data.frame(Rdm)
Kdm <- as.data.frame(Kdm)
Rddm <- as.data.frame(Rddm)

fwrite(Rm,'~/Research/K_cluster_series/de/Rm.csv')
fwrite(Km,'~/Research/K_cluster_series/de/Km.csv')
fwrite(Rdm,'~/Research/K_cluster_series/de/Rdm.csv')
fwrite(Kdm,'~/Research/K_cluster_series/de/Kdm.csv')
fwrite(Rddm,'~/Research/K_cluster_series/de/Rddm.csv')

##fwrite(Rm,'~/scratch/Rcode/de/Rm.csv')
##fwrite(Km,'~/scratch/Rcode/de/Km.csv')
##fwrite(Rdm,'~/scratch/Rcode/de/Rdm.csv')
##fwrite(Kdm,'~/scratch/Rcode/de/Kdm.csv')
##fwrite(Rddm,'~/scratch/Rcode/de/Rddm.csv')

maxnc <- sapply(csep, function(x){max(sapply(x, length))})
for(i in 1:n){
  csepexp <- matrix(NA,p,maxnc[[i]])
  for(j in 1:p){
    csepexp[j,(1:length(csep[[i]][[j]]))] <- csep[[i]][[j]]
  }
  csepexp <- as.data.frame(csepexp)
  fwrite(csepexp,paste('~/Research/K_cluster_series/de/csep',toString(i),'.csv',sep = ''))
  
  ##fwrite(csepexp,paste('~/scratch/Rcode/den/csep',toString(i),'.csv',sep = ''))
}



# Cluster Radius ----------------------------------------------------------

# Cluster radius series
rm(list=ls())
# see how the K-function changes
under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

set.seed(10)

cluster.size <- c(60,60,60)

# detector efficiency series values
r <- list(1, 2, 3, 4, 5, 8)
n <- length(r)

Rm <- matrix(0,101,n)
Km <- matrix(0,101,n)
Rdm <- matrix(0,101,n)
Kdm <- matrix(0,101,n)
Rddm <- matrix(0,101,n)

csep <- list()
for(i in 1:n){
  csep[[i]] <- list()
}

toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)
##toSub <- fread('~/scratch/Rcode/cubetoSub.csv',drop=1)
##env.r <- fread('~/scratch/Rcode/cube.csv',select=2)

p <- length(under.nums)

#loop
for(i in(1:p)){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  #test on different radius values below here
  for(j in 1:n){
    #make cluster
    cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=r[[j]],toPlot = TRUE)

    #test
    result <- anomK3est(cluster[[1]],toSub,max(env.r),nrow(env.r),correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    #plot(rvals,tvals,type = "n",xlab = "", ylab = "")
    rvals.new <- rvals[13:length(rvals)]
    tvals.new <- tvals[13:length(rvals)]
    
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new)
    
    Km[i,j] <- metrics[[1]]
    Rm[i,j] <- metrics[[2]]
    Rdm[i,j] <- metrics[[3]]
    Rddm[i,j] <- metrics[[4]]
    Kdm[i,j] <- metrics[[5]]
    
    csep[[j]][[i]] <- cluster[[4]]
    
    #output
    #fwrite(peak,paste('~/Research/K_cluster_series/radius/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    ##fwrite(peak,paste('~/scratch/Rcode/rad/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    
    
    print(paste(toString(i),"_",toString(j)))
    rm(cluster, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
}


Rm <- as.data.frame(Rm)
Km <- as.data.frame(Km)
Rdm <- as.data.frame(Rdm)
Kdm <- as.data.frame(Kdm)
Rddm <- as.data.frame(Rddm)

fwrite(Rm,'~/Research/K_cluster_series/rad/Rm.csv')
fwrite(Km,'~/Research/K_cluster_series/rad/Km.csv')
fwrite(Rdm,'~/Research/K_cluster_series/rad/Rdm.csv')
fwrite(Kdm,'~/Research/K_cluster_series/rad/Kdm.csv')
fwrite(Rddm,'~/Research/K_cluster_series/rad/Rddm.csv')
##fwrite(Rm,'~/scratch/Rcode/rad/Rm.csv')
##fwrite(Km,'~/scratch/Rcode/rad/Km.csv')
##fwrite(Rdm,'~/scratch/Rcode/rad/Rdm.csv')
##fwrite(Kdm,'~/scratch/Rcode/rad/Kdm.csv')
##fwrite(Rddm,'~/scratch/Rcode/rad/Rddm.csv')

maxnc <- sapply(csep, function(x){max(sapply(x, length))})
for(i in 1:n){
  csepexp <- matrix(NA,p,maxnc[[i]])
  for(j in 1:p){
    csepexp[j,(1:length(csep[[i]][[j]]))] <- csep[[i]][[j]]
  }
  csepexp <- as.data.frame(csepexp)
  fwrite(csepexp,paste('~/Research/K_cluster_series/rad/csep',toString(i),'.csv',sep = ''))
  
  ##fwrite(csepexp,paste('~/scratch/Rcode/rad/csep',toString(i),'.csv',sep = ''))
}


# Cluster Density ---------------------------------------------------------

# Cluster density series 
rm(list=ls())
gc()
# see how the K-function changes
under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

set.seed(10)

cluster.size <- c(60,60,60)

r <- list(1, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15, 0.1)
n <- length(r)

Rm <- matrix(0,101,n)
Km <- matrix(0,101,n)
Rdm <- matrix(0,101,n)
Kdm <- matrix(0,101,n)
Rddm <- matrix(0,101,n)

csep <- list()
for(i in 1:n){
  csep[[i]] <- list()
}

toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)
##toSub <- fread('~/scratch/Rcode/cubetoSub.csv',drop=1)
##env.r <- fread('~/scratch/Rcode/cube.csv',select=2)

p <- length(under.nums)

#loop
for(i in(1:p)){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  #test on different densities
  for(j in 1:n){
    #make cluster
    cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=2,den = r[[j]])
    
    #test
    result <- anomK3est(cluster[[1]],toSub,max(env.r),nrow(env.r),correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    #plot(rvals,tvals,type = "n",xlab = "", ylab = "")
    rvals.new <- rvals[13:length(rvals)]
    tvals.new <- tvals[13:length(rvals)]
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new)
    
    Km[i,j] <- metrics[[1]]
    Rm[i,j] <- metrics[[2]]
    Rdm[i,j] <- metrics[[3]]
    Rddm[i,j] <- metrics[[4]]
    Kdm[i,j] <- metrics[[5]]
    
    csep[[j]][[i]] <- cluster[[4]]
    
    #output
    #fwrite(peak,paste('~/Research/K_cluster_series/den/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    ##fwrite(peak,paste('~/scratch/Rcode/den/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
  
    print(paste(toString(i),"_",toString(j)))
    rm(cluster, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
}

Rm <- as.data.frame(Rm)
Km <- as.data.frame(Km)
Rdm <- as.data.frame(Rdm)
Kdm <- as.data.frame(Kdm)
Rddm <- as.data.frame(Rddm)

fwrite(Rm,'~/Research/K_cluster_series/den/Rm.csv')
fwrite(Km,'~/Research/K_cluster_series/den/Km.csv')
fwrite(Rdm,'~/Research/K_cluster_series/den/Rdm.csv')
fwrite(Kdm,'~/Research/K_cluster_series/den/Kdm.csv')
fwrite(Rddm,'~/Research/K_cluster_series/den/Rddm.csv')
##fwrite(Rm,'~/scratch/Rcode/den/Rm.csv')
##fwrite(Km,'~/scratch/Rcode/den/Km.csv')
##fwrite(Rdm,'~/scratch/Rcode/den/Rdm.csv')
##fwrite(Kdm,'~/scratch/Rcode/den/Kdm.csv')
##fwrite(Rddm,'~/scratch/Rcode/den/Rddm.csv')

maxnc <- sapply(csep, function(x){max(sapply(x, length))})
for(i in 1:n){
  csepexp <- matrix(NA,p,maxnc[[i]])
  for(j in 1:p){
    csepexp[j,(1:length(csep[[i]][[j]]))] <- csep[[i]][[j]]
  }
  csepexp <- as.data.frame(csepexp)
  fwrite(csepexp,paste('~/Research/K_cluster_series/den/csep',toString(i),'.csv',sep = ''))
  
  ##fwrite(csepexp,paste('~/scratch/Rcode/den/csep',toString(i),'.csv',sep = ''))
}



# Den + Rad ---------------------------------------------------------------

# Cluster density + radius size series - [1, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15]
rm(list=ls())
# see how the K-function changes
under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

set.seed(10)

cluster.size <- c(60,60,60)

den <- list(1, 0.6, 0.4, 0.3, 0.2, 0.15)
r <- list(2, 3 ,4, 5, 6, 7)

m <- length(den)
n <- length(r)

Rm <- matrix(0,101,n*m)
Km <- matrix(0,101,n*m)
Rdm <- matrix(0,101,n*m)
Kdm <- matrix(0,101,n*m)
Rddm <- matrix(0,101,n*m)

csep <- list()
for(i in 1:n*m){
  csep[[i]] <- list()
}

toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)
##toSub <- fread('~/scratch/Rcode/cubetoSub.csv',drop=1)
##env.r <- fread('~/scratch/Rcode/cube.csv',select=2)

p <- length(under.nums)

#loop
for(i in(1:p)){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  #test on different densities
  count <- 1
  for(j in 1:n){
    for(k in 1:m){
      #make cluster
      cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=r[[j]],den = den[[k]])
      
      #test
      result <- anomK3est(cluster[[1]],toSub,max(env.r),nrow(env.r),correction="trans")
      rvals <- result$r
      tvals <- result$trans
      
      # get out that peak info son
      #plot(rvals,tvals)
      rvals.new <- rvals[13:length(rvals)]
      tvals.new <- tvals[13:length(rvals)]
      
      #get those metrics out
      metrics <- k3metrics(rvals.new, tvals.new)
      
      Km[i,count] <- metrics[[1]]
      Rm[i,count] <- metrics[[2]]
      Rdm[i,count] <- metrics[[3]]
      Rddm[i,count] <- metrics[[4]]
      Kdm[i,count] <- metrics[[5]]
      
      csep[[count]][[i]] <- cluster[[4]]
      count <- count + 1
      
      #output
      #fwrite(peak,paste('~/Research/K_cluster_series/denr/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
      ##fwrite(peak,paste('~/scratch/Rcode/denr/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
      
      print(paste(toString(i),"_",toString(j),"_",toString(k)))
      rm(cluster, result, rvals, tvals, rvals.new, tvals.new)
      gc()
    }
  }
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
}

Rm <- as.data.frame(Rm)
Km <- as.data.frame(Km)
Rdm <- as.data.frame(Rdm)
Kdm <- as.data.frame(Kdm)
Rddm <- as.data.frame(Rddm)

fwrite(Rm,'~/Research/K_cluster_series/denr/Rm.csv')
fwrite(Km,'~/Research/K_cluster_series/denr/Km.csv')
fwrite(Rdm,'~/Research/K_cluster_series/denr/Rdm.csv')
fwrite(Kdm,'~/Research/K_cluster_series/denr/Kdm.csv')
fwrite(Rddm, '~/Research/K_cluster_series/denr/Rddm.csv')
##fwrite(Rm,'~/scratch/Rcode/denr/Rm.csv')
##fwrite(Km,'~/scratch/Rcode/denr/Km.csv')
##fwrite(Rdm,'~/scratch/Rcode/denr/Rdm.csv')
##fwrite(Kdm,'~/scratch/Rcode/denr/Kdm.csv')
##fwrite(Rddm,'~/scratch/Rcode/denr/Rddm.csv')

maxnc <- sapply(csep, function(x){max(sapply(x, length))})
for(i in 1:n){
  csepexp <- matrix(NA,p,maxnc[[i]])
  for(j in 1:p){
    csepexp[j,(1:length(csep[[i]][[j]]))] <- csep[[i]][[j]]
  }
  csepexp <- as.data.frame(csepexp)
  fwrite(csepexp,paste('~/Research/K_cluster_series/denr/csep',toString(i),'.csv',sep = ''))
  
  ##fwrite(csepexp,paste('~/scratch/Rcode/denr/csep',toString(i),'.csv',sep = ''))
}


# Position Blur -----------------------------------------------------------

# Gaussian Blur of Cluster centers Series
rm(list=ls())
# see how the K-function changes
under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

set.seed(10)

cluster.size <- c(60,60,60)

# percent of seperation distance for blur series
gbperc <- list(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
n <- length(gbperc)

Rm <- matrix(0,101,n)
Km <- matrix(0,101,n)
Rdm <- matrix(0,101,n)
Kdm <- matrix(0,101,n)
Rddm <- matrix(0,101,n)

csep <- list()

for(i in 1:n){
  csep[[i]] <- list()
}

toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)
##toSub <- fread('~/scratch/Rcode/cubetoSub.csv',drop=1)
##env.r <- fread('~/scratch/Rcode/cube.csv',select=2)

p <- length(under.nums)

sep <- 0
cnt <- 0

#loop
for(i in(1:p)){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  if(i == 1){
    #figure out what the cluster separation is
    a <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=3)
    sep <- mean(a[[4]])
    rm(a)
  }
  
  #test on different gaussian blurs
  for(j in 1:n){
    cnt <- cnt + 1
    #make cluster
    cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=3,gb=TRUE,gbp=c(0,sep*gbperc[[j]]),gbmethod = 2,s=cnt)
    
    #test
    result <- anomK3est(cluster[[1]],toSub,max(env.r),nrow(env.r),correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    #plot(rvals,tvals,type = "n",xlab = "", ylab = "")
    rvals.new <- rvals[13:length(rvals)]
    tvals.new <- tvals[13:length(rvals)]
    
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new)
    
    Km[i,j] <- metrics[[1]]
    Rm[i,j] <- metrics[[2]]
    Rdm[i,j] <- metrics[[3]]
    Rddm[i,j] <- metrics[[4]]
    Kdm[i,j] <- metrics[[5]]
    
    csep[[j]][[i]] <- cluster[[4]]
    
    #output
    #fwrite(peak,paste('~/Research/K_cluster_series/gb/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    ##fwrite(peak,paste('~/scratch/Rcode/rad/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    
    
    print(paste(toString(i),"_",toString(j)))
    rm(cluster, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
}

Rm <- as.data.frame(Rm)
Km <- as.data.frame(Km)
Rdm <- as.data.frame(Rdm)
Kdm <- as.data.frame(Kdm)
Rddm <- as.data.frame(Rddm)

fwrite(Rm,'~/Research/K_cluster_series/gb/Rm.csv')
fwrite(Km,'~/Research/K_cluster_series/gb/Km.csv')
fwrite(Rdm,'~/Research/K_cluster_series/gb/Rdm.csv')
fwrite(Kdm,'~/Research/K_cluster_series/gb/Kdm.csv')
fwrite(Rddm,'~/Research/K_cluster_series/gb/Rddm.csv')
##fwrite(Rm,'~/scratch/Rcode/gb/Rm.csv')
##fwrite(Km,'~/scratch/Rcode/gb/Km.csv')
##fwrite(Rdm,'~/scratch/Rcode/gb/Rdm.csv')
##fwrite(Kdm,'~/scratch/Rcode/gb/Kdm.csv')
##fwrite(Rddm,'~/scratch/Rcode/gb/Rddm.csv')

maxnc <- sapply(csep, function(x){max(sapply(x, length))})
for(i in 1:n){
  csepexp <- matrix(NA,p,maxnc[[i]])
  for(j in 1:p){
    csepexp[j,(1:length(csep[[i]][[j]]))] <- csep[[i]][[j]]
  }
  csepexp <- as.data.frame(csepexp)
  fwrite(csepexp,paste('~/Research/K_cluster_series/gb/csep',toString(i),'.csv',sep = ''))
  
  ##fwrite(csepexp,paste('~/scratch/Rcode/gb/csep',toString(i),'.csv',sep = ''))
}



# Radius Blur -------------------------------------------------------------
# Gaussian Blur of cluster radius series
rm(list=ls())
# see how the K-function changes
under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

set.seed(10)

cluster.size <- c(60,60,60)

# percent of cluster radius to set blur sds
rbperc <- list(0, 0.05, 0.1, 0.2, 0.3, 0.5,0.6)
n <- length(rbperc)

Rm <- matrix(0,101,n)
Km <- matrix(0,101,n)
Rdm <- matrix(0,101,n)
Kdm <- matrix(0,101,n)
Rddm <- matrix(0,101,n)

csep <- list()
crs <- list()

for(i in 1:n){
  csep[[i]] <- list()
  crs[[i]] <- list()
}

toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)
##toSub <- fread('~/scratch/Rcode/cubetoSub.csv',drop=1)
##env.r <- fread('~/scratch/Rcode/cube.csv',select=2)

p <- length(under.nums)
#p <- 2
cnt <- 0
tot <- p*n

#loop
for(i in(1:p)){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  #test on different radius gaussian blurs
  for(j in 1:n){
    cnt <- cnt + 1
    #make cluster
    cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=3,rb=TRUE,rbp=3*rbperc[[j]],s=cnt)
    
    #test
    result <- anomK3est(cluster[[1]],toSub,max(env.r),nrow(env.r),correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son

    #plot(rvals,tvals,type = "n",xlab = "", ylab = "", ylim = c(-10,25))
    rvals.new <- rvals[13:length(rvals)]
    tvals.new <- tvals[13:length(rvals)]
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new)
    
    Km[i,j] <- metrics[[1]]
    Rm[i,j] <- metrics[[2]]
    Rdm[i,j] <- metrics[[3]]
    Rddm[i,j] <- metrics[[4]]
    Kdm[i,j] <- metrics[[5]]
    
    csep[[j]][[i]] <- cluster[[4]]
    crs[[j]][[i]] <- cluster[[5]]
    
    #output
    #fwrite(peak,paste('~/Research/K_cluster_series/gb/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    ##fwrite(peak,paste('~/scratch/Rcode/rad/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    
    perc_done <- round((cnt/tot)*100)
    print(paste(toString(i),"_",toString(j)," - ",toString(perc_done),"%"))
    rm(cluster, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
}


Rm <- as.data.frame(Rm)
Km <- as.data.frame(Km)
Rdm <- as.data.frame(Rdm)
Kdm <- as.data.frame(Kdm)
Rddm <- as.data.frame(Rddm)

fwrite(Rm,'~/Research/K_cluster_series/rb/Rm.csv')
fwrite(Km,'~/Research/K_cluster_series/rb/Km.csv')
fwrite(Rdm,'~/Research/K_cluster_series/rb/Rdm.csv')
fwrite(Kdm,'~/Research/K_cluster_series/rb/Kdm.csv')
fwrite(Rddm,'~/Research/K_cluster_series/rb/Rddm.csv')
##fwrite(Rm,'~/scratch/Rcode/rb/Rm.csv')
##fwrite(Km,'~/scratch/Rcode/rb/Km.csv')
##fwrite(Rdm,'~/scratch/Rcode/rb/Rdm.csv')
##fwrite(Kdm,'~/scratch/Rcode/rb/Kdm.csv')
##fwrite(Rddm,'~/scratch/Rcode/rb/Rddm.csv')

maxnc <- sapply(crs, function(x){max(sapply(x, length))})
for(i in 1:n){
  crsexp <- matrix(NA,p,maxnc[[i]])
  csepexp <- matrix(NA,p,maxnc[[i]])
  for(j in 1:p){
    crsexp[j,(1:length(crs[[i]][[j]]))] <- crs[[i]][[j]]
    csepexp[j,(1:length(csep[[i]][[j]]))] <- csep[[i]][[j]]
  }
  crsexp <- as.data.frame(crsexp)
  csepexp <- as.data.frame(csepexp)
  fwrite(crsexp,paste('~/Research/K_cluster_series/rb/crs',toString(i),'.csv',sep = ''))
  fwrite(csepexp,paste('~/Research/K_cluster_series/rb/csep',toString(i),'.csv',sep = ''))
  
  ##fwrite(crsexp,paste('~/scratch/Rcode/rb/crs',toString(i),'.csv',sep = ''))
  ##fwrite(csepexp,paste('~/scratch/Rcode/rb/csep',toString(i),'.csv',sep = ''))
}

# Binomial Radius Blur ----------------------------------------------------

# Two different radii with different proportions series
rm(list=ls())
# see how the K-function changes
under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

set.seed(10)

cluster.size <- c(60,60,60)

# percent of cluster radius to set blur sds
bp <- list(0, 0.2, 0.4, 0.6, 0.8, 1)
n <- length(bp)

Rm <- matrix(0,101,n)
Km <- matrix(0,101,n)
Rdm <- matrix(0,101,n)
Kdm <- matrix(0,101,n)
Rddm <- matrix(0,101,n)

csep <- list()
crs <- list()

for(i in 1:n){
  csep[[i]] <- list()
  crs[[i]] <- list()
}

toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)
##toSub <- fread('~/scratch/Rcode/cubetoSub.csv',drop=1)
##env.r <- fread('~/scratch/Rcode/cube.csv',select=2)

p <- length(under.nums)
#p <- 2
cnt <- 0
tot <- p*n

#loop
for(i in(1:p)){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  #test on different radius gaussian blurs
  for(j in 1:n){
    cnt <- cnt + 1
    #make cluster
    cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",rb=TRUE,rbmethod = 2,rbp=c(bp[[j]],2,5),s=cnt)
    
    #test
    result <- anomK3est(cluster[[1]],toSub,max(env.r),nrow(env.r),correction="trans")
    rvals <- result$r
    tvals <- result$trans
    
    # get out that peak info son
    rvals.new <- rvals[13:length(rvals)]
    tvals.new <- tvals[13:length(rvals)]
    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new, TRUE)
    
    Km[i,j] <- metrics[[1]]
    Rm[i,j] <- metrics[[2]]
    Rdm[i,j] <- metrics[[3]]
    Rddm[i,j] <- metrics[[4]]
    Kdm[i,j] <- metrics[[5]]
    
    csep[[j]][[i]] <- cluster[[4]]
    crs[[j]][[i]] <- cluster[[5]]
    
    #output
    #fwrite(peak,paste('~/Research/K_cluster_series/gb/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    ##fwrite(peak,paste('~/scratch/Rcode/rad/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
    
    perc_done <- round((cnt/tot)*100)
    print(paste(toString(i),"_",toString(j)," - ",toString(perc_done),"%"))
    rm(cluster, result, rvals, tvals, rvals.new, tvals.new)
    gc()
  }
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
}

Rm <- as.data.frame(Rm)
Km <- as.data.frame(Km)
Rdm <- as.data.frame(Rdm)
Kdm <- as.data.frame(Kdm)
Rddm <- as.data.frame(Rddm)

fwrite(Rm,'~/Research/K_cluster_series/rbbinom/Rm.csv')
fwrite(Km,'~/Research/K_cluster_series/rbbinom/Km.csv')
fwrite(Rdm,'~/Research/K_cluster_series/rbbinom/Rdm.csv')
fwrite(Kdm,'~/Research/K_cluster_series/rbbinom/Kdm.csv')
fwrite(Rddm,'~/Research/K_cluster_series/rbbinom/Rddm.csv')
##fwrite(Rm,'~/scratch/Rcode/rbbinom/Rm.csv')
##fwrite(Km,'~/scratch/Rcode/rbbinom/Km.csv')
##fwrite(Rdm,'~/scratch/Rcode/rbbinom/Rdm.csv')
##fwrite(Kdm,'~/scratch/Rcode/rbbinom/Kdm.csv')
##fwrite(Rddm,'~/scratch/Rcode/rbbinom/Rddm.csv')

maxnc <- sapply(crs, function(x){max(sapply(x, length))})
for(i in 1:n){
  crsexp <- matrix(NA,p,maxnc[[i]])
  csepexp <- matrix(NA,p,maxnc[[i]])
  for(j in 1:p){
    crsexp[j,(1:length(crs[[i]][[j]]))] <- crs[[i]][[j]]
    csepexp[j,(1:length(csep[[i]][[j]]))] <- csep[[i]][[j]]
  }
  crsexp <- as.data.frame(crsexp)
  csepexp <- as.data.frame(csepexp)
  fwrite(crsexp,paste('~/Research/K_cluster_series/rbbinom/crs',toString(i),'.csv',sep = ''))
  fwrite(csepexp,paste('~/Research/K_cluster_series/rbbinom/csep',toString(i),'.csv',sep = ''))
  
  ##fwrite(crsexp,paste('~/scratch/Rcode/rbbinom/crs',toString(i),'.csv',sep = ''))
  ##fwrite(csepexp,paste('~/scratch/Rcode/rbbinom/csep',toString(i),'.csv',sep = ''))
}


# Den + Rad + Rad Blur ----------------------------------------------------

# Density + Radius + Radius Blur
# Cluster density + radius size + radius blur series 

rm(list=ls())
gc()
# see how the K-function changes
under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

cluster.size <- c(60,60,60)

den <- list(1, 0.4, 0.2, 0.15)
r <- list(2, 3, 4, 5)
rbperc <- list(0, 0.2, 0.4, 0.6)

m <- length(den)
n <- length(r)
o <- length(rbperc)

Rm <- matrix(0,101,n*m*o)
Km <- matrix(0,101,n*m*o)
Rdm <- matrix(0,101,n*m*o)
Kdm <- matrix(0,101,n*m*o)
Rddm <- matrix(0,101,n*m*o)

csep <- list()
crs <- list()
for(i in 1:(m*n*o)){
  csep[[i]] <- list()
  crs[[i]] <- list()
}

toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)
##toSub <- fread('~/scratch/Rcode/cubetoSub.csv',drop=1)
##env.r <- fread('~/scratch/Rcode/cube.csv',select=2)

p <- length(under.nums)
cnt <- 0
tot <- p*m*n*o

#loop
for(i in(1:p)){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  #test on different densities
  count <- 1
  for(l in 1:o){
    for(j in 1:n){
      for(k in 1:m){
        cnt <- cnt + 1
        #make cluster
        cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=r[[j]],den = den[[k]],rb = TRUE, rbp = r[[j]]*rbperc[[l]], s = cnt)
        
        #test
        result <- anomK3est(cluster[[1]],toSub,max(env.r),nrow(env.r),correction="trans")
        rvals <- result$r
        tvals <- result$trans
        
        # get out that peak info son
        #plot(rvals,tvals, type = "n")
        #points(rvals,tvals)
        rvals.new <- rvals[13:length(rvals)]
        tvals.new <- tvals[13:length(rvals)]
        
        #get those metrics out
        metrics <- k3metrics(rvals.new, tvals.new)
        
        Km[i,count] <- metrics[[1]]
        Rm[i,count] <- metrics[[2]]
        Rdm[i,count] <- metrics[[3]]
        Rddm[i,count] <- metrics[[4]]
        Kdm[i,count] <- metrics[[5]]
        
        csep[[count]][[i]] <- cluster[[4]]
        crs[[count]][[i]] <- cluster[[5]]
        
        count <- count + 1
        
        #output
        #fwrite(peak,paste('~/Research/K_cluster_series/denr/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
        ##fwrite(peak,paste('~/scratch/Rcode/denr/result',toString(i),'_',toString(j),'.csv',sep=""),sep=',')
        
        perc_done <- round((cnt/tot)*100)
        print(paste(toString(i),"_",toString(l),"_",toString(j),"_",toString(k)," - ",toString(perc_done),"%"))
        rm(cluster, result, rvals, tvals, rvals.new, tvals.new, perc_done)
        gc()
      }
    }
  }
  #clear memory
  rm(under,over,under.big,over.big)
  gc()
  #repeat
}

Rm <- as.data.frame(Rm)
Km <- as.data.frame(Km)
Rdm <- as.data.frame(Rdm)
Kdm <- as.data.frame(Kdm)
Rddm <- as.data.frame(Rddm)

fwrite(Rm,'~/Research/K_cluster_series/denrb/Rm.csv')
fwrite(Km,'~/Research/K_cluster_series/denrb/Km.csv')
fwrite(Rdm,'~/Research/K_cluster_series/denrb/Rdm.csv')
fwrite(Kdm,'~/Research/K_cluster_series/denrb/Kdm.csv')
fwrite(Rddm, '~/Research/K_cluster_series/denrb/Rddm.csv')
##fwrite(Rm,'~/scratch/Rcode/denrb/Rm.csv')
##fwrite(Km,'~/scratch/Rcode/denrb/Km.csv')
##fwrite(Rdm,'~/scratch/Rcode/denrb/Rdm.csv')
##fwrite(Kdm,'~/scratch/Rcode/denrb/Kdm.csv')
##fwrite(Rddm,'~/scratch/Rcode/denrb/Rddm.csv')

maxnc <- sapply(crs, function(x){max(sapply(x, length))})
for(i in 1:(n*m*o)){
  crsexp <- matrix(NA,p,maxnc[[i]])
  csepexp <- matrix(NA,p,maxnc[[i]])
  for(j in 1:p){
    print(paste(toString(i),"_",toString(j)," - length: ",toString(length(crs[[i]][[j]]))))
    crsexp[j,(1:length(crs[[i]][[j]]))] <- crs[[i]][[j]]
    csepexp[j,(1:length(csep[[i]][[j]]))] <- csep[[i]][[j]]
  }
  crsexp <- as.data.frame(crsexp)
  csepexp <- as.data.frame(csepexp)
  fwrite(crsexp,paste('~/Research/K_cluster_series/denrb/crs',toString(i),'.csv',sep = ''))
  fwrite(csepexp,paste('~/Research/K_cluster_series/denrb/csep',toString(i),'.csv',sep = ''))
  
  ##fwrite(crsexp,paste('~/scratch/Rcode/denrb/crs',toString(i),'.csv',sep = ''))
  ##fwrite(csepexp,paste('~/scratch/Rcode/denrb/csep',toString(i),'.csv',sep = ''))
}


