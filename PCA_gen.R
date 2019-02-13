# Script to output random data points for PCA

library(data.table)
library(zoo)
library(rapt)
library(parallel)

rm(list=ls())
gc()

under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)
p <- length(under.nums)

cluster.size <- c(60,60,60)

# Number of random simulations desired
nrand <- 10

# Set up a 5x5x5x5 grid of simulation conditions
grdsize <- 2
r <- seq(1, 7, len = grdsize)
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
set.seed(10)
rr <- runif(nrand, min = 1, max = 7)
denr <- runif(nrand, min = 0.15, max = 1) 
rbpr <- runif(nrand, min = 0, max = 0.6)
gbpr <- runif(nrand, min = 0, max = 0.6)
for(i in 1:nrand){
  tot[[cnt]] <- c(rr[i], denr[i], rbpr[i], gbpr[i])
  cnt <- cnt + 1
}

#Set up the final matrix
outfin <- matrix(NA, nrow = length(tot), ncol = 14)
for(i in 1:length(tot)){
  outfin[i,1:4] <- tot[[i]]  
}

outtemp <- list()
for(i in 1:p){
  outtemp[[i]] <- matrix(NA, nrow = length(tot), ncol = 10)
}

# Upload cube RRL files
toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)

# Max r value and number of r values to go to in the k tests
maxr <- max(env.r)
nr <- nrow(env.r)

cores2use <- detectCores() - 1 
cl <- makePSOCKcluster(cores2use)
clusterEvalQ(cl, library(rapt))
cluserExport(cl, "kseries", "p", "maxr", "nr", "toSub")


for(i in 1:p){
  #upload
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  ##over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[i]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  #clusterExport(cl, "under.big", "over.big")
  
  out <- matrix(unlist(parLapply(cl,tot, kseries, maxr, nr, under.big, over.big, toSub)),nrow = length(tot), ncol = 5, byrow = TRUE)
  
  
}

stopCluster(cl)


##############################
kseries <- function(props, maxr, nr, under, over, toSub){
  #make cluster
  sd <- props[1]*props[3]
  # how to get the seed in here?
  a <- rnorm(1, 0, 3)
  
  # cluster <- makecluster(under,over,0.5,0.5,type="cr",speed="superfast",cr=props[1],den = props[2],rb = TRUE, rbp = sd, gb = TRUE, gbp = c(0,props[4]), s = seed)
  # 
  # #test
  # result <- anomK3est(cluster[[1]],toSub,maxr,nr,correction="trans")
  # rvals <- result$r
  # tvals <- result$trans
  # 
  # # get out that peak info son
  # rvals.new <- rvals[13:length(rvals)]
  # tvals.new <- tvals[13:length(rvals)]
  # 
  # #get those metrics out
  # metrics <- k3metrics(rvals.new, tvals.new, FALSE))
  
  #Km[i,count] <- metrics[[1]]
  #Rm[i,count] <- metrics[[2]]
  #Rdm[i,count] <- metrics[[3]]
  #Rddm[i,count] <- metrics[[4]]
  #Kdm[i,count] <- metrics[[5]]
  
  return(c(props[1],props[2],props[3],props[4],a)) #c(metrics[[1]],metrics[[2]],metrics[[3]],metrics[[4]],metrics[[5]]))
}


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