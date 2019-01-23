# Effect of Box Size on Function Results
# Written to be tested locally, run on HPC for full analysis
library(rapt)
library(parallel)
library(data.table)

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

  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
#make cluster
  cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=2)
  
#test on each of the different cluster sizes below here
  
  for(j in(1:10)){
    #cut it down to size
    cluster.cut <- subSquare(cluster[[1]],cluster.sizes[[j]])
    
    #test
    result <- anomK3est(cluster.cut,toSub[[j]],max(env.r[[j]]),200,correction="trans")
    
    #output
    fwrite(result,paste('~/Research/box_size_effect/density_1/result',toString(j),'_',toString(i),'.csv',sep=""),sep=',')
    print(paste(toString(i),"_",toString(j)))
    rm(cluster.cut)
    gc()
  }
#clear memory
  rm(under,over,under.big,over.big,cluster,result)
  gc()
#repeat
}

#### 50% CLUSTER DENSITY
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
