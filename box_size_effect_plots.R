# box size effect plots
library(data.table)
library(rapt)
library(RColorBrewer)

percentile <- .999
plot.big <- matrix(NaN,200,10)
plot.small <- matrix(NaN,200,10)
rvals <- matrix(NaN,200,10)

#### RRL Envelopes only ####
for(i in (1:10)){
  rtemp <- fread(paste('~/Research/box_size_effect/RRL',toString(i),'.csv',sep=""),select=2) # r values
  tvals <- fread(paste('~/Research/box_size_effect/RRL',toString(i),'.csv',sep=""),drop=c(1,2)) # results 
  nTests <- ncol(tvals) # number of tests done
  prange <- percentile*nTests # get the range of indeces for which each percentile spans
  
  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  
  ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  ind.small <- round(nTests/2)-floor(prange/2) # do the same for the low end
  
  plot.big[,i] <- sortedtVals[,ind.big]
  plot.small[,i] <- sortedtVals[,ind.small]
  rvals[,i] <- as.numeric(rtemp$V1)
  
  rm(rtemp,tvals,nTests,prange,sortedtVals,ind.big,ind.small)
  gc()
}


color = brewer.pal(10, "Paired")
for(i in 1:2){
  temp <- color[2*i]
  color[2*i] <- color[2*i+5]
  color[2*i+5] <- temp
}
temp <- color[1:5]
color[1:5] <- color[6:10]
color[6:10] <- temp

xlim = c(0,max(rvals))
ylim = c(min(plot.small),max(plot.big))
  
par(mgp = c(2,1,0),mar = c(3.5,3,3.5,2.5))
plot(rvals[,1],plot.big[,1],type="n",main="Box Size Comparison - 99.9% AI",xlab="r",ylab=expression(sqrt('K'[3]*'(r)')*'  Anomaly'),ylim=ylim,xlim=xlim)
for(i in 1:10){
  lines(rvals[,i],plot.big[,i],col=color[i],lwd = 2)
  lines(rvals[,i],plot.small[,i],col=color[i],lwd = 2)
}
abline(h=0,lty=2,lwd=1,col="black")
legend(0,11, legend=c("15x15x15",
                           "20x20x20",
                           "30x30x30",
                           "40x40x40",
                           "60x60x60"),
       col=color[1:5],lty=rep(1,10),lwd=rep(2,10),bty="n")
legend(0,-4.2, legend=c("10x60x60",
                           "15x60x60",
                           "20x60x60",
                           "30x60x60",
                           "40x60x60"),
       col=color[6:10],lty=rep(1,10),lwd=rep(2,10),bty="n")

#### Cluster only #### 
percentile <- .95
plot.big <- matrix(NaN,200,10)
plot.small <- matrix(NaN,200,10)
plot.toSub <- matrix(NaN,200,10)

rvals <- matrix(NaN,200,10)

for(i in (1:10)){
  rtemp <- fread(paste('~/Research/box_size_effect/result',toString(i),'_',toString(1),'.csv',sep=""),select=1) # r values
  tvals <- matrix(NaN,200,40)
  for(j in (1:40)){
    tvals[,j] <- fread(paste('~/Research/box_size_effect/result',toString(i),'_',toString(j),'.csv',sep=""),select=2)$trans # results 
  }
  nTests <- ncol(tvals) # number of tests done
  prange <- percentile*nTests # get the range of indeces for which each percentile spans
  
  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  
  ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  ind.small <- round(nTests/2)-floor(prange/2) # do the same for the low end
  
  plot.big[,i] <- sortedtVals[,ind.big]
  plot.small[,i] <- sortedtVals[,ind.small]
  plot.toSub[,i] <- (sortedtVals[,20]+sortedtVals[,21])/2
  rvals[,i] <- as.numeric(rtemp$r)
  
  rm(rtemp,tvals,nTests,prange,sortedtVals,ind.big,ind.small)
  gc()
}

#Uncomment for cluster anom 
#plot.big <- plot.big-plot.toSub
#plot.small <- plot.small-plot.toSub


color = brewer.pal(10, "Paired")
for(i in 1:2){
  temp <- color[2*i]
  color[2*i] <- color[2*i+5]
  color[2*i+5] <- temp
}
temp <- color[1:5]
color[1:5] <- color[6:10]
color[6:10] <- temp

xlim = c(0,max(rvals))
ylim = c(min(plot.small),max(plot.big))

par(mgp = c(2,1,0),mar = c(3.5,3.5,3.5,2.5))
plot(rvals[,1],plot.big[,1],type="n",main="Box Size Comparison - Clusters - 100% density - 95% AI",xlab="r",ylab=expression(sqrt('K'[3]*'(r)')*'  Anomaly'),ylim=ylim,xlim=xlim)
for(i in 1:10){
  lines(rvals[,i],plot.big[,i],col=color[i],lwd = 2)
  lines(rvals[,i],plot.small[,i],col=color[i],lwd = 2)
}
abline(h=0,lty=2,lwd=1,col="black")
legend(-3,16.5, legend=c("15x15x15",
                      "20x20x20",
                      "30x30x30",
                      "40x40x40",
                      "60x60x60"),
       col=color[1:5],lty=rep(1,10),lwd=rep(2,10),bty="n")
legend(-3,-9, legend=c("10x60x60",
                        "15x60x60",
                        "20x60x60",
                        "30x60x60",
                        "40x60x60"),
       col=color[6:10],lty=rep(1,10),lwd=rep(2,10),bty="n")

#### Plot each cluster against its envelope ####
percentile_env <- .999
percentile_clust <- .95
env.plot.big <- matrix(NaN,200,10)
env.plot.small <- matrix(NaN,200,10)
clust.plot.big <- matrix(NaN,200,10)
clust.plot.small <- matrix(NaN,200,10)
rvals <- matrix(NaN,200,10)

for(i in (1:10)){
  rtemp <- fread(paste('~/Research/box_size_effect/RRL',toString(i),'.csv',sep=""),select=2) # r values
  tvals <- fread(paste('~/Research/box_size_effect/RRL',toString(i),'.csv',sep=""),drop=c(1,2)) # results 
  nTests <- ncol(tvals) # number of tests done
  prange <- percentile_env*nTests # get the range of indeces for which each percentile spans
  
  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  
  ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  ind.small <- round(nTests/2)-floor(prange/2) # do the same for the low end
  
  env.plot.big[,i] <- sortedtVals[,ind.big]
  env.plot.small[,i] <- sortedtVals[,ind.small]
  rvals[,i] <- as.numeric(rtemp$V1)
  rm(rtemp,tvals,nTests,prange,sortedtVals,ind.big,ind.small)
  ##
  tvals <- matrix(NaN,200,40)
  for(j in (1:40)){
    tvals[,j] <- fread(paste('~/Research/box_size_effect/result',toString(i),'_',toString(j),'.csv',sep=""),select=2)$trans # results 
  }
  nTests <- ncol(tvals) # number of tests done
  prange <- percentile_clust*nTests # get the range of indeces for which each percentile spans
  
  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  
  ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  ind.small <- round(nTests/2)-floor(prange/2) # do the same for the low end
  
  clust.plot.big[,i] <- sortedtVals[,ind.big]
  clust.plot.small[,i] <- sortedtVals[,ind.small]
  
  rm(tvals,nTests,prange,sortedtVals,ind.big,ind.small)
  gc()

}

color = brewer.pal(10, "Paired")
for(i in 1:2){
  temp <- color[2*i]
  color[2*i] <- color[2*i+5]
  color[2*i+5] <- temp
}
temp <- color[1:5]
color[1:5] <- color[6:10]
color[6:10] <- temp

xlim = c(0,max(rvals))
ylim = c(min(clust.plot.small),max(clust.plot.big))
titles = c("15x15x15",
           "20x20x20",
           "30x30x30",
           "40x40x40",
           "60x60x60",
           "10x60x60",
           "15x60x60",
           "20x60x60",
           "30x60x60",
           "40x60x60")

par(mfrow=c(3,3),mar=c(3,3.5,2,1))#,mgp = c(2,1,0),mar = c(3.5,3.5,3.5,2.5))
for (i in 2:10){
  plot(rvals[,1],env.plot.big[,1],type="n",main=titles[i],xlab="r",ylab=expression(sqrt('K'[3]*'(r)')*'  Anomaly'),ylim=ylim,xlim=xlim)
  polygon(c(rvals[,i],rev(rvals[,i])),c(env.plot.big[,i],rev(env.plot.small[,i])),col=color[i])
  lines(rvals[,i],clust.plot.big[,i],col="black",lwd=1.5)
  lines(rvals[,i],clust.plot.small[,i],col="black",lwd=1.5)
  abline(h=0,lty=2,lwd=1,col="black")
  if(i ==2){
    text(0,-12.2,"Random - 99.9% AI",pos=4)
    text(0,-15,"Cluster - 95% AI",pos=4)
  }
}
