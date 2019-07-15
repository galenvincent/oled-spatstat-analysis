# File to check if the homewritten border correction is working as expected...
# checking by comparing to the 2D version
library(spatstat)

a <- rpoispp(1000)

b <- Kest(a,rmax=.25,correction="border")
c <- bK2est(a,rmax = .25,nrval=513)

d <- b$border-c$bord

plot(b$r,d/b$border)

####################################3
bK2est <- function(X,rmax=NULL,nrval=128){
  
  verifyclass(X,"ppp")
  
  bi <- bdist.points2(X)
  n <- npoints(X)
  lambda <- n/area(domain(X))
  
  if(is.null(rmax)){
    rmax <- max(bi)/2
  }else if(rmax > max(bi)){
    print("rmax is too large for data set")
    return()
  }
  
  cp <- closepairs(X,rmax,twice=FALSE,what="indices")
  cpm <- cbind(cp[[1]],cp[[2]])
  cpm<-cpm[order(cpm[,1]),]
  distmat <- as.matrix(dist(coords(X)))
  cpmdist <- rep(0,nrow(cpm))
  for(i in 1:nrow(cpm)){
    temp <- sort(cpm[i,])
    cpmdist[i] <- distmat[temp[2],temp[1]]
  }
  
  rlist <- seq(from=0,to=rmax,length.out=nrval)
  Kb <- rep(0,nrval)
  
  np <- 0
  for(i in 1:n){
    if(bi[i] >= rmax){
      np <- np + 1
    }
  }
  
  for(j in 1:length(rlist)){
    t <- 0
    r <- rlist[j]
    for(i in 1:nrow(cpm)){
      if(cpmdist[i] <= r){
        if((bi[cpm[i,1]] >= rmax) & (bi[cpm[i,2]] >= rmax)){
          t <- t + 2
        }else if((bi[cpm[i,1]] < rmax) & (bi[cpm[i,2]] < rmax)){
        }else{
          t <- t + 1
        }
      }
    }
    Kb[j] <- t/(lambda*np)
  }
  
  K <- as.data.frame(cbind(rlist,Kb))
  colnames(K)<-c("r","bord")
  
  return(K)
}


bdist.points2 <- function (X) {
  
  verifyclass(X, "ppp")
  
  x <- X$x
  y <- X$y
  #z <- X$data$z
  d <- X$window
  
  xmin <- min(d$xrange)
  xmax <- max(d$xrange)
  ymin <- min(d$yrange)
  ymax <- max(d$yrange)
  #zmin <- min(d$zrange)
  #zmax <- max(d$zrange)
  result <- pmin.int(x - xmin, xmax - x, y - ymin, ymax - y)# , z - zmin , zmax - z)
  
  return(result)
}