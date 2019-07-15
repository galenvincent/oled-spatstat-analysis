# Home Made Matern Proccess
library(spatstat)
source("R/rapt-file.R")

rm(list = ls())

modMatern <- function(inhibitionR = 0.043,nOver = 2,pull = 0,percMove = .5){

  n <- nOver
  percentMove <- percMove

  oglist <- vector("list", n)
  fin <- vector("list",n)
  radiusCheck <- inhibitionR

  for (i in 1:n){
    oglist[[i]] <- rpoispp3(3000, domain = box3(c(0,1),c(0,1),c(0,1)), nsim = 1, drop = TRUE)
  }

  nnInd <- nnwhich(oglist[[1]], k = 1)
  nnDist <- nndist(oglist[[1]], k = 1)

  #######################
  # Original random poisson pattern generation and removal

  print("All set up.")

  i=1
  while(i <= npoints(oglist[[1]])) {
    if(nnDist[i] < radiusCheck){
      oglist[[1]] <- oglist[[1]][-nnInd[i]]
      nnInd <- nnwhich(oglist[[1]], k = 1)
      nnDist <- nndist(oglist[[1]], k = 1)
    }else {
      direction <- (coords(oglist[[1]])[nnInd[i],]-coords(oglist[[1]])[i,])/nnDist[i]
      coords(oglist[[1]])[nnInd[i],] <- coords(oglist[[1]])[nnInd[i],]-(nnDist[i]-radiusCheck)*direction*percentMove
      nnInd <- nnwhich(oglist[[1]], k = 1)
      nnDist <- nndist(oglist[[1]], k = 1)
    }
    i = i+1
  }

  fin[[1]] <- oglist[[1]]

  print("Entering the loop.")

  ##########################################
  # Loop through the put on top patterns, sorting, and pulling closer

  if(pull==0){
    for(j in 2:n){
      print(cbind("Starting Loop # ",toString(j)))
      nn <- nncross(fin[[j-1]],oglist[[j]])

      i<-1
      while(i <= npoints(oglist[[j]])) {
        if(nn$dist[i] < radiusCheck){
          oglist[[j]] <- oglist[[j]][-nn$which[i]]
          nn <- nncross(fin[[j-1]],oglist[[j]])
        }
        i <- i+1
      }

      nnInd <- nnwhich(oglist[[j]], k = 1)
      nnDist <- nndist(oglist[[j]], k = 1)
      i<-1
      while(i <= npoints(oglist[[j]])) {
        if(nnDist[i] < radiusCheck){
          oglist[[j]] <- oglist[[j]][-nnInd[i]]
          nnInd <- nnwhich(oglist[[j]], k = 1)
          nnDist <- nndist(oglist[[j]], k = 1)
        }
        i <- i+1
      }

      binded <- rbind(coords(fin[[j-1]]),coords(oglist[[j]]))
      fin[[j]] <- pp3(binded$x,binded$y,binded$z,as.box3(c(0,1,c(0,1),c(0,1))))

      nnInd <- nnwhich(fin[[j]], k = 1)
      nnDist <- nndist(fin[[j]], k = 1)
      i<-1
      while(i <= npoints(fin[[j]])) {
        if(nnDist[i] < radiusCheck){
          fin[[j]] <- fin[[j]][-nnInd[i]]
          nnInd <- nnwhich(fin[[j]], k = 1)
          nnDist <- nndist(fin[[j]], k = 1)
        }
        i <- i+1
      }
    }
  } else {
    for(j in 2:n){
      print(cbind("Starting Loop # ",toString(j)))
      nn <- nncross(fin[[j-1]],oglist[[j]])

      i<-1
      while(i <= npoints(oglist[[j]])) {
        if(nn$dist[i] < radiusCheck){
          oglist[[j]] <- oglist[[j]][-nn$which[i]]
          nn <- nncross(fin[[j-1]],oglist[[j]])
        }
        i <- i+1
      }

      nnInd <- nnwhich(oglist[[j]], k = 1)
      nnDist <- nndist(oglist[[j]], k = 1)
      i<-1
      while(i <= npoints(oglist[[j]])) {
        if(nnDist[i] < radiusCheck){
          oglist[[j]] <- oglist[[j]][-nnInd[i]]
          nnInd <- nnwhich(oglist[[j]], k = 1)
          nnDist <- nndist(oglist[[j]], k = 1)
        }
        i <- i+1
      }

      binded <- rbind(coords(fin[[j-1]]),coords(oglist[[j]]))
      fin[[j]] <- pp3(binded$x,binded$y,binded$z,as.box3(c(0,1,c(0,1),c(0,1))))

      nnInd <- nnwhich(fin[[j]], k = 1)
      nnDist <- nndist(fin[[j]], k = 1)
      i<-1
      while(i <= npoints(fin[[j]])) {
        if(nnDist[i] < radiusCheck){
          fin[[j]] <- fin[[j]][-nnInd[i]]
          nnInd <- nnwhich(fin[[j]], k = 1)
          nnDist <- nndist(fin[[j]], k = 1)
        }else if (nnDist[i] > radiusCheck*1.05) {
          direction <- (coords(fin[[j]])[nnInd[i],]-coords(fin[[j]])[i,])/nnDist[i]
          coords(fin[[j]])[nnInd[i],] <- coords(fin[[j]])[nnInd[i],]-(nnDist[i]-radiusCheck)*direction*percentMove
          nnInd <- nnwhich(fin[[j]], k = 1)
          nnDist <- nndist(fin[[j]], k = 1)
        }
        i <- i+1
      }
    }

  }

  return(fin)
}

#########################
# Check out the results! (make return the full array to play with this)

mat1 <- rMaternI(3000,radiusCheck,win = box3(c(0,1),c(0,1),c(0,1)),stationary = TRUE,nsim=1,drop=TRUE)
mat2 <- rMaternII(3000,radiusCheck,win = box3(c(0,1),c(0,1),c(0,1)),stationary = TRUE,nsim=1,drop=TRUE)


g3est <- list(G3est(mat1,rmax=.2,nrval=100),G3est(fin[[2]],rmax=.2,nrval=100),G3est(fin[[3]],rmax=.2,nrval=100),G3est(fin[[10]],rmax=.2,nrval=100))
plot(g3est[[1]]$r,g3est[[1]]$theo,xlim=c(0,.1),type = "l",lwd = 2.5,col="black",main="G3est of Matern vs Modified Matern",ylab="G(r)",xlab="r")
lines(g3est[[1]]$r,g3est[[1]]$rs,lwd=2,col="red")
#lines(g3est[[2]]$r,g3est[[2]]$rs,lwd=2,col="orange")
lines(g3est[[3]]$r,g3est[[3]]$rs,lwd=2,col="blue")
lines(g3est[[4]]$r,g3est[[4]]$rs,lwd=2,col="forestgreen")
legend(.001,.8,c("Theoretical (Poisson)","Matern","Modified Matern (1)","Modified Matern (10)"),lty=c(1,1,1,1),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","forestgreen"),box.lty=0)

p3est <- list(pcf3est(mat1,rmax=.2,nrval=100),pcf3est(fin[[1]],rmax=.2,nrval=100),pcf3est(fin[[2]],rmax=.2,nrval=100),pcf3est(fin[[10]],rmax=.2,nrval=100))
plot(p3est[[1]]$r,p3est[[1]]$theo,type = "l",lwd = 2.5,col="black",xlim = c(0,.15),ylim = c(0,1.2),main="PCF of Matern vs Modified Matern",ylab="PCF(r)",xlab="r")
lines(p3est[[1]]$r,p3est[[1]]$iso,lwd=2,col="red")
#lines(p3est[[2]]$r,p3est[[2]]$iso,lwd=2,col="orange")
lines(p3est[[3]]$r,p3est[[3]]$iso,lwd=2,col="blue")
lines(p3est[[4]]$r,p3est[[4]]$iso,lwd=2,col="forestgreen")
legend(.075,.5,c("Theoretical (Poisson)","Matern","Modified Matern (1)","Modified Matern (10)"),lty=c(1,1,1,1),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","forestgreen"),box.lty=0)

f3est <- list(F3est(mat1),F3est(fin[[1]]),F3est(fin[[2]]),F3est(fin[[10]]))
plot(f3est[[1]]$r,f3est[[1]]$theo,type = "l",lwd = 2.5,col="black",xlim = c(0,.15),ylim = c(0,1.2),main="F3est of Matern vs Modified Matern",ylab="F(r)",xlab="r")
lines(f3est[[1]]$r,f3est[[1]]$rs,lwd=2,col="red")
#lines(f3est[[2]]$r,f3est[[2]]$rs,lwd=2,col="orange")
lines(f3est[[3]]$r,f3est[[3]]$rs,lwd=2,col="blue")
lines(f3est[[4]]$r,f3est[[4]]$rs,lwd=2,col="forestgreen")
legend(.075,.5,c("Theoretical (Poisson)","Matern","Modified Matern (1)","Modified Matern (10)"),lty=c(1,1,1,1),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","forestgreen"),box.lty=0)

k3est <- list(K3est(mat1,rmax=.2,nrval=100),K3est(fin[[1]],rmax=.2,nrval=100),K3est(fin[[2]],rmax=.2,nrval=100),K3est(fin[[10]],rmax=.2,nrval=100))
plot(k3est[[1]]$r,k3est[[1]]$theo,type = "l",lwd = 2.5,col="black",xlim = c(0,.1),ylim = c(0,.005),main="K3est of Matern vs Modified Matern",ylab="K(r)",xlab="r")
lines(k3est[[1]]$r,k3est[[1]]$iso,lwd=2,col="red")
#lines(k3est[[2]]$r,k3est[[2]]$iso,lwd=2,col="orange")
lines(k3est[[3]]$r,k3est[[3]]$iso,lwd=2,col="blue")
lines(k3est[[4]]$r,k3est[[4]]$iso,lwd=2,col="forestgreen")
legend(.0,.004,c("Theoretical (Poisson)","Matern","Modified Matern (1)","Modified Matern (10)"),lty=c(1,1,1,1),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","forestgreen"),box.lty=0)


#2D analysis
mat <- rMaternII(80,0.1,win = owin(c(0,1),c(0,1)),stationary = TRUE,nsim=1,drop=TRUE)

plot(begin,type = "n",cex = 2)
lines(og,type = "p",pch=19)
lines(fin,type = "p",pch=10)
lines(fin2,type = "p",pch=19,cex=2)
lines(mat,type = "p",cex=2)

#####################
# Density plot

densityPlot <- function(pp3List,inhibition=0.043){
  n <- length(pp3List)
  d <- oglist <- vector("numeric", n)
  fillPerc <- vector("numeric", n)
  for(i in 1:n){
    d[i] <- npoints(pp3List[[i]])
  }
  fillPerc <- d*(4/3)*pi*(inhibition/2)^3

  plot(1:n,fillPerc,type="p",pch=19,col="red",xlab="Number of Runs",ylab = "Packing Density Percent Volume",main = "% Volume Versus # Trials")

  return(fillPerc)
}
