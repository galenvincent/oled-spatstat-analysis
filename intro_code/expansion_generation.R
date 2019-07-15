#Experimentations in pattern generation
library(spatstat)

source("R/rapt-file.R")

og <- rpoispp3(3000, domain = box3(c(0,1),c(0,1),c(0,1)), nsim = 1, drop = TRUE)

nnInd <- nnwhich(og, k = 1)
nnDist <- nndist(og, k = 1)

radiusCheck <- 0.03

while(any(nnDist < radiusCheck)) {
  for(i in 1:npoints(og)) {
    if(nnDist[i] < radiusCheck) {
      direction <- (coords(og)[nnInd[i],]-coords(og)[i,])/nnDist[i]
      coords(og)[i,] <- coords(og)[i,]+(nnDist[i]-radiusCheck)*direction
    }
  }
  nnInd <- nnwhich(og, k = 1)
  nnDist <- nndist(og, k = 1)
}

og_scaled = createSpat((.5/.0360755)*coords(og))

#plot data

ktest <- list(K3est(matern_1_scaled),K3est(matern_2_scaled),K3est(og_scaled))

plot(ktest[[1]]$r,ktest[[1]]$theo,type="l",lwd=2.5,col="black",main = "K3est on Experimental vs Matern Model",xlab = "r",ylab = "K3(r)",xlim=c(0,2.5),ylim=c(0,70))
lines(ktest[[1]]$r,ktest[[1]]$iso,lwd=2,col="red")
lines(ktest[[2]]$r,ktest[[2]]$iso,lwd=2,col="green",lty="dashed")
lines(ktest[[3]]$r,ktest[[3]]$iso,lwd=2,col="blue")
legend(0,65,c("Theoretical (Poisson)","Matern 1","Matern 2","Experimental"),lty=c(1,1,2,1),lwd=c(1.5,1.5,1.5,1.5),col=c("black","red","green","blue"),box.lty=0)


gtest <- list(G3est(matern_1_scaled),G3est(matern_2_scaled),G3est(og_scaled))

plot(gtest[[1]]$r,gtest[[1]]$theo,type="l",lwd=2.5,col="black",main = "G3est on Experimental vs Matern Model",xlab = "r",ylab = "G3(r)",xlim=c(0,1.5))
lines(gtest[[1]]$r,gtest[[1]]$rs,lwd=2,col="red")
lines(gtest[[2]]$r,gtest[[2]]$rs,lwd=2,col="green",lty="dashed")
lines(gtest[[3]]$r,gtest[[3]]$rs,lwd=2,col="blue")
legend(0,1,c("Theoretical (Poisson)","Matern 1","Matern 2","Experimental"),lty=c(1,1,2,1),lwd=c(1.5,1.5,1.5,1.5),col=c("black","red","green","blue"),box.lty=0)


ftest <- list(F3est(matern_1_scaled),F3est(matern_2_scaled),F3est(og_scaled))

plot(ftest[[1]]$r,ftest[[1]]$theo,type="l",lwd=2.5,col="black",main = "F3est on Experimental vs Matern Model",xlab = "r",ylab = "F3(r)",xlim=c(0,1.75))
lines(ftest[[1]]$r,ftest[[1]]$rs,lwd=2,col="red")
lines(ftest[[2]]$r,ftest[[2]]$rs,lwd=2,col="green",lty="dashed")
lines(ftest[[3]]$r,ftest[[3]]$rs,lwd=2,col="blue")
legend(1,.5,c("Theoretical (Poisson)","Matern 1","Matern 2","Experimental"),lty=c(1,1,2,1),lwd=c(1.5,1.5,1.5,1.5),col=c("black","red","green","blue"),box.lty=0)


ptest <- list(pcf3est(matern_1_scaled),pcf3est(matern_2_scaled),pcf3est(og_scaled))

plot(ptest[[1]]$r,ptest[[1]]$theo,type="l",lwd=2.5,col="black",main = "PCF on Experimental vs Matern Model",xlab = "r",ylab = "F3(r)",xlim=c(0,2),ylim=c(0,1.5))
lines(ptest[[1]]$r,ptest[[1]]$iso,lwd=2,col="red")
lines(ptest[[2]]$r,ptest[[2]]$iso,lwd=2,col="green",lty="dashed")
lines(ptest[[3]]$r,ptest[[3]]$iso,lwd=2,col="blue")
legend(3,.75,c("Theoretical (Poisson)","Matern 1","Matern 2","Experimental"),lty=c(1,1,2,1),lwd=c(1.5,1.5,1.5,1.5),col=c("black","red","green","blue"),box.lty=0)

