# Looking into the strauss model vs matern and rcp data in 2D
library(spatstat)
source("R/rapt-file.R")

# Matern data:
matern_1 <- rMaternII(3000,.043,win = owin(c(0,2),c(0,2)),stationary=FALSE,nsim=1,drop=TRUE)
matern_2 <- rMaternII(3000,.043,win = owin(c(0,2),c(0,2)),stationary=FALSE,nsim=1,drop=TRUE)
npoints(matern_1)
plot(matern_1)

#Strauss data
strauss_1 <- rStrauss(500, gamma = .1, R = .043, W = owin(c(0,2),c(0,2)), expand=TRUE, nsim=1, drop=TRUE)
strauss_2 <- rStrauss(500, gamma = .1, R = .043, W = owin(c(0,2),c(0,2)), expand=TRUE, nsim=1, drop=TRUE)
npoints(strauss_1)
plot(strauss_1)

#Analysis
#Kest
ktest <- list(Kest(matern_1),Kest(matern_2),Kest(strauss_1),Kest(strauss_2))

plot(ktest[[1]]$r,ktest[[1]]$theo,type="l",lwd=2.5,col="black",main = "Kest on Strauss vs Matern Model",xlab = "r",ylab = "K(r)",xlim=c(0,.2),ylim=c(0,.1))
lines(ktest[[1]]$r,ktest[[1]]$iso,lwd=2,col="red")
lines(ktest[[2]]$r,ktest[[2]]$iso,lwd=2,col="green",lty="dashed")
lines(ktest[[3]]$r,ktest[[3]]$iso,lwd=2,col="blue")
lines(ktest[[4]]$r,ktest[[4]]$iso,lwd=2,col="orange",lty="dashed")
legend(0,0.08,c("Theoretical (Poisson)","Matern 1","Matern 2","Strauss 1","Strauss 2"),lty=c(1,1,2,1,2),lwd=c(1.5,1.5,1.5,1.5,1.5),col=c("black","red","green","blue","orange"),box.lty=0)


#Gest
gtest <- list(Gest(matern_1),Gest(matern_2),Gest(strauss_1),Gest(strauss_2))

plot(gtest[[1]]$r,gtest[[1]]$theo,type="l",lwd=2.5,col="black",main = "Gest on Strauss vs Matern Model",xlab = "r",ylab = "G(r)")#,xlim=c(0,2.5),ylim=c(0,70))
lines(gtest[[1]]$r,gtest[[1]]$rs,lwd=2,col="red")
lines(gtest[[2]]$r,gtest[[2]]$rs,lwd=2,col="green",lty="dashed")
lines(gtest[[3]]$r,gtest[[3]]$rs,lwd=2,col="blue")
lines(gtest[[4]]$r,gtest[[4]]$rs,lwd=2,col="orange",lty="dashed")
legend(0,1,c("Theoretical (Poisson)","Matern 1","Matern 2","Strauss 1","Strauss 2"),lty=c(1,1,2,1,2),lwd=c(1.5,1.5,1.5,1.5,1.5),col=c("black","red","green","blue","orange"),box.lty=0)


#Fest
ftest <- list(Fest(matern_1),Fest(matern_2),Fest(strauss_1),Fest(strauss_2))

plot(ftest[[1]]$r,ftest[[1]]$theo,type="l",lwd=2.5,col="black",main = "Fest on Strauss vs Matern Model",xlab = "r",ylab = "F(r)")#,xlim=c(0,2.5),ylim=c(0,70))
lines(ftest[[1]]$r,ftest[[1]]$rs,lwd=2,col="red")
lines(ftest[[2]]$r,ftest[[2]]$rs,lwd=2,col="green",lty="dashed")
lines(ftest[[3]]$r,ftest[[3]]$rs,lwd=2,col="blue")
lines(ftest[[4]]$r,ftest[[4]]$rs,lwd=2,col="orange",lty="dashed")
legend(.1,.4,c("Theoretical (Poisson)","Matern 1","Matern 2","Strauss 1","Strauss 2"),lty=c(1,1,2,1,2),lwd=c(1.5,1.5,1.5,1.5,1.5),col=c("black","red","green","blue","orange"),box.lty=0)


#PCF
ptest <- list(pcf(matern_1),pcf(matern_2),pcf(strauss_1),pcf(strauss_2))

plot(ptest[[1]]$r,ptest[[1]]$theo,type="l",lwd=2.5,col="black",main = "PCF on Strauss vs Matern Model",xlab = "r",ylab = "PCF(r)",xlim=c(0,.2),ylim=c(0,1.3))
lines(ptest[[1]]$r,ptest[[1]]$iso,lwd=2,col="red")
lines(ptest[[2]]$r,ptest[[2]]$iso,lwd=2,col="green",lty="dashed")
lines(ptest[[3]]$r,ptest[[3]]$iso,lwd=2,col="blue")
lines(ptest[[4]]$r,ptest[[4]]$iso,lwd=2,col="orange",lty="dashed")
legend(.125,.5,c("Theoretical (Poisson)","Matern 1","Matern 2","Strauss 1","Strauss 2"),lty=c(1,1,2,1,2),lwd=c(1.5,1.5,1.5,1.5,1.5),col=c("black","red","green","blue","orange"),box.lty=0)

#not a favorable process. Runs for a long time with any high intensity. Non-physical with some very close points together
