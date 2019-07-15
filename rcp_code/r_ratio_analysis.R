library(spatstat)

source("R/rapt-file.R")
source("R/rcp_functions.R")

# Import the 3 different radius ratio files
rcp_1_upload <- read.table("C:/Users/galen/Documents/Research/point_patterns/041718rcp_1_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
rcp_11 <- scaleRCP(createSpat(rcp_1_upload[,c("x","y","z")]),newRadius = 0.5,oldRadius = 0.0259528)

rcp_2_upload <- read.table("C:/Users/galen/Documents/Research/point_patterns/041718rcp_2_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
rcp_12 <- scaleRCP(createSpat(rcp_2_upload[,c("x","y","z")]),newRadius = 0.5,oldRadius = 0.0252329)

rcp_3_upload <- read.table("C:/Users/galen/Documents/Research/point_patterns/041718rcp_3_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
rcp_13 <- scaleRCP(createSpat(rcp_3_upload[,c("x","y","z")]),newRadius = 0.5,oldRadius = 0.0244723)

# Caculate function on the different patterns

a <- list(K3est(rcp_11,rmax=5,nrval=300),K3est(rcp_12,rmax=5,nrval=300),K3est(rcp_13,rmax=5,nrval=300))
b <- list(G3est(rcp_11,rmax=1.5,nrval=1000),G3est(rcp_12,rmax=1.5,nrval=1000),G3est(rcp_13,rmax=1.5,nrval=1000))
c <- list(F3est(rcp_11,rmax=1.5,nrval=1000),F3est(rcp_12,rmax=1.5,nrval=1000),F3est(rcp_13,rmax=1.5,nrval=1000))
d <- list(pcf3est(rcp_11,rmax=5,nrval=300),pcf3est(rcp_12,rmax=5,nrval=300),pcf3est(rcp_13,rmax=5,nrval=300))

# Plot the different functions

plot(a[[1]]$r,a[[1]]$theo,type="l",lwd=2,col="black",main = "K3est on r Ratio RCPs",xlab = "r",ylab = "K3(r)",xlim = c(0,2.5),ylim=c(0,70))
lines(a[[1]]$r,a[[1]]$iso,lwd=2,col="red")
lines(a[[2]]$r,a[[2]]$iso,lwd=2,col="blue")
lines(a[[3]]$r,a[[3]]$iso,lwd=2,col="forestgreen")
abline(v=1,lwd = 2,lty="dashed")
legend(0,60,c("Theoretical (Poisson)",expression(paste(sigma," = 1.1")),expression(paste(sigma," = 1.2")),expression(paste(sigma," = 1.3")),"Expected Inhibition"),lty=c(1,1,1,1,2),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","forestgreen","black"),box.lty=0)

plot(b[[1]]$r,b[[1]]$theo,type="l",lwd=2,col="black",main = "G3est on r Ratio RCPs",xlab = "r",ylab = "G3(r)",xlim = c(0,1.5),ylim=c(0,1.1))
lines(b[[1]]$r,b[[1]]$rs,lwd=2,col="red")
lines(b[[2]]$r,b[[2]]$rs,lwd=2,col="blue")
lines(b[[3]]$r,b[[3]]$rs,lwd=2,col="forestgreen")
abline(v=1,lwd = 2,lty="dashed")
legend(0,1,c("Theoretical (Poisson)",expression(paste(sigma," = 1.1")),expression(paste(sigma," = 1.2")),expression(paste(sigma," = 1.3")),"Expected Inhibition"),lty=c(1,1,1,1,2),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","forestgreen","black"),box.lty=0)

plot(c[[1]]$r,c[[1]]$theo,type="l",lwd=2,col="black",main = "F3est on r Ratio RCPs",xlab = "r",ylab = "F3(r)",xlim = c(0,1.5),ylim=c(0,1.1))
lines(c[[1]]$r,c[[1]]$rs,lwd=2,col="red")
lines(c[[2]]$r,c[[2]]$rs,lwd=2,col="blue")
lines(c[[3]]$r,c[[3]]$rs,lwd=2,col="forestgreen")
abline(v=1,lwd = 2,lty="dashed")
legend(0,1,c("Theoretical (Poisson)",expression(paste(sigma," = 1.1")),expression(paste(sigma," = 1.2")),expression(paste(sigma," = 1.3")),"Expected Inhibition"),lty=c(1,1,1,1,2),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","forestgreen","black"),box.lty=0)

plot(d[[1]]$r,d[[1]]$theo,type="l",lwd=2,col="black",main = "PCF on r Ratio RCPs",xlab = "r",ylab = "PCF(r)",xlim = c(0,5),ylim=c(0,1.9))
lines(d[[1]]$r,d[[1]]$iso,lwd=2,col="red")
lines(d[[2]]$r,d[[2]]$iso,lwd=2,col="blue")
lines(d[[3]]$r,d[[3]]$iso,lwd=2,col="forestgreen")
abline(v=1,lwd = 2,lty="dashed")
legend(3,.75,c("Theoretical (Poisson)",expression(paste(sigma," = 1.1")),expression(paste(sigma," = 1.2")),expression(paste(sigma," = 1.3")),"Expected Inhibition"),lty=c(1,1,1,1,2),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","forestgreen","black"),box.lty=0)
