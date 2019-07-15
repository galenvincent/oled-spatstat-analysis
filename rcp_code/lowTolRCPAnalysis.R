library(spatstat)

source("R/rapt-file.R")

# import rcp generated points
rcp_1 <- read.table("C:/Users/galen/Documents/Research/point_patterns/022818rcp_1_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
rcp_1_xyz <- rcp_1[,c("x","y","z")]
rcp_1_pp3 <- createSpat(rcp_1_xyz)
# expand to be .5 unit radius (.5nm for a particle)
rcp_1_scaled <- createSpat((.5/0.0360893)*rcp_1_xyz)
plot(rcp_1_scaled)

rcp_2 <- read.table("C:/Users/galen/Documents/Research/point_patterns/040318rcp_1_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
rcp_2_xyz <- rcp_2[,c("x","y","z")]
rcp_2_pp3 <- createSpat(rcp_2_xyz)
# expand to be .5 unit radius (.5nm for a particle)
rcp_2_scaled <- createSpat((.5/0.02594)*rcp_2_xyz)
plot(rcp_2_scaled)

############################################################################

# Calculate

a <- list(K3est(rcp_1_scaled,rmax=5,nrval=200),K3est(rcp_2_scaled,rmax=5,nrval=200))
b <- list(G3est(rcp_1_scaled,rmax=1.25,nrval=500),G3est(rcp_2_scaled,rmax=1.25,nrval=500))
c <- list(F3est(rcp_1_scaled,rmax=1.1,nrval=500),F3est(rcp_2_scaled,rmax=1.1,nrval=500))
d <- list(pcf3est(rcp_1_scaled,rmax=5,nrval=200),pcf3est(rcp_2_scaled,rmax=5,nrval=200))
#d2 <- list(pcf3est(rcp_1_scaled,rmax=5,nrval=200,delta=.25),pcf3est(rcp_2_scaled,rmax=5,nrval=200,delta=.25))

e <- cbind(a[[1]]$r,(a[[1]]$iso-a[[2]]$iso)*100/70)
f <- cbind(b[[1]]$r,(b[[1]]$rs-b[[2]]$rs)*100/1)
g <- cbind(c[[1]]$r,(c[[1]]$rs-c[[2]]$rs)*100/1)
h <- cbind(d[[1]]$r,(d[[1]]$iso-d[[2]]$iso)*100/1.91)



############################

# Plot the residuals

plot(e[,1],e[,2],type="l",lwd=2,main = "K3est Difference: High Tol RCP - Low Tol RCP",xlab = "r",ylab = "K3(r) Difference (% of Max on Interval)",xlim = c(0,5),ylim=c(-.7,.7))
abline(v=1,lwd = 2,lty="dashed")
grid()

plot(f[,1],f[,2],type="l",lwd=2,main = "G3est Difference: High Tol RCP - Low Tol RCP",xlab = "r",ylab = "G3(r) Difference (% of Max on Interval)",xlim = c(.9,1.1),ylim=c(-.1,.7))
abline(v=1,lwd = 2,lty="dashed")
grid()

plot(g[,1],g[,2],type="l",lwd=2,main = "F3est Difference: High Tol RCP - Low Tol RCP",xlab = "r",ylab = "F3(r) Difference (% of Max on Interval)",xlim = c(0,1),ylim=c(-5,5))
abline(v=1,lwd = 2,lty="dashed")
grid()

plot(h[,1],h[,2],type="l",lwd=2,main = "PCF Difference: High Tol RCP - Low Tol RCP",xlab = "r",ylab = "PCF(r) Difference (% of Max on Interval)",xlim = c(0,5),ylim=c(-1.5,1.5))
abline(v=1,lwd = 2,lty="dashed")
grid()

#################

# Plot the straight up values


plot(a[[1]]$r,a[[1]]$theo,type="l",lwd=2,col="black",main = "K3est on Low Tolerance RCP",xlab = "r",ylab = "K3(r)",xlim = c(0,2.5),ylim=c(0,70))
lines(a[[1]]$r,a[[1]]$iso,lwd=2,col="red")
lines(a[[2]]$r,a[[2]]$iso,lwd=2,lty="dashed",col="blue")
abline(v=1,lwd = 2,lty="dashed")
legend(0,60,c("Theoretical (Poisson)","Low Tol RCP","Regular RCP","Expected Inhibition"),lty=c(1,1,2,2),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","black"),box.lty=0)

plot(b[[1]]$r,b[[1]]$theo,type="l",lwd=2,col="black",main = "G3est on Low Tolerance RCP",xlab = "r",ylab = "G3(r)",xlim = c(0,2.5),ylim=c(0,1.1))
lines(b[[1]]$r,b[[1]]$rs,lwd=2,col="red")
lines(b[[2]]$r,b[[2]]$rs,lwd=2,col="blue",lty="dashed")
abline(v=1,lwd = 2,lty="dashed")
legend(1.5,.5,c("Theoretical (Poisson)","Low Tol RCP","Regular RCP","Expected Inhibition"),lty=c(1,1,2,2),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","black"),box.lty=0)

plot(c[[1]]$r,c[[1]]$theo,type="l",lwd=2,col="black",main = "F3est on Low Tolerance RCP",xlab = "r",ylab = "F3(r)",xlim = c(0,2.5),ylim=c(0,1.1))
lines(c[[1]]$r,c[[1]]$rs,lwd=2,col="red")
lines(c[[2]]$r,c[[2]]$rs,lwd=2,col="blue",lty="dashed")
abline(v=1,lwd = 2,lty="dashed")
legend(1.25,.4,c("Theoretical (Poisson)","Low Tol RCP","Regular RCP","Expected Inhibition"),lty=c(1,1,2,2),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","black"),box.lty=0)

plot(d[[1]]$r,d[[1]]$theo,type="l",lwd=2,col="black",main = "PCF on Low Tolerance RCP",xlab = "r",ylab = "PCF(r)",xlim = c(0,2.5),ylim=c(0,1.9))
lines(d[[1]]$r,d[[1]]$iso,lwd=2,col="red")
lines(d[[2]]$r,d[[2]]$iso,lwd=2,col="blue",lty="dashed")
abline(v=1,lwd = 2,lty="dashed")
legend(0,1.75,c("Theoretical (Poisson)","Low Tol RCP","Regular RCP","Expected Inhibition"),lty=c(1,1,2,2),lwd=c(2.5,2.5,2.5,2.5),col=c("black","red","blue","black"),box.lty=0)
