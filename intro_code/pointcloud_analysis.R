# Import Libraries and Functions
library(spatstat)
library(parallel)

source("R/rapt-file.R")
source("R/galen_baseline.R")


#envelopes for the mattern process
matG3est <- function(nEvals){

  # nEvals = number of times to run K3est

  #Create list with different matern model simulations
  toTest <- rMaternII(4540,.0360755,win = box3(c(0,1),c(0,1),c(0,1)),stationary=TRUE,nsim=nEvals)

  #Scale each entry
  for(i in 1:nEvals){
    toTest[[i]] = createSpat((1/.0360755)*coords(toTest[[i]]))
  }

  # run the test once to fill the r values
  tst <- G3est(toTest[[1]])
  tests <- tst$r

  #find cores and initialize the cluster
  cores2use <- detectCores()-1
  cl <- makeCluster(cores2use)

  # apply K3est function to each of the pp3 patterns in parallel
  result <- parLapply(cl,toTest,G3est)

  # convert the results into the matrix tests
  for(i in 1:length(result)){
    tests <- cbind(tests,result[[i]]$rs)
  }

  # stop the cluster and revert computer to normal
  stopCluster(cl)

  return(tests)
}

matK3est <- function(nEvals){
  # nEvals = number of times to run K3est

  #Create list with different matern model simulations
  toTest <- rMaternII(4540,.0360755,win = box3(c(0,1),c(0,1),c(0,1)),stationary=TRUE,nsim=nEvals)

  #Scale each entry
  for(i in 1:nEvals){
    toTest[[i]] = createSpat((1/.0360755)*coords(toTest[[i]]))
  }

  # run the test once to fill the r values
  tst <- K3est(toTest[[1]])
  tests <- tst$r

  #find cores and initialize the cluster
  cores2use <- detectCores()-1
  cl <- makeCluster(cores2use)

  # apply K3est function to each of the pp3 patterns in parallel
  result <- parLapply(cl,toTest,K3est)

  # convert the results into the matrix tests
  for(i in 1:length(result)){
    tests <- cbind(tests,result[[i]]$iso)
  }

  # stop the cluster and revert computer to normal
  stopCluster(cl)

  return(tests)
}


  # import rcp generated points
  rcp_1 <- read.table("C:/Users/galen/Documents/Research/point_patterns/022718rcp_1_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
  rcp_1_xyz <- rcp_1[,c("x","y","z")]
  rcp_1_pp3 <- createSpat(rcp_1_xyz)
  # expand to be .5 unit radius (.5nm for a particle)
  rcp_1_scaled <- createSpat((.5/.0360755)*rcp_1_xyz)
  plot(rcp_1_scaled)

  rcp_2 <- read.table("C:/Users/galen/Documents/Research/point_patterns/022818rcp_1_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
  rcp_2_xyz <- rcp_2[,c("x","y","z")]
  rcp_2_pp3 <- createSpat(rcp_2_xyz)
  # expand to be .5 unit radius (.5nm for a particle)
  rcp_2_scaled <- createSpat((.5/.0360755)*rcp_2_xyz)
  plot(rcp_2_scaled)

  #import different parameter rcp pattern
  rcp_3 <- read.table("C:/Users/galen/Documents/Research/point_patterns/030818rcp_1_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
  rcp_3_xyz <- rcp_3[,c("x","y","z")]
  rcp_3_pp3 <- createSpat(rcp_3_xyz)
  # expand to be .5 unit radius (.5nm for a particle)
  rcp_3_scaled <- createSpat((.5/.0360755)*rcp_3_xyz)
  plot(rcp_3_scaled)

  #create two matern patterns
  matern_1 <- rMaternII(13000,.043,win = box3(c(0,1),c(0,1),c(0,1)),stationary=FALSE,nsim=1,drop=TRUE)
  npoints(matern_1)
  matern_1_scaled <- createSpat((.5/.0360755)*coords(matern_1))

  matern_2 <- rMaternII(13000,.043,win = box3(c(0,1),c(0,1),c(0,1)),stationary=TRUE,nsim=1,drop=TRUE)
  matern_2_scaled <- createSpat((.5/.0360755)*coords(matern_1))

####################################################################################

# Run tests on all of the point patterns and display comparison plots

#K3est
b1 <- K3est(rcp_1_scaled)
b2 <- K3est(rcp_2_scaled)
b3 <- K3est(matern_1_scaled)
b4 <- K3est(matern_2_scaled)
b5 <- K3est(rcp_3_scaled)

plot(b1$r,b1$theo,type="l",lwd=2,col="black",main = "K3est on RCP vs Matern Model Point Cloud",xlab = "r",ylab = "K3(r)")
lines(b1$r,b1$iso,lwd=2,col="red")
lines(b2$r,b2$iso,lwd=2,col="green",lty="dashed")
lines(b3$r,b3$iso,lwd=2,col="blue")
lines(b4$r,b4$iso,lwd=2,col="orange",lty="dashed")
lines(b5$r,b5$iso,lwd=2,col="pink")
legend(0,6500,c("Theoretical (Poisson)","RCP 1","RCP 2","Matern 1","Matern 2"),lty=c(1,1,2,1,2),lwd=c(1.5,1.5,1.5,1.5,1.5),col=c("black","red","green","blue","orange"),box.lty=0)
#zoomed in
plot(b1$r,b1$theo,type="l",lwd=2,col="black",main = "K3est on RCP vs Matern Model Point Cloud",xlab = "r",ylab = "K3(r)",xlim=c(0,2.5),ylim=c(0,70))
lines(b1$r,b1$iso,lwd=2,col="red")
lines(b2$r,b2$iso,lwd=2,col="green",lty="dashed")
lines(b3$r,b3$iso,lwd=2,col="blue")
lines(b4$r,b4$iso,lwd=2,col="orange",lty="dashed")
lines(b5$r,b5$iso,lwd=2,col="purple")
legend(0,65,c("Theoretical (Poisson)","RCP 1","RCP 2","Matern 1","Matern 2"),lty=c(1,1,2,1,2),lwd=c(1.5,1.5,1.5,1.5,1.5),col=c("black","red","green","blue","orange"),box.lty=0)


#G3est
a1 <- G3est(rcp_1_scaled)
a2 <- G3est(rcp_2_scaled)
a3 <- G3est(matern_1_scaled)
a4 <- G3est(matern_2_scaled)
a5 <- G3est(rcp_3_scaled)

plot(a1$r,a1$theo,type="l",lwd=2.5,col="black",main = "G3est on RCP vs Matern Model Point Cloud",xlab = "r",ylab = "G3(r)",xlim=c(0,1.75))
lines(a1$r,a1$rs,lwd=2,col="red")
lines(a2$r,a2$rs,lwd=2,col="green",lty="dashed")
lines(a3$r,a3$rs,lwd=2,col="blue")
lines(a4$r,a4$rs,lwd=2,col="orange",lty="dashed")
lines(a5$r,a5$rs,lwd=2,col="pink")
legend(0,1,c("Theoretical (Poisson)","RCP 1","RCP 2","Matern 1","Matern 2"),lty=c(1,1,2,1,2),lwd=c(1.5,1.5,1.5,1.5,1.5),col=c("black","red","green","blue","orange"),box.lty=0)


#F3est
c1 <- F3est(rcp_1_scaled)
c2 <- F3est(rcp_2_scaled)
c3 <- F3est(matern_1_scaled)
c4 <- F3est(matern_2_scaled)
c5 <- F3est(rcp_3_scaled)

plot(c1$r,c1$theo,type="l",lwd=2.5,col="black",main = "F3est on RCP vs Matern Model Point Cloud",xlab = "r",ylab = "F3(r)",xlim=c(0,1.75))
lines(c1$r,c1$rs,lwd=2,col="red")
lines(c2$r,c2$rs,lwd=2,col="green",lty="dashed")
lines(c3$r,c3$rs,lwd=2,col="blue")
lines(c4$r,c4$rs,lwd=2,col="orange",lty="dashed")
lines(c5$r,c5$rs,lwd=2,col="pink")
legend(1,.5,c("Theoretical (Poisson)","RCP 1","RCP 2","Matern 1","Matern 2"),lty=c(1,1,2,1,2),lwd=c(1.5,1.5,1.5,1.5,1.5),col=c("black","red","green","blue","orange"),box.lty=0)


#Pair Correlation Function
d1 <- pcf3est(rcp_1_scaled)
d2 <- pcf3est(rcp_2_scaled)
d3 <- pcf3est(matern_1_scaled)
d4 <- pcf3est(matern_2_scaled)
d5 <- pcf3est(rcp_3_scaled)

plot(d1$r,d1$theo,type="l",lwd=2.5,col="black",main = "PCF on RCP vs Matern Model Point Cloud",xlab = "r",ylab = "PCF(r)",xlim=c(0,5),ylim=c(0,2))
lines(d1$r,d1$iso,lwd=2,col="red")
lines(d2$r,d2$iso,lwd=2,col="green",lty="dashed")
lines(d3$r,d3$iso,lwd=2,col="blue")
lines(d4$r,d4$iso,lwd=2,col="orange",lty="dashed")
lines(d5$r,d5$iso,lwd=2,col="purple")
legend(3,.75,c("Theoretical (Poisson)","RCP 1","RCP 2","Matern 1","Matern 2"),lty=c(1,1,2,1,2),lwd=c(1.5,1.5,1.5,1.5,1.5),col=c("black","red","green","blue","orange"),box.lty=0)

###################################################################################

# Last RCP to First RCP Comparisons
#K3et
plot(b1$r,b1$theo,type="l",lwd=2,col="black",main = "K3est on RCP 10% vs 25% r difference",xlab = "r",ylab = "K3(r)",xlim=c(0,2.5),ylim=c(0,70))
lines(b1$r,b1$iso,lwd=2,col="red")
lines(b2$r,b2$iso,lwd=2,col="green",lty="dashed")
lines(b5$r,b5$iso,lwd=2,col="purple")
legend(0,65,c("Theoretical (Poisson)","RCP 10% r diff.","RCP 10% r diff","RCP 25% r diff"),lty=c(1,1,2,1),lwd=c(1.5,1.5,1.5,1.5),col=c("black","red","green","purple"),box.lty=0)

#G3est
plot(a1$r,a1$theo,type="l",lwd=2.5,col="black",main = "G3est on RCP 10% vs 25% r difference",xlab = "r",ylab = "G3(r)",xlim=c(0,1.75))
lines(a1$r,a1$rs,lwd=2,col="red")
lines(a2$r,a2$rs,lwd=2,col="green",lty="dashed")
lines(a5$r,a5$rs,lwd=2,col="purple")
legend(0,1,c("Theoretical (Poisson)","RCP 10% r diff.","RCP 10% r diff","RCP 25% r diff"),lty=c(1,1,2,1),lwd=c(1.5,1.5,1.5,1.5),col=c("black","red","green","purple"),box.lty=0)

#F3est
plot(c1$r,c1$theo,type="l",lwd=2.5,col="black",main = "F3est on RCP 10% vs 25% r difference",xlab = "r",ylab = "F3(r)",xlim=c(0,1.75))
lines(c1$r,c1$rs,lwd=2,col="red")
lines(c2$r,c2$rs,lwd=2,col="green",lty="dashed")
lines(c5$r,c5$rs,lwd=2,col="purple")
legend(1,.5,c("Theoretical (Poisson)","RCP 10% r diff.","RCP 10% r diff","RCP 25% r diff"),lty=c(1,1,2,1),lwd=c(1.5,1.5,1.5,1.5),col=c("black","red","green","purple"),box.lty=0)

#PCF
plot(d1$r,d1$theo,type="l",lwd=2.5,col="black",main = "PCF on RCP 10% vs 25% r difference",xlab = "r",ylab = "PCF(r)",xlim=c(0,5),ylim=c(0,2))
lines(d1$r,d1$iso,lwd=2,col="red")
lines(d2$r,d2$iso,lwd=2,col="green",lty="dashed")
lines(d5$r,d5$iso,lwd=2,col="purple")
legend(3,.75,c("Theoretical (Poisson)","RCP 10% r diff.","RCP 10% r diff","RCP 25% r diff"),lty=c(1,1,2,1),lwd=c(1.5,1.5,1.5,1.5),col=c("black","red","green","purple"),box.lty=0)
