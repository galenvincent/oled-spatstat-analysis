library(spatstat)

source("R/rapt-file.R")
source("R/rcp_functions.R")
source("R/cluster_functions.R")
r1 <- 0.0259528 #1.1 sigma
r2 <-0.0252329 #1.2 sigma
r3 <-0.0244723 #1.3 sigma
r4 <-0.0217032 #1.4 sigma 50%
r5 <-0.0237094 #1.4 sigma 75%

rcp_1_upload <- read.table("C:/Users/galen/Documents/Research/point_patterns/042918rcp_1_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
rcp_1 <- scale(createSpat(rcp_1_upload[,c("x","y","z")]),newRadius = 0.5,oldRadius = r5)

a <- pcf3est(rcp_1,rmax=2.5,nrval=200)
plot(a,xlim=c(0,2.5),main="PCF of RCP Pattern")

a$r
b<-a[a$r>1 & a$r <2]
c<-which.min(b$trans)

check.rad<-b$r[c]
#check.rad <- 1.5

nn <- nnR(rcp_1,check.rad)
nnd <- nn[[1]]
nnw <- nn[[2]]

rcp_xyz <- coords(rcp_1)

npts <- vector("numeric",npoints(rcp_1))
for(i in 1:8000){
  npts[i] <- length(nnd[,i][!is.na(nnd[,i])])
}
totalpts <- sum(npts)
xvec <- vector("numeric",totalpts)
yvec <- vector("numeric",totalpts)
zvec <- vector("numeric",totalpts)
nnwlst <- vector("numeric",totalpts)

k = 1
for(i in 1:8000){
  for(j in 1:npts[i]){
    xvec[k] <- (rcp_xyz$x[nnw[j,i]]-rcp_xyz$x[i])/nnd[j,i]
    yvec[k] <- (rcp_xyz$y[nnw[j,i]]-rcp_xyz$y[i])/nnd[j,i]
    zvec[k] <- (rcp_xyz$z[nnw[j,i]]-rcp_xyz$z[i])/nnd[j,i]
    nnwlst[k] <- nnw[j,i]
    k = k+1
  }
}


write.table(xvec, file="C:/Users/galen/Desktop/xvec5.csv",sep=",",row.names=F)
write.table(yvec, file="C:/Users/galen/Desktop/yvec5.csv",sep=",",row.names=F)
write.table(zvec, file="C:/Users/galen/Desktop/zvec5.csv",sep=",",row.names=F)
write.table(nnwlst, file="C:/Users/galen/Desktop/nnw5.csv",sep=",",row.names=F)
write.table(npts, file="C:/Users/galen/Desktop/npt5.csv",sep=",",row.names=F)
#################################
# Pass to mathematica
##################################
# Pass crystal points back in
cry_upload <- read.table("C:/Users/galen/Documents/Research/crystal_checking/5crycoords.csv",sep = ",")
rcp_cry <- createSpat(rcp_xyz[cry_upload[,1],])

plot3d.pp3(rcp_1,col="lightgray")
plot3d.pp3(rcp_cry,col="red",size=5,add=TRUE)
