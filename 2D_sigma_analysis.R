## Cool 2D RCP Picture

library(spatstat)
library(plotrix)
library(shape)
source("R/rapt-file.R")

rcp_1 <- read.table("C:/Users/galen/Documents/Research/point_patterns/2D/042618_2D_1.2_FinalConfig",sep = " ",col.names = c("x","y","type"))
r1s <- 0.00982629
r1b <- 0.0117916

scale <- 4
rcp_2 <- rcp_1
rcp_2$x <- rcp_1$x*scale
rcp_2$y <- rcp_1$y*scale-.2

x=.35
y=.65

plot.box <- owin(c(x,y),c(x,y))
rcp_1.inboxtf <- inside.owin(rcp_1$x,rcp_1$y,plot.box)
rcp_2.inboxtf<- inside.owin(rcp_2$x,rcp_2$y,plot.box)
rcp_1.inbox <- rcp_1[rcp_1.inboxtf,]
rcp_2.inbox <- rcp_2[rcp_2.inboxtf,]

par(pty="s")
plot(rcp_1.inbox$x,rcp_1.inbox$y,type="n",xlim=c(x,y),ylim=c(x,y),xlab="",ylab="",main="Point and Cluster Locations",xaxt='n',yaxt='n')
#for(i in 1:length(rcp_1.inbox$x)){
#  if(rcp_1.inbox$type[i]==1){
#    draw.circle(rcp_1.inbox$x[i],rcp_1.inbox$y[i],r1b,border="black",col=NA,lty=1,lwd=1.75)
#  }
#  if(rcp_1.inbox$type[i]==2){
#    draw.circle(rcp_1.inbox$x[i],rcp_1.inbox$y[i],r1s,border="black",col=NA,lty=1,lwd=1.75)
#  }
#}
#for(i in 1:length(rcp_2.inbox$x)){
#  if(rcp_2.inbox$type[i]==1){
#    draw.circle(rcp_2.inbox$x[i],rcp_2.inbox$y[i],radius=r1b*scale,border="red",col=NA,lty=1,lwd=2.25)
#  }
#  if(rcp_2.inbox$type[i]==2){
#    draw.circle(rcp_2.inbox$x[i],rcp_2.inbox$y[i],radius=r1s*scale,border="red",col=NA,lty=1,lwd=2.25)
#  }
#}
points(rcp_1.inbox$x,rcp_1.inbox$y,pch=20,cex=1,xlim=c(x,y),ylim=c(x,y))
points(rcp_2.inbox$x,rcp_2.inbox$y,pch=17,cex=1.25,xlim=c(x,y),ylim=c(x,y),col='red')

#Arrows(rcp_2.inbox$x[3],rcp_2.inbox$y[3],rcp_2.inbox$x[2],rcp_2.inbox$y[2],arr.type="triangle",code=3,col="blue",arr.adj = 1,lwd=2,arr.length = .25)

a <- as.ppp(rcp_1.inbox,W=plot.box)
b<- as.ppp(rcp_2.inbox,W=plot.box)

#nn <- nncross(b,a,what="which",k=1:5)
#for(i in 1:npoints(b)){
#  for(j in 1:5){
#    points(rcp_1.inbox$x[nn[i,j]],rcp_1.inbox$y[nn[i,j]],pch=15,cex=1,col="red")
#  }
#}

#for(i in 1:length(rcp_2.inbox$x)){
#  draw.circle(rcp_2.inbox$x[i],rcp_2.inbox$y[i],radius=.03,border="red",col=NA,lty=2,lwd=2.25)
#}

c <- crosspairs(b,a,.03,what="indices")

#for(i in 1:length(c[[1]])){
#  points(rcp_1.inbox$x[c[[2]]],rcp_1.inbox$y[c[[2]]],pch=15,cex=1,col="red")
#}
