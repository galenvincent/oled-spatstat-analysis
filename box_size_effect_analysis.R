# box size effect plots
library(data.table)
library(rapt)
library(RColorBrewer)
library(pracma)
library(scales)

percentile <- .999
plot.big <- matrix(NaN,200,10)
plot.small <- matrix(NaN,200,10)
rvals <- matrix(NaN,200,10)

#### RRL Envelopes only ####
for(i in (1:10)){
  rtemp <- fread(paste('Z:/Galen/Box\ Size\ RRLs/RRL',toString(i),'.csv',sep=""),select=2) # r values
  tvals <- fread(paste('Z:/Galen/Box\ Size\ RRLs/RRL',toString(i),'.csv',sep=""),drop=c(1,2)) # results 
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
ylim = c(min(plot.small[,5:10]),max(plot.big[,5:10]))
  
par(mgp = c(2,1,0),mar = c(3.5,4,3.5,2.5))
plot(rvals[,1],plot.big[,1],type="n",main="99.9% AI Envelopes vs Box Size",
     xlab="r",ylab=expression(sqrt('K'[3]*'(r)')*'  Anomaly'),
     ylim=ylim,xlim=xlim,
     cex.lab = 1.75,cex.main = 1.75,cex.axis = 1.25)
for(i in 6:10){
  lines(rvals[,i],plot.big[,i],col=color[i],lwd = 2)
  lines(rvals[,i],plot.small[,i],col=color[i],lwd = 2)
}
abline(h=0,lty=2,lwd=1,col="black")
# legend(34,3, legend=c("15x15x15",
#                            "20x20x20",
#                            "30x30x30",
#                            "40x40x40",
#                            "60x60x60"),
#        col=color[1:5],lty=rep(1,10),lwd=rep(2,10),bty="n",cex = 1.25)
legend(-1,11, legend=c("10x60x60",
                       "15x60x60",
                       "20x60x60",
                       "30x60x60",
                       "40x60x60",
                       "60x60x60"),
       col=c(color[6:10],color[5]),lty=rep(1,6),lwd=rep(2,6),bty="n",cex = 1.25)


#make plots for paper - make sure to run the stuff above
# Get one cluster RRL on there
n <- 523
tvals.clust <- matrix(NaN,200,n)
percentile <- 0.90
for(j in (1:n)){
  tvals.clust[,j] <- fread(paste('~/Research/box_size_effect/190624_density_0.5/result',toString(j),'_5.csv',sep=""),select=2)$trans # results 
}
prange <- percentile*n# get the range of indeces for which each percentile spans

sortedtVals <- t(apply(tvals.clust,1,sort)) # sort the results at each r value from lowest to highest

ind.big <- round(n/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
ind.small <- round(n/2)-floor(prange/2) # do the same for the low end

plot.big.clust <- sortedtVals[,ind.big]
plot.small.clust <- sortedtVals[,ind.small]


#envelopes
newcolors <- c(color[1], color[2], color[3], color[4], color[1], "red")

par(mgp = c(2.5,1,0),mar = c(4,4.5,3.5,2.5))
plot(rvals[,1], plot.big[,1], type="n", main="99.9% AI Envelopes vs Box Size",
     xlab="r (arb.)", ylab= 'T(r)',
     ylim=c(-8,8), xlim=xlim,
     cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.25)
for(i in c(2, 4, 5)){
  a <- c(rvals[,i], rev(rvals[,i]))
  b <- c(plot.big[,i], rev(plot.small[,i]))
  polygon(a, b, col=alpha(newcolors[i], 0.85), border="black", lwd=1.5)
}
lines(rvals[,9], plot.big[,9], col = newcolors[6], lwd = 2, lty = 2)
lines(rvals[,9], plot.small[,9], col = newcolors[6], lwd = 2, lty = 2)

lines(rvals[,5], plot.big.clust, col = 'black', lwd = 2)
lines(rvals[,5], plot.small.clust, col = 'black', lwd = 2)

abline(h=0,lty=2,lwd=1,col="black")
abline(v=7,lty=2,lwd=2,col="black")

legend(10,8.5, legend = c("20x20x20",
                      "40x40x40",
                      "60x60x60"),
       col=c(newcolors[c(2, 4, 5)]), lty=c(1, 1, 1), lwd = c(15, 15, 15), bty="n", cex = 1.25)
legend(10,-5, legend = "30x60x60", col = newcolors[6], lty = 2, lwd = 2, bty = "n", cex = 1.25)
legend(10,-6.5, legend = "60x60x60 clusters", col = 'black', lty = 1, lwd = 2, bty = 'n', cex = 1.25)

#envelope width as f_n of box size
rs <- c(7)
widths <- matrix(NaN,length(rs),ncol(plot.big))

for(j in (1:length(rs))){
  for(i in (1:ncol(plot.big))){
    rind <- which(abs(rs[j]-rvals[,i]) == min(abs(rs[j]-rvals[,i])))
    widths[j,i] <- plot.big[rind,i] #- plot.small[rind,i]
  }
}

x1 <- c(15,20,30,40,60)

par(mfrow = c(1,1),mar = c(4,4,2,2), mgp = c(2.7,1,0))
plot(x1,widths[,1:5], ylim = c(0,max(widths[,1:5])),
     pch = 16, col="blue",
     xlab = 'Cube Side Length (arb.)',ylab = "Envelope Width",
     main = "Cubic Volumes", cex = 1.5, cex.lab = 1.5, cex.axis = 1.5)

mdl <- lm(widths[,1:5] ~ I(x1^(-3/2)))

x <- seq(0, 60^3, 1)
y <- predict(mdl, newdata = data.frame(x1 = x))

lines(x, y, col = "black", lty = 2, lwd = 2)


#### Cluster only #### 
rm(list=ls())
gc()

percentile <- .90
plot.big <- matrix(NaN,200,10)
plot.small <- matrix(NaN,200,10)
plot.toSub <- matrix(NaN,200,10)

rvals <- matrix(NaN,200,10)

n <- 523

for(i in (1:10)){
  rtemp <- fread(paste('~/Research/box_size_effect/190624_density_0.5/result',toString(1),'_',toString(i),'.csv',sep=""),select=1) # r values
  tvals <- matrix(NaN,200,n)
  for(j in (1:n)){
    tvals[,j] <- fread(paste('~/Research/box_size_effect/190624_density_0.5/result',toString(j),'_',toString(i),'.csv',sep=""),select=2)$trans # results 
  }
  nTests <- ncol(tvals) # number of tests done
  prange <- percentile*nTests # get the range of indeces for which each percentile spans
  
  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  
  ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  ind.small <- round(nTests/2)-floor(prange/2) # do the same for the low end
  
  plot.big[,i] <- sortedtVals[,ind.big]
  plot.small[,i] <- sortedtVals[,ind.small]
  plot.toSub[,i] <- sortedtVals[,round(nTests/2)]
  rvals[,i] <- as.numeric(rtemp$r)
  
  rm(rtemp, tvals, nTests, prange, sortedtVals, ind.big, ind.small)
  gc()
}

#Uncomment for cluster anom 
plot.big <- plot.big-plot.toSub
plot.small <- plot.small-plot.toSub


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
ylim = c(min(plot.small[,1:5]),max(plot.big[,1:5]))

par(mgp = c(2.25,1,0),mar = c(3.5,4.25,3.5,2.5))
plot(rvals[,1],plot.big[,1],type="n",main="95% Cluster AI Envelopes vs Box Size",
     xlab="r",ylab=expression(sqrt('K'[3]*'(r)')*'  Anomaly'),
     ylim=ylim,xlim=xlim,
     cex.main = 1.75, cex.lab = 1.75, cex.axis = 1.25)
for(i in 1:5){
  lines(rvals[,i],plot.big[,i],col=color[i],lwd = 2)
  lines(rvals[,i],plot.small[,i],col=color[i],lwd = 2)
}
abline(h=0,lty=2,lwd=1,col="black")
legend(35,-2.75, legend=c("15x15x15",
                      "20x20x20",
                      "30x30x30",
                      "40x40x40",
                      "60x60x60"),
       col=color[1:5],lty=rep(1,10),lwd=rep(2,10),cex = 1.3)
legend(-3,-9, legend=c("10x60x60",
                        "15x60x60",
                        "20x60x60",
                        "30x60x60",
                        "40x60x60"),
       col=color[6:10],lty=rep(1,10),lwd=rep(2,10),bty="n")


#make plots for paper - make sure to run the stuff above
#envelopes
newcolors <- c(color[1], color[2], color[3], color[4], color[1], "red")

plot.inds <- c(2, 4, 5)

par(mgp = c(2.5,1,0),mar = c(4,4.5,3.5,2.5))
plot(rvals[,1], plot.big[,1], type="n", main="90% AI Envelopes vs Box Size",
     xlab="r (arb.)", ylab=expression(sqrt('K'[3]*'(r)')*'  Anomaly'),
     #ylim=c(-max(abs(c(min(plot.small[,plot.inds]), max(plot.big[,plot.inds])))) - 1, max(abs(c(min(plot.small[,plot.inds]), max(plot.big[,plot.inds])))) + 1), 
     ylim = c(-8,8),
     xlim=xlim,
     cex.lab = 1.75, cex.main = 1.75, cex.axis = 1.25)
for(i in c(2, 4, 5)){
  a <- c(rvals[,i], rev(rvals[,i]))
  b <- c(plot.big[,i], rev(plot.small[,i]))
  polygon(a, b, col=alpha(newcolors[i], 0.85), border="black", lwd=1.5)
}
lines(rvals[,9], plot.big[,9], col = newcolors[6], lwd = 2, lty = 2)
lines(rvals[,9], plot.small[,9], col = newcolors[6], lwd = 2, lty = 2)

abline(h=0,lty=2,lwd=1,col="black")
abline(v=7,lty=2, lwd=2, col = "black")

legend(-1,9, legend = c("20x20x20",
                          "40x40x40",
                          "60x60x60"),
       col=c(newcolors[c(2, 4, 5)]), lty=c(1, 1, 1), lwd = c(15, 15, 15), bty="n", cex = 1.25, horiz = FALSE)
legend(-1,-5, legend = "30x60x60", col = newcolors[6], lty = 2, lwd = 2, bty = "n", cex = 1.25)

#envelope width as f_n of box size
rs <- c(7)
widths <- matrix(NaN,length(rs),ncol(plot.big))

for(j in (1:length(rs))){
  for(i in (1:ncol(plot.big))){
    rind <- which(abs(rs[j]-rvals[,i]) == min(abs(rs[j]-rvals[,i])))
    widths[j,i] <- plot.big[rind,i] #- plot.small[rind,i]
  }
}

x1 <- c(20,30,40,60)

par(mfrow = c(1,1),mar = c(4,4,2,2), mgp = c(2.7,1,0))
plot(x1,widths[,2:5], ylim = c(0,max(widths[,2:5])),
     pch = 16, col="blue",
     xlab = 'Cube Side Length (arb.)' ,ylab = "Envelope Width",
     main = "Cubic Volumes", cex = 1.5, cex.lab = 1.5, cex.axis = 1.5)

mdl <- lm(widths[,2:5] ~ I(x1^(-3/2)))
#mdl <- nls(widths[,2:5] ~ a*(x1 - b)^(-3/2), start = list(a = 191, b = 0))

x <- seq(0, 60^3, 1)
y <- predict(mdl, newdata = data.frame(x1 = x))

lines(x, y, col = "black", lty = 2, lwd = 2)


#### Plot each cluster against its envelope ####
rm(list=ls())
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
  prange <- percentile_env*nTests # get the range of indices for which each percentile spans
  
  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  
  ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  ind.small <- round(nTests/2)-floor(prange/2) # do the same for the low end
  
  env.plot.big[,i] <- sortedtVals[,ind.big]
  env.plot.small[,i] <- sortedtVals[,ind.small]
  rvals[,i] <- as.numeric(rtemp$V1)
  rm(rtemp,tvals,nTests,prange,sortedtVals,ind.big,ind.small)
  ##
  tvals <- matrix(NaN,200,101)
  for(j in (1:101)){
    tvals[,j] <- fread(paste('~/Research/box_size_effect/density_1/result',toString(i),'_',toString(j),'.csv',sep=""),select=2)$trans # results 
    print(paste(toString(i)," - ",toString(j)))
    }
  nTests <- ncol(tvals) # number of tests done
  prange <- percentile_clust*nTests # get the range of indices for which each percentile spans
  
  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  
  ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  ind.small <- round(nTests/2)-floor(prange/2) # do the same for the low end
  
  clust.plot.big[,i] <- sortedtVals[,ind.big]
  clust.plot.small[,i] <- sortedtVals[,ind.small]
  
  rm(tvals,nTests,prange,sortedtVals,ind.big,ind.small)
  gc()

}

#find peak values
peak.big.inds <- matrix(NaN,2,10)
peak.big.vals <- matrix(NaN,2,10)
peak.small.inds <- matrix(NaN,2,10)
peak.small.vals <- matrix(NaN,2,10)
peak.size <- matrix(NaN,2,10)
peak.pos <- matrix(NaN,2,10)

nud <- 1
for(i in (1:10)){
  a <- findpeaks(clust.plot.big[,i],nups=nud,ndowns=nud,npeaks = 1)
  b <- findpeaks(-1*clust.plot.big[,i],nups=nud,ndowns=nud,npeaks = 1)
  c <- findpeaks(clust.plot.small[,i],nups=nud,ndowns=nud,npeaks = 1)
  d <- findpeaks(-1*clust.plot.small[,i],nups=nud,ndowns=nud,npeaks = 1)
  
  peak.big.inds[,i] <- c(a[2],b[2])
  peak.big.vals[,i] <- c(a[1],b[1])
  peak.small.inds[,i] <- c(c[2],d[2])
  peak.small.vals[,i] <- c(c[1],d[1])
  
  peak.size[,i] <- c(abs(a[1]-c[1]),abs(b[1]-d[1]))
  peak.pos[,i] <- c((rvals[a[2],i] + rvals[c[2],i])/2, (rvals[b[2],i] + rvals[d[2],i])/2)
  
  plot(clust.plot.small[,i],col="blue",ylim = c(min(clust.plot.small[,i]),max(clust.plot.big[,i])),main = toString(i))
  points(clust.plot.big[,i],col="blue")
  points(a[2],a[1],col="red",pch=5,cex = 1)
  points(b[2],-1*b[1],col="red",pch=5,cex = 1)
  points(c[2],c[1],col="red",pch=5,cex = 1)
  points(d[2],-1*d[1],col="red",pch=5,cex = 1)
  text(a[2],a[1],paste("r = ",toString(round(peak.pos[1,i],2))," & size = ",toString(round(peak.size[1,i],2))),pos=4)
  text(b[2],-1*b[1],paste("r = ",toString(round(peak.pos[2,i],2))," & size = ",toString(round(peak.size[2,i],2))),pos=4)
  
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
ylim = c(-15,15)
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

cb <- brewer.pal(6,"Paired")
color[2] <- cb[2]
color[4] <- cb[4]
color[5] <- cb[6]
color[8] <- cb[1]
color[9] <- cb[3]
color[10] <- cb[5]

par(mfcol=c(3,2),mar=c(3.5,4.25,2,1),mgp = c(2,1,0))#,mar = c(3.5,3.5,3.5,2.5))
for (i in c(2, 4, 5, 8, 10)){
  plot(rvals[,1],env.plot.big[,1],type="n",main=titles[i],
       xlab="r",ylab=expression(sqrt('K'[3]*'(r)')*'  Anomaly'),
       ylim=ylim,xlim=xlim,
       cex.main = 1.75, cex.lab = 1.75, cex.axis = 1.25)
  
  polygon(c(rvals[,i],rev(rvals[,i])),c(env.plot.big[,i],rev(env.plot.small[,i])),col=color[i])
  lines(rvals[,i],clust.plot.big[,i],col="black",lwd=2)
  lines(rvals[,i],clust.plot.small[,i],col="black",lwd=2)
  abline(h=0,lty=2,lwd=1,col="black")
  # text(2,-11,paste("Peak: (",toString(round(peak.pos[1,i],2)),
  #                    ", ",toString(round(peak.big.vals[1,i],2)),
  #                    ", ",toString(round(peak.size[1,i],2)),")"),pos=4)
  # text(-1,-14.5,paste("Trough: (",toString(round(peak.pos[2,i],2)),
  #                   ", ",toString(round(peak.small.vals[2,i],2)),
  #                   ", ",toString(round(peak.size[2,i],2)),")"),pos=4)
    if(i == 2){
    # text(15,14,"Random: 99.9% AI",pos=4,cex = 1.25)
    # text(15,10.5,"Clusters: 95% AI",pos=4,cex = 1.25)
    # text(30,7,"100% Den",pos=4,cex = 1.25)
    # text(30,4,"R = 3",pos=4,cex = 1.25)
    #text(15,-7,"(r, amplitude, spread)",pos=4)
    
    #legend(20,14,c("Random - 99.9% AI", "Clusters - 95% AI"), col = c(color[2],"black"), lty = c(1,1,1), lwd = c(10,2), bty = "n")
    }
  legend(-2,-6,c("Random - 99.9% AI", "Clusters - 95% AI"), col = c(color[i],"black"), lty = c(1,1,1), lwd = c(10,2), bty = "n", cex = 1.1, y.intersp = 0.75)
  # if(i == 7){
  #   text(-0.5,15,"Random: 99.9% AI",pos=4)
  #   text(0,-10,"Clusters: 95% AI",pos=4)
  #   text(12,-13,"100% Den",pos=4)
  #   text(12,-15.5,"R = 3",pos=4)
  #   #text(15,-7,"(r, amplitude, spread)",pos=4)
  # }
}


#### radii at fixed envelope width
widths <- c(5)
rs <- matrix(NaN,length(widths),ncol(env.plot.big))

env.widths <- env.plot.big - env.plot.small

for(j in (1:length(widths))){
  for(i in (1:ncol(env.plot.big))){
    rind <- which(abs(widths[j]-env.widths[,i]) == min(abs(widths[j]-env.widths[,i])))
    rs[j,i] <- rvals[rind,i]
  }
}

par(mfrow = c(1,2),mar = c(4.5,3.5,3,1),mgp = c(2.25,1,0))
x1 <- c(15,20,30,40,60)
x2 <- c(10,15,20,30,40)

plot(x1,rs[1,1:5],pch = 1, col="blue",xlab = "Cube Side Length [nm]",ylab = "Radius (Env width = 1)",main = "Cubic Volumes")
#points(x1,rs[2,1:5],pch = 1, col="blue",xlab = "Cube Side Length [nm]",ylab = "Radius (Env width = 2)",main = "Cubic Volumes")
plot(x2,rs[1,6:10],pch = 6, col="red",xlab = "Box Height (h) [nm]", ylab = "Radius (Env width = 1)", main = "Non-Cubic Volumes - 60x60xh")
#points(x2,rs[2,6:10],pch = 6, col="red",xlab = "Box Height (h) [nm]", ylab = "Radius (Env width = 2)", main = "Non-Cubic Volumes - 60x60xh")




#### Metric SD vs box size ####
load('~/Research/box_size_effect/...')
load('~/Research/box_size_effect/...')
