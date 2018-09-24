library(rapt)

# Do child envelopes from different parent RCP sets have different behavior?
# Check for radius similarities:
rs <- rep(0,96)
for(i in 1:96){
  a <- read.table(paste('~/Research/point_patterns/Final/system',toString(i),sep=""))
  rs[i]<-as.numeric(levels(a$V1)[2])
}
#plot(rs)
#hist(rs)

neval <- 2000
cor <- "iso"
nrval <- 128
rmax <- 15

#r <- 0.02523824
p <- vector("list",96)
temp_upload <- read.table(paste('~/Research/point_patterns/Final/FinalConfig',toString(1),sep=""),sep=" ",col.names=c("x","y","z","type"))
temp <- scaleRCP(createSpat(temp_upload[,c("x","y","z")]),newRadius = 0.5,oldRadius = rs[1])
base <- pK3est(0.06,temp,neval,correction = cor,nrval=nrval,rmax=rmax,anom=TRUE)
p[[1]] <- base[[1]]

for(i in 2:5){
  t1 <- Sys.time()
  temp_upload <- read.table(paste('~/Research/point_patterns/Final/FinalConfig',toString(i),sep=""),sep=" ",col.names=c("x","y","z","type"))
  temp <- scaleRCP(createSpat(temp_upload[,c("x","y","z")]),newRadius = 0.5,oldRadius = rs[i])
  temp2  <- pK3est(0.06,temp,neval,correction = cor,anom=TRUE,toSub=base[[2]],nrval=nrval,rmax=rmax)
  p[[i]] <- temp2[[1]]
  t2 <- Sys.time()
  print(i)
  print((t2-t1))
  }

#toPlot <- list(p[[3]],p[[4]])

#envPlotAdj(toPlot,percentiles = c(.97))


percentOutside <- rep(0,sum(1:95))
cnt <- 1
for(j in 1:4){
  for(k in (j+1):5){
    ai <- .97
    p1 <- p[[j]]
    p2 <- p[[k]]
    airange <- ai*neval
    aidist <- neval-(round(neval/2)+floor(airange/2)) # select the indexes based on being 1/2 of the percentile span from the middle of the tests
    piBig <- neval-aidist
    piSmall <- 1+aidist
    p1h <- p1[,piBig]
    p1l <- p1[,piSmall]
    p2h <- p2[,piBig]
    p2l <- p2[,piSmall]
    di <- rep(0,nrval)
    d1 <- p1h-p1l
    for(i in 1:nrval){
      if((p2l[i] > p1h[i])|(p1l[i] > p2h[i])){
        di[i] <- d1[i]
      }else if((p2h[i] > p1h[i]) & (p2l[i] < p1l[i])){
        di[i] <- 0
      }else if((p2h[i] < p1h[i]) & (p2l[i] > p1l[i])){
        di[i] <- p1h[i]-p2h[i]+p2l[i]-p1l[i]
      }else if((p2h[i] > p1h[i]) & (p2l[i] > p1l[i])){
        di[i] <- p2l[i]-p1l[i]
      }else if((p2h[i] < p1h[i]) & (p2l[i] < p1l[i])){
        di[i] <- p1h[i]-p2h[i]
      }else if(p2h[i]==0 & p2l[i]==0 & p1h[i]==0 & p1l[i]==0){
        di[i] <- 0
      }else{
        browser()
        print("What went wrong?")
      }
    }
    percentOutside[cnt] <- (sum(di)/sum(d1))*100
    print(cnt)
    cnt <- cnt + 1
  }
}

#####
# FROM THE HPC
sameRCP <- read.table('C:/Users/galen/OneDrive/Documents/Research/parent_child/percentOut_sameRCP_oneradius_20nmbox.csv',sep=",")
sameRCP <- as.numeric(sameRCP$V1)
plot(hist(sameRCP),col=rgb(1,0,0,.5))

diffRCP <- read.table('C:/Users/galen/OneDrive/Documents/Research/parent_child/percentOut_diffRCP_oneradius_20nmbox.csv',sep=",")
diffRCP <- as.numeric(diffRCP$V1)

p1 <- hist(sameRCP)
p2 <- hist(diffRCP)

plot(p1,col=rgb(1,0,0,.5),xlim=c(0,8),main='% Non-Overlap Area of Random Relabeled RCPs',xlab='Percent Non-Overlapping Area',ylab='Counts')
plot(p2,col=rgb(0,0,1,.5),xlim=c(0,8),add=T)
legend("topright",c('Same RCP 96 times','96 different RCPs'),col=c(rgb(1,0,0),rgb(0,0,1)),lwd=10)
text(8,650,'6% random relabeling of',pos=2)
text(8,600,'8000 point RCP patterns',pos=2)
text(8,500, 'Used one radius for all patterns',pos=2)

######################################
envPlotAdj <- function(tests,percentiles=c(.999),ylim=c(-3,3),xlim=c(0,ceiling(max(tests[[1]][,1])))){
  # tests = lists of array of values returned from the rrK3est function above
  # percentiles = vector including the different percentiles you would like to see on the plot (0-1)
  # do these in descending order please

  color <- palette(rainbow(length(tests)))

  # break up data into r values and test results
  rvals <- tests[[1]][,1]
  tvals <- tests[[1]][,2:ncol(tests[[1]])]

  nTests <- ncol(tvals) # number of tests done
  prange <- percentiles*nTests # get the range of indeces for which each percentile spans

  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  percentileIndicesBig <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  percentileIndicesSmall <- round(nTests/2)-floor(prange/2) # do the same for the low end

  # grab out the columns from the sorted test results that we will plot
  toPlotBigs <- matrix(0,nrow=nrow(tvals),ncol=length(percentiles))
  toPlotSmalls <- matrix(0,nrow=nrow(tvals),ncol=length(percentiles))
  for(i in 1:length(percentiles)){
    toPlotBigs[,i] <- sortedtVals[,percentileIndicesBig[i]]
    toPlotSmalls[,i] <- sortedtVals[,percentileIndicesSmall[i]]
  }

  # plot the envelopes from the percentile data
  par(oma = c(0, 2, 0, 0))
  plot(rvals,tvals[,1],type="n",main="Envelopes for K Function",xlab="r",ylab="",ylim=ylim,xlim=xlim)
  mtext(text=expression(sqrt('K'[3]*'(r)')*'  Anomaly'),side=2,line=0,outer=TRUE)
  axis(1,at=0:xlim[2],labels=FALSE)
  axis(1,at=seq(0,xlim[2],by=5))
  polygon(c(rvals,rev(rvals)),c(toPlotBigs[,1],rev(toPlotSmalls[,1])),col=color[1])#,border=color[i],lwd=2)

  abline(h=0,lty=2,lwd=1,col="black")
  #legend(0, ylim[2], legend=c(paste(toString(percentiles[1]*100),"% AI"), paste(toString(percentiles[2]*100),"% AI"),paste(toString(percentiles[3]*100),"% AI")),col=c(color[1],color[2],color[3]), lty=c(1,1,1), lwd=c(10,10,10))

  for(i in 2:length(tests)){
    # break up data into r values and test results
    rvals <- tests[[i]][,1]
    tvals <- tests[[i]][,2:ncol(tests[[i]])]

    nTests <- ncol(tvals) # number of tests done
    prange <- percentiles*nTests # get the range of indeces for which each percentile spans

    sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
    percentileIndicesBig <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
    percentileIndicesSmall <- round(nTests/2)-floor(prange/2) # do the same for the low end

    # grab out the columns from the sorted test results that we will plot
    toPlotBigs <- matrix(0,nrow=nrow(tvals),ncol=length(percentiles))
    toPlotSmalls <- matrix(0,nrow=nrow(tvals),ncol=length(percentiles))

    toPlotBigs[,1] <- sortedtVals[,percentileIndicesBig[1]]
    toPlotSmalls[,1] <- sortedtVals[,percentileIndicesSmall[1]]

    polygon(c(rvals,rev(rvals)),c(toPlotBigs[,1],rev(toPlotSmalls[,1])),col=color[i])#,border=color[i],lwd=2)
  }
  text(1,2.75,paste(toString(percentiles[1]*100),"% AI"))

}
