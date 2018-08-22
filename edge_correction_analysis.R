# Comparing border vs trans & iso edge corrections
library(rapt)

# Upload an RCP pattern
r1 <-0.0252429 #1.2 sigma
rcp_1_upload <- read.table("C:/Users/galen/Documents/Research/point_patterns/Final/FinalConfig1",sep = " ",col.names = c("x","y","z","type"))
rcp_1 <- scale(createSpat(rcp_1_upload[,c("x","y","z")]),newRadius = 0.5,oldRadius = r1)

# Run random relabeling envelopes with border correction on the full sample up to about 5nm radius
maxR <- 5
bordRRL <- panomK3est(0.06,rcp_1,1000,rmax=maxR,correction="bord")

# Run random relabeling envelopes with trans on smaller inner box sample
inner.size <- domain(rcp_1)$xrange[2]-maxR*2
inner.box <- subSquare(rcp_1,rep(inner.size,3))

transRRL <- panomK3est(0.06,inner.box,1000,rmax=maxR,correction="trans",toSub=bordRRL[[2]])

isoRRL <- panomK3est(0.06,inner.box,1000,rmax=maxR,correction="iso",toSub=bordRRL[[2]])

envPlotAdj(list(transRRL[[1]],bordRRL[[1]],isoRRL[[1]]),labels = c("Trans","Bord","Iso"),percentiles = c(.97))

ai <- .97
nTests <- ncol(bordRRL[[1]])-1
irange <- ai*nTests
piBig <- round(nTests/2)+floor(irange/2)
piSmall <- round(nTests/2)-floor(irange/2)
bordh <- bordRRL[[1]][,piBig+1]
bordl <- bordRRL[[1]][,piSmall+1]

nTests <- ncol(transRRL[[1]])-1
irange <- ai*nTests
piBig <- round(nTests/2)+floor(irange/2)
piSmall <- round(nTests/2)-floor(irange/2)
transh <- transRRL[[1]][,piBig+1]
transl <- transRRL[[1]][,piSmall+1]

nTests <- ncol(isoRRL[[1]])-1
irange <- ai*nTests
piBig <- round(nTests/2)+floor(irange/2)
piSmall <- round(nTests/2)-floor(irange/2)
isoh <- isoRRL[[1]][,piBig+1]
isol <- isoRRL[[1]][,piSmall+1]

bwidth <- bordh-bordl
twidth <- transh-transl
iwidth <- isoh-isol

btdiff <- cbind(bordRRL[[1]][,1],(twidth-bwidth))
bidiff <- cbind(bordRRL[[1]][,1],(iwidth-bwidth))

plot(btdiff,type="l",lty=1,lwd=2,col="black",ylim=c(-1,2),main="Envelope Width Difference From Border Method",xlab="r",ylab="Difference")
lines(bidiff,lty=1,lwd=2,col="red")
abline(h=0,lty=2,lwd=1,col="gray")
legend(0,2.0,legend=c("Translation","Isotropic","Border"),lty=c(1,1,2),lwd=c(2,2,2),col=c("black","red","gray"))
text(3,1.9,"Envelopes on 10nm cube, particle radius = 0.5nm")

#####################
# Adjusted envelope plot for this analysis
envPlotAdj <- function(tests,percentiles=c(.999),ylim=c(-3,3),xlim=c(0,ceiling(max(tests[[1]][,1]))),labels=c("Translation","Isotropic","Border")){
  # tests = lists of array of values returned from the rrK3est function above
  # percentiles = vector including the different percentiles you would like to see on the plot (0-1)
  # do these in descending order please

  color <- c("lightskyblue","mediumpurple","lightpink")

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
  axis(1,at=seq(0,xlim[2],by=2))
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
  legend(0,ylim[2],legend=labels,col=c(color[1],color[2],color[3]), lty=c(1,1,1),lwd=c(10,10,10))
  text(1.7,2.75,paste(toString(percentiles[1]*100),"% AI"))

  }

