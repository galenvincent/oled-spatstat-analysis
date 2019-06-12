# Analyze the results from cluster_series_analysis
library(data.table)
library(rapt)
library(RColorBrewer)


# Metric Example ----------------------------------------------------------

# Things to make nice metric example plot
under <- read.rcp('~/Research/point_patterns/Final/FinalConfig1','~/Research/point_patterns/Final/system1',scaleUp = TRUE,newRadius = 0.5)
over <- read.rcp('~/Research/point_patterns/Final/FinalConfig2','~/Research/point_patterns/Final/system2',scaleUp = TRUE,newRadius = 0.5)
under.big <- stitch.size(under, boxSize = c(60,60,60))
over.big <- stitch.size(over, boxSize = c(60,60,60))
toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)
env <- fread('~/Research/K_cluster_series/cube.csv',drop = 1)

cluster <- list()
result <- list()
cluster[[1]] <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=1)
cluster[[2]] <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=3)
cluster[[3]] <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=5)
cluster[[4]] <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",cr=7)

result[[1]] <- anomK3est(cluster[[1]][[1]],toSub,max(env.r),nrow(env.r),correction="trans")
result[[2]] <- anomK3est(cluster[[2]][[1]],toSub,max(env.r),nrow(env.r),correction="trans")
result[[3]] <- anomK3est(cluster[[3]][[1]],toSub,max(env.r),nrow(env.r),correction="trans")
result[[4]] <- anomK3est(cluster[[4]][[1]],toSub,max(env.r),nrow(env.r),correction="trans")


cols <- c('black','red','blue','purple')
plot(result[[1]]$r,result[[1]]$trans, col = "black", lwd = 2, type = "l", ylim = c(-70,70),
     xlab = "r", ylab = expression(sqrt('K'[3]*'(r)')*'  Anomaly'), main = "K vs Cluster Radius",
     cex.main = 1.75, cex.lab = 1.75, cex.axis = 1.25)
for(i in 2:4){
  lines(result[[i]]$r, result[[i]]$trans, col = cols[i], lwd = 2)
}
legend(0,-10,legend = c("R = 1", "R = 3", "R = 5", "R = 7"), col = cols, lwd = rep(2,4), cex = 1.3)

#test
result <- anomK3est(cluster[[2]][[1]],toSub,max(env.r),nrow(env.r),correction="trans")
rvals <- result$r
tvals <- result$trans

par(mgp = c(2.25,1,0),mar = c(3.5,4.25,3,2))
envPlot(env, percentiles = c(0.999,0.99,0.95),xlim = c(0,12),ylim = c(-0.75,0.75),
        leg = FALSE, main = "K-Function Example", cex.main = 1.75, cex.lab = 1.75, cex.axis = 1.5)
lines(rvals,tvals,lwd = 2)
color <- c("lightskyblue","mediumpurple","lightpink")
legend(0, -0.15, legend=c("Observed","Med. Sim.", paste(toString(99.9),"% AI"),
                                   paste(toString(99.0),"% AI"),
                                   paste(toString(95.0),"% AI")),
              col=c("black", "black", color[1],color[2],color[3]), lty=c(1,2,1,1,1), lwd=c(2,2,10,10,10), cex = 1.4)

# get out that peak info son
rvals.new <- rvals[13:length(rvals)]
tvals.new <- tvals[13:length(rvals)]
peak.info <- argmax(rvals.new, tvals.new, w = 3, span = 0.08)

#points(peak.info$x,tvals.new[peak.info$i],pch = 6, cex = 2, col="red")
span <- (peak.info$x[1]/7)*(0.3)
peak.info <- argmax(rvals.new, tvals.new, w = 3, span = span)

peak.info$deriv <- finite_deriv(rvals.new, peak.info$y.hat)
#points(rvals.new,peak.info$deriv)
peak.info$derivsm <- argmax(rvals.new, -1*peak.info$deriv, w = 3, span = span)

#points(peak.info$derivsm$x,peak.info$deriv[peak.info$derivsm$i], col="blue", cex = 2)

peak.info$dderiv <- finite_deriv(rvals.new, -1*peak.info$derivsm$y.hat)
peak.info$dderivsm <- argmax(rvals.new, peak.info$dderiv, w = 3, span = span)
#points(rvals.new,peak.info$dderiv)

#points(peak.info$dderivsm$x,peak.info$dderiv[peak.info$dderivsm$i], col="green", cex = 2, pch = 19)

peak.info$ddderiv <- finite_deriv(rvals.new, peak.info$dderivsm$y.hat)
peak.info$ddderivsm <- argmax(rvals.new, peak.info$ddderiv, w = 3, span = span)
#points(rvals.new, peak.info$ddderiv)

Rddm_ind <- peak.info$ddderivsm$i[peak.info$ddderivsm$i > peak.info$i[1] & peak.info$ddderivsm$i < peak.info$derivsm$i[1]]
#points(rvals.new[Rddm_ind],peak.info$ddderiv[Rddm_ind], col="green", cex = 2, pch = 19)
par(mfrow = c(1,1), mar=c(3.5,4,2,1), mgp=c(2,1,0))
plot(rvals,tvals,type = 'n',xlim = c(1,11),xlab = "r",
     ylab = expression(sqrt('K'[3]*'(r)')*'  Anomaly'),
     cex.lab = 1.75, cex.axis = 1.25)
lines(rvals.new, peak.info$y.hat,col="black",lwd = 2)
lines(rvals.new, -peak.info$derivsm$y.hat, col = "red", lwd = 2)
lines(rvals.new, peak.info$dderivsm$y.hat, col = "purple", lwd = 2)
lines(rvals.new, peak.info$ddderivsm$y.hat, col = "blue", lwd = 2)
legend(8.25,23, c("K","K'","K''","K'''"), col = c("black","red","purple","blue"), lty = 1, lwd = 2, bty = "n", cex = 1.5, y.intersp = 1.35)
  

# Detector Efficiency -----------------------------------------------------

# Detector efficiency analysis
rm(under, over, under.big, over.big, env.r, tvals, rvals, cluster, result, toSub)
#tvals <- rep(list(matrix(0,101,200)),6)
#for(i in 1:101){
#  for(j in 1:6){
#    tvals[[j]][i,] <- fread(paste("~/Research/K_cluster_series/181029_de/result",toString(i),"_",toString(j),".csv",sep=""),select=3)$K
#  }
#}
#rvals <- fread("~/Research/K_cluster_series/de/result1_1.csv",select=1)$r

de <- list()
de$vals <- list()
de$mean <- list()
de$sd <- list()
de$eff <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4)

de$vals$Km <- fread("~/Research/K_cluster_series/181029_de/Km.csv")
de$vals$Rm <- fread("~/Research/K_cluster_series/181029_de/Rm.csv")
#de$vals$Kdm <- fread("~/Research/K_cluster_series/181029_de/Kdm.csv")
de$vals$Rdm <- fread("~/Research/K_cluster_series/181029_de/Rdm.csv")
de$vals$Rddm <- fread("~/Research/K_cluster_series/181029_de/Rddm.csv")

# Lets look at the mean and SDs of the matrices
de$mean <- lapply(de$vals, apply, 2, mean)
de$sd <- lapply(de$vals, apply, 2, sd)

#Plots!
ylabs <- c(expression('K'['max']),
           expression('R'['max']),
           expression('Rd'['min']),
           expression('Rd'['max']^3))
titles <- c("Km vs Detector Efficiency",
            "Rm vs Detector Efficiency",
            "Rdm vs Detector Efficiency",
            "Rddm vs Detector Efficiency")
par(mfrow = c(2,2),mgp = c(2.5,1,0),mar = c(3.5,5,1.5,0.75))
for(i in 1:4){
  plot(de$eff, de$mean[[i]],
       ylim=range(c(de$mean[[i]]-de$sd[[i]], de$mean[[i]]+de$sd[[i]])),
       pch=19, xlab="Detector Efficiency", ylab=ylabs[i],
      cex.main = 1.75, cex.lab = 1.75, cex.axis = 1.25
  )
  # hack: we draw arrows but with very special "arrowheads"
  arrows(de$eff[de$sd[[i]] != 0], de$mean[[i]][de$sd[[i]] != 0]-de$sd[[i]][de$sd[[i]] != 0], de$eff[de$sd[[i]] != 0], 
         de$mean[[i]][de$sd[[i]] != 0]+de$sd[[i]][de$sd[[i]] != 0], length=0.05, angle=90, code=3)
}


# Cluster Radius ----------------------------------------------------------

# radius analysis
rm(list=ls())
gc()

de <- list()
de$vals <- list()
de$mean <- list()
de$sd <- list()
de$eff <- c(1,2,3,4,5,8)

de$vals$Km <- fread("~/Research/K_cluster_series/190211_rad/Km.csv")
de$vals$Rm <- fread("~/Research/K_cluster_series/190211_rad/Rm.csv")
#de$vals$Kdm <- fread("~/Research/K_cluster_series/190211_rad/Kdm.csv")
de$vals$Rdm <- fread("~/Research/K_cluster_series/190211_rad/Rdm.csv")
de$vals$Rddm <- fread("~/Research/K_cluster_series/190211_rad/Rddm.csv")

# Lets look at the mean and SDs of the matrices
de$mean <- lapply(de$vals, apply, 2, mean, na.rm=TRUE)
de$sd <- lapply(de$vals, apply, 2, sd, na.rm = TRUE)

#Upload the envelopes to see if we fall outside of them
env <- fread("~/Research/K_cluster_series/cube.csv",drop=1)
percentile <- .999

rvals <- env[,1]
tvals <- env[,2:ncol(env)]
nTests <- ncol(tvals) # number of tests done
prange <- percentile*nTests # get the range of indeces for which each percentile spans

sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
plot.big <- sortedtVals[,ind.big]

env.val <- vector("numeric",length(de$eff))
for(i in (1:length(de$eff))){
  rind <- which(abs(rvals-de$eff[i]) == min(abs(rvals-de$eff[i])))
  env.val[i] <- plot.big[rind]
}

#Plots!
ylabs <- c(expression('K'['max']),
           expression('R'['max']),
           expression('Rd'['min']),
           expression('Rd'['max']^3))
titles <- c("Km vs Cluster Radius",
            "Rm vs Cluster Radius",
            "Rdm vs Cluster Radius",
            "Rddm vs Cluster Radius")

#fit some lines
# Full linear models 
a <- lm(de$mean$Rdm ~ de$eff)
summary(a)
b <- lm(de$mean$Rm ~ de$eff)
summary(b)
c <- lm(de$mean$Rddm ~ de$eff)
summary(c)

#no intercept models
an <- lm(de$mean$Rdm ~ de$eff -1)
summary(an)
bn <- lm(de$mean$Rm ~ de$eff -1)
summary(bn)
cn <- lm(de$mean$Rddm ~ de$eff -1)
summary(cn)

#info criterion
AIC(a)
AIC(an)
BIC(a)
BIC(an)

AIC(b)
AIC(bn)
BIC(b)
BIC(bn)

AIC(c)
AIC(cn)
BIC(c)
BIC(cn)

par(mfrow = c(2,2),mgp = c(2.5,1,0),mar = c(3.5,5,1.5,0.75))
for(i in 1:4){
  plot(de$eff, de$mean[[i]],
       ylim=range(c(de$mean[[i]]-de$sd[[i]], de$mean[[i]]+de$sd[[i]]),0,na.rm=TRUE),
       pch=19, xlab="Cluster Radius", ylab=ylabs[i],
       cex.lab = 1.75, cex.axis= 1.25
  )
  # hack: we draw arrows but with very special "arrowheads"
  arrows(de$eff[de$sd[[i]] != 0], de$mean[[i]][de$sd[[i]] != 0]-de$sd[[i]][de$sd[[i]] != 0], de$eff[de$sd[[i]] != 0], 
         de$mean[[i]][de$sd[[i]] != 0]+de$sd[[i]][de$sd[[i]] != 0], length=0.05, angle=90, code=3)
  
  if(i == 2){
    abline(bn,lty=2, col= "red")
    text(1,10,paste("Rm =",toString(round(b$coefficients[1],3)),"+",toString(round(b$coefficients[2],3)),"r"),pos=4, cex = 1.25)
  }else if(i == 3){
    abline(an,lty = 2,col="red")
    text(1,17,paste("Rdm =",toString(round(a$coefficients[1],3)),"+",toString(round(a$coefficients[2],3)),"r"),pos=4, cex = 1.25)
  }else if (i==1){
    lines(de$eff, env.val, col = "blue",lwd = 2)
    text(8.5,9,"99.9% AI RRL Envelope width",pos=2, cex = 1.1)
  }else if(i == 4){
    abline(cn,lty = 2,col="red")
    text(1,14,paste("Rddm =",toString(round(c$coefficients[1],3)),"+",toString(round(c$coefficients[2],3)),"r"),pos=4, cex = 1.25)
  }

}


# Cluster Density ---------------------------------------------------------

# cluster density series
rm(list=ls())
#tvals <- rep(list(matrix(0,101,200)),8)
#for(i in 1:101){
#  for(j in 1:8){
#    tvals[[j]][i,] <- fread(paste("~/Research/K_cluster_series/181022_den/result",toString(i),"_",toString(j),".csv",sep=""),select=3)$K
#  }
#}

rvals <- fread("~/Research/K_cluster_series/181022_den/result1_1.csv",select=1)$r

de <- list()
de$vals <- list()
de$mean <- list()
de$sd <- list()
de$eff <- c(1, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15)

de$vals$Km <- fread("~/Research/K_cluster_series/181029_den/Km.csv")
de$vals$Rm <- fread("~/Research/K_cluster_series/181029_den/Rm.csv")
#de$vals$Kdm <- fread("~/Research/K_cluster_series/181029_den/Kdm.csv")
de$vals$Rdm <- fread("~/Research/K_cluster_series/181029_den/Rdm.csv")
de$vals$Rddm <- fread("~/Research/K_cluster_series/181029_den/Rddm.csv")

# Lets look at the mean and SDs of the matrices
# Get rid of ridiculous radius values
de$vals$Rm[de$vals$Rm > 5] <- NaN

de$mean <- lapply(de$vals, apply, 2, mean, na.rm=TRUE)
de$sd <- lapply(de$vals, apply, 2, sd, na.rm = TRUE)

#Upload the envelopes to see if we fall outside of them
env <- fread("~/Research/K_cluster_series/cube.csv",drop=1)
percentile <- .999

rvals <- env[,1]
tvals <- env[,2:ncol(env)]
nTests <- ncol(tvals) # number of tests done
prange <- percentile*nTests # get the range of indeces for which each percentile spans

sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
plot.big <- sortedtVals[,ind.big]

env.val <- vector("numeric",length(de$eff))
for(i in (1:length(de$eff))){
  rind <- which(abs(rvals-2) == min(abs(rvals-2)))
  env.val[i] <- plot.big[rind]
}

#Plots!
ylabs <- c(expression('K'['max']),
           expression('R'['max']),
           expression('Rd'['min']),
           expression('Rd'['max']^3))
titles <- c("Km vs Cluster Density",
            "Rm vs Cluster Density",
            "Rdm vs Cluster Density",
            "Rddm vs Cluster Density")
par(mfrow = c(2,2),mgp = c(2.5,1,0), mar = c(3.5,5,2,0.75))
for(i in 1:4){
  plot(de$eff[1:7], de$mean[[i]][1:7],
       ylim=range(c(de$mean[[i]][1:7]-de$sd[[i]][1:7], de$mean[[i]][1:7]+de$sd[[i]][1:7])),
       pch=19, xlab="Cluster Density", ylab=ylabs[i],
       cex.lab = 1.75, cex.axis = 1.25
  )
  # hack: we draw arrows but with very special "arrowheads"
  arrows(de$eff[de$sd[[i]] != 0], de$mean[[i]][de$sd[[i]] != 0]-de$sd[[i]][de$sd[[i]] != 0], de$eff[de$sd[[i]] != 0], 
         de$mean[[i]][de$sd[[i]] != 0]+de$sd[[i]][de$sd[[i]] != 0], length=0.05, angle=90, code=3)
 if (i==1){
  lines(de$eff[1:7], env.val[1:7], col = "blue",lwd = 2)
  text(1.05,1,"99.9% AI RRL Envelope width",pos=2, cex = 1.1)
}
}


# Den + Rad ---------------------------------------------------------------

# cluster density + radius series
rm(list=ls())

de <- list()
de$vals <- list()
de$mean <- list()
de$sd <- list()
de$den <- c(1, 0.6, 0.4, 0.3, 0.2, 0.15)
de$r <- c(2, 3 ,4, 5, 6, 7)


de$vals$Km <- fread("~/Research/K_cluster_series/181105_denr/Km.csv")
de$vals$Rm <- fread("~/Research/K_cluster_series/181105_denr/Rm.csv")
de$vals$Rdm <- fread("~/Research/K_cluster_series/181105_denr/Rdm.csv")
de$vals$Rddm <- fread("~/Research/K_cluster_series/181105_denr/Rddm.csv")
de$vals$csep <- fread("~/Research/K_cluster_series/181105_denr/csep.csv")

for(i in 1:5){
  de$mean <- lapply(de$vals, apply, 2, mean, na.rm=TRUE)
  de$sd <- lapply(de$vals, apply, 2, sd, na.rm = TRUE)
}


de$mean$mat <- list()
de$sd$mat <- list()
for(i in 1:5){
  de$mean$mat[[i]] <- matrix(de$mean[[i]],6,6)
  de$sd$mat[[i]] <- matrix(de$sd[[i]],6,6)
  
  }
names(de$mean$mat) <- c("Km", "Rm", "Rdm", "Rddm", "csep")
names(de$sd$mat) <- c("Km", "Rm", "Rdm", "Rddm", "csep")

#Upload the envelopes to see if we fall outside of them
env <- fread("~/Research/K_cluster_series/cube.csv",drop=1)
percentile <- .999

rvals <- env[,1]
tvals <- env[,2:ncol(env)]
nTests <- ncol(tvals) # number of tests done
prange <- percentile*nTests # get the range of indeces for which each percentile spans

sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
plot.big <- sortedtVals[,ind.big]

env.val <- vector("numeric",length(de$r))
for(i in (1:length(de$r))){
  rind <- which(abs(rvals-2) == min(abs(rvals-2)))
  env.val[i] <- plot.big[rind]
}

#Plots!
ylabs <- c(expression('K'['max']),
           expression('R'['max']),
           expression('Rd'['min']),
           expression('Rd'['max']^3))
titles <- c("Km",
            "Rm",
            "Rdm",
            "Rddm")
color = c("red","darkorange","darkgreen","blue","purple","black")
par(mfrow = c(2,2),mgp = c(2.5,1,0),mar = c(3.5,5,2,0.75))
for(i in 1:4){
  plot(de$r, de$mean$mat$Km[,1],
       ylim=range(c(min(de$mean$mat[[i]]-de$sd$mat[[i]]), max(de$mean$mat[[i]]+de$sd$mat[[i]]))),
       pch=19, xlab="Cluster Radius", ylab=ylabs[i],type = "n",
       cex.lab = 1.75, cex.axis = 1.25)
       for(j in c(1,2,3,4,5,6)){
         points(de$r, de$mean$mat[[i]][j,],pch= 19,col=color[j])
         lines(de$r, de$mean$mat[[i]][j,],col=color[j])
         arrows(de$r[de$sd$mat[[i]][j,] != 0], de$mean$mat[[i]][j,][de$sd$mat[[i]][j,] != 0]-de$sd$mat[[i]][j,][de$sd$mat[[i]][j,] != 0], 
                de$r[de$sd$mat[[i]][j,] != 0], de$mean$mat[[i]][j,][de$sd$mat[[i]][j,] != 0]+de$sd$mat[[i]][j,][de$sd$mat[[i]][j,] != 0], 
                length=0.05, angle=90, code=3,col=color[j])
       }
  if (i==1){
    lines(de$r, env.val, col = "blue",lwd = 2)
    #text(1,1,"99.9% AI RRL Envelope width",pos=2)
    c(1, 0.6, 0.4, 0.3, 0.2, 0.15)
    legend(2,77,legend = c("1","0.6","0.4","0.3","0.2","0.15"), col = color, pch = 19, lwd = 1.3,bty = "n",cex = 1.1,y.intersp= 0.9)
    text(1.8,77,"Cluster Density",pos = 4,cex = 1.1)
  }
}

par(mfrow = c(1,1),mgp = c(2.5,1,0),mar = c(3.5,5,2,0.75))
i <- 1
plot(de$r, de$mean$mat$Km[,1],
     ylim=range(c(min(de$mean$mat[[i]]-de$sd$mat[[i]]), max(de$mean$mat[[i]]+de$sd$mat[[i]]))),
     pch=19, xlab="Cluster Radius", ylab=ylabs[i],type = "n",log = "y",
     cex.lab = 1.75, cex.axis = 1.25)
lines(de$r, env.val, col = "blue",lwd = 2,lty = 2)
text(7,0.11,"99.9% AI RRL Envelope width",pos=2, cex = 1.25)
for(j in c(1,2,3,4,5,6)){
  points(de$r, de$mean$mat[[i]][j,],pch= 19,col=color[j])
  lines(de$r, de$mean$mat[[i]][j,],col=color[j])
  arrows(de$r[de$sd$mat[[i]][j,] != 0], de$mean$mat[[i]][j,][de$sd$mat[[i]][j,] != 0]-de$sd$mat[[i]][j,][de$sd$mat[[i]][j,] != 0], 
         de$r[de$sd$mat[[i]][j,] != 0], de$mean$mat[[i]][j,][de$sd$mat[[i]][j,] != 0]+de$sd$mat[[i]][j,][de$sd$mat[[i]][j,] != 0], 
         length=0.05, angle=90, code=3,col=color[j])
}


# Position Blur -----------------------------------------------------------

# gaussian blur series
rm(list=ls())

de <- list()
de$vals <- list()
de$mean <- list()
de$sd <- list()
de$eff <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5 ,0.6)

de$vals$Km <- fread("~/Research/K_cluster_series/190114_gblur_method2/Km.csv")
de$vals$Rm <- fread("~/Research/K_cluster_series/190114_gblur_method2/Rm.csv")
#de$vals$Kdm <- fread("~/Research/K_cluster_series/181112_gblur/Kdm.csv")
de$vals$Rdm <- fread("~/Research/K_cluster_series/190114_gblur_method2/Rdm.csv")
de$vals$Rddm <- fread("~/Research/K_cluster_series/190114_gblur_method2/Rddm.csv")

# Lets look at the mean and SDs of the matrices
# Get rid of ridiculous radius values
#de$vals$Rm[de$vals$Rm > 5] <- NaN

de$mean <- lapply(de$vals, apply, 2, mean, na.rm=TRUE)
de$sd <- lapply(de$vals, apply, 2, sd, na.rm = TRUE)

#Upload the envelopes to see if we fall outside of them
env <- fread("~/Research/K_cluster_series/cube.csv",drop=1)
percentile <- .999

rvals <- env[,1]
tvals <- env[,2:ncol(env)]
nTests <- ncol(tvals) # number of tests done
prange <- percentile*nTests # get the range of indeces for which each percentile spans

sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
plot.big <- sortedtVals[,ind.big]

env.val <- vector("numeric",length(de$eff))
for(i in (1:length(de$eff))){
  rind <- which(abs(rvals-2) == min(abs(rvals-2)))
  env.val[i] <- plot.big[rind]
}

#Plots!
ylabs <- c(expression('K'['max']),
           expression('R'['max']),
           expression('Rd'['min']),
           expression('Rd'['max']^3))
titles <- c("Km vs Cluster Position Blur",
            "Rm vs Cluster Position Blur",
            "Rdm vs Cluster Position Blur",
            "Rddm vs Cluster Position Blur")
par(mfrow = c(2,2),mgp = c(2.5,1,0),mar = c(3.5,5,1.5,0.75))
for(i in 1:4){
  plot(de$eff, de$mean[[i]],
       ylim=range(c(de$mean[[i]]-de$sd[[i]], de$mean[[i]]+de$sd[[i]]),0),
       pch=19, xlab="Blur SD % of Separation", ylab=ylabs[i],
       cex.lab = 1.75, cex.axis = 1.25
  )
  # hack: we draw arrows but with very special "arrowheads"
  arrows(de$eff[de$sd[[i]] != 0], de$mean[[i]][de$sd[[i]] != 0]-de$sd[[i]][de$sd[[i]] != 0], de$eff[de$sd[[i]] != 0], 
         de$mean[[i]][de$sd[[i]] != 0]+de$sd[[i]][de$sd[[i]] != 0], length=0.05, angle=90, code=3)
  if (i==1){
    lines(de$eff, env.val, col = "blue",lwd = 2)
    text(0.5,1.75,"99.9% AI RRL Envelope width",pos=2, cex = 1.1)
    text(-0.025,18,"Cluster R = 3",pos=4, cex = 1.25)
  }
}

#check out the cluster separation data
graphics.off()

csep1 <- list()
csep1$means <- vector("numeric", 8)
csep1$sd <- vector("numeric", 8)
csep2 <- list()
csep2$means <- vector("numeric", 8)
csep2$sd <- vector("numeric", 8)
for(i in 1:8){
  csep1[[i]] <- fread(paste("~/Research/K_cluster_series/190114_gblur_method1/csep",toString(i),".csv",sep = ""))
  csep2[[i]] <- fread(paste("~/Research/K_cluster_series/190114_gblur_method2/csep",toString(i),".csv",sep = ""))
  
  csep1$means[i] <- mean(as.matrix(csep1[[i]]), na.rm = TRUE)
  csep1$sd[i] <- sd(as.matrix(csep1[[i]]), na.rm = TRUE)
  csep2$means[i] <- mean(as.matrix(csep2[[i]]), na.rm = TRUE)
  csep2$sd[i] <- sd(as.matrix(csep2[[i]]), na.rm = TRUE)
  
  par(mfrow = c(1,2))
  hist(as.matrix(csep1[[i]]),xlim = c(6,20),ylim = c(0,8000),breaks = 0:50,main = paste("Blur % of Separation: ",toString(de$eff[i]),sep = ""),xlab = "Method 1")
  hist(as.matrix(csep2[[i]]),xlim = c(6,20),ylim = c(0,8000),breaks = 0:50,main = paste("Blur % of Separation: ",toString(de$eff[i]),sep = ""),xlab = "Method 2")
}



# Radius Blur -------------------------------------------------------------

# radius blur series
rm(list=ls())
gc()

de <- list()
de$vals <- list()
de$mean <- list()
de$sd <- list()
de$eff <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6)

cr <- 3
de$rw <- (cr^4 + 6*cr^2*(de$eff*cr)^2  + 3*(de$eff*cr)^4)/(cr^3 + 3*cr*(de$eff*cr)^2) 

de$vals$Km <- fread("~/Research/K_cluster_series/190116_rblur/Km.csv")
de$vals$Rm <- fread("~/Research/K_cluster_series/190116_rblur/Rm.csv")
#de$vals$Kdm <- fread("~/Research/K_cluster_series/190116_rblur/Kdm.csv")
de$vals$Rdm <- fread("~/Research/K_cluster_series/190116_rblur/Rdm.csv")
de$vals$Rddm <- fread("~/Research/K_cluster_series/190116_rblur/Rddm.csv")

# Lets look at the mean and SDs of the matrices
# Get rid of ridiculous radius values
#de$vals$Rm[de$vals$Rm > 5] <- NaN

de$mean <- lapply(de$vals, apply, 2, mean, na.rm=TRUE)
de$sd <- lapply(de$vals, apply, 2, sd, na.rm = TRUE)

#Upload the envelopes to see if we fall outside of them
env <- fread("~/Research/K_cluster_series/cube.csv",drop=1)
percentile <- .999

rvals <- env[,1]
tvals <- env[,2:ncol(env)]
nTests <- ncol(tvals) # number of tests done
prange <- percentile*nTests # get the range of indeces for which each percentile spans

sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
plot.big <- sortedtVals[,ind.big]

env.val <- vector("numeric",length(de$eff))
for(i in (1:length(de$eff))){
  rind <- which(abs(rvals-3) == min(abs(rvals-3)))
  env.val[i] <- plot.big[rind]
}

#Plots!
ylabs <- c(expression('K'['max']),
           expression('R'['max']),
           expression('Rd'['min']),
           expression('Rd'['max']^3))
titles <- c("Km vs Radius Blur",
            "Rm vs Radius Blur",
            "Rdm vs Radius Blur",
            "Rddm vs Radius Blur")
par(mfrow = c(2,2),mgp = c(2.5,1,0),mar = c(3.5,5,1.5,0.75))
for(i in 1:4){
  plot(de$rw, de$mean[[i]],
       ylim=range(c(de$mean[[i]]-de$sd[[i]], de$mean[[i]]+de$sd[[i]])),
       pch=19, xlab="Weighted Radius", ylab=ylabs[i],
       cex.lab = 1.75, cex.axis = 1.25
  )
  # hack: we draw arrows but with very special "arrowheads"
  arrows(de$rw[de$sd[[i]] != 0], de$mean[[i]][de$sd[[i]] != 0]-de$sd[[i]][de$sd[[i]] != 0], de$rw[de$sd[[i]] != 0], 
         de$mean[[i]][de$sd[[i]] != 0]+de$sd[[i]][de$sd[[i]] != 0], length=0.05, angle=90, code=3)
  if (i==1){
    lines(de$rw, env.val, col = "blue",lwd = 2)
    text(0.5,1.75,"99.9% AI RRL Envelope width",pos=2, cex = 1.1)
    text(-0.025,18,"Cluster R = 3",pos=4, cex = 1.25)
  }
}

# Binomial Radius Blur ----------------------------------------------------
rm(list=ls())
gc()

de <- list()
de$vals <- list()
de$mean <- list()
de$sd <- list()
de$eff <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

# Calculate average volume for each value above
r1 <- 2
r2 <- 5
de$vol <- ((r1^3)*de$eff + (r2^3)*(1-de$eff))^(1/3)

# Calculate the volume percent of the small size for each percentage
de$perc <- de$eff*r1^3/(de$eff*r1^3 + (1-de$eff)*r2^3)

#weird k metric
de$new <- (r1*r1^3*de$eff + r2*r2^3*(1-de$eff))/(r1^3*de$eff + r2^3*(1-de$eff))


de$vals$Km <- fread("~/Research/K_cluster_series/190130_rbbinom/Km.csv")
de$vals$Rm <- fread("~/Research/K_cluster_series/190130_rbbinom/Rm.csv")
#de$vals$Kdm <- fread("~/Research/K_cluster_series/190124_rbbinom/Kdm.csv")
de$vals$Rdm <- fread("~/Research/K_cluster_series/190130_rbbinom/Rdm.csv")
de$vals$Rddm <- fread("~/Research/K_cluster_series/190130_rbbinom/Rddm.csv")

# Lets look at the mean and SDs of the matrices
# Get rid of ridiculous radius values
#de$vals$Rm[de$vals$Rm > 5] <- NaN

de$mean <- lapply(de$vals, apply, 2, mean, na.rm=TRUE)
de$sd <- lapply(de$vals, apply, 2, sd, na.rm = TRUE)

#Upload the envelopes to see if we fall outside of them
env <- fread("~/Research/K_cluster_series/cube.csv",drop=1)
percentile <- .999

rvals <- env[,1]
tvals <- env[,2:ncol(env)]
nTests <- ncol(tvals) # number of tests done
prange <- percentile*nTests # get the range of indeces for which each percentile spans

sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
plot.big <- sortedtVals[,ind.big]

env.val <- vector("numeric",length(de$eff))
for(i in (1:length(de$eff))){
  rind <- which(abs(rvals-de$vol[i]) == min(abs(rvals-de$vol[i])))
  env.val[i] <- plot.big[rind]
}

#Plots!
ylabs <- c(expression('K'['max']),
           expression('R'['max']),
           expression('Rd'['min']),
           expression('Rd'['max']^3))

par(mfrow = c(2,2),mgp = c(2.5,1,0),mar = c(3.5,5,1.5,0.75))
for(i in 1:4){
  plot(de$new, de$mean[[i]],
       ylim=range(c(de$mean[[i]]-de$sd[[i]], de$mean[[i]]+de$sd[[i]])),
       pch=19, xlab="Fraction R1", ylab=ylabs[i],
       cex.lab = 1.75, cex.axis = 1.25
  )
  # hack: we draw arrows but with very special "arrowheads"
  arrows(de$new[de$sd[[i]] != 0], de$mean[[i]][de$sd[[i]] != 0]-de$sd[[i]][de$sd[[i]] != 0], de$new[de$sd[[i]] != 0], 
         de$mean[[i]][de$sd[[i]] != 0]+de$sd[[i]][de$sd[[i]] != 0], length=0.05, angle=90, code=3)
  if (i==1){
    lines(de$new, env.val, col = "blue",lwd = 2)
    text(0.5,1.75,"99.9% AI RRL Envelope width",pos=2, cex = 1.1)
    text(-0.025,18,"Cluster R1 = 2, R2 = 5",pos=4, cex = 1.25)
  }
}

a <- lm(de$mean[[1]] ~ de$perc)
summary(a)

graphics.off()
plot(de$mean[[1]] ~ de$perc)
abline(a)

# Density + Radius + Radius Blur ------------------------------------------

#denrb
rm(list=ls())
gc()

de <- list()
de$vals <- list()
de$mean <- list()
de$sd <- list()

de$den <- c(1, 0.4, 0.2, 0.1)
de$r <- c(2, 3, 4, 5)
de$rbperc <- c(0, 0.2, 0.4, 0.6)

avgvol.r <- matrix(0,4,4)
# i = radius, j = blur amount
for(i in 1:4){
  for(j in 1:4){
    #avgvol[i,j] <- (4/3)*pi*de$r[i]*(de$r[i]^2 + 3*(de$rbperc[j]*de$r[i])^2)
    avgvol.r[i,j] <- (de$r[i]^4 + 6*(de$r[i]^2)*(de$rbperc[j]*de$r[i])^2 + 3*(de$rbperc[j]*de$r[i])^4)/(de$r[i]*(de$r[i]^2 + 3*(de$rbperc[j]*de$r[i])^2))
  }
}

de$vals$Km <- fread("~/Research/K_cluster_series/190206_denrb/Km.csv")
de$vals$Rm <- fread("~/Research/K_cluster_series/190206_denrb/Rm.csv")
de$vals$Rdm <- fread("~/Research/K_cluster_series/190206_denrb/Rdm.csv")
de$vals$Rddm <- fread("~/Research/K_cluster_series/190206_denrb/Rddm.csv")

#de$vals$csep <- fread("~/Research/K_cluster_series/181105_denr/csep.csv")


for(i in 1:4){
  de$mean <- lapply(de$vals, apply, 2, mean, na.rm=TRUE)
  de$sd <- lapply(de$vals, apply, 2, sd, na.rm = TRUE)
}
#names(de$mean) <- c("Km", "Rm", "Rdm", "Rddm", "csep")
#names(de$sd) <- c("Km", "Rm", "Rdm", "Rddm", "csep")

de$mean$mat <- list()
de$sd$mat <- list()
seps <- matrix(c(1,16,17,32,33,48,49,64),nrow = 4, ncol = 2, byrow = TRUE)
for(i in 1:4){
  de$mean$mat[[i]] <- list() 
  de$sd$mat[[i]] <- list()
  for(j in 1:4){
    de$mean$mat[[i]][[j]] <- matrix(de$mean[[i]][seps[j,1]:seps[j,2]],4,4)
    de$sd$mat[[i]][[j]] <- matrix(de$sd[[i]][seps[j,1]:seps[j,2]],4,4)
  }
}
#NOTE: j indicies run through rblur changes, i indices run through metrics,
# within the matrix, each column represents a different radius while each row represents a different density

names(de$mean$mat) <- c("Km", "Rm", "Rdm", "Rddm")
names(de$sd$mat) <- c("Km", "Rm", "Rdm", "Rddm")

#Upload the envelopes to see if we fall outside of them
env <- fread("~/Research/K_cluster_series/cube.csv",drop=1)
percentile <- .999

rvals <- env[,1]
tvals <- env[,2:ncol(env)]
nTests <- ncol(tvals) # number of tests done
prange <- percentile*nTests # get the range of indeces for which each percentile spans

sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
plot.big <- sortedtVals[,ind.big]

env.val <- vector("numeric",length(de$r))
for(i in (1:length(de$r))){
  rind <- which(abs(rvals-de$r[i]) == min(abs(rvals-de$r[i])))
  env.val[i] <- plot.big[rind]
}

#Plots!
ylabs <- c(expression('K'['max']),
           expression('R'['max']),
           expression('Rd'['min']),
           expression('Rd'['max']^3))
rblurvals <- c("0",
               "0.2",
               "0.4",
               "0.6")
color = c("red","darkorange","darkgreen","blue")
par(mfrow = c(2,2),mgp = c(2.5,1,0),mar = c(3.5,5,2,0.75))
for(k in 1:4){
  for(i in 1:4){
    plot(de$r, de$mean$mat$Km[[1]][1,],
         ylim=range(c(min(unlist(de$mean$mat[[i]])-unlist(de$sd$mat[[i]])), max(unlist(de$mean$mat[[i]])+unlist(de$sd$mat[[i]])))),
         pch=19, xlab="Cluster Radius", ylab=ylabs[i],type = "n",
         cex.lab = 1.75, cex.axis = 1.25)
    for(j in c(1,2,3,4)){
      points(de$r, de$mean$mat[[i]][[k]][j,],pch= 19,col=color[j])
      lines(de$r, de$mean$mat[[i]][[k]][j,],col=color[j])
      arrows(de$r[de$sd$mat[[i]][[k]][j,] != 0], de$mean$mat[[i]][[k]][j,][de$sd$mat[[i]][[k]][j,] != 0]-de$sd$mat[[i]][[k]][j,][de$sd$mat[[i]][[k]][j,] != 0], 
             de$r[de$sd$mat[[i]][[k]][j,] != 0], de$mean$mat[[i]][[k]][j,][de$sd$mat[[i]][[k]][j,] != 0]+de$sd$mat[[i]][[k]][j,][de$sd$mat[[i]][[k]][j,] != 0], 
             length=0.05, angle=90, code=3,col=color[j])
    }
    if (i==1){
      lines(de$r, env.val, col = "black",lwd = 2)
      #text(1,1,"99.9% AI RRL Envelope width",pos=2)
      legend(2,60,legend = c("1","0.4","0.2","0.15"), col = color, pch = 19, lwd = 1.3,bty = "n",cex = 1.1,y.intersp= 0.9)
      text(1.8,60.5,"Cluster Density",pos = 4,cex = 1.1)
    }else if(i == 4){
      mtext(paste("Radius Blur SD Fraction: ", rblurvals[[k]]), side = 3, outer = TRUE, line = -1.5, cex = 1.25)
    }
  }
}

# Average each of the radius blurs and plot on top of one another
de$mean$avgs <- list()
de$sd$avgs <- list()
for(i in 1:4){
  de$mean$avgs[[i]] <- lapply(de$mean$mat[[i]], apply, 2, mean)
  de$sd$avgs <- lapply(de$sd$mat[[i]], apply, 2, function(x){sqrt(x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2)})
}

#plot them 
par(mfrow = c(2,2),mgp = c(2.5,1,0),mar = c(3.5,5,2,0.75))
for(i in 1:4){
  plot(de$r, de$mean$avgs[[1]][[1]],
       ylim=range(c(min(unlist(de$mean$avgs[[i]])-unlist(de$sd$avgs[[i]])), max(unlist(de$mean$avgs[[i]])+unlist(de$sd$avgs[[i]])))),
       pch=19, xlab="Cluster Radius", ylab=ylabs[i],type = "n",
       cex.lab = 1.75, cex.axis = 1.25)
  if(i == 1){
    legend(2.75,17,legend = c("0","0.2","0.4","0.6"), col = color, pch = 19, lwd = 1.3,bty = "n",cex = 1.4,y.intersp= 0.9)
    text(2,18,"Radius Blur SD Fraction",pos = 4,cex = 1.4)
    next
  }
  for(j in c(1,2,3,4)){
    points(de$r, de$mean$avgs[[i]][[j]],pch= 19,col=color[j])
    lines(de$r, de$mean$avgs[[i]][[j]],col=color[j])
    # arrows(de$r[de$sd$avgs[[i]][[j]] != 0], de$mean$avgs[[i]][[j]][de$sd$avgs[[i]][[j]] != 0]-de$sd$avgs[[i]][[j]][de$sd$avgs[[i]][[j]] != 0],
    #        de$r[de$sd$avgs[[i]][[j]] != 0], de$mean$avgs[[i]][[j]][de$sd$avgs[[i]][[j]] != 0]+de$sd$avgs[[i]][[j]][de$sd$avgs[[i]][[j]] != 0],
    #        length=0.05, angle=90, code=3,col=color[j])
  }
  if(i == 4){
    mtext("Averaged Over Cluster Density", side = 3, outer = TRUE, line = -1.5, cex = 1.25)
  }
}



# Plot by WEIGHTED RADIUS
x <- rep(0,4^3)
y <- rep(0,4^3)
cnt <- 1
color = c("red","darkorange","darkgreen","blue")
par(mfrow = c(2,2),mgp = c(2.5,1,0),mar = c(3.5,5,2,0.75))
for(k in 1:4){
  for(i in 1:4){
    plot(diag(avgvol.r), de$mean$mat$Km[[1]][1,],
         ylim=range(c(min(unlist(de$mean$mat[[i]])-unlist(de$sd$mat[[i]])), max(unlist(de$mean$mat[[i]])+unlist(de$sd$mat[[i]])))),
         pch=19, xlab="Weighted Radius", ylab=ylabs[i],type = "n",
         cex.lab = 1.75, cex.axis = 1.25)
    
    for(j in c(1,2,3,4)){
      if(i == 2){
        x[cnt:(cnt+3)] <- avgvol.r[,k]
        y[cnt:(cnt+3)] <- de$mean$mat[[i]][[k]][j,]
        cnt <- cnt + 4
      }
      points(avgvol.r[,k], de$mean$mat[[i]][[k]][j,],pch= 19,col=color[j])
      lines(avgvol.r[,k], de$mean$mat[[i]][[k]][j,],col=color[j])
      arrows(avgvol.r[,k][de$sd$mat[[i]][[k]][j,] != 0], de$mean$mat[[i]][[k]][j,][de$sd$mat[[i]][[k]][j,] != 0]-de$sd$mat[[i]][[k]][j,][de$sd$mat[[i]][[k]][j,] != 0], 
             avgvol.r[,k][de$sd$mat[[i]][[k]][j,] != 0], de$mean$mat[[i]][[k]][j,][de$sd$mat[[i]][[k]][j,] != 0]+de$sd$mat[[i]][[k]][j,][de$sd$mat[[i]][[k]][j,] != 0], 
             length=0.05, angle=90, code=3,col=color[j])
    }
    if (i==1){
      lines(diag(avgvol.r), env.val, col = "black",lwd = 2)
      #text(1,1,"99.9% AI RRL Envelope width",pos=2)
      legend(2,60,legend = c("1","0.4","0.2","0.15"), col = color, pch = 19, lwd = 1.3,bty = "n",cex = 1.1,y.intersp= 0.9)
      text(1.8,60.5,"Cluster Density",pos = 4,cex = 1.1)
    }else if(i == 4){
      mtext(paste("Radius Blur SD Fraction: ", rblurvals[[k]]), side = 3, outer = TRUE, line = -1.5, cex = 1.25)
    }
  }
}
# Fit linear relationship in Rm:
rmlm <- lm(y~x)

#Plot averaged across density vs average radius
par(mfrow = c(2,2),mgp = c(2.5,1,0),mar = c(3.5,5,2,0.75))
for(i in 1:4){
  plot(diag(avgvol.r), de$mean$avgs[[1]][[1]],
       ylim=range(c(min(unlist(de$mean$avgs[[i]])-unlist(de$sd$avgs[[i]])), max(unlist(de$mean$avgs[[i]])+unlist(de$sd$avgs[[i]])))),
       pch=19, xlab="Weighted Radius", ylab=ylabs[i],type = "n",
       cex.lab = 1.75, cex.axis = 1.25)
  if(i == 1){
    legend(2.75,24,legend = c("0","0.2","0.4","0.6"), col = color, pch = 19, lwd = 1.3,bty = "n",cex = 1.4,y.intersp= 0.9)
    text(2,25,"Radius Blur SD Fraction",pos = 4,cex = 1.4)
    next
  }
  for(j in c(1,2,3,4)){
    points(avgvol.r[,j], de$mean$avgs[[i]][[j]],pch= 19,col=color[j])
    lines(avgvol.r[,j], de$mean$avgs[[i]][[j]],col=color[j])
    arrows(avgvol.r[,j], de$mean$avgs[[i]][[j]] - de$sd$avgs[[i]][[j]], 
           avgvol.r[,j], de$mean$avgs[[i]][[j]] + de$sd$avgs[[i]][[j]], 
           length=0.05, angle=90, code=3,col=color[j])
  }
  if(i == 4){
    mtext("Averaged Over Cluster Density", side = 3, outer = TRUE, line = -1.5, cex = 1.25)
  }
}


#Now focus on just Kmax (to find density)
par(mfrow = c(1,1), mar = c(4,5,2,1))
#for(i in 1:4){
plot(de$r, de$mean$mat$Km[[1]][1,],
     ylim=range(c(min(unlist(de$mean$mat$Km)), max(unlist(de$mean$mat$Km)))),
     pch=19, xlab="Cluster Radius", ylab=ylabs[1],type = "n", log= "y",
     cex.lab = 1.75, cex.axis = 1.25)
for(i in 1:3){
  for(j in c(1,2,3,4)){
    points(de$r, de$mean$mat$Km[[i]][j,],pch= 19,col=color[j])
    lines(de$r, de$mean$mat$Km[[i]][j,],col=color[j])
    arrows(de$r[de$sd$mat$Km[[i]][j,] != 0], de$mean$mat$Km[[i]][j,][de$sd$mat$Km[[i]][j,] != 0]-de$sd$mat$Km[[i]][j,][de$sd$mat$Km[[i]][j,] != 0], 
           de$r[de$sd$mat$Km[[i]][j,] != 0], de$mean$mat$Km[[i]][j,][de$sd$mat$Km[[i]][j,] != 0]+de$sd$mat$Km[[i]][j,][de$sd$mat$Km[[i]][j,] != 0], 
           length=0.05, angle=90, code=3,col=color[j])
  }
}
lines(de$r, env.val, col = "black",lwd = 2)
legend(4.25,0.5,legend = c("1","0.4","0.2","0.15"), col = color, pch = 19, lwd = 1.3,bty = "n",cex = 1.4,y.intersp= 0.9)
text(3.3,0.2,"Cluster Density",pos = 4,cex = 1.4)
#}
