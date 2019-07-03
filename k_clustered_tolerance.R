# What is the lowest density of clustering which is indicated with the K function?
library(rapt)
library(data.table)
library(parallel)

# Generate Data -----------------------------------------------------------

under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(1),sep=""),paste('~/Research/point_patterns/Final/system',toString(1),sep=""),scaleUp = TRUE,newRadius = 0.5)
over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(2),sep=""),paste('~/Research/point_patterns/Final/system',toString(2),sep=""),scaleUp = TRUE,newRadius = 0.5)
under.big <- stitch.size(under, boxSize = c(60,60,60))
over.big <- stitch.size(over, boxSize = c(60,60,60))

r <- 3
pic <- seq(0.001, 0.2, by = 0.001)
den <- seq(0.1, 1, by = 0.1)

#upload rrl envelopes
toSub <- fread("C:/Users/galen/Documents/Research/K_cluster_series/cubetoSub.csv")
tests <- fread("C:/Users/galen/Documents/Research/K_cluster_series/cube.csv", drop = 1)
env.r <- tests$V1

#find upper edge of envelopes
percentiles <- c(0.999, 0.99, 0.95)

rvals <- tests[,1]
tvals <- tests[,2:ncol(tests)]

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


library(parallel)
cl <- makePSOCKcluster(detectCores())
clusterExport(cl, c("under.big","over.big","den","r","toSub","toPlotBigs"), envir = environment())
clusterEvalQ(cl,library(spatstat))
clusterEvalQ(cl,library(rapt))

pic <- seq(0.001, 0.2, by = 0.001)


res <- parLapply(cl, pic, function(x){
  
  dic_temp <- matrix(FALSE, nrow = 300, ncol = 3)
  dic_fin <- matrix(FALSE, nrow = length(den), ncol = 3)
  
  for(i in 1:length(den)){
    clust <- makecluster(under.big, over.big, 0.5, 0.5, type = "cr", cr = r, speed = "superfast", pcp = 0.1, den = den[i], pic = x)
    a <- anomK3est(clust[[1]], toSub = toSub[,2], rmax = 25, nrval = 300, correction = "trans")
    for(j in 1:300){
      dic_temp[j,] <- (a$trans[j] > toPlotBigs[j,])
    }
    for(j in 1:3){
      if(any(dic_temp[,j] == TRUE)){
        dic_fin[i,j] <- TRUE
      }
    }
    #cnt <- cnt + 1
    #print(paste(toString(round(cnt/tot, 5)*100),"%", sep = ""))
    b <- 1
    write.csv(b, file = paste("C:/Users/galen/Documents/Research/oled-spatstat-analysis/",toString(x),"_",toString(i),".csv",sep = ""))
  }
  colnames(dic_fin) <- c("99.9%","99%","95%")
  return(dic_fin)})
stopCluster(cl)

save(res, file = "ktoleracedata.RData")

# Analyze Data ------------------------------------------------------------

rm(list = ls())
gc()

library(ggplot2)

load("C:/Users/galen/Documents/Research/K_cluster_series/ktolerancedata_big.RData")
pic <- seq(0.001, 0.5, by = 0.001)
den <- seq(0.1, 1, by = 0.1)

result <- list("99.9%" = matrix(NA, nrow = length(pic), ncol= length(den)),
               "99%" = matrix(NA, nrow = length(pic), ncol= length(den)),
               "95%" = matrix(NA, nrow = length(pic), ncol= length(den)))

for(i in 1:3){
  for(j in 1:500){
    for(k in 1:10){
      result[[i]][j,k] <- res[[j]][k,i] 
    }
  }
}
rm(res,i,j,k)

first_above <- matrix(NA, nrow = 3, ncol = 10)
rownames(first_above) <- c("99.9", "99", "95")

for(i in 1:3){
  for(j in 1:10){
    first_above[i,j] <- match(TRUE, result[[i]][,j])
  }
}

first_above.pic <- matrix(pic[first_above], nrow = 3, byrow = F)
fa.df <- data.frame(ai = c(rep(99.9,10),rep(99,10),rep(95,10)), pic = as.vector(t(first_above.pic))*100, den = rep(den,3))
fa.df$den <- as.factor(fa.df$den)

ggplot(fa.df, aes(x = as.factor(ai), y = pic)) +
  geom_point(aes(col = den), size = 3) +
  xlab("Acceptance Interval %") +
  ylab("Percent of Points in Clusters Needed for Significance") + 
  scale_color_discrete(name="Cluster Density") +
  theme_minimal()




# p-value analysis with DCLF test -----------------------------------------

#upload rrl envelopes
toSub <- fread("C:/Users/galen/Documents/Research/K_cluster_series/cubetoSub_unsorted.csv", select = 2)
tests <- fread("C:/Users/galen/Documents/Research/K_cluster_series/cube_unsorted.csv", drop = 1)
toSub <- fread("~/scratch/Rcode/cubetoSub_unsorted.csv", select = 2)
tests <- fread("~/scratch/Rcode/cube_unsorted.csv", drop = 1)

rvals <- tests[,1]
tvals <- tests[,2:ncol(tests)]

r.cutoff <- 12
r.ind <- which(min(abs(rvals - r.cutoff)) == abs(rvals - r.cutoff))
tvals.chopped <- tvals[1:r.ind,]

diffsum <- apply(tvals.chopped, 2, function(x){sum(x[x > 0])}) #DCLF
mad <- apply(tvals.chopped, 2, max, na.rm = TRUE)

# Calculate vector of integrated deviations for DCLF test

#diffsum <- apply(tvals, 2, function(x){sum(x[x > 0])}) 
#plot(diffsum)
#hist(diffsum)

diffsum.sorted <- sort(diffsum)
mad.sorted <- sort(mad)

r <- 3
pic <- seq(0.01, 0.4, by = 0.01)
den <- seq(0.1, 1, by = 0.1)

pic <- seq(0.01, 0.2, by = 0.075)
den <- c(0.5, 0.8)

n <- 3

under.nums <- seq(2,(n+1),1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,n,1)

cl <- makePSOCKcluster(detectCores())
clusterExport(cl, c("r.ind","pic","den","r","toSub","diffsum.sorted","mad.sorted","under.nums","over.nums"), envir = environment())
clusterEvalQ(cl,library(rapt))

ktol <- parLapply(cl, 1:n, function(j){
  under <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(under.nums[j]),sep=""),paste('~/scratch/Rcode/systems/system',toString(under.nums[j]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/scratch/Rcode/FinalConfigs/FinalConfig',toString(over.nums[j]),sep=""),paste('~/scratch/Rcode/systems/system',toString(over.nums[j]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[j]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[j]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[j]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[j]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  res <- lapply(1:length(pic), function(q){
    x <- pic[q]
    dclf.ps <- rep(NA, length(den))
    mad.ps <- rep(NA, length(den))
    
    for(i in 1:length(den)){
      clust <- makecluster(under.big, over.big, 0.5, 0.5, type = "cr", cr = r, speed = "superfast", pcp = 0.1, den = den[i], pic = x)
      a <- anomK3est(clust[[1]], toSub = toSub$x, rmax = 25, nrval = 300, correction = "trans")
      obs <- a$trans[1:r.ind]
      diffsum.obs <- sum(obs[obs > 0])
      mad.obs <- max(obs)
      
      dclf.ind.obs <- as.numeric(which(min(abs(diffsum.obs-diffsum.sorted)) == abs(diffsum.obs-diffsum.sorted)))
      mad.ind.obs <- as.numeric(which(min(abs(mad.obs-mad.sorted)) == abs(mad.obs-mad.sorted)))
        
      dclf.ps[i] <- 1 - (dclf.ind.obs/length(diffsum.sorted))
      mad.ps[i] <- 1 -  (mad.ind.obs/length(mad.sorted))
    }
    
    b <- 1
    #write.csv(b, file = paste("C:/Users/galen/Documents/Research/oled-spatstat-analysis/",toString(x), '_', toString(j), '.csv',sep = ''))
    if(q > 1){
      #file.remove(paste('~/scratch/Rcode/junk/', toString(q-1), '_', toString(j), '.csv', sep = ''))
      file.remove(paste('~/Research/junk/', toString(q-1), '_', toString(j), '.csv', sep = ''))
    }
    #write.csv(b, file = paste('~/scratch/Rcode/junk/', toString(q), '_', toString(j), '.csv',sep = ''))
    write.csv(b, file = paste('~/Research/junk/', toString(q), '_', toString(j), '.csv',sep = ''))
    
    rm(clust, a, obs, diffsum.obs, mad.obs, mad.ind.obs, dclf.ind.obs)
    gc()
    return(list(dclf = dclf.ps, mad = mad.ps))})
  
  dclf <- matrix(NA, nrow = length(pic), ncol = length(den))
  mad <- matrix(NA, nrow = length(pic), ncol = length(den))
  for(i in 1:length(res)){
    dclf[i,] <- res[[i]]$dclf
    mad[i,] <- res[[i]]$mad
  }
  
  names <- sapply(1:length(den), function(x){paste('den.',toString(x),sep = '')})
  colnames(dclf) <- names
  colnames(mad) <- names
  
  return(list(dclf = dclf, mad = mad))
})
stopCluster(cl)

ktol$pic <- pic
ktol$den <- den

save(ktol, file = "ktolerancedata.RData")


# Analyze p-value data ----------------------------------------------------
load("~/Research/K_cluster_series/190703_ktolerancedata.RData")

pic <- ktol[[length(ktol)-1]]
den <- ktol[[length(ktol)]]

p.dclf <- list()
p.mad <- list()
for(i in 1:(length(ktol)-2)){
  p.dclf[[i]] <- ktol[[i]]$dclf
  p.mad[[i]] <- ktol[[i]]$mad
}

p.vals <- c(0.05, 0.01)

p.toplot.dclf <- list()
p.toplot.mad <- list()
for(i in 1:length(p.vals)){
  p.toplot.dclf[[i]] <- matrix(NA, nrow = length(den), ncol = length(p.dclf))
  p.toplot.mad[[i]] <- matrix(NA, nrow = length(den), ncol = length(p.mad))
}


for(j in 1:length(p.dclf)){
  for(i in 1:length(p.vals)){
    p.check <- p.vals[i]
    
    p.tf.dclf <- p.dclf[[j]] < p.check
    p.tf.mad <- p.mad[[j]] < p.check
    
    p.ind.dclf <- apply(p.tf.dclf, 2, function(x){
      if(!any(x == TRUE) | !any(x == FALSE)){
        return(NA)
      }else{
        return(match(TRUE, x))
      }})
    
    p.ind.mad <- apply(p.tf.mad, 2, function(x){
      if(!any(x == TRUE) | !any(x == FALSE)){
        return(NA)
      }else{
        return(match(TRUE, x))
      }})
    print(paste(toString(i),'_',toString(j)))
    if(!all(is.na(p.ind.mad))){
      p.toplot.mad[[i]][,j] <- pic[p.ind.mad]
    }
    if(!all(is.na(p.ind.dclf))){
      p.toplot.dclf[[i]][,j] <- pic[p.ind.dclf]
    }
    
  }
}

#TO PLOT
p.toplot <- p.toplot.dclf

p.toplot.mean <- lapply(p.toplot, function(x){apply(x, 1, mean, na.rm = TRUE)})
p.toplot.sd <- lapply(p.toplot, function(x){apply(x, 1, function(x){sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))})})

# OR, other averaging method:
library(abind)
all.matrix <- abind(p.dclf, along = 3)
p.mean <- apply(all.matrix, c(1,2), mean, na.rm = TRUE)
p.sd <- apply(all.matrix, c(1,2), sd, na.rm = TRUE)

p.toplot <- list()
for(i in 1:length(p.vals)){
  p.check <- p.vals[i]
  p.tf <- p.mean < p.check
  p.ind <- apply(p.tf, 2, function(x){
    if(!any(x == TRUE) | !any(x == FALSE)){
      return(NA)
    }else{
      return(match(TRUE, x))
    }})
  
  p.toplot[[i]] <- pic[p.ind]
}


plotpts <- 2:length(den)
cols <- c('red', 'blue')
plot(100*den[plotpts], 100*p.toplot.mean[[1]][plotpts], col = cols[1], pch = 19, xlab = 'Intra-cluster concentration', ylab = '% of points in clusters',
     ylim = c(0,20))
lines(100*den[plotpts], 100*p.toplot.mean[[1]][plotpts], col = cols[1], lwd = 2)
arrows(100*den[plotpts], 100*(p.toplot.mean[[1]][plotpts] - p.toplot.sd[[1]][plotpts]), 100*den[plotpts], 100*(p.toplot.mean[[1]][plotpts] + p.toplot.sd[[1]][plotpts]), 
       length=0.05, angle=90, code=3, col = cols[1])
for(i in 2:length(p.toplot.mean)){
  points(100*den[plotpts], 100*p.toplot.mean[[i]][plotpts], col = cols[i], pch = 19)
  lines(100*den[plotpts], 100*p.toplot.mean[[i]][plotpts], col = cols[i], lwd = 2)
  arrows(100*den[plotpts], 100*(p.toplot.mean[[i]][plotpts] - p.toplot.sd[[i]][plotpts]), 100*den[plotpts], 100*(p.toplot.mean[[i]][plotpts] + p.toplot.sd[[i]][plotpts]), 
         length=0.05, angle=90, code=3, col = cols[i])
}
legend(80, 20, legend = sapply(p.vals, toString), col = cols, lwd = 2, title = "p-value", text.width = 10)


