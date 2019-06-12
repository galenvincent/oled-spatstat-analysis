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

