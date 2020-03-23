# Cluster envelopes
library(parallel)
library(data.table)
library(rapt)


#### Generate Data ####
r.params <- data.frame('r' = c(4, 5), 'rb' = c(0.3476223, 0.1043725))

toSub <- fread('RCP_RRL_toSub_5.1percent.csv')

cluster.envs <- list()

for(j in 1:2){
  
  cl <- makePSOCKcluster(8)
  clusterEvalQ(cl, library(rapt))
  clusterExport(cl, c('toSub', 'r.params', 'j'))
  
  #for(i in 1:500){
  n <- 5
  res <- parLapply(cl, 1:n, function(i){
    rcp_path <- 'C:/Users/galen/Documents/Research/point_patterns/Final/'
    #rcp_path <- '~/scratch/Rcode/RCP/'
    under <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i), sep = ''),
                      paste(rcp_path, 'system', toString(i), sep = ''),
                      scaleUp = TRUE,newRadius = 0.5)
    over <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i+1), sep = ''),
                     paste(rcp_path, 'system', toString(i+1), sep = ''),
                     scaleUp = TRUE,newRadius = 0.5)
    under.big <- stitch.size(under, boxSize = c(60,60,60))
    over.big <- stitch.size(over, boxSize = c(60,60,60))
    
    clust <- clustersim(under.big, over.big, 0.5,
                        pcp = 0.05114235,
                        cr = r.params$r[j],
                        rho1 = 0.2118379,
                        rho2 = 0.02494626,
                        rb = r.params$rb[j],
                        pb = 0,
                        tol = 0.005,
                        s = i)
    if(is.numeric(clust)){
      return(NA)
    }
    
    result <- anomK3est(clust[[1]], toSub, 35, 400)
    return(result)
  })
    
  res.mat <- matrix(NA, nrow = 400, ncol = n+1)
  res.mat[,1] <- res[[1]]$r
  for(k in 1:n){
    res.mat[,k+1] <- res[[k]]$trans
  }
  
  cluster.envs[[j]] <- res.mat
}

save(cluster.envs, file = 'cluster.envelopes.RData')
#### Test Data ####
cl <- makePSOCKcluster(8)
clusterEvalQ(cl, library(rapt))
clusterExport(cl, c('toSub'))

n <- 32
t1 <- Sys.time()
res <- parLapply(cl, 1:n, function(i){
  rcp_path <- 'C:/Users/galen/Documents/Research/point_patterns/Final/'
  #rcp_path <- '~/scratch/Rcode/RCP/'
  under <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i), sep = ''),
                    paste(rcp_path, 'system', toString(i), sep = ''),
                    scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, 'FinalConfig', toString(i+1), sep = ''),
                   paste(rcp_path, 'system', toString(i+1), sep = ''),
                   scaleUp = TRUE,newRadius = 0.5)
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  clust <- clustersim(under.big, over.big, 0.5,
                      pcp = 0.05114235,
                      cr = 4.5,
                      rho1 = 0.17,
                      rho2 = 0.034,
                      rb = 0.2384064,
                      pb = 0.2,
                      tol = 0.005,
                      s = i)
  if(is.numeric(clust)){
    return(NA)
  }
  
  result <- anomK3est(clust[[1]], toSub, 35, 400)
  return(result)
})
stopCluster(cl)
t2 <- Sys.time()
print(t2 - t1)

res.mat <- matrix(NA, nrow = 400, ncol = n+1)
res.mat[,1] <- res[[1]]$r
for(k in 1:n){
  res.mat[,k+1] <- res[[k]]$trans
}


#### Analyze Data ####
load('Z:/Galen/Stephan\ APT\ Analysis/cluster.envelopes.RData')

# Compare to the k-function measured on stephan's data. 
# You have two different cr/rb combinations that give the correct weighted radius; try them both.
r.params <- data.frame('r' = c(4, 5), 'rb' = c(0.3476223, 0.1043725))

i <- 1
tests <- cluster.envs[[i]]
percentiles = c(0.99, 0.95, 0.90)

#tests <- res.mat
#percentiles = c(0.95, 0.92, 0.90)

ylab = expression(sqrt('K'[3]*'(r)')*' Anomaly')
xlab = 'r'
#color <- c("lightpink", "mediumpurple", "lightskyblue")
color <- c('cornflowerblue','palegreen3', 'firebrick1')
# break up data into r values and test results
rvals <- tests[,1]
tvals <- tests[,2:ncol(tests)]

# number of tests done
nTests <- ncol(tvals)
# get the range of indices for which each percentile spans
prange <- percentiles * nTests

# sort the results at each r value from lowest to highest
sortedtVals <- t(apply(tvals, 1, sort))
# select the high end indexes based on being 1/2 of the percentile span
# from the middle of the tests
percentileIndicesSmall <- round(nTests/2) - floor(prange/2)
# do the same for the low end
percentileIndicesBig <- round(nTests/2) + floor(prange/2)

# grab out the columns from the sorted test results that we will plot
toPlotBigs <- matrix(0, nrow = nrow(tvals), ncol = length(percentiles))
toPlotSmalls <- matrix(0, nrow = nrow(tvals), ncol = length(percentiles))
for(i in 1:length(percentiles)) {
  toPlotBigs[,i] <- sortedtVals[,percentileIndicesBig[i]]
  toPlotSmalls[,i] <- sortedtVals[,percentileIndicesSmall[i]]
}

# plot the envelopes from the percentile data
#par(oma = c(0, 2, 0, 0))
xlim <- c(0, 20)
ylim <- c(-2.5, 4.5)

toplt <- data.frame(rvals,tvals[,1])
par(mar = c(4, 4, 2, 2), mgp = c(2, 1, 0))
plot(toplt, type = "n",
     ylab = ylab, xlab = xlab,
     ylim = ylim, xlim = xlim, xaxt = "n")
#axis(1, at = 0:xlim[2], labels=FALSE)
axis(1, at = seq(0,xlim[2],by=5))
a <- c(rvals, rev(rvals))
for(i in 1:length(percentiles)) {
  polygon(a, c(toPlotBigs[,i], rev(toPlotSmalls[,i])), col = color[i])
  #,border=color[i],lwd=2)
}
abline(h = 0, lty = 2, lwd = 1, col="black")

legend(0, -0.5, legend = c(paste(toString(percentiles[1]*100), "% AI"),
                              paste(toString(percentiles[2]*100), "% AI"),
                              paste(toString(percentiles[3]*100), "% AI")),
       col=c(color[1], color[2], color[3]),
       lty = c(1,1,1), lwd = c(10,10,10))
legend(6, -1.75, legend = 'Reconstruction', col = 'black', lwd = 2, lty = 1)

# Add stephan's data - run 'stephan_apt_analysis.R' code first to get out t.agg.
lines(t.agg, col = 'black', lwd = 2, lty = 1)
