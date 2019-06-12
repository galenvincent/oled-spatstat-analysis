# Exploring K and G function 
library(data.table)
library(rapt)
library(moments)
library(parallel)

under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(1),sep=""),paste('~/Research/point_patterns/Final/system',toString(1),sep=""),scaleUp = TRUE,newRadius = 0.5)
over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(2),sep=""),paste('~/Research/point_patterns/Final/system',toString(2),sep=""),scaleUp = TRUE,newRadius = 0.5)
under.big <- stitch.size(under, boxSize = c(60,60,60))
over.big <- stitch.size(over, boxSize = c(60,60,60))

clust <- makecluster(under.big, over.big, 0.5, 0.5, type= "ppc", ppc = 5, pic = 0.25, toPlot = TRUE, pcp = 0.1)


nn <- nndist(clust[[1]], k = 1:10)
h <- apply(nn, 2, function(x){hist(x, breaks = seq(0,max(x),len = 300))})

g <- lapply(h, function(x){wid <- x$breaks[2]-x$breaks[1] 
                           return(cumsum(wid*x$density))})

plot(h[[1]]$mids, g[[1]], type = "l", xlim = c(0,6))
for(i in 2:10){
  lines(h[[i]]$mids, g[[i]])
}

means <- apply(nn, 2, mean)
plot(1:10,means, xlim = c(0,11), ylim = c(0,4), main = "mean")

sds <- apply(nn, 2, sd)
plot(1:10, sds, xlim = c(0,11), main = "sd")

skew <- apply(nn, 2, skewness)
plot(1:10,skew, xlim = c(0,11), main = "skew")

kurt <- apply(nn, 2, kurtosis)
plot(1:10, kurt, xlim = c(0,11), main = "kurtosis")

nsim <- 2000
nk <- 10
maxr <- 20
step <- 0.1
pselect <- 0.1
plist <- seq(1,nsim)
cl <- makePSOCKcluster(detectCores())
clusterExport(cl,c("percentSelect"))
clusterExport(cl, c("under.big","nk","step","maxr","pselect"))
clusterEvalQ(cl,library(spatstat))
hapl <- parLapply(cl, plist, function(x){
  X <- percentSelect(pselect, under.big, x)
  nn <- nndist(X, k = 1:nk)
  h <- apply(nn, 2, function(x){hist(x, breaks = seq(0,maxr, step), plot = FALSE)})
  g <- sapply(h, function(x){cumsum((x$breaks[2]-x$breaks[1])*x$density)})
  return(g)
})
stopCluster(cl)

nns <- list()
for(i in 1:nk){
  nns[[i]] <- matrix(0, nrow = nrow(hapl[[1]]), ncol = nsim)
  for(j in 1:nsim){
    nns[[i]][,j] <- hapl[[j]][,i]
  }
}

nnss <- lapply(nns, function(x){t(apply(x,1,sort))})
mids <- lapply(nnss, function(x){x[,nsim/2]})
nnssanom <- lapply(nnss, function(x){mid <- x[,nsim/2]
                                    apply(x, 2, function(y, mid){y - mid}, mid)})

obs <- nndist(clust[[1]], k = 1:nk)
hobs <- apply(obs, 2, function(x){hist(x, breaks = seq(0,maxr, step), plot = FALSE)})
gobs <- sapply(hobs, function(x){cumsum((x$breaks[2]-x$breaks[1])*x$density)})

obsanom <- matrix(NA, ncol = ncol(gobs), nrow = nrow(gobs))
for(i in 1:length(mids)){
  obsanom[,i] <- gobs[,i]- mids[[i]]
}

for(j in 1:nk){
  envp <- c(0.95, 0.9, 0.8)
  envibig <- (1-((1-envp)/2))*ncol(nnssanom[[j]])
  envismall <- ((1-envp)/2)*ncol(nnssanom[[j]])
  color <- c("lightskyblue","mediumpurple","lightpink")
  a <- c(seq(0,maxr-step,step), seq(maxr-step,0,-step))
  plot(seq(0,maxr-step,step),obsanom[,j], type = "l", 
       xlim = c(0,8), ylim = c(-max(obsanom), max(obsanom)), 
       main = paste("nn",toString(j), sep = ""), xlab = "r", ylab = "G anomaly")
  for(i in 1:3){
    polygon(a,c(nnssanom[[j]][,envibig[i]],rev(nnssanom[[j]][,envismall[i]])),col=color[i])#,border=color[i],lwd=2)
  }
}





