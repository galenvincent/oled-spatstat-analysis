# Maximum separation algorithm code
library(rapt)
library(dplyr)
library(caret)
library(parallel)

## Set up ML stuff
# Upload training data to get PCA transformation:
load('Z:/Galen/Stephan\ APT\ Analysis/stephan_params.RData')
load('Z:/Galen/Stephan\ APT\ Analysis/stephan_sim.res.RData')

data.unlisted <- matrix(NA, nrow = length(sim.res)*nrow(sim.res[[1]]), ncol = 9)

cnt <- 1
for(i in 1:length(sim.res)){
  for(j in 1:nrow(sim.res[[1]])){
    data.unlisted[cnt, 1:4] <- params[[i]]
    data.unlisted[cnt, 5:9] <- sim.res[[i]][j,]
    cnt <- cnt + 1
  }
}

data.unlisted <- as.data.frame(data.unlisted) 
names(data.unlisted) <- c("r","den","rb","gb","Km","Rm","Rdm","Rddm","Kdm")
data.unlisted$rw <- (data.unlisted$r^4 + 6* (data.unlisted$r^2) *(data.unlisted$rb * data.unlisted$r)^2 + 3 *(data.unlisted$rb * data.unlisted$r)^4)/(data.unlisted$r^3 + 3*data.unlisted$r*(data.unlisted$rb * data.unlisted$r)^2)
data.unlisted$sigma <- data.unlisted$rb * data.unlisted$r

#drop NA rows
data.unlisted <- data.unlisted[complete.cases(data.unlisted),]

mets <- c("Km","Rm","Rdm","Rddm","Kdm")
pp.train <- preProcess(data.unlisted[,mets], method = c("scale", "center", "pca", "BoxCox"), thresh = 1)

# Upload RRLs
load('Z:/Galen/Stephan\ APT\ Analysis/stephan_RRL_data_r35.RData')
rrls <- stephan.RRL.data
rm(stephan.RRL.data)
toSub <- rrls[[2]]
rmax <- rrls[[3]]
nr <- rrls[[4]]
rm(rrls)

#Upload models
load('Z:/Galen/Stephan\ APT\ Analysis/ml.models_rw.RData')
rw.model <- models[[3]]
rm(models)
gc()

load('Z:/Galen/Stephan\ APT\ Analysis/ml.models_den.RData')
den.model <- models[[3]]
rm(models)
gc()

#### MSA on simulated data ####
rw <- function(mu, rb){
  sig <- mu*rb
  return((mu^4 + 6*mu^2*sig^2 + 3*sig^4)/(mu^3 + 3*mu*sig^2))
}

rw.test.set <- data.frame('r' = seq(3, 4.5, by = 0.25), 'rb' = c(0.496, 0.423, 0.356, 0.294, 0.232, 0.167, 0.081))
res <- matrix(0, nrow = nrow(rw.test.set), ncol = 3)

nmins <- c(5, 10, 15, 20)
dmaxs <- seq(1, 2, 0.1)

nmins <- c(10, 11)
dmaxs <- c(1.3)
params.to.sweep <- expand.grid(nmins, dmaxs)

msa.all.all <- list()

  
for(i in 1:nrow(rw.test.set)){
  print(i)
  msa.all.all[[i]] <- list()
  
  cl <- makePSOCKcluster(8)
  clusterEvalQ(cl, library(rapt))
  clusterExport(cl, c('msa', 'params.to.sweep', 'rw.test.set', 'i', 'rmax', 'nr', 'toSub', 'rw.model', 'den.model', 'pp.train'))
  
  num.to.do <- 2
  
  msa.all.all[[i]] <- parLapply(cl, 1:num.to.do, function(q){
    print(q)
    under <- read.rcp(paste('C:/Users/galen/Documents/Research/point_patterns/Final/FinalConfig', toString(q), sep = ''),
                      paste('C:/Users/galen/Documents/Research/point_patterns/Final/system', toString(q), sep = ''),
                      scaleUp = TRUE,newRadius = 0.5)
    over <- read.rcp(paste('C:/Users/galen/Documents/Research/point_patterns/Final/FinalConfig', toString(q+1), sep = ''),
                     paste('C:/Users/galen/Documents/Research/point_patterns/Final/system', toString(q+1), sep = ''),
                     scaleUp = TRUE,newRadius = 0.5)
    under.big <- stitch.size(under, boxSize = c(60,60,60))
    over.big <- stitch.size(over, boxSize = c(60,60,60))
    
    clust <- makecluster(under.big, over.big, 0.5, 0.5, pcp = 0.051, den = 0.267, cr = rw.test.set$r[i],
                         rb = TRUE, rbp = rw.test.set$r[i]*rw.test.set$rb[i],
                         toPlot = FALSE, s = q)
    
    X <- clust[[6]]
    
    res <- lapply(1:nrow(params.to.sweep), function(x){
      return(msa(X, params.to.sweep[x,2], params.to.sweep[x,1], params.to.sweep[x,2], params.to.sweep[x,2]))
    })
    
    msa.res <- lapply(res, function(x){
      if(is.na(x)){return(NA)}
      return(data.frame('rad' = x$radius, 'den' = x$den))
    })
    
    msa.summary <- lapply(msa.res, function(x){
      if(is.na(x)){return(data.frame('rc' = NA, 'rwT' = NA, 'rwE' = NA, 'rho' = NA))}
      sigma <- sd(x$rad)
      Rc <- mean(x$rad)
      Rw_true <- mean(x$rad*(4/3)*pi*(x$rad)^3)/mean((4/3)*pi*(x$rad)^3)
      Rw_est <- (Rc^4 + 6*Rc^2*sigma^2 + 3*sigma^4)/(Rc^3 + 3*Rc*sigma^2)
      rho <- mean(x$den)
      return(data.frame('rc' = Rc, 'rwT' = Rw_true, 'rwE' = Rw_est, 'rho' = rho, 'sig'= sigma))
    })
    
    msa.all <- matrix(NA, nrow = length(msa.summary), ncol = 9)
    for(j in 1:length(msa.summary)){
      if(any(is.na(msa.summary[[j]]))){msa.all[j,] <- rep(NA, 9)}
      else {msa.all[j,] <- as.numeric(c(rw.test.set[i,], params.to.sweep[j,], msa.summary[[j]]))}
    }
    msa.all <- as.data.frame(msa.all)
    names(msa.all) <- c('r','rb','Nmin', 'dmax', 'rc', 'rwT', 'rwE', 'rho', 'sig')
    
    #Do k stuff
    k.sim <- K3est(clust[[1]], rmax = rmax, nrval = nr, correction = 'translation')
    t.sim <- data.frame('r' = k.sim$r, 't' = sqrt(k.sim$trans) - toSub)
    met.sim <- as.data.frame(t(unlist(k3metrics(t.sim$r, t.sim$t, toplot = FALSE))))
    names(met.sim) <- c("Km","Rm","Rdm","Rddm","Kdm")
    
    test.pca <- predict(pp.train, met.sim)
    test <- data.frame(met.sim, test.pca)
    
    predictors <- sapply(1:5,function(x){paste('PC',toString(x),sep = '')})
    
    preds.rw <- predict(rw.model, newdata = test[,predictors])
    preds.den <- predict(den.model, newdata = test[,predictors])
    
    msa.all$ml_rw <- rep(preds.rw, nrow(msa.all))
    msa.all$ml_den <- rep(preds.den, nrow(msa.all))
    
    
    return(msa.all)})
  
  stopCluster(cl)
  
}

load('Z:/Galen/MSA/msa.all.all_with_k.RData')

load('Z:/Galen/MSA/msa.all.all.RData')
msa.all.all <- lapply(msa.all.all, function(x){lapply(x, as.matrix)})

msa.avgs <- lapply(msa.all.all, function(x){
  a <- do.call(cbind, x)
  dim(a) <- c(44, 11, 200)
  avgs <- as.data.frame(apply(a, 1:2, mean, na.rm = TRUE))
  names(avgs) <- c('r','rb','Nmin', 'dmax', 'rc', 'rwT', 'rwE', 'rho', 'sig', 'ml_rw', 'ml_den')
  sds <- as.data.frame(apply(a, 1:2, sd, na.rm = TRUE))
  names(sds) <- c('r','rb','Nmin', 'dmax', 'rc', 'rwT', 'rwE', 'rho', 'sig', 'ml_rw', 'ml_den')
  return(list('mean' = avgs, 'sd' = sds))
})

CI.90 <- lapply(msa.all.all, function(x){
  a <- do.call(cbind, x)
  dim(a) <- c(44, 11, 200)
  
  rws <- a[1, 10, ]
  dens <- a[1, 11, ]
  
  n <- 200
  inds <- c(round((n*0.1)/2), round(200 - (n*0.1)/2))
  
  rw.sorted <- sort(rws)
  den.sorted <- sort(dens)
  
  return(data.frame('rw' = rw.sorted[inds], 'den' = den.sorted[inds]))
})

par(mar = c(4, 4, 2, 2), mgp = c(2, 1, 0))

i <- 6
plot(msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rwT, col = msa.avgs[[i]]$mean$Nmin/5, pch = 16,
     xlab = 'dmax', ylab = 'Weighted Radius', ylim = c(1.75, 6))
arrows(msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rwT-msa.avgs[[i]]$sd$rwT, 
       msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rwT+msa.avgs[[i]]$sd$rwT, 
       length=0.05, angle=90, code=3, col = msa.avgs[[i]]$mean$Nmin/5)
abline(h = 4.586, col = 'red', lwd = 2) # True value
# Add mean ML estimate and 90% CI
abline(h = msa.avgs[[i]]$mean$ml_rw[1], col = 'black', lwd = 2, lty = 1)
abline(h = CI.90[[i]]$rw, col = 'black', lwd = 2, lty = 2)
legend(1.5, 3.5, legend = c('Nmin = 5', 'Nmin = 10', 'Nmin = 15', 'Nmin= 20'), col = c(1, 2, 3, 4), pch = 16, bty = 'n')
legend(1.5, 2, legend = c('True Value'), col = 'red', lwd = 2, lty = c(1, 2), bty = 'n')

plot(msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rho, col = msa.avgs[[i]]$mean$Nmin/5, pch = 16,
     xlab = 'dmax', ylab = 'Rho',
     ylim = c(0.2, 1))
arrows(msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rho-msa.avgs[[i]]$sd$rho, 
       msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rho+msa.avgs[[i]]$sd$rho, 
       length=0.05, angle=90, code=3, col = msa.avgs[[i]]$mean$Nmin/5)
abline(h = 0.267, col = 'red', lwd = 2) # True value
# Add mean ML estimate and 90% CI
abline(h = msa.avgs[[i]]$mean$ml_den[1], col = 'black', lwd = 2, lty = 1)
abline(h = CI.90[[i]]$den, col = 'black', lwd = 2, lty = 2)
legend(1, 0.7, legend = c('Nmin = 5', 'Nmin = 10', 'Nmin = 15', 'Nmin= 20'), col = c(1, 2, 3, 4), pch = 16, bty = 'n')
legend(1.35, 1, legend = c('True Value'), col = 'red', lwd = 2, lty = c(1, 2), bty = 'n')

#### Acutal MSA ####

#Data for testing 
under <- read.rcp('C:/Users/galen/Documents/Research/point_patterns/Final/FinalConfig1','C:/Users/galen/Documents/Research/point_patterns/Final/system1',scaleUp = TRUE,newRadius = 0.5)
over <- read.rcp('C:/Users/galen/Documents/Research/point_patterns/Final/FinalConfig2','C:/Users/galen/Documents/Research/point_patterns/Final/system2',scaleUp = TRUE,newRadius = 0.5)
under.big <- stitch.size(under, boxSize = c(60,60,60))
over.big <- stitch.size(over, boxSize = c(60,60,60))

clust <- clustersim(under.big, over.big, 0.5,
                    pcp = 0.0511, 
                    cr = 3.2, 
                    rb = 0.436, 
                    rho1 = 0.267,
                    rho2 = 0.01,
                    pb = 0.1,
                    s = 103,
                    toplot = TRUE)

t.sim <- anomK3est(clust[[1]], toSub, maxr, nr)

plot(t.agg$r, t.agg$T, type = 'l', lwd = 2 ,col = 'red', ylim = c(-1, 4.5), xlim = c(0, 25),
     xlab = 'r (arb.)', ylab = 'T(r)')

nsim <- 500

env.mat <- matrix(0, nrow = nr, ncol = nsim+1)

cl <- makePSOCKcluster(8)
clusterEvalQ(cl, library(rapt))
clusterExport(cl, c('under.big', 'toSub', 'maxr', 'nr'))


env.list <- parLapply(cl, 1:nsim, function(i){

#for(i in 1:nsim){
  over <- read.rcp(paste('C:/Users/galen/Documents/Research/point_patterns/Final/FinalConfig', toString(i+1), sep = ''),
                   paste('C:/Users/galen/Documents/Research/point_patterns/Final/system', toString(i+1), sep = ''),
                   scaleUp = TRUE,newRadius = 0.5)
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  clust <- makecluster(under.big, over.big, 0.5, 0.5, pcp = 0.0511, 
                       rb = TRUE,
                       cr = 4, rbp = 4*0.232, den = 0.267,
                       gb = FALSE, gbp = c(0, 3), 
                       s = i)
  
  t.sim <- anomK3est(clust[[1]], toSub, maxr, nr)
  
  #lines(t.sim$r, t.sim$trans)
  
  return(t.sim)
  
  #if(i == 1){
  #  env.mat[,1] <- t.sim$r
  #}
  
  #env.mat[,i+1] <- t.sim$trans
  #print(i)
})

stopCluster(cl)

env.mat[,1] <- env.list[[1]]$r

for(i in 1:nsim){
  env.mat[,i+1] <- env.list[[i]]$trans
}

par(mar = c(3.5, 3.5, 1, 1), mgp= c(2, 1, 0))
envPlot(env.mat, ylim = c(-4,6), xlim = c(0, 25), leg = FALSE, percentiles = c(0.99, 0.95, 0.9))
lines(t.agg$r, t.agg$T, col = 'black', lwd = 2)

legend(9, 4.5, legend = c('Scaled Reconstruction', 'Simulated'), lwd = 2, lty = 1, col = c('red', 'black'), bty = 'n')


msa <- function(X, dmax, Nmin, denv, der, clust.mark = 'A'){
    X.A <- X[marks(X) == clust.mark]
    X.B <- X[!(marks(X) == clust.mark)]
    
    marks(X.A) <- which(marks(X) == clust.mark)
    marks(X.B) <- which(marks(X) != clust.mark)
    
    #find the nns within dmax of all type A points:
    cp <- closepairs(X.A, rmax = dmax, twice = TRUE, what = 'indices')
    cp <- data.frame('i' = cp$i[order(cp$i)], 'j' = cp$j[order(cp$i)])
    
    #change results to list (nns): for each i, create an entry in a list that contains a vector of its nns
    nns <- list()
    nns.gb <- dplyr::group_by(cp, i)
    nns.labs <- attr(nns.gb, 'labels')$i
    nns.inds <- attr(nns.gb, 'indices')
    nns.inds.no <- (1:npoints(X.A))[-unlist(nns.labs)]
    
    for(k in 1:length(nns.labs)){nns[[nns.labs[k]]] <- cp$j[nns.inds[[k]]+1]}
    nns[nns.inds.no] <- c(0)
    
    # Find individual clusers with more points than Nmax
    diveDeep <- function(is, inds){
      inds.all <- c()
      for(j in is){
        if(nns[[j]][1] == 0){next}
        inds.all <- append(inds.all, nns[[j]])
      }
      inds.all <- unique(inds.all)
      inds.new <- inds.all[!inds.all%in%inds]
      
      #base case
      if(length(inds.new) == 0){return(inds)}
      
      #recursive case
      inds <- append(inds, inds.new)
      return(diveDeep(inds.new, inds))
    }
    
    clusters <- list()
    #ind.list <- 1:npoints(X.A)
    ind.list <- nns.labs
    cnt <- 1
    
    while(length(ind.list) > 0){
      
      to.do <- ind.list[1]
      clusters[[cnt]] <- diveDeep(to.do, to.do) 
      
      ind.list <- ind.list[!ind.list%in%clusters[[cnt]]]
      
      cnt <- cnt + 1
    }
    
    cluster.sizes <- sapply(clusters, length)
    clusters.pass <- clusters[cluster.sizes >= Nmin]
    if(length(clusters.pass) == 0){
      print('No clusters found')
      return(NA)  
    }
    
    # Get background points in clusters
    X.clusters.A <- X.A[unlist(clusters.pass)]
    mks <- rep(1:length(clusters.pass), sapply(clusters.pass, length))
    marks(X.clusters.A) <- mks
    
    B.in.clusters <- list()
    
    cp.AB <- crosspairs(X.clusters.A, X.B, rmax = denv, what = 'indices')
    cp.AB <- data.frame('i' = cp.AB$i[order(cp.AB$i)], 'j' = cp.AB$j[order(cp.AB$i)])
    
    AB.gb <- dplyr::group_by(cp.AB, i)
    AB.labs <- attr(AB.gb, 'labels')$i
    AB.inds <- attr(AB.gb, 'indices')
    
    for(k in 1:length(clusters.pass)){
      clust.inds <- which(marks(X.clusters.A) == k)
      clust.inds <- clust.inds[clust.inds %in% AB.labs]
      data.inds <- which((AB.labs %in% clust.inds) == TRUE)
      to.pull <- unlist(lapply(data.inds, function(l){AB.inds[[l]]+1}))
      B.in.clusters[[k]] <- unique(cp.AB$j[to.pull])
    }

    #crosspairs to erode
    matrix.bgnd <- X.B[-unique(unlist(B.in.clusters))]
    matrix.A <- X.A[-unique(unlist(clusters.pass))]
    #combine into full matrix
    coo.mat <- rbind(coords(matrix.bgnd), coords(matrix.A))
    matrix.all <- pp3(coo.mat$x, coo.mat$y, coo.mat$z, domain(X))
    
    #make pp3s of the cluster stuff so far:
    marks(X.clusters.A) <- unlist(clusters.pass)
    
    X.clusters.B <- X.B[unlist(B.in.clusters)]
    marks(X.clusters.B) <- unlist(B.in.clusters)
    
    # Check these
    cp.AM <- crosspairs(X.clusters.A, matrix.all, rmax = der, what = 'indices')
    cp.BM <- crosspairs(X.clusters.B, matrix.all, rmax = der, what = 'indices')
    
    A.remove <- marks(X.clusters.A[unique(cp.AM$i)])
    B.remove <- marks(X.clusters.B[unique(cp.BM$i)])
    
    
    A.clusters.eroded <- lapply(clusters.pass, function(x){
      return(x[!x%in%A.remove])
    })
    
    B.clusters.eroded <- lapply(B.in.clusters, function(x){
      return(x[!x%in%B.remove])
    })
    
    nclusters <- length(A.clusters.eroded)
    
    cluster.den <- sapply(1:nclusters, function(x){
      length(A.clusters.eroded[[x]])/(length(B.clusters.eroded[[x]]) + length(A.clusters.eroded[[x]]))
    })
    
    # Guinier Radius
    cluster.Rg <- sapply(A.clusters.eroded, function(x){
      coo <- coords(X.A[x])
      com <- apply(coo, 2, mean)
      rs2 <- apply(t(t(coo) - com), 1, function(y){sum(y^2)})
      Rg <- sqrt(sum(rs2)/nrow(coo))
      Dg <- 2*sqrt(5/3)*Rg
      return(Dg/2)
    })
    
    A.cluster.inds.orig <- lapply(A.clusters.eroded, function(x){marks(X.A)[x]})
    B.cluster.inds.orig <- lapply(B.clusters.eroded, function(x){marks(X.B)[x]})
    
    return(list('radius'=cluster.Rg, 'den'=cluster.den, 'A' = A.cluster.inds.orig, 'B' = B.cluster.inds.orig))
}



