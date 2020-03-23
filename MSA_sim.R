# MSA simulations
library(rapt)
library(dplyr)
library(caret)
library(parallel)
library(data.table)

#Load MSA function
msa <- function(X, dmax, Nmin, denv, der, clust.mark = c('A', 'B')){
  X.A <- X[marks(X) %in% clust.mark]
  X.B <- X[!(marks(X) %in% clust.mark)]
  
  marks(X.A) <- which(marks(X) %in% clust.mark)
  marks(X.B) <- which(!(marks(X) %in% clust.mark))
  
  #find the nns within dmax of all type A points:
  cp <- closepairs(X.A, rmax = dmax, twice = TRUE, what = 'indices')
  cp <- data.frame('i' = cp$i[order(cp$i)], 'j' = cp$j[order(cp$i)])
  
  #change results to list (nns): for each i, create an entry in a list that contains a vector of its nns
  nns <- list()
  nns.gb <- dplyr::group_by(cp, i)
  nns.labs <- attr(nns.gb, 'groups')$i
  nns.inds <- attr(nns.gb, 'groups')$.rows
  nns.inds.no <- (1:npoints(X.A))[-unlist(nns.labs)]
  
  for(k in 1:length(nns.labs)){nns[[nns.labs[k]]] <- cp$j[nns.inds[[k]]]}
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
  AB.labs <- attr(AB.gb, 'groups')$i
  AB.inds <- attr(AB.gb, 'groups')$.rows
  
  for(k in 1:length(clusters.pass)){
    clust.inds <- which(marks(X.clusters.A) == k)
    clust.inds <- clust.inds[clust.inds %in% AB.labs]
    data.inds <- which((AB.labs %in% clust.inds) == TRUE)
    to.pull <- unlist(lapply(data.inds, function(l){AB.inds[[l]]}))
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
  
  # Intra-cluster density
  cluster.den <- sapply(1:nclusters, function(x){
    length(A.clusters.eroded[[x]])/(length(B.clusters.eroded[[x]]) + length(A.clusters.eroded[[x]]))
  })
  
  #background density
  bgnd.total <- npoints(X) - length(unlist(A.clusters.eroded)) - length(unlist(B.clusters.eroded))
  bgnd.A <- npoints(X.A) - length(unlist(A.clusters.eroded))
  bgnd.den <- bgnd.A/bgnd.total
  
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
  
  return(list('radius'=cluster.Rg, 'den'=cluster.den, 'bgnd.den' = bgnd.den, 'A' = A.cluster.inds.orig, 'B' = B.cluster.inds.orig))
}

## Set up ML stuff
# Upload training data to get PCA transformation:
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/200131_params.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/200131_sim.res.RData')

data.unlisted <- matrix(NA, nrow = length(sim.res)*nrow(sim.res[[1]]), ncol = 10)

cnt <- 1
for(i in 1:length(sim.res)){
  for(j in 1:nrow(sim.res[[1]])){
    data.unlisted[cnt, 1:5] <- params[[i]]
    data.unlisted[cnt, 6:10] <- sim.res[[i]][j,]
    cnt <- cnt + 1
  }
}

data.unlisted <- as.data.frame(data.unlisted) 
names(data.unlisted) <- c("cr","rho1","rho2","rb","pb","Km","Rm","Rdm","Rddm","Kdm")
data.unlisted$sigma <- data.unlisted$rb * data.unlisted$cr
data.unlisted$rw <- (data.unlisted$cr^4 + 6* data.unlisted$cr^2 * data.unlisted$sigma^2 + 3*data.unlisted$sigma^4)/
  (data.unlisted$cr^3 + 3*data.unlisted$cr*data.unlisted$sigma^2)

#drop NA rows
data.unlisted <- data.unlisted[complete.cases(data.unlisted),]

mets <- c("Km","Rm","Rdm","Rddm","Kdm")
pp.train <- preProcess(data.unlisted[,mets], method = c("scale", "center", "pca", "BoxCox"), thresh = 1)

# Upload RRLs
load('C:/Users/galen/Documents/Research/K_cluster_series/RCP_RRL_5.1percent.RData')
toSub <- fread('RCP_RRL_toSub_5.1percent.csv')$`rcp_rrl[[2]]`
toSub <- rcp_rrl[[2]]
rmax <- rcp_rrl[[3]]
nr <- rcp_rrl[[4]]
rm(rcp_rrl)

#Upload models
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_rw.RData')
rw.model <- models[[3]]
rm(models)
gc()

load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_rho1.RData')
den.model <- models[[3]]
rm(models)
gc()

load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_rho2.RData')
bgnd.model <- models[[3]]
rm(models)
gc()

#### MSA on simulated data ####

# Extract neccesary rbs
# ESTIMATED RW: 5.192789
# ESTIMATED RHO1: 0.2118379
# ESTIMATED RHO2: 0.02494626
rw.solve <- function(rb, mu, rw){
  sig <- mu*rb
  diff <- abs((mu^4 + 6*mu^2*sig^2 + 3*sig^4)/(mu^3 + 3*mu*sig^2) - rw)
  return(diff)
}
mu2 = 5
rw2 = 5.192789
op <- optimize(rw.solve, c(0,1), mu2 = mu2, rw2 = rw2)
op$minimum

r.params <- matrix(c(4, 5, 0.3476223, 0.1043725), nrow = 2, ncol = 2, byrow = FALSE)
clust.param.set <- matrix(NA, nrow = 6, ncol = 3)
clust.param.set[,1:2] <- rbind(r.params, r.params, r.params)
clust.param.set[,3] <- c(0.01, 0.01, 0.02, 0.02, 0.03, 0.03)
clust.param.set <- as.data.frame(clust.param.set)
names(clust.param.set) <- c('r', 'rb', 'rho2')

nmins <- c(5, 10, 15, 20)
dmaxs <- seq(1, 2, 0.1)

#nmins <- c(10, 11)
#dmaxs <- c(1.3)
params.to.sweep <- expand.grid(nmins, dmaxs)

msa.all.all <- list()


for(i in 1:nrow(clust.param.set)){
  print(i)
  msa.all.all[[i]] <- list()
  
  cl <- makePSOCKcluster(8)
  clusterEvalQ(cl, library(rapt))
  clusterExport(cl, c('msa', 'params.to.sweep', 'clust.param.set', 'i', 'rmax', 'nr', 
                      'toSub', 'rw.model', 'den.model', 'bgnd.model', 'pp.train'))
  
  num.to.do <- 2
  
  msa.all.all[[i]] <- parLapply(cl, 1:num.to.do, function(q){
    print(q)
    rcp_path <- 'C:/Users/galen/Documents/Research/point_patterns/Final/'
    #rcp_path <- '~/scratch/Rcode/RCP/'
    under <- read.rcp(paste(rcp_path, 'FinalConfig', toString(q), sep = ''),
                      paste(rcp_path, 'system', toString(q), sep = ''),
                      scaleUp = TRUE,newRadius = 0.5)
    over <- read.rcp(paste(rcp_path, 'FinalConfig', toString(q+1), sep = ''),
                     paste(rcp_path, 'system', toString(q+1), sep = ''),
                     scaleUp = TRUE,newRadius = 0.5)
    under.big <- stitch.size(under, boxSize = c(60,60,60))
    over.big <- stitch.size(over, boxSize = c(60,60,60))
    
    clust <- clustersim(under.big, over.big, 0.5,
                        pcp = 0.05114235,
                        cr = clust.param.set$r[i],
                        rho1 = 0.2118379,
                        rho2 = clust.param.set$rho2[i],
                        rb = clust.param.set$rb[i],
                        pb = 0,
                        tol = 0.005,
                        s = q)
    if(is.numeric(clust)){
      return(NA)
    }
    X <- clust[[2]]
    
    res <- lapply(1:nrow(params.to.sweep), function(x){
      return(msa(X, params.to.sweep[x,2], params.to.sweep[x,1], params.to.sweep[x,2], params.to.sweep[x,2]))
    })
    
    msa.res <- lapply(res, function(x){
      if(is.na(x)){return(NA)}
      return(data.frame('rad' = x$radius, 'den' = x$den, 'bgnd.den' = x$bgnd.den))
    })
    
    msa.summary <- lapply(msa.res, function(x){
      if(is.na(x)){return(data.frame('rc' = NA, 'rwT' = NA, 'rwE' = NA, 'rho1' = NA, 'rho2' = NA))}
      sigma <- sd(x$rad)
      Rc <- mean(x$rad)
      Rw_true <- mean(x$rad*(4/3)*pi*(x$rad)^3)/mean((4/3)*pi*(x$rad)^3)
      Rw_est <- (Rc^4 + 6*Rc^2*sigma^2 + 3*sigma^4)/(Rc^3 + 3*Rc*sigma^2)
      rho1 <- mean(x$den)
      rho2 <- mean(x$bgnd.den)
      return(data.frame('rc' = Rc, 'rwT' = Rw_true, 'rwE' = Rw_est, 'rho1' = rho1, 'rho2' = rho2, 'sig'= sigma))
    })
    
    msa.all <- matrix(NA, nrow = length(msa.summary), ncol = 11)
    for(j in 1:length(msa.summary)){
      if(any(is.na(msa.summary[[j]]))){msa.all[j,] <- rep(NA, 11)}
      else {msa.all[j,] <- as.numeric(c(clust.param.set[i,], params.to.sweep[j,], msa.summary[[j]]))}
    }
    msa.all <- as.data.frame(msa.all)
    names(msa.all) <- c('r','rb','rho2_T','min', 'dmax', 'rc', 'rwT', 'rwE', 'rho1', 'rho2', 'sig')
    
    #Do k stuff
    result <- anomK3est(clust[[1]], toSub, rmax, nr)
    rvals <- result$r
    tvals <- result$trans
    rvals.new <- rvals[15:length(rvals)]
    tvals.new <- tvals[15:length(rvals)]
    
    met.sim <- as.data.frame(t(unlist(k3metrics(rvals.new, tvals.new, toplot = FALSE))))
    names(met.sim) <- c("Km","Rm","Rdm","Rddm","Kdm")
    
    test.pca <- predict(pp.train, met.sim)
    test <- data.frame(met.sim, test.pca)
    
    predictors <- sapply(1:5,function(x){paste('PC',toString(x),sep = '')})
    
    preds.rw <- predict(rw.model, newdata = test[,predictors])
    preds.rho1 <- predict(den.model, newdata = test[,predictors])
    preds.rho2 <- predict(bgnd.model, newdata = test[,predictors])
    
    msa.all$ml_rw <- rep(preds.rw, nrow(msa.all))
    msa.all$ml_rho1 <- rep(preds.rho1, nrow(msa.all))
    msa.all$ml_rho2 <- rep(preds.rho2, nrow(msa.all))
    
    
    return(msa.all)})
  
  stopCluster(cl)
}

#### PLOTS ####
load('Z:/Galen/MSA/msa.sim.results.RData')

msa.all.all <- lapply(msa.all.all, function(x){lapply(x, as.matrix)})
n <- 500

msa.avgs <- lapply(msa.all.all, function(x){
  ibad <- which(sapply(x, nrow) == 1)
  if(!is.empty(ibad)){
    x <- x[-ibad]
    n.new <- n - length(ibad)
  }else{
    n.new <- n
  }

  a <- do.call(cbind, x)
  dim(a) <- c(44, 14, n.new)
  avgs <- as.data.frame(apply(a, 1:2, mean, na.rm = TRUE))
  names(avgs) <- c('r','rb', 'rho2_T', 'Nmin', 'dmax', 'rc', 'rwT', 'rwE', 'rho1', 'rho2', 'sig', 'ml_rw', 'ml_rho1', 'ml_rho2')
  sds <- as.data.frame(apply(a, 1:2, sd, na.rm = TRUE))
  names(sds) <- c('r','rb','rho2_T', 'Nmin', 'dmax', 'rc', 'rwT', 'rwE', 'rho1', 'rho2', 'sig', 'ml_rw', 'ml_rho1', 'ml_rho2')
  return(list('mean' = avgs, 'sd' = sds))
})

CI.90 <- lapply(msa.all.all, function(x){
  ibad <- which(sapply(x, nrow) == 1)
  if(!is.empty(ibad)){
    x <- x[-ibad]
    n.new <- n - length(ibad)
  }else{
    n.new <- n
  }
  
  a <- do.call(cbind, x)
  dim(a) <- c(44, 14, n.new)
  
  rws <- a[1, 12, ]
  rho1 <- a[1, 13, ]
  rho2 <- a[1, 14, ]
  
  inds <- c(round((n*0.1)/2), round(n - (n*0.1)/2))
  
  rw.sorted <- sort(rws)
  rho1.sorted <- sort(rho1)
  rho2.sorted <- sort(rho2)
  
  return(data.frame('rw' = rw.sorted[inds], 'rho1' = rho1.sorted[inds], 'rho2' = rho2.sorted[inds]))
})

par(mar = c(4, 4, 2, 2), mgp = c(2, 1, 0))

i <- 5 # Used i = 5 for paper plots

plot(msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rwT, col = msa.avgs[[i]]$mean$Nmin/5, pch = 16,
     xlab = 'dmax', ylab = 'Weighted Radius', ylim = c(1.5, 7.5))
arrows(msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rwT-msa.avgs[[i]]$sd$rwT, 
       msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rwT+msa.avgs[[i]]$sd$rwT, 
       length=0.05, angle=90, code=3, col = msa.avgs[[i]]$mean$Nmin/5)
# Add mean ML estimate and 90% CI
abline(h = msa.avgs[[i]]$mean$ml_rw[1], col = 'black', lwd = 2, lty = 1)
abline(h = CI.90[[i]]$rw, col = 'black', lwd = 2, lty = 2)
abline(h = 5.192789, col = 'red', lwd = 2, lty = 2) # True value
legend(1.5, 3.5, legend = c('Nmin = 5', 'Nmin = 10', 'Nmin = 15', 'Nmin= 20'), col = c(1, 2, 3, 4), pch = 16, bty = 'n')
legend(1.5, 2, legend = c('True Value'), col = 'red', lwd = 2, lty = 2, bty = 'n')
legend(1, 7, legend = c('Mean ML Estimate', 'ML Estimate 90% CI'), col = 'black', lwd = 2, lty = c(1,2), bty = 'n')


plot(msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rho1, col = msa.avgs[[i]]$mean$Nmin/5, pch = 16,
     xlab = 'dmax', ylab = 'Rho1',
     ylim = c(0.1, 1))
arrows(msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rho1-msa.avgs[[i]]$sd$rho1, 
       msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rho1+msa.avgs[[i]]$sd$rho1, 
       length=0.05, angle=90, code=3, col = msa.avgs[[i]]$mean$Nmin/5)
abline(h = 0.2118379, col = 'red', lwd = 2) # True value
# Add mean ML estimate and 90% CI
abline(h = msa.avgs[[i]]$mean$ml_rho1[1], col = 'black', lwd = 2, lty = 1)
abline(h = CI.90[[i]]$rho1, col = 'black', lwd = 2, lty = 2)
legend(1, 0.7, legend = c('Nmin = 5', 'Nmin = 10', 'Nmin = 15', 'Nmin= 20'), col = c(1, 2, 3, 4), pch = 16, bty = 'n')
legend(1.35, 1, legend = c('True Value'), col = 'red', lwd = 2, lty = 1, bty = 'n')
legend(1.35, 0.9, legend = c('Mean ML Estimate', 'ML Estimate 90% CI'), col = 'black', lwd = 2, lty = c(1,2), bty = 'n')


plot(msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rho2, col = msa.avgs[[i]]$mean$Nmin/5, pch = 16,
     xlab = 'dmax', ylab = 'Rho2',
     ylim = c(0, 0.06))
arrows(msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rho2-msa.avgs[[i]]$sd$rho2, 
       msa.avgs[[i]]$mean$dmax, msa.avgs[[i]]$mean$rho2+msa.avgs[[i]]$sd$rho2, 
       length=0.05, angle=90, code=3, col = msa.avgs[[i]]$mean$Nmin/5)
# Add mean ML estimate and 90% CI
abline(h = msa.avgs[[i]]$mean$ml_rho2[1], col = 'black', lwd = 2, lty = 1)
abline(h = CI.90[[i]]$rho2, col = 'black', lwd = 2, lty = 2)
#abline(h = 0.01, col = 'red', lwd = 2) # True value
#abline(h = 0.02, col = 'red', lwd = 2) # True value
abline(h = 0.03, col = 'red', lwd = 2, lty = 2) # True value
legend(1, 0.02, legend = c('Nmin = 5', 'Nmin = 10', 'Nmin = 15', 'Nmin= 20'), col = c(1, 2, 3, 4), pch = 16, bty = 'n')
legend(1.35, 0.02, legend = c('True Value'), col = 'red', lwd = 2, lty = 2, bty = 'n')
legend(1.35, 0.01, legend = c('Mean ML Estimate', 'ML Estimate 90% CI'), col = 'black', lwd = 2, lty = c(1,2), bty = 'n')



