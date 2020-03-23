# Stephan APT analysis for paper

library(rapt)
library(data.table)
library(parallel)
library(dplyr)
library(caret)

# Upload APT data
pos <- readPOS('Z:/Galen/Stephan\ APT\ Analysis/R34_06365-v01-CentralCubeExclPole.pos')
rng.full <- readRRNG('Z:/Galen/Stephan\ APT\ Analysis/6365_Alalloy_rngs.RRNG')
rng.simple <- readRRNG('Z:/Galen/Stephan\ APT\ Analysis/6365_Alalloy_SimpleRngs.RRNG')

ms <- createSpec(pos, res = 0.05)
ms.log <- transformIntensity(ms, method = 'log10')

prettyPlot(ms.log, rng = rng.full, xlim = c(0, 75), main= 'Full RRNG')
prettyPlot(ms.log, rng = rng.simple, xlim = c(0, 75), main = 'Simple RRNG')

pos.r <- rngPOS(pos, rng.simple)
cnts <- rngCount(pos.r, rng.simple)
cnts.totals <- group_by(cnts, name) %>%
               summarize(counts = sum(counts), fraction = sum(fraction))
cnts.totals

frac <- sum(cnts.totals$fraction[2:3])

win.r <- box3(c(min(pos.r$x), max(pos.r$x)), 
            c(min(pos.r$y), max(pos.r$y)), 
            c(min(pos.r$z), max(pos.r$z)))

# Create pp3 of different elements
pp3.full <- createSpat(pos.r, win = win.r)
marks(pp3.full) <- pos.r$mark

pp3.aggregates <- subset(pp3.full, marks == 'Mg1' | marks == 'Zn1' )
pp3.mg <- subset(pp3.full, marks == 'Mg1')
pp3.zn <- subset(pp3.full, marks == 'Zn1')

#plot3d.pp3(pp3.mg)
#plot3d.pp3(pp3.zn)

# Scale up pp3's to simulation intensity (for consistency)

# Get original simulation scaling:
# rcp_path = 'C:/Users/galen/Documents/Research/point_patterns/Final'
# lambda <- rep(NA, 509)
# for(i in 1:509){
#   under <- read.rcp(paste(rcp_path, '/FinalConfig', toString(i), sep=''),
#                   paste(rcp_path, '/system', toString(i), sep=''),
#                   scaleUp = TRUE,newRadius = 0.5)
#   under.big <- stitch.size(under, boxSize = c(60,60,60))
# 
#   lambda[i] <- npoints(under.big)/volume(domain(under.big))
#   print(i)
#   print(lambda[i])
# }
# lambda.mu <- mean(lambda[lambda > 1.02])

lambda.mu <- 1.02896826

vol.new <- npoints(pp3.full)/lambda.mu
vol.scale.factor <- vol.new/volume(domain(pp3.full))
xyz.scale.factor <- (vol.scale.factor)^(1/3)

coo <- coords(pp3.full)
coo.scaled <- coo*xyz.scale.factor
#shift to center @ zero
shifts <- c((domain(pp3.full)$xrange*xyz.scale.factor)[1], 
          (domain(pp3.full)$yrange*xyz.scale.factor)[1], 
          (domain(pp3.full)$zrange*xyz.scale.factor)[1])

coo.scaled$x <- coo.scaled$x - shifts[1]
coo.scaled$y <- coo.scaled$y - shifts[2]
coo.scaled$z <- coo.scaled$z - shifts[3]

win.scaled <- box3((domain(pp3.full)$xrange*xyz.scale.factor) - shifts[1], 
                   (domain(pp3.full)$yrange*xyz.scale.factor) - shifts[2], 
                   (domain(pp3.full)$zrange*xyz.scale.factor) - shifts[3])

pp3.full.scaled <- pp3(coo.scaled$x, coo.scaled$y, coo.scaled$z, win.scaled)
marks(pp3.full.scaled) <- pos.r$mark

pp3.mg.scaled <- subset(pp3.full.scaled, marks == 'Mg1')
pp3.zn.scaled <- subset(pp3.full.scaled, marks == 'Zn1')
pp3.agg.scaled <- subset(pp3.full.scaled, marks == 'Zn1' | marks == 'Mg1')

#Set up for msa:
X <- pp3.full.scaled
mks.ind <- which(marks(pp3.full.scaled) == 'Mg1' | marks(pp3.full.scaled) == 'Zn1')
mks <- rep('B', npoints(pp3.full.scaled))
mks[mks.ind] <- 'A'
marks(X) <- mks

# Get K-functions for all aggregate combinations
to.run <- list(pp3.mg.scaled, pp3.zn.scaled, pp3.agg.scaled)
cl <- makePSOCKcluster(3)
clusterEvalQ(cl, library(spatstat))

t1 <- Sys.time()
res <- parLapply(cl, to.run, function(x){
  return(K3est(x, rmax = 35, nrval = 400, correction = 'translation'))
})
stopCluster(cl)
t2 <- Sys.time()
print(t2 - t1)

k.mg <- res[[1]]
k.zn <- res[[2]]
k.agg <- res[[3]]

# Upload 10,000 RRL data & test models on it:
load('Z:/Galen/Stephan\ APT\ Analysis/stephan_RRL_data_r35.RData')
rrls <- stephan.RRL.data
rm(stephan.RRL.data)
gc()

# Transform aggregate k-functions to T(r)
t.mg <- data.frame('r' = k.mg$r, 'T' = sqrt(k.mg$trans) - rrls[[2]])
t.zn <- data.frame('r' = k.zn$r, 'T' = sqrt(k.zn$trans) - rrls[[2]])
t.agg <- data.frame('r' = k.agg$r, 'T' = sqrt(k.agg$trans) - rrls[[2]])

# Get the metrics from each
met.mg <- as.data.frame(t(unlist(k3metrics(t.mg$r, t.mg$T, toplot = TRUE))))
names(met.mg) <- c("Km","Rm","Rdm","Rddm","Kdm")
met.zn <- as.data.frame(t(unlist(k3metrics(t.zn$r, t.zn$T, toplot = TRUE))))
names(met.zn) <- c("Km","Rm","Rdm","Rddm","Kdm")
met.agg <- as.data.frame(t(unlist(k3metrics(t.agg$r, t.agg$T, toplot = TRUE))))
names(met.agg) <- c("Km","Rm","Rdm","Rddm","Kdm")


# Upload training data to get PCA transformation:
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/200131_params.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/200131_sim.res.RData')

# unpack data
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

#test set of interest (needs to at least have the metrics)
#test.df <- test.out #from simulations
test.df <- met.agg
test.pca <- predict(pp.train, test.df[,mets])

test <- data.frame(test.df, test.pca)

#OR use test data
load('Z:/Galen/Machine\ Learning\ Files/Test\ Data/test.results.RData')
test.set <- as.data.frame(test.set)
names(test.set) <- c("cr","rho1","rho2","rb","pb","Km","Rm","Rdm","Rddm","Kdm")
test.set <- test.set[complete.cases(test.set),]
test.set$sigma <- test.set$cr * test.set$rb
test.set$rw <- (test.set$cr^4 + 6* test.set$cr^2 * test.set$sigma^2 + 3*test.set$sigma^4)/
  (test.set$cr^3 + 3*test.set$cr*test.set$sigma^2)

test.pca <- predict(pp.train, test.set[,mets])
test <- data.frame(test.set, test.pca)


# Load Model
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_rw.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_rho1.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_rho2.RData')

#Test model on simulated (or real) data:
predictors <- sapply(1:length(mets),function(x){paste('PC',toString(x),sep = '')})
out <- 'rw'

# Test models
toRound <- 6
model.res <- lapply(models, function(x){
  a <- predict(x, newdata = test[,predictors])
  #RMSE <- sqrt(mean((a - test[,out])^2))
  #plot(test[,out], a, main = x$method, xlab = "Input Simulation Value", ylab = "Estimated Value")
  #abline(0,1, col = "red", lwd = 1.25)
  print(x$method)
  print(a)
  #print(paste('RMSE: ', toString(round(RMSE,toRound)), sep = ''))
})

rm(models)
gc()

## Get confidence intervals from the simulated data
preds.rw <- predict(models[[3]], newdata = test[,predictors])
preds.rho1 <- predict(models[[3]], newdata = test[,predictors])
preds.rho2 <- predict(models[[3]], newdata = test[,predictors])

preds.all <- data.frame('rw' = test$rw, 'rho1' = test$rho1, 'rho2' = test$rho2, 
                        'rw.pred' = preds.rw, 'rho1.pred' = preds.rho1, 'rho2.pred' = preds.rho2)

perc.buffer <- 0.075
rw.buffer <- diff(range(preds.all$rw))*perc.buffer
rho1.buffer <- diff(range(preds.all$rho1))*perc.buffer
rho2.buffer <- diff(range(preds.all$rho2))*perc.buffer

preds.cut <- preds.all[abs(preds.all$rw.pred - 5.192789) < rw.buffer &
                         abs(preds.all$rho1.pred - 0.2118379) < rho1.buffer &
                         abs(preds.all$rho2.pred - 0.02494626) < rho2.buffer, ]

preds.cut$rw.diff.perc <- (preds.cut$rw.pred - preds.cut$rw)/preds.cut$rw
preds.cut$rho1.diff.perc <- (preds.cut$rho1.pred - preds.cut$rho1)/preds.cut$rho1
preds.cut$rho2.diff.perc <- (preds.cut$rho2.pred - preds.cut$rho2)/preds.cut$rho2

pd.sorted <- data.frame('rw' = sort(preds.cut$rw.diff.perc), 
                        'rho1' = sort(preds.cut$rho1.diff.perc), 
                        'rho2' = sort(preds.cut$rho2.diff.perc))

level <- 0.9
percentiles <- c((1-level)/2, level + (1-level)/2)

nobs <- nrow(pd.sorted)
inds <- round(percentiles*nobs)

CI.90 <- data.frame('rw' = 5.192789*(1+pd.sorted$rw[inds]), 
                    'rho1' = 0.2118379*(1+pd.sorted$rho1[inds]),
                    'rho2' = 0.02494626*(1+pd.sorted$rho2[inds]))
CI.90

#### Run msa algoritm with varied parameter space ####
nmins <- c(5, 10, 15, 20)
dmaxs <- seq(1, 2, 0.1)
params.to.sweep <- expand.grid(nmins, dmaxs)

msa.sweep <- lapply(1:nrow(params.to.sweep), function(i){
  print(i)
  return(msa(X, params.to.sweep[i,2], params.to.sweep[i,1], params.to.sweep[i,2], params.to.sweep[i,2], clust.mark = c('A')))
})

msa.sweep.reduced <- lapply(msa.sweep, function(x){
  if(is.na(x)){return(NA)}
  return(data.frame('rad' = x$radius, 'den' = x$den, 'bgnd.den' = x$bgnd.den))
})

save(msa.sweep.reduced, file = 'msa.sweep.RData')
save(params.to.sweep, file = 'msa.sweep.params.RData')

## Analyze the above:
load('Z:/Galen/MSA/stephan_msa.RData')
load('Z:/Galen/MSA/stephan_msa.params.RData')

msa.params.cut <- params.to.sweep
msa.cut <- msa.sweep.reduced

msa.cut.summary <- lapply(msa.cut, function(x){
  if(is.na(x)){return(data.frame('rc' = NA, 'rwT' = NA, 'rwE' = NA, 'rho' = NA))}
  sigma <- sd(x$rad)
  Rc <- mean(x$rad)
  Rw_true <- mean(x$rad*(4/3)*pi*(x$rad)^3)/mean((4/3)*pi*(x$rad)^3)
  Rw_est <- (Rc^4 + 6*Rc^2*sigma^2 + 3*sigma^4)/(Rc^3 + 3*Rc*sigma^2)
  rho1 <- mean(x$den)
  rho2 <- mean(x$bgnd.den)
  return(data.frame('rc' = Rc, 'rwT' = Rw_true, 'rwE' = Rw_est, 'rho1' = rho1, 'rho2' = rho2, 'sig'= sigma))
})

#re-arrange data for plotting
msa.all <- matrix(NA, nrow = length(msa.cut), ncol = 8)
for(i in 1:length(msa.cut)){
  if(any(is.na(msa.cut.summary[[i]]))){msa.all[i,] <- rep(NA, 8)}
  else {msa.all[i,] <- as.numeric(c(msa.params.cut[i,], msa.cut.summary[[i]]))}
}
msa.all <- as.data.frame(msa.all)
names(msa.all) <- c('Nmin', 'dmax', 'rc', 'rwT', 'rwE', 'rho1', 'rho2', 'sig')

plot(msa.all$dmax, msa.all$rwT, col = msa.all$Nmin/5, pch = 16,
     xlab = 'dmax', ylab = 'Weighted Radius')
abline(h = 5.192, col = 'black', lwd = 2)
abline(h = 5.700, col = 'black', lwd = 2, lty = 2)
abline(h = 4.915, col = 'black', lwd = 2, lty = 2)
legend(0.95, 13, legend = c('Nmin = 5', 'Nmin = 10', 'Nmin = 15', 'Nmin= 20'), col = c(1, 2, 3, 4), pch = 16, bty = 'n')
legend(1.3, 13, legend = c('ML Model Prediction', 'ML Model 90% CI'), col = 'black', lwd = 2, lty = c(1, 2), bty = 'n')

plot(msa.all$dmax, msa.all$rho1, col = msa.all$Nmin/5, pch = 16,
     xlab = 'dmax', ylab = 'Rho1',
     ylim = c(0.15, 0.7))
abline(h = 0.212, col = 'black', lwd = 2)
abline(h = 0.171, col = 'black', lwd = 2, lty = 2)
abline(h = 0.255, col = 'black', lwd = 2, lty = 2)
legend(1, 0.50, legend = c('Nmin = 5', 'Nmin = 10', 'Nmin = 15', 'Nmin= 20'), col = c(1, 2, 3, 4), pch = 16, bty = 'n')
legend(1.35, 0.73, legend = c('ML Model Prediction', 'ML Model 90% CI'), col = 'black', lwd = 2, lty = c(1, 2), bty = 'n')

plot(msa.all$dmax, msa.all$rho2, col = msa.all$Nmin/5, pch = 16,
     xlab = 'dmax', ylab = 'Rho2',
     ylim = c(0, 0.06))
abline(h = 0.0249, col = 'black', lwd = 2)
abline(h = 0.0212, col = 'black', lwd = 2, lty = 2)
abline(h = 0.0307, col = 'black', lwd = 2, lty = 2)
legend(1, 0.02, legend = c('Nmin = 5', 'Nmin = 10', 'Nmin = 15', 'Nmin= 20'), col = c(1, 2, 3, 4), pch = 16, bty = 'n')
legend(1.4, 0.02, legend = c('ML Model Prediction', 'ML Model 90% CI'), col = 'black', lwd = 2, lty = c(1, 2), bty = 'n')

plot(msa.all$dmax, msa.all$rc, col = msa.all$Nmin/5, pch = 16,
     xlab = 'dmax', ylab = 'Mean Cluster Radius (arb.)',
     ylim = c(0, 8))
legend(1, 0.5, legend = c('Nmin = 5', 'Nmin = 10', 'Nmin = 15', 'Nmin= 20'), col = c(1, 2, 3, 4), pch = 16, bty = 'n')



#### Simulate clusters on stephan's background and test the models: ####
# Note that this isn't really a very good thing to test... inhomogeneity in
# stephan's UPP could make simulating clusters on top of it strange... we should
# trust the simulated models and that the RRL will take out those potential
# inhomogeneous effects (Andrew showed that this should work)

nrand <- 8
rcp.nums <- 1:nrand
s <- 1345
set.seed(s)
test.params <- list()
r <- runif(nrand, min = 2, max = 6.5)
den <- runif(nrand, min = 0.15, max = 1) 
rbp <- runif(nrand, min = 0, max = 0.5)
gbp <- runif(nrand, min = 0, max = 0.6)
for(i in 1:nrand){
  test.params[[i]] <- c(r[i], den[i], rbp[i], gbp[i])
}
rm(r, den, rbp, gbp)


seeds <- round(runif(nrand, 1, 1e6))
run.over <- 1:nrand

cl <- makePSOCKcluster(8)
clusterEvalQ(cl, library(rapt))
clusterExport(cl, c('win.scaled', 'test.params' ,'pp3.full.scaled', 'seeds', 'rrls','frac'))

results <- parLapply(cl, run.over, function(i){
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(i),sep=""),
                   paste('~/Research/point_patterns/Final/system',toString(i),sep=""),
                   scaleUp = TRUE,newRadius = 0.5)
  over.big <- stitch.size(over, domain(over), c(115, 115, 115))
  over.big <- over.big[inside.boxx(over.big, w = win.scaled)]
  over.big$domain <- win.scaled
  
  t.p <- test.params[[i]] 
  clust <- makecluster(pp3.full.scaled, over.big, 0.5, 0.5, cr = t.p[1], 
                            pcp = frac, 
                            den = t.p[2], 
                            gb = FALSE, 
                            gbp = c(0,t.p[4]), 
                            rb = TRUE, 
                            rbp = t.p[3]*t.p[1], 
                            s = seeds[i])
  
  t.test <- anomK3est(clust[[1]], rrls[[2]], rmax = rrls[[3]], nrval = rrls[[4]])
  #plot(t.test, type = 'l', col = 'red', lwd = 2)
  
  mets.test <- as.data.frame(t(unlist(k3metrics(t.test$r, t.test$trans, toplot = FALSE))))
  names(mets.test) <- c("Km","Rm","Rdm","Rddm","Kdm")

  t.p <- as.data.frame(t(t.p))
  names(t.p) <- c('r', 'den', 'rb', 'gb')
  t.p$rw <- (t.p$r^4 + 6* (t.p$r^2) *(t.p$rb * t.p$r)^2 + 3 *(t.p$rb * t.p$r)^4)/(t.p$r^3 + 3*t.p$r*(t.p$rb * t.p$r)^2)
  
  out.full <- data.frame(t.p, mets.test)
  return(out.full)
})

stopCluster(cl)

test.out <- as.data.frame(matrix(unlist(results), nrow = nrand, byrow = TRUE))
names(test.out) <- names(results[[1]])

#from the hpc
load('Z:/Galen/Stephan\ APT\ Analysis/stephan_sims.RData')
test.out <- test.out[complete.cases(test.out),]
nrow(test.out)

# Now go up to the model loading section and test the models on these results

#
#### Exporations ####
library(ks)
contourplot <- function(estimate, kind, new = FALSE, ...){
  est.toplot <- estimate
  cont <- quantile(est.toplot$lambda.est[,kind], probs = seq(0.01, 0.99, by = 0.01))
  cont <- rev(cont)
  grid.size <- round(length(est.toplot$estimate.coords$x)^(1/3), 0)
  imap3d.est <- list(x = est.toplot$x,
                     eval.points = list(sort(unique(est.toplot$estimate.coords$x)), sort(unique(est.toplot$estimate.coords$y)), sort(unique(est.toplot$estimate.coords$z))),
                     estimate = array(est.toplot$lambda.est[,kind], dim = c(grid.size,grid.size,grid.size)),
                     H = NULL,
                     gridtype = 'linear',
                     gridded = TRUE,
                     binned = TRUE,
                     names = c("x","y","z"),
                     w = rep(1, nrow(est.toplot$x)),
                     cont = cont)
  class(imap3d.est) <- 'kde'
  if(new == TRUE){
    open3d()
  }
  plot(imap3d.est, ...)
}

a <- nndensity(pp3.temp, k = c(10,20), nx = 50, ny = 50, nz = 50, dz = 0.2, par = FALSE)

cont <- c(10, 20, 30)
k.ind <- 2
contourplot(a, 1, cont = c(cont), alphavec = c(0.3), colors = c('blue'))


# do a kde estimate (kernel estimate)
est.kde <- kde(coords(pp3.temp), gridsize = 50)
plot(est.kde, cont = c(cont), alphavec = c(0.3, 0.5, 1), colors = c('blue', 'red', 'green'))

est.t <- kde(coords(b), gridsize = 50)
plot(est.t, cont = c(cont), alphavec = c(0.3, 0.5, 1), colors = c('blue', 'red', 'green'))
