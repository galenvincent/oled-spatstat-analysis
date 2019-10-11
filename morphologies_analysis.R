library(rapt)
library(parallel)


# Run K-rrl on each of them
sim <- list(1, 2, 3, 4)

sim.res <- lapply(sim, function(x){
  print('starting...')
  nsim <- 8*3
  nrank <- 1
  sim <- switch(x, expression(morph_lamellar(216000, 0.1, c(1, 0, 0), 4)),
                expression(morph_rods(216000, 0.1, 'z', c(3, 3), 'hexagonal',
                                      rcp.path = "C:/Users/galen/Documents/Research/point_patterns/2D/Final")),
                expression(morph_gb(216000, 0.1, rcp.rad = 0.2, 
                                    rcp.path = "C:/Users/galen/Documents/Research/point_patterns/Final",
                                    rcp.number = 'rand')),
                expression(morph_gyroid(216000, 0.1, 1)))
  cl <- makePSOCKcluster(detectCores())
  clusterEvalQ(cl, library(rapt))
  X <- rpoispp3(1000, box3(c(0,1), c(0,1), c(0,1))) # unimportant
  clusterExport(cl, c('X','nrank','sim'), envir = environment())
  
  pNsim <- rle(cut(seq_len(nsim), length(cl), labels = FALSE))$lengths
  env <- parLapply(cl, pNsim, function(n) {
    envelope(X, fun=K3est, nsim=n, nrank=nrank, funargs = list('rmax' = 1, nrval = 200, correction = 'translation'),
             simulate=sim, savefuns=TRUE, verbose=FALSE)
  })
  
  stopCluster(cl)
  
  po <- do.call(pool, c(env, savefuns=TRUE))
  dat <- envelope(po, nrank=nrank, savefuns=TRUE)
  
  sims <- as.data.frame(attr(dat, 'simfuns'))
  sims[,2:length(sims)] <- sims[,2:length(sims)] - (4/3)*pi*sims[,1]^3
  print('finished...')
  return(sims)
})

names(sim.res) <- c('lamellar', 'rods', 'gb', 'gyroid')

save(sim.res, file = 'sim.res.RData')


# From HPC
load('sim.res.gyroid.RData')
lam <- sim.res$lamellar
rod <- sim.res$rods
gb <- sim.res$gb
gyr <- sim.res[[1]]

par(mar = c(4, 4.5, 2, 2), mgp = c(2.25, 1, 0))
envPlot(lam, ylim = c(-0.15, 0.15), percentiles = c(0.99, 0.90, 0.75))
envPlot(rod, ylim = c(-0.15, 0.15), percentiles = c(0.99, 0.90, 0.75))
envPlot(gb, ylim = c(-0.15, 0.15), percentiles = c(0.99, 0.90, 0.75))
envPlot(gyr, ylim = c(-0.15, 0.15), percentiles = c(0.99, 0.90, 0.75), leg = F)
