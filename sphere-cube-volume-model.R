# We're going to do some modeling for sphere-cube interactions
library(rapt)
library(caret)

# Create domain (60x60x60 to match simulations)
sl <- 60
box.domain <- box3(c(0,sl), c(0,sl), c(0,sl))

# Generate a bunch of points in the region of interest to estimate volume:
vp.x <- seq(0, 60, length.out = 200)
vp.y <- seq(0, 60, length.out = 200)
vp.z <- seq(0, 60, length.out = 200)
vol.coo <- expand.grid(vp.x, vp.y, vp.z)
names(vol.coo) <- c('x', 'y', 'z')
vol.pp3 <- pp3(vol.coo$x, vol.coo$y, vol.coo$z, box.domain)
vol.lambda <- npoints(vol.pp3)/(60^3)

# Create point pattern around the outside of the box
a <- seq(0, 60, length.out = 400)
b <- seq(0, 60, length.out = 400)
d <- expand.grid(a, b)
possibles <- list(c(T, T, F), c(T, F, T), c(F, T, T))
surfs <- list()
for(i in 1:3){
  q <- nrow(d)
  surfs[[i]] <- matrix(NA, nrow = q*2, ncol = 3)
  surfs[[i]][1:q, which(possibles[[i]])] <- as.matrix(d)
  surfs[[i]][1:q, -which(possibles[[i]])] <- rep(0, q)
  
  surfs[[i]][(q + 1):(2*q), which(possibles[[i]])] <- as.matrix(d)
  surfs[[i]][(q + 1):(2*q), -which(possibles[[i]])] <- rep(60, q)
}

surfs.pp3s <- lapply(surfs, function(x){pp3(x[,1], x[,2], x[,3], box.domain)})

rm(a, b, d, possibles, surfs)
gc()
n <- 250

for(k in 20:40){
  print(k)
  # Generate n random cluster radii
  cr <- runif(n, 2, 12)
  
  # Generate n random cluster centers within 
  # striking distance of the edge of the domain
  centers <- matrix(NA, nrow = n, ncol = 3)
  for(i in 1:n){
    centers[i,] <- runif(3, -cr[i], sl+cr[i])
  }
  centers.df <- as.data.frame(centers)
  names(centers.df) <- c('x', 'y', 'z')
  
  bdist <- bdist.complex(centers.df, box.domain)
  good.inds <- !apply(bdist > cr, 1, all)
  bdist <- bdist[good.inds,]
  centers.df <- centers.df[good.inds,]
  cr <- cr[good.inds]
  
  maxes <- apply(centers.df, 2, max)
  mins <- apply(centers.df, 2, min)
  centers.box <- box3(c(mins[1], maxes[1]), c(mins[2], maxes[2]), c(mins[3], maxes[3]))
  centers.pp3 <- pp3(centers.df$x, centers.df$y, centers.df$z, centers.box)
  bdist.min <- apply(bdist, 1, function(x){min(abs(x))})
  io <- as.numeric(inside.boxx(centers.pp3, w = box.domain))
  
  # Approx volume around each cluster center with given radius
  nnR <- crosspairs.pp3(centers.pp3, vol.pp3, rmax = max(cr), what = 'ijd', neat = FALSE, distinct = TRUE, twice = FALSE)
  nnR.split <- list()
  nnR.split$d <- split(nnR$d, nnR$i, drop=FALSE)
  nnR.split$j <- split(nnR$j, nnR$i, drop=FALSE)
  nnR.split$i <- as.numeric(attr(nnR.split$d, 'name'))
  
  vol.totals <- sapply(1:length(nnR.split$i), function(f){length(nnR.split$j[[f]][nnR.split$d[[f]] < cr[nnR.split$i[[f]]]])/vol.lambda})
  
  # Put it all together
  sph.data <- matrix(NA, nrow = length(cr), ncol = 7)
  colnames(sph.data) <- c('dx', 'dy', 'dz', 'dmin', 'r', 'vol', 'io')
  sph.data[,1:3] <- as.matrix(bdist)
  sph.data[,4] <- bdist.min
  sph.data[,5] <- cr
  sph.data[nnR.split$i,6] <- vol.totals
  sph.data[-nnR.split$i, 6] <- 0
  sph.data[, 7] <- io
  
  if(k == 1){smoosh <- sph.data}else{smoosh <- rbind(smoosh, sph.data)}
}

#### Try something new
for(k in 1:20){
  print(k)
  # Generate n random cluster radii
  cr <- runif(n, 2, 12)
  
  # Generate n random cluster centers within 
  # striking distance of the edge of the domain
  centers <- matrix(NA, nrow = n, ncol = 3)
  for(i in 1:n){
    centers[i,1:2] <- runif(2, -cr[i], sl+cr[i])
    centers[i,3] <- runif(1, sl-cr[i], sl+cr[i])
  }
  centers.df <- as.data.frame(centers)
  names(centers.df) <- c('x', 'y', 'z')
  maxes <- apply(centers.df, 2, max)
  mins <- apply(centers.df, 2, min)
  centers.box <- box3(c(mins[1], maxes[1]), c(mins[2], maxes[2]), c(mins[3], maxes[3]))
  centers.pp3 <- pp3(centers.df$x, centers.df$y, centers.df$z, centers.box)
  nnc <- lapply(surfs.pp3s, function(x){nncross(centers.pp3, x, what = 'dist')})
  nnc.mat <- matrix(unlist(nnc), ncol = 3, byrow = F)
  colnames(nnc.mat) <- c('dxy', 'dxz', 'dyz')
  
  nnc.xyzio <- t(apply(centers.df, 1, function(x){
    xyzio <- c(-1, -1, -1)
    if(x[1] > 0 & x[1] < 60){xyzio[1] <- 1}
    if(x[2] > 0 & x[2] < 60){xyzio[2] <- 1}
    if(x[3] > 0 & x[3] < 60){xyzio[3] <- 1}
    return(xyzio)
  }))
  nnc.mat <- nnc.mat*nnc.xyzio[,c(3, 2, 1)]
  
  nnc.mins <- unlist(apply(nnc.mat, 1, function(x){min(x[which(abs(x) == min(abs(x)))])}))
  #nnc.mins <- unlist(apply(nnc.mat, 1, function(x){min(abs(x))}))
  nnc.numtouch <- apply(abs(nnc.mat) < cr, 1, sum)
  
  good.inds <- which(nnc.numtouch > 1)
  centers.df <- centers.df[good.inds,]
  cr <- cr[good.inds]
  nnc.mat <- nnc.mat[good.inds,]
  nnc.mins <- nnc.mins[good.inds]
  nnc.numtouch <- nnc.numtouch[good.inds]
  rm(nnc)
  
  maxes <- apply(centers.df, 2, max)
  mins <- apply(centers.df, 2, min)
  centers.box <- box3(c(mins[1], maxes[1]), c(mins[2], maxes[2]), c(mins[3], maxes[3]))
  centers.pp3 <- pp3(centers.df$x, centers.df$y, centers.df$z, centers.box)
  
  
  # Approx volume around each cluster center with given radius
  nnR <- crosspairs.pp3(centers.pp3, vol.pp3, rmax = max(cr), what = 'ijd', neat = FALSE, distinct = TRUE, twice = FALSE)
  nnR.split <- list()
  nnR.split$d <- split(nnR$d, nnR$i, drop=FALSE)
  nnR.split$j <- split(nnR$j, nnR$i, drop=FALSE)
  nnR.split$i <- as.numeric(attr(nnR.split$d, 'name'))
  
  vol.totals <- sapply(1:length(nnR.split$i), function(f){length(nnR.split$j[[f]][nnR.split$d[[f]] < cr[nnR.split$i[[f]]]])/vol.lambda})
  
  io <- as.numeric(inside.boxx(centers.pp3, w = box.domain))

  # Put it all together
  sph.data <- matrix(NA, nrow = length(cr), ncol = 8)
  colnames(sph.data) <- c('dxy', 'dxz', 'dyz', 'dmin', 'r', 'vol', 'nts', 'io')
  sph.data[,1:3] <- nnc.mat
  sph.data[,4] <- nnc.mins
  sph.data[,5] <- cr
  sph.data[nnR.split$i,6] <- vol.totals
  sph.data[-nnR.split$i, 6] <- 0
  sph.data[, 7] <- nnc.numtouch
  sph.data[, 8] <- io

  if(k == 1){smoosh <- sph.data}else{smoosh <- rbind(smoosh, sph.data)}
  
}

sph.data.df <- as.data.frame(smoosh)
rm(nnR, nnR.split)
gc()

# Plot as function of dmin and r
rgl::plot3d(sph.data.df$dmin, sph.data.df$r, sph.data.df$vol)

# Fit a ML model to the data
predictors <- c('dx', 'dy', 'dz', 'dmin', 'r', 'io')
model <- train(sph.data.df[,predictors], sph.data.df$vol, 
                     method = "brnn", 
                     trControl = trainControl("cv", number = 10),
               tuneGrid = data.frame('neurons' = c(9)))

library(parallel)
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
RFmodel <- train(sph.data.df[,predictors], sph.data.df$vol,
                 method = "parRF", trControl = trainControl("cv", number = 10))
stopCluster(cl)
registerDoSEQ()

BRNNmodel <- model

save(RFmodel, file = 'sphere-cube-RFmodel.RData')
save(BRNNmodel, file = 'sphere-cube-BRNNmodel.RData')

# Helper function
bdist.complex <- function(X.df, win){
  
  x <- X.df$x
  y <- X.df$y
  z <- X.df$z
  
  xmin <- min(win$xrange)
  xmax <- max(win$xrange)
  ymin <- min(win$yrange)
  ymax <- max(win$yrange)
  zmin <- min(win$zrange)
  zmax <- max(win$zrange)
  
  result <- data.frame(x = pmin.int(x - xmin, xmax - x), y =  pmin.int(y - ymin, ymax - y), z = pmin.int(z - zmin, zmax - z))
  
  return(result)
}
