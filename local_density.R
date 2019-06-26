rm(list = ls())
gc()

# Local Density Estimates
library(rapt)
library(ks)
library(parallel)
library(plot3Drgl)
library(plot3D)
source("C:/Users/galen/Documents/Research/oled-spatstat-analysis/rpoint3_extended.R")


# Different data sets to test on
load("imap3d.diff_cold.RData")
cold.map <- imap3d.diff
load("imap3d.diff_warm.RData")
warm.map <- imap3d.diff
rm(imap3d.diff)
load("apt.pp3.60cube.irppy.cold.RData")
load("apt.pp3.60cube.irppy.warm.RData")
ir.cold <- apt.pp3.60cube.irppy.cold
ir.warm <- apt.pp3.60cube.irppy.warm
rm(apt.pp3.60cube.irppy.cold, apt.pp3.60cube.irppy.warm)
load("apt.pp3.60cube.sup.cold.RData")
load("apt.pp3.60cube.sup.warm.RData")
sup.cold <- apt.pp3.60cube.sup.cold
sup.warm <- apt.pp3.60cube.sup.warm
rm(apt.pp3.60cube.sup.cold, apt.pp3.60cube.sup.warm)

X.cold <- rpoint3(n = npoints(apt.pp3.60cube.irppy.cold), f = cold.map, im = TRUE, w = box3(c(-300,300),c(-300,300),c(100,700)))
X.warm <- rpoint3(n = npoints(apt.pp3.60cube.irppy.warm), f = warm.map, im = TRUE, w = box3(c(-300,300),c(-300,300),c(100,700)))

# Test on some simulations
under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(1),sep=""),paste('~/Research/point_patterns/Final/system',toString(1),sep=""),scaleUp = TRUE,newRadius = 0.5)
over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(2),sep=""),paste('~/Research/point_patterns/Final/system',toString(2),sep=""),scaleUp = TRUE,newRadius = 0.5)
under.big <- stitch.size(under, boxSize = c(60,60,60))
over.big <- stitch.size(over, boxSize = c(60,60,60))
clust.sim <- makecluster(under.big, over.big, 0.5, 0.5, cr = 4, pcp = 0.1, den = 1, toPlot = TRUE)

X.sim <- clust.sim[[1]]

grid.size <- 30
dz <- 0.5
kby <- 5
k <- seq(20,400, kby)

# intensity map cold
t1 <- Sys.time()
est.cold <- nndensity.pp3(X.cold, k, grid.size, grid.size, grid.size, dz, par = TRUE, cores = 8)
t2 <- Sys.time()
print(t2-t1)

#intensity map warm
t1 <- Sys.time()
est.warm <- nndensity.pp3(X.warm, k, grid.size, grid.size, grid.size, dz, par = TRUE, cores = 8)
t2 <- Sys.time()
print(t2-t1)

dz <- 5
#actual apt data cold
t1 <- Sys.time()
est.cold.real <- nndensity.pp3(ir.cold, k, grid.size, grid.size, grid.size, dz, par = TRUE, cores = 8)
t2 <- Sys.time()
print(t2-t1)

dz <- 5
#actual apt data warm
t1 <- Sys.time()
est.warm.real <- nndensity.pp3(ir.warm, k, grid.size, grid.size, grid.size, dz, par = TRUE, cores = 8)
t2 <- Sys.time()
print(t2-t1)

#simulated cluster
t1 <- Sys.time()
est.sim <- nndensity.pp3(X.sim, k, grid.size, grid.size, grid.size, dz, par = TRUE, cores = 8)
t2 <- Sys.time()
print(t2-t1)

# normalize to kde for comparison purposes (only do if measuring from kde intensity map data)
normfactor.c <- mean(apply(est.cold$lambda.est, 2, sum)[2:ncol(est.cold$lambda.est)])/sum(cold.map$estimate)
cold.map.scaled <- cold.map
cold.map.scaled$estimate <- cold.map$estimate * normfactor.c

normfactor.w <- mean(apply(est.warm$lambda.est, 2, sum)[2:ncol(est.warm$lambda.est)])/sum(warm.map$estimate)
warm.map.scaled <- warm.map
warm.map.scaled$estimate <- warm.map$estimate * normfactor.w

# for visualizations, create KDE objects (note, only works on grid)

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

kind <- 20
cont <- 15
contourplot(est.cold, kind, cont = c(cont), alphavec = c(0.3), colors = c('blue'), main = paste("nn",toString(k[kind]),sep = ""))
contourplot(est.warm, kind, new = TRUE, cont = c(cont), alphavec = c(0.3), colors = c('red'), main = paste("nn",toString(k[kind]),sep = ""))

# now we need to go through and check the performance at different values of k
actual.cold <- as.vector(cold.map.scaled$estimate)
actual.warm <- as.vector(warm.map.scaled$estimate)

check.cold <- apply(est.cold$lambda.est, 2, function(x){sum((x-actual.cold)^2)})
check.warm <- apply(est.warm$lambda.est, 2, function(x){sum((x-actual.warm)^2)})

# Plot on log scale
range <- 2:length(check.cold)
plot(k[range],log(check.cold[range]), 
     main = "Differnce KDE: NN Density Estimates", xlab = "NN #", ylab = "log(MISE)", pch = 19, col = "blue", ylim = c(-12,-10))
points(k[range],log(check.warm[range]),
     pch = 19, col = "red")
legend(300,-10.25, legend = c("cold dep.", "warm dep."), col = c("blue", "red"), pch = 19, bty = "n")

# Plot on linear scale
plot(k[range],check.cold[range], 
     main = "Differnce KDE: NN Density Estimates", xlab = "NN #", ylab = "MISE", pch = 19, col = "blue", ylim = c(0,4*10^(-5)))
points(k[range],check.warm[range],
       pch = 19, col = "red")
legend(300,4*10^(-5), legend = c("cold dep.", "warm dep."), col = c("blue", "red"), pch = 19, bty = "n")


### Histogram analysis ###
# Normalize actual data cold and warm so that they have the same overall concentration
den_warm <- npoints(ir.warm)/volume(domain(ir.warm))
den_cold <- npoints(ir.cold)/volume(domain(ir.cold))

est.warm.real$lambda.est <- est.warm.real$lambda.est * den_cold/den_warm

est.warm$lambda.est <- est.warm$lambda.est * den_cold/den_warm

# choose pattern to analyze
est.hist <- est.warm

#linear scale
e.h <- data.frame(est = unlist(est.hist$lambda.est), k = rep(k, each = nrow(est.hist$lambda.est)))
max(e.h$est)
min(e.h$est)
est_min <- 0
est_max <- 0.0007
est_by <- (est_max-est_min)/100
est_c <- cut(e.h$est, seq(est_min,est_max, est_by))
k_c <- cut(e.h$k, seq(min(k)-0.5,max(k)+0.5+kby, kby))
z <- table(est_c, k_c)
z.df <- as.data.frame(z)

#for normalization
area <- est_by*kby
tot <- area * sum(z.df$Freq)

#log on the color scale if you want (don't normalize if you do this)
z <- log(z)
z[z == -Inf] = 0

xb <- seq(est_min + 0.5*est_by, est_max - 0.5*est_by, est_by)
yb <- k

par(mar = c(5,4.5,2,4))
image2D(x = xb, y = yb, z = z, xlab = "estimate", ylab = "k", xlim = c(est_min, est_max), clab = "log(freq)",zlim = c(0,9))

open3d()
hist3Drgl(x = xb, y = yb, z = z/tot, xlab = "estimate", ylab = "k", zlab = "frequency", axes = TRUE, label = TRUE, nticks = 5, ticktype = "detailed")

# Log scale
e.h <- data.frame(est = log(unlist(est.hist$lambda.est)), k = rep(k, each = nrow(est.hist$lambda.est)))
e.h$est[is.nan(e.h$est)] <- mean(e.h$est, na.rm = TRUE)
max(e.h$est)
min(e.h$est)
est_min <- -9.5
est_max <- -8
est_by <- (est_max-est_min)/100
est_c <- cut(e.h$est, seq(est_min,est_max, est_by))
k_c <- cut(e.h$k, seq(min(k)-0.5,max(k)+0.5+kby, kby))
z <- table(est_c, k_c)
z.df <- as.data.frame(z)

#for normalization
area <- est_by*kby
tot <- area * sum(z.df$Freq)

xb <- seq(est_min + 0.5*est_by, est_max - 0.5*est_by, est_by)
yb <- k

par(mar = c(5,4.5,2,4))
image2D(x = xb, y = yb, z = z/tot, xlab = "log(estimate)", ylab = "k", xlim = c(est_min,est_max), zlim = c(0,0.02))

open3d()
hist3Drgl(x = xb, y = yb, z = z/tot, 
          xlab = "log(estimate)", ylab = "k", zlab = "frequency", border = "black", axes = TRUE, label = TRUE, nticks = 5, ticktype = "detailed")


### global intensity removal for the actual data ###
grid.size <- 30
kby <- 5
k <- seq(5,15, kby)
est.cold <- nncrossden.pp3(ir.cold, sup.cold, k, grid.size, grid.size, grid.size, at.points = FALSE, nsplit = 8*5, os = "windows")
est.warm <- nncrossden.pp3(ir.warm, sup.warm, k = 10, grid.size, grid.size, grid.size, at.points = FALSE, os = "windows")

hist(est.cold[[1]]$nn10, breaks = seq(0,0.0005,length.out = 100), ylim = c(0,70))
hist(est.warm[[1]]$nn10, breaks = seq(0,0.0005,length.out = 100), ylim = c(0,70))

# visualize this^ data run on the HPC
load("C:/Users/galen/Documents/Research/local_intensity_estimate/est.cold.full.RData")
load("C:/Users/galen/Documents/Research/local_intensity_estimate/est.warm.full.RData")
k <- est.cold$k
kby <- k[2]-k[1]



