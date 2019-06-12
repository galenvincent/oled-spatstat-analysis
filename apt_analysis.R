# RRL for APT data
library(rapt)
library(parallel)
library(data.table)

# R44_03122-v02 - COLD -----------------------------------------------------------
# Volume selection

oled.ato <- readATO('Z:/LEAP\ Reconstructions/R44_03122-v02.ato')
oled.pp3 <- createSpat(oled.ato)

oled.irppy <- with(oled.ato, oled.pp3[(mass > 652.016 & mass < 660.269) |
                                             (mass > 326.311 & mass < 328.715),])
oled.mcbp <- with(oled.ato, oled.pp3[(mass > 482.545 & mass < 491.163) |
                                            (mass > 241.366 & mass < 245.787),])

oled.in.irppy <- list()
oled.in.mcbp <- list()
oled.in.sup <- list()

# pick a domain
#xrange <- list(domain(oled.pp3)$xrange)
#yrange <- list(domain(oled.pp3)$yrange)
#zrange <- list(domain(oled.pp3)$zrange)

#xrange <- list(c(-300,300))
#yrange <- list(c(-300,300))
#zrange <- list(c(100,700))

xrange <- list(domain(oled.pp3)$xrange,c(-500,500),c(-300,300))
yrange <- list(domain(oled.pp3)$yrange,c(-500,500),c(-300,300))
zrange <- list(domain(oled.pp3)$zrange,c(0,1000),c(100,700))

#xrange = list(c(-300,0), c(-300,0), c(-300,0), c(-300,0), c(0,300), c(0,300), c(0,300), c(0,300))
#yrange = list(c(-300,0), c(0,300), c(-300,0), c(0,300), c(-300,0), c(0,300), c(-300,0), c(0,300))
#zrange = list(c(100,400), c(100,400), c(400,700), c(400,700), c(100,400), c(100,400), c(400,700), c(400,700))


for(i in 1:length(xrange)){
  oled.in.irppy[[i]] <- oled.irppy[inside.pp3(oled.irppy,
                                                   box3(xrange = xrange[[i]], yrange = yrange[[i]], zrange = zrange[[i]]))]
  oled.in.mcbp[[i]] <- oled.mcbp[inside.pp3(oled.mcbp,
                                                 box3(xrange = xrange[[i]], yrange = yrange[[i]], zrange = zrange[[i]]))]
  
  oled.in.irppy[[i]]$domain <- box3(xrange = xrange[[i]], yrange = yrange[[i]], zrange = zrange[[i]])
  oled.in.mcbp[[i]]$domain <- box3(xrange = xrange[[i]], yrange = yrange[[i]], zrange = zrange[[i]])
  oled.in.sup[[i]] <- superimpose.pp3(oled.in.mcbp[[i]], oled.in.irppy[[i]])
  marks(oled.in.sup[[i]]) <- c(rep('mCBP', npoints(oled.in.mcbp[[i]])),
                               rep('Ir(ppy)3', npoints(oled.in.irppy[[i]])))
  class(oled.in.sup[[i]]) <- c('pp3','ppx')
  oled.in.sup[[i]]$domain <- box3(xrange = xrange[[i]], yrange = yrange[[i]], zrange = zrange[[i]])
}


#visualize selection
plot3d.pp3(percentSelect(0.025, oled.irppy))
for(i in 1:length(xrange)){
  plot3d.pp3(oled.in.irppy[[i]], col = i, add = TRUE)
}

perc <- list()
for(i in 1:length(xrange)){
  perc[[i]] <- npoints(oled.in.irppy[[i]])/npoints(oled.in.sup[[i]])
}

a <- list()
t1 <- Sys.time()
for(i in 1:length(xrange)){
  a[[i]] <- pK3est(perc[[i]], oled.in.sup[[i]], 200, nrval=200, correction="trans",anom=TRUE, sorted = FALSE)
}
t2 <- Sys.time()
print('Time to complete: ')
print(t2-t1)


b <- list()
for(i in 1:length(xrange)){
  b[[i]] <- anomK3est(oled.in.irppy[[i]], a[[i]][[2]], rmax = a[[i]][[3]], nrval = a[[i]][[4]], correction = "trans")
  b[[i]]$r <- b[[i]]$r/10
  a[[i]][[1]][,1] <- a[[i]][[1]][,1]/10 
}

par(mar = c(4,4.5,2,2))
envPlot(a[[6]][[1]], percentiles = c(0.9, 0.8, 0.7), ylim = c(-70,70))
for(i in 1:length(xrange)){
  lines(b[[i]]$r, b[[i]]$trans, lwd = 2, col = i)
}
legend(0,40, c("1","2","3","4","5","6"), col = i, lwd = 2, lty = 1)


# R44_03151-v02 - WARM -----------------------------------------------------------
rm(list = ls())
gc()

oled.ato <- readATO('Z:/LEAP\ Reconstructions/R44_03151-v02.ato')
oled.pp3 <- createSpat(oled.ato)
oled.irppy <- with(oled.ato, oled.pp3[(mass > 652.016 & mass < 660.269) |
                                        (mass > 326.311 & mass < 328.715),])
oled.mcbp <- with(oled.ato, oled.pp3[(mass > 481.545 & mass < 491.163) |
                                       (mass > 241.366 & mass < 245.787),])

oled.in.irppy <- list()
oled.in.mcbp <- list()
oled.in.sup <- list()

#volume selection
#xrange <- list(domain(oled.pp3)$xrange)
#yrange <- list(domain(oled.pp3)$yrange)
#zrange <- list(domain(oled.pp3)$zrange)

#xrange <- list(c(-300,300))
#yrange <- list(c(-300,300))
#zrange <- list(c(100,700))

xrange <- list(domain(oled.pp3)$xrange,c(-500,500),c(-300,300))
yrange <- list(domain(oled.pp3)$yrange,c(-500,500),c(-300,300))
zrange <- list(domain(oled.pp3)$zrange,c(0,1000),c(100,700))

#xrange = list(c(-300,0), c(-300,0), c(-300,0), c(-300,0), c(0,300), c(0,300), c(0,300), c(0,300))
#yrange = list(c(-300,0), c(0,300), c(-300,0), c(0,300), c(-300,0), c(0,300), c(-300,0), c(0,300))
#zrange = list(c(100,400), c(100,400), c(400,700), c(400,700), c(100,400), c(100,400), c(400,700), c(400,700))

for(i in 1:length(xrange)){
  oled.in.irppy[[i]] <- oled.irppy[inside.pp3(oled.irppy,
                                              box3(xrange = xrange[[i]], yrange = yrange[[i]], zrange = zrange[[i]]))]
  oled.in.mcbp[[i]] <- oled.mcbp[inside.pp3(oled.mcbp,
                                            box3(xrange = xrange[[i]], yrange = yrange[[i]], zrange = zrange[[i]]))]
  
  oled.in.irppy[[i]]$domain <- box3(xrange = xrange[[i]], yrange = yrange[[i]], zrange = zrange[[i]])
  oled.in.mcbp[[i]]$domain <- box3(xrange = xrange[[i]], yrange = yrange[[i]], zrange = zrange[[i]])
  oled.in.sup[[i]] <- superimpose.pp3(oled.in.mcbp[[i]], oled.in.irppy[[i]])
  marks(oled.in.sup[[i]]) <- c(rep('mCBP', npoints(oled.in.mcbp[[i]])),
                               rep('Ir(ppy)3', npoints(oled.in.irppy[[i]])))
  class(oled.in.sup[[i]]) <- c('pp3','ppx')
  oled.in.sup[[i]]$domain <- box3(xrange = xrange[[i]], yrange = yrange[[i]], zrange = zrange[[i]])
}

#visualize selection
plot3d.pp3(percentSelect(0.025, oled.irppy))
for(i in 1:length(xrange)){
  plot3d.pp3(oled.in.irppy[[i]], col = i, add = TRUE)
}

perc <- list()
for(i in 1:length(xrange)){
  perc[[i]] <- npoints(oled.in.irppy[[i]])/npoints(oled.in.sup[[i]])
}

a <- list()
t1 <- Sys.time()
for(i in 1:length(xrange)){
  a[[i]] <- pK3est(perc[[i]], oled.in.sup[[i]], 200, nrval=200, correction="trans",anom=TRUE, sorted = FALSE)
}
t2 <- Sys.time()
print('Time to complete: ')
print(t2-t1)

b <- list()
for(i in 1:length(xrange)){
  b[[i]] <- anomK3est(oled.in.irppy[[i]], a[[i]][[2]], rmax = a[[i]][[3]], nrval = a[[i]][[4]], correction = "trans")
  b[[i]]$r <- b[[i]]$r/10
  a[[i]][[1]][,1] <- a[[i]][[1]][,1]/10 
}

par(mar = c(4,4.5,2,2))
envPlot(a[[6]][[1]], percentiles = c(0.9, 0.8, 0.7), ylim = c(-70,70))
for(i in 1:length(xrange)){
  lines(b[[i]]$r, b[[i]]$trans, lwd = 2, col = i)
}
legend(0,40, c("1","2","3","4","5","6"), col = 1:6, lwd = 2, lty = 1)


# 60x60x60 analysis -------------------------------------------------------

rm(list = ls())
gc()

# R44_03151-v02
data.r <- fread("C:/Users/galen/Documents/Research/apt_data_analysis/R44_03151-v02.csv", select = 1)
data.rrl <- fread("C:/Users/galen/Documents/Research/apt_data_analysis/R44_03151-v02.csv")
data.toSub <- fread("C:/Users/galen/Documents/Research/apt_data_analysis/R44_03151-v02_toSub.csv")
data.res <- fread("C:/Users/galen/Documents/Research/apt_data_analysis/R44_03151-v02_result.csv")


# R44_03122-v02
data.r <- fread("C:/Users/galen/Documents/Research/apt_data_analysis/R44_03122-v02.csv", select = 1)
data.rrl <- fread("C:/Users/galen/Documents/Research/apt_data_analysis/R44_03122-v02.csv")
data.toSub <- fread("C:/Users/galen/Documents/Research/apt_data_analysis/R44_03122-v02_toSub.csv")
data.res <- fread("C:/Users/galen/Documents/Research/apt_data_analysis/R44_03122-v02_result.csv")


#rescale r to  from angstrom
data.r <- data.r/10
data.rrl[,1] <- data.rrl[,1]/10
data.res[,1] <- data.res[,1]/10


envPlot(data.rrl, c(0.999, 0.99, 0.97), ylim = c(-165,165), xlim = c(0,40))
lines(data.res, col = "black", lwd = 2)
lines(b$r/10, b$trans, col = "red", lwd = 2, lty = 2)


# Envelope smoothness analysis --------------------------------------------

rm(list = ls())
gc()
ex.rrl <- as.matrix(fread("C:/Users/galen/Documents/Research/apt_data_analysis/R44_03122-v02.csv"))
ex.rrl[,1] <- ex.rrl[,1]/10

nums <- c(2500, 5000, 10000, 50000)

par(mfrow = c(2,2), mar = c(2.5, 2, 1.5, 2))
for(i in 1:length(nums)){
  samp <- as.numeric(sample(2:50001, nums[i]))
  temp.rrl <- cbind(ex.rrl[,1], ex.rrl[,samp])
  envPlot(temp.rrl, ylim = c(-100,100), leg = FALSE, main = paste(toString(nums[i])))
}


cols <- c(rgb(1,0,0,alpha = 0.4), rgb(0,0,1,alpha = 0.4), rgb(0,1,0,alpha = 0.4), rgb(1,0,1,alpha = 0.4))
par(mfrow = c(1,1))
for(i in 1:4){
  samp <- as.numeric(sample(2:50001, nums[i]))
  temp.rrl <- cbind(ex.rrl[,1], ex.rrl[,samp])
  if(i != 1){
    par(new = TRUE)
  }  
  envPlot(temp.rrl, percentiles = c(0.999), ylim = c(-100,100), leg = FALSE, col = cols[i])
  if(i == length(nums)){
    par(new = FALSE)
    break
  }
  n <- readline(prompt="Enter to continue: ")
}



# Density Mapping ---------------------------------------------------------
### 2D Density Mapping
coo.irppy <- coords(oled.in.irppy[[1]])
coo.sup <- coords(oled.in.sup[[1]])

xyplane.irppy <- which(abs(coo.irppy$y) < 10)
xyplane.sup <- which(abs(coo.sup$y) < 10)

oled.in.irppy.csec <- ppp(coo.irppy$x[xyplane.irppy], coo.irppy$z[xyplane.irppy], 
                          window = owin(c(min(coo.irppy$x[xyplane.irppy]), max(coo.irppy$x[xyplane.irppy])), 
                                        c(min(coo.irppy$z[xyplane.irppy]), max(coo.irppy$z[xyplane.irppy]))))
oled.in.sup.csec <- ppp(coo.sup$x[xyplane.sup], coo.sup$z[xyplane.sup], 
                          window = owin(c(min(coo.sup$x[xyplane.sup]), max(coo.sup$x[xyplane.sup])), 
                                        c(min(coo.sup$z[xyplane.sup]), max(coo.sup$z[xyplane.sup]))))



plot(oled.in.irppy.csec$x, oled.in.irppy.csec$y, asp = 1)
plot(oled.in.sup.csec$x, oled.in.sup.csec$y, asp = 1)

#dimyx parameter controls grid spacing
imap.irppy <- density.ppp(oled.in.irppy.csec, sigma = 20)
imap.sup <- density.ppp(oled.in.sup.csec, sigma = 20)

par(mfrow = c(2,1), mar = c(1,1,1,2))
plot(imap.irppy)
plot(imap.sup)

#subtract the two density maps
# first, normalize the sets
max.irppy <- max(imap.irppy$v)
max.sup <- max(imap.sup$v)

imap.n.irppy <- imap.irppy$v/max.irppy
imap.n.sup <- imap.sup$v/max.sup

imap.diff <- as.im(imap.n.sup - imap.n.irppy)
graphics.off()
plot(imap.diff)



### 3D Density Mapping
library(ks)
imap3d.irppy <- kde(coords(oled.in.irppy[[2]]), gridsize = 50, xmin = c(-500,-500,0), xmax = c(500,500,1000))
imap3d.sup <- kde(coords(oled.in.sup[[2]]), gridsize = 50, xmin = c(-500,-500,0), xmax = c(500,500,1000))

#plot(imap3d.irppy)
#for(i in 1:8){
#  plot3d.pp3(percentSelect(0.1,oled.in.irppy[[i]]), col = i, add = TRUE)
#}

# look at the difference in densities in 3D
imap3d.n.irppy <- imap3d.irppy$estimate/max(imap3d.irppy$estimate)
imap3d.n.sup <- imap3d.sup$estimate/max(imap3d.sup$estimate)
imap3d.diff.raw <- imap3d.n.irppy - imap3d.n.sup
imap3d.diff.shifted <- imap3d.diff.raw + abs(min(imap3d.diff.raw))
imap3d.diff.scaled <- imap3d.diff.shifted*(max(imap3d.irppy$estimate)/max(imap3d.diff.shifted))

# Here, we create a subsetted kde object which only keeps estimates within the 60x60x60 cube (to avoid edge effects)
good <- imap3d.irppy$x$x > -300 & imap3d.irppy$x$x < 300 &
        imap3d.irppy$x$y > -300 & imap3d.irppy$x$y < 300 &
        imap3d.irppy$x$z > 100 & imap3d.irppy$x$z < 700
x.new <- imap3d.irppy$x[good,]

eval.points.new <- list()
eval.points.new[[1]] <- imap3d.irppy$eval.points[[1]][imap3d.irppy$eval.points[[1]] > -300 & imap3d.irppy$eval.points[[1]] < 300]
eval.points.new[[2]] <- imap3d.irppy$eval.points[[2]][imap3d.irppy$eval.points[[2]] > -300 & imap3d.irppy$eval.points[[2]] < 300]
eval.points.new[[3]] <- imap3d.irppy$eval.points[[3]][imap3d.irppy$eval.points[[3]] > 100 & imap3d.irppy$eval.points[[3]] < 700]

estimate.new <- imap3d.diff.scaled[11:40,11:40,11:40]

cont <- quantile(estimate.new, probs = seq(0.01, 0.99, by = 0.01))
cont <- rev(cont)
imap3d.diff <- list(x = x.new,
                    eval.points = eval.points.new,
                    estimate = estimate.new,
                    H = imap3d.irppy$H,
                    gridtype = 'linear',
                    gridded = TRUE,
                    binned = TRUE,
                    names = c("x","y","z"),
                    w = rep(1, nrow(x.new)),
                    cont = cont)

class(imap3d.diff) <- 'kde'
#
cont <- quantile(imap3d.diff.scaled, probs = seq(0.01, 0.99, by = 0.01))
cont <- rev(cont)
imap3d.diff <- list(x = imap3d.irppy$x,
                    eval.points = imap3d.irppy$eval.points,
                    estimate = imap3d.diff.scaled,
                    H = imap3d.irppy$H,
                    gridtype = 'linear',
                    gridded = TRUE,
                    binned = TRUE,
                    names = c("x","y","z"),
                    w = rep(1, nrow(imap3d.irppy$x)),
                    cont = cont)
class(imap3d.diff) <- 'kde'


# plot(imap3d.diff, xlim = c(-300,300), ylim = c(-300,300), zlim = c(100,700))

# printBox <- function(x,y,z,x1,y1,z1, ...) {
#   mycube <- scale3d(cube3d(...),x1,y1,z1)
#   return(translate3d(mycube,x,y,z))
# }
# c3d <- printBox(0,0,400,300,300,300, color="blue", alpha=0.3)  # nothing happens yet
# c3d   # Look at structure

#note that cont = c(20) means the upper 20% density range
save(imap3d.diff, file = "imap3d.diff_warm.RData")

plot(imap3d.diff, cont = c(20), alphavec = c(0.25))
rm(imap3d.irppy, imap3d.sup)
gc()


load("imap3d.diff_cold.RData")
im_cold <- imap3d.diff
plot(imap3d.diff, cont = c(20), alphavec = c(0.25))
load("imap3d.diff_warm.RData")
im_warm <- imap3d.diff
open3d()
plot(imap3d.diff, cont = c(20), alphavec = c(0.25))

# Compare the two kernel densties with hypthesis test
kde.test(im_cold$x, im_warm$x, binned = TRUE)


source("C:/Users/galen/Documents/Research/oled-spatstat-analysis/rpoint3_extended.R")
# sim.diff <- rpoint3(n = npoints(oled.in.irppy[[3]]), f = imap3d.diff, im = TRUE, w = box3(c(-300,300),c(-300,300),c(100,700)))
# sim.rand <- rpoint3(n = npoints(oled.in.sup[[3]]), win = domain(sim.diff))
# 
# perc <- npoints(sim.diff)/npoints(sim.rand)
# 
# t1 <- Sys.time()
# a <- pK3est(perc, sim.rand, nEvals = 25 ,correction = "trans", rmax = 300, nrval = 200, anom = TRUE)
# t2 <- Sys.time()
# print(t2-t1)
# 
# obs <- anomK3est(sim.diff, a[[2]], rmax = 300, nrval = 200)
# #
sim.diff <- rpoint3(n = npoints(oled.in.irppy[[2]]), f = imap3d.diff, im = TRUE, w = box3(c(-500,500),c(-500,500),c(0,1000)))

sim.rand <- rpoint3(n = npoints(oled.in.sup[[2]]), win = domain(sim.diff))

win.fin <- box3(c(-300,300),c(-300,300),c(100,700))
sim.diff.full <- sim.diff[inside.boxx(sim.diff, w = win.fin)]$data
sim.rand.full <- sim.rand[inside.boxx(sim.rand, w = win.fin)]$data

sim.diff.cut <- pp3(sim.diff.full$x, sim.diff.full$y, sim.diff.full$z, win = win.fin) 
sim.rand.cut <- pp3(sim.rand.full$x, sim.rand.full$y, sim.rand.full$z, win = win.fin) 

perc <- npoints(sim.diff.cut)/npoints(sim.rand.cut)

t1 <- Sys.time()
a <- pK3est(perc, sim.rand.cut, nEvals = 100 ,correction = "trans", rmax = 300, nrval = 200, anom = TRUE)
t2 <- Sys.time()
print(t2-t1)

obs <- anomK3est(sim.diff.cut, a[[2]], rmax = 300, nrval = 200)

#fix r values to nm
a[[1]][,1] <- a[[1]][,1]/10
obs$r <- obs$r/10

envPlot(a[[1]], percentiles = c(0.95), ylim = c(-125,125), leg = FALSE, main = "Warm Dep.", colors = "red")
lines(obs$r, obs$trans, lwd = 2, col = "black")
legend(0,125, c("95% RRL AI", "Observed"), col = c("red", "black"), lwd = c(10, 2), lty = 1, bty = "n")

# cross section density analysis
coo.a <- coords(sim.diff.cut)
xyplane.a <- which(abs(coo.a$y) < 20)
a.csec <- ppp(coo.a$x[xyplane.a], coo.a$z[xyplane.a], 
                          window = owin(c(min(coo.a$x[xyplane.a]), max(coo.a$x[xyplane.a])), 
                                        c(min(coo.a$z[xyplane.a]), max(coo.a$z[xyplane.a]))))
plot(a.csec$x, a.csec$y, asp = 1)
imap.a <- density.ppp(a.csec, sigma = 20)
par(mfrow = c(2,1), mar = c(1,1,1,2))
plot(imap.irppy)
plot(imap.a)


# Hypothesis Testing ------------------------------------------------------

# Load Irppy 60x60x60 apt data
load("apt.pp3.60cube.irppy.cold.RData")
load("apt.pp3.60cube.irppy.warm.RData")
apt.pp3.ir.c <- apt.pp3.60cube.irppy.cold
apt.pp3.ir.w <- apt.pp3.60cube.irppy.warm
rm(apt.pp3.60cube.irppy.cold, apt.pp3.60cube.irppy.warm)

# Load superposition 60x60x60 apt data
load("apt.pp3.60cube.sup.cold.RData")
load("apt.pp3.60cube.sup.warm.RData")
apt.pp3.s.c <- apt.pp3.60cube.sup.cold
apt.pp3.s.w <- apt.pp3.60cube.sup.warm
rm(apt.pp3.60cube.sup.cold, apt.pp3.60cube.sup.warm)

# Load difference 60x60x60 data (these are kde objects)
load("imap3d.diff_warm.RData")
im3d.diff.warm <- imap3d.diff
rm(imap3d.diff)
load("imap3d.diff_cold.RData")
im3d.diff.cold <- imap3d.diff
rm(imap3d.diff)

#Generate realization of difference intensity map:
source("C:/Users/galen/Documents/Research/oled-spatstat-analysis/rpoint3_extended.R")

f <- im3d.diff.warm
diff.pp3 <- rpoint3(n = length(f$w), f = f, im = TRUE, w = box3(c(-300,300),c(-300,300),c(100,700)))
plot3d.pp3(diff.pp3)

# Poisson hypothesis test: H0 = poisson process, H1 = some other process (see section 10.4.2 of SPP:MAR)
nx <- 12
ny <- 12
nz <- 12
X <- diff.pp3

qX <- quadratcount.pp3(X, nx = nx, ny = ny, nz = nz) #observed counts
mu <- mean(qX$count) # null hypothesis poission mean
tX <- table(factor(as.numeric(qX$count), levels = 0:max(qX$count))) # counts of number of quadrants where each value was observed
eX <- nx*ny*nz * dpois(0:max(qX$count), mu) # expected number of quadrants where each value is observed
eX
lower <- 5
upper <- 26

fX <- factor(as.numeric(qX$count), levels = 0:max(qX$count))
fX <- mergeLevels(fX, "lower.tail" = 0:lower) #combine lower and upper ends to avoid low probabiltiy problems
fX <- mergeLevels(fX, "upper.tail" = upper:max(qX$count)) 
tX <- table(fX) # observations for each count
tX

eX.mid <- dpois((lower+1):(upper-1), mu) #find expected probabilities for the upper and lower summed sections
eX.upper <- 1 - sum(dpois(0:(upper-1), mu))
eX.lower <- 1 - sum(eX.mid) - eX.upper
sum(eX.mid,eX.upper,eX.lower) == 1
eX <- nx*ny*nz*c(eX.lower, eX.mid, eX.upper) # expectations for each count
eX

X2 <- sum((tX - eX)^2/eX) #chi-square stat
pval <- pchisq(X2, df = length(tX)-1, lower.tail = FALSE)
pval # reject or confirm that null hypothesis

