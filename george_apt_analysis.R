# George APT analysis
library(rapt)
library(ks)
library(parallel)

file <- 'R45_03277-v02'
file <- 'R45_03279-v02'

ato.1 <- readATO(paste('X:/LEAP\ Reconstructions/From\ George\ for\ Galen\ Paper/', file,'.ato', sep = ''))
rng.1 <- readRRNG(paste('X:/LEAP\ Reconstructions/From\ George\ for\ Galen\ Paper/', file, '.rrng', sep = ''))
#ms <- createSpec(ato.1, res = 0.05)
#ms.log <- transformIntensity(ms, method = 'log10')
#prettyPlot(ms.log, rng = rng.1, xlim = c(0,600))

if(file == 'R45_03277-v02'){
  rng.sup <- rng.1[1:33,]
  rng.as <- rng.1[21:28,]
  xlit <- c(-250,250)
  ylit <- xlit
  zlit <- c(4000, 5500)
}else if(file == 'R45_03279-v02'){
  rng.sup <- rng.1[1:35,]
  rng.as <- rng.1[21:28,]
  xlit <- c(-250,250)
  ylit <- xlit
  zlit <- c(6000, 7500)
}

ato.sup <- rngPOS(ato.1, rng.sup)
rm(ato.1)
gc()


ato.as <- rngPOS(ato.sup, rng.as)
as.pp3 <- createSpat(ato.as, win = box3(xlit, ylit, zlit)) 
marks(as.pp3) <-ato.as$mark
as.pp3 <- as.pp3[inside.boxx(as.pp3, w = domain(as.pp3))]

sup.pp3 <- createSpat(ato.sup, win = box3(xlit, ylit, zlit))
marks(sup.pp3) <- ato.sup$mark
sup.pp3 <- sup.pp3[inside.boxx(sup.pp3, w = domain(sup.pp3))]

rm(ato.sup, ato.as, rng.1, rng.sup, rng.as, file)
gc()

ratio <- npoints(as.pp3)/npoints(sup.pp3)

plot3d.pp3(as.pp3, aspect = FALSE, col = "red")
plot3d.pp3(percentSelect(0.005, sup.pp3), aspect = FALSE, col = 'black')

as.map <- kde(coords(as.pp3), gridsize = 75)
plot(as.map, cont = c(5, 10, 50))

a <- pK3est(ratio, sup.pp3, 100, anom = TRUE, sorted = FALSE, rmax = 250)
b <- anomK3est(as.pp3, a[[2]], a[[3]], a[[4]])

envPlot(a[[1]], ylim = c(-5000,5000))
lines(b$r, b$trans, lwd = 2)
b <- K3est(as.pp3)
