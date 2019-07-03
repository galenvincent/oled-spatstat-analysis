# George APT analysis
library(rapt)

ato.1 <- readATO('Z:/LEAP\ Reconstructions/From\ George\ for\ Galen\ Paper/R45_03277-v02.ato')
rng.1 <- readRRNG('Z:/LEAP\ Reconstructions/From\ George\ for\ Galen\ Paper/R45_03277-v02.rrng')
#ms <- createSpec(ato.1, res = 0.05)
#ms.log <- transformIntensity(ms, method = 'log10')
#prettyPlot(ms.log, rng = rng.1, xlim = c(0,600))

rng.fin <- rng.1[1:33,]
ato.fin <- rngPOS(ato.1, rng.fin)



ato.2 <- readATO('Z:/LEAP\ Reconstructions/From\ George\ for\ Galen\ Paper/R45_03279-v02.ato')
