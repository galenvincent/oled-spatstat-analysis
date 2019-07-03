library(rapt)
library(parallel)
library(data.table)

# Do RRL on RCP, BCC, and Random

nr <- 300
rmax <- 15

bcc.pp3 <- bcc.gen(216000, box3(c(0, 60), c(0, 60), c(0, 60)))
bcc.res <- pK3est(perc= 0.1, pattern = bcc.pp3, nEvals = 5, rmax = rmax, nrval = nr, anom = FALSE, sorted = FALSE)
save(bcc.res, file = "bcc.res.RData")

rand.pp3 <- rpoint3(216000, win = box3(c(0, 60), c(0, 60), c(0, 60)))
rand.res <- pK3est(0.1, rand.pp3, 5, rmax, nr, anom = FALSE)
save(rand.res, file = "rand.res.RData")

under <- read.rcp('~/Research/point_patterns/Final/FinalConfig1','~/Research/point_patterns/Final/system1',scaleUp = TRUE,newRadius = 0.5)
rcp.pp3 <- stitch.size(under, boxSize = c(60,60,60))
rcp.res <- pK3est(0.1, rcp.pp3, 50000, rmax, nr, anom = FALSE, sorted = FALSE)
save(rcp.res, file = "rcp.res.RData")




