library(data.table)
library(rapt)

morph <- fread('C:/Users/galen/Desktop/morphology_0.txt', skip = 12 ,col.names = c('x', 'y', 'z', 'mark'))
morph.box <- box3(c(0, 60), c(0, 60), c(0, 60))

morph.pp3 <- pp3(morph$x, morph$y, morph$z, morph.box, marks = morph$mark)
morph.pp3.1 <- morph.pp3[marks(morph.pp3)==1]
morph.pp3.2 <- morph.pp3[marks(morph.pp3)==2]

plot3d.pp3(morph.pp3.1, col = 'blue')
plot3d.pp3(morph.pp3.2, col = 'red', add = TRUE)

#re-scale morphology
scale.factors <- c(1/60, 1/60, 1/60)
morph.pp3.scaled <- pp3(morph$x*scale.factors[1], 
                        morph$y*scale.factors[2], 
                        morph$z*scale.factors[3], 
                        box3(c(0,1), c(0,1), c(0,1)), marks = morph$mark)

# re-assign to poisson points
lambda <- 50000
pois <- rpoispp3(lambda, box3(c(0,1), c(0,1), c(0,1)))

nnw <- nncross(pois, morph.pp3.scaled, what = 'which')
pois.marks <- marks(morph.pp3.scaled[nnw])

marks(pois) <- pois.marks


plot3d.pp3(pois[marks(pois) == 1], col = 'blue')
plot3d.pp3(pois[marks(pois) == 2], col = 'red', add = TRUE)

open3d()
plot3d.pp3(morph.pp3.scaled[marks(morph.pp3.scaled) == 1], col = 'blue')
plot3d.pp3(morph.pp3.scaled[marks(morph.pp3.scaled) == 2], col = 'red', add = TRUE)

# check out the cross sections
# lnt <- 0.1
# scng <- 0.01
# slab1 <- seq(0, 1 - lnt, by = scng)
# slab2 <- seq(lnt, 1, by = scng)
# coo <- coords(pois)
# for(i in 1:length(slab1)){
#   coo.cut <- coo[coo$z > slab1[i] & coo$z < slab2[i],]
#   marks.cut <- pois.marks[coo$z > slab1[i] & coo$z < slab2[i]]
#   
#   plot(coo.cut[marks.cut == 1,1], coo.cut[marks.cut == 1,2], col = 'red', pch = 16, cex = 0.75)
#   points(coo.cut[marks.cut == 2,1], coo.cut[marks.cut == 2,2], col = 'blue', pch = 16, cex = 0.75)
# }




