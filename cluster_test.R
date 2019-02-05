# testing the makecluster function 
# 1/28/19

library(spatstat)
library(rapt)
library(data.table)

under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

p <- 100

more <- rep(0,p)
tot <- rep(0,p)
nc <- rep(0,p)
ncreal <- rep(0,p)
ncpred <- rep(0,p)

# more.o <- rep(0,p)
# over.rf <- rep(0,p)
# over.rf.o <- rep(0,p)

cnt <- 1

for(i in 1:p){
  under <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste('~/Research/point_patterns/Final/FinalConfig',toString(over.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(over.nums[i]),sep=""),scaleUp = TRUE,newRadius = 0.5)
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  a <- helpme(under.big, over.big, 0.5, 0.5, cr = 4, rb = TRUE, rbp = 1, rbmethod = 1, s = cnt)
  
  more[i] <- a[[3]]
  tot[i] <- npoints(a[[1]])
  nc[i] <- npoints(a[[2]])
  ncreal[i] <- sum(inside.boxx(a[[2]], w = domain(a[[2]])))
  ncpred[i] <- a[[6]]
  
  # more.o[i] <- a[[9]]
  # over.rf[i] <- a[[7]]
  # over.rf.o[i] <- a[[8]]
  
  print(toString(i))
  cnt <- cnt + 1
  
  rm(over, under, over.big, under.big, a)
}

