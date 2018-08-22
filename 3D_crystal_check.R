# 3D crystal check
library(spatstat)

source("R/rapt-file.R")
source("R/rcp_functions.R")

nperfectnn <- function(rcp_1,rs,rb,ncutoff = 8){
  cutoffs <- c(2*rs+rs*.03,rb+rs+((rb+rs)*.03)/2,2*rb+rb*.03)

  nnd <- nndist(rcp_1,k=1:12)
  nnw <- nnwhich(rcp_1,k=1:12)
  mrks <- marks(rcp_1)

  cry <- vector('numeric',npoints(rcp_1))
  for(i in 1:npoints(rcp_1)){
    ptsum <- 0
    for(j in 1:12){
      if(mrks[i]==mrks[nnw[i,j]]){
        if(mrks[i]==2){
          if(nnd[i,j] < cutoffs[1]){
            ptsum <- ptsum + 1
          }
        }else{
          if(nnd[i,j] < cutoffs[3]){
            ptsum <- ptsum + 1
          }
        }
      }else{
        if(nnd[i,j] < cutoffs[2]){
          ptsum <- ptsum + 1
        }
      }
    }
    cry[i] <- ptsum
  }

  hist(cry,freq=F,breaks = 12,xlim=c(0,12),ylim=c(0,.25))

  perc <- sum(cry >= ncutoff)/length(cry)
  return(perc)
}
