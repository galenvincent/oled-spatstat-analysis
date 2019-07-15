# Quick crystal check - 2D

library(spatstat)

source("R/rapt-file.R")
source("R/rcp_functions.R")

rcp_1_upload <- read.table("C:/Users/galen/Desktop/042418_2D_1.1_FinalConfig",sep = " ",col.names = c("x","y","type"))
rcp_1 <- ppp(rcp_1_upload[,'x'],rcp_1_upload[,'y'], marks = rcp_1_upload[,'type'])
marks(rcp_1) <- as.factor(marks(rcp_1))

nn <- nndist(rcp_1,k=1:6)
cutoff <- 2*0.0101351

cry <- apply(nn,1,function(d){sum(d < cutoff + 0.001)})

plot(rcp_1,pch=c(1,19))
#points(rcp_1[cry==6],pch=19,col='red')
points(rcp_1[cry >= 5],pch=19,col='red')

sum(cry >= 5)/length(cry)
