library(spatstat)

source("R/rapt-file.R")
source("R/rcp_functions.R")
source("R/3D_crystal_check.R")

rcp_1_upload <- read.table("C:/Users/galen/Documents/Research/point_patterns/041718rcp_1_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
rcp_1 <- createSpat(rcp_1_upload[,c("x","y","z")])
marks(rcp_1) <- as.factor(rcp_1_upload$type)
r1s <- 0.0259528
r1b <- 0.0285481

rcp_2_upload <- read.table("C:/Users/galen/Documents/Research/point_patterns/041718rcp_2_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
rcp_2 <- createSpat(rcp_2_upload[,c("x","y","z")])
marks(rcp_2) <- as.factor(rcp_2_upload$type)
r2s <-0.0252329
r2b <-0.0302794

rcp_3_upload <- read.table("C:/Users/galen/Documents/Research/point_patterns/041718rcp_3_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
rcp_3 <- createSpat(rcp_3_upload[,c("x","y","z")])
marks(rcp_3) <- as.factor(rcp_3_upload$type)
r3s <- 0.0244723
r3b <- 0.031814

rcp_4_upload <- read.table("C:/Users/galen/Documents/Research/point_patterns/042818rcp_1_FinalConfig",sep = " ",col.names = c("x","y","z","type"))
rcp_4 <- createSpat(rcp_4_upload[,c("x","y","z")])
marks(rcp_4) <- as.factor(rcp_4_upload$type)
r4s <- 0.0217032
r4b <- 0.0303845

perc1 <- nperfectnn(rcp_1,r1s,r1b)
perc2 <- nperfectnn(rcp_2,r2s,r2b)
perc3 <- nperfectnn(rcp_3,r3s,r3b)
perc4 <- nperfectnn(rcp_4,r4s,r4b)
