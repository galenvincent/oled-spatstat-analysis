# Script to run in on the HPC for envelopes
library(spatstat)
library(parallel)

source("R/rapt-file.R")
source("R/rcp_functions.R")
source("R/galen_functions.R")

# Upload the RCP data set(un-stitched)
r1 <-0.0252429 #1.2 sigma
rcp_1_upload <- read.table("C:/Users/galen/Documents/Research/point_patterns/Final/FinalConfig1",sep = " ",col.names = c("x","y","z","type"))
rcp_1 <- scaleRCP(createSpat(rcp_1_upload[,c("x","y","z")]),newRadius = 0.5,oldRadius = r1)

#Stitch the set together
rcp <- stitch(rcp_1)

#run envelopes on it
t1 <- Sys.time()
a <- pK3est(0.06,rcp_1,50,nrval=20,anom=TRUE)
t2 <- Sys.time()

print(t2-t1)

#export the data from the envelope calculations
write.table(a[[1]], file = "C:/Users/galen/Desktop/data.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")

