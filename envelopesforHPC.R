# Script to run in on the HPC for envelopes
library(spatstat)
library(parallel)

source("/u/cb/dm/galenvincent/bins/rtest/rapt-file.R")
source("/u/cb/dm/galenvincent/bins/rtest/rcp_functions.R")
source("/u/cb/dm/galenvincent/bins/rtest/galen_functions.R")

# Upload the RCP data set(un-stitched)
r1 <-0.0252429 #1.2 sigma
rcp_1_upload <- read.table("/u/cb/dm/galenvincent/bins/rtest/FinalConfig1",sep = " ",col.names = c("x","y","z","type"))
rcp_1 <- scaleRCP(createSpat(rcp_1_upload[,c("x","y","z")]),newRadius = 0.5,oldRadius = r1)

#Stitch the set together
rcp <- stitch(rcp_1)

#run envelopes on it
a <- panomK3est(0.06,rcp,50,nrval=20)

#export the data from the envelope calculations
write.table(a[[1]], file = "/u/cb/dm/galenvincent/bins/rtest/data.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")

