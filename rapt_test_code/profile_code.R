# Script for code profiling and figuring out what to make faster.

# Source files from the correct folders. Note that they need to be sourced for profiling purposes.
source('C:/Users/galen/Documents/Research/rapt/R/cluster_functions.R')
source('C:/Users/galen/Documents/Research/rapt/R/rcp_functions.R')
source('C:/Users/galen/Documents/Research/rapt/R/envelope_functions.R')
source('C:/Users/galen/Documents/Research/rapt/R/rapt-file.R')
source('C:/Users/galen/Documents/Research/rapt/R/rapt-extend.R')
library(spatstat)
library(rgl)
library(profvis)
library(parallel)
#makecluster radius gaussian blur improvement - 11/26/18
profvis({
  under <- read.rcp('~/Research/point_patterns/Final/FinalConfig1','~/Research/point_patterns/Final/system1',scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp('~/Research/point_patterns/Final/FinalConfig2','~/Research/point_patterns/Final/system2',scaleUp = TRUE,newRadius = 0.5)
  #under.big <- stitch.size(under, boxSize = c(60,60,60))
  #over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  a <- makecluster(under, over, 0.5,0.5, type = "cr", speed = "superfast", cr = 3, rb = TRUE, rbp = 1)
  c <- makecluster(under, over, 0.5,0.5, type = "cr", speed = "superfast", cr = 3)
  
})


#makecluster function improvement 2 - 11/8/18
profvis({
  under <- read.rcp('~/Research/point_patterns/Final/FinalConfig1','~/Research/point_patterns/Final/system1',scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp('~/Research/point_patterns/Final/FinalConfig3','~/Research/point_patterns/Final/system3',scaleUp = TRUE,newRadius = 0.5)
  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))
  
  a <- makecluster(under.big, over.big, 0.5,0.5, type = "cr", speed = "superfast", cr = 3)
  b <- helpme(under.big, over.big, 0.5,0.5, type = "cr", speed = "superfast", cr = 3)
  c <- makecluster(under.big, over.big, 0.5,0.5, type = "cr", speed = "superfast", cr = 3, gb = TRUE)
  
})

# makecluster function improvement
profvis({
  
  # Upload RCP patterns
  r1 <- as.numeric(levels(read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/system_test_1')$V1)[2])
  r2 <- as.numeric(levels(read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/system_test_2')$V1)[2])
  rcp1 <- read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/rcp_test_1',sep=" ",col.names=c("x","y","z","type"))
  rcp1 <- scaleRCP(createSpat(rcp1[,c("x","y","z")]),newRadius = 0.5,oldRadius = r1)
  rcp2 <- read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/rcp_test_2',sep=" ",col.names=c("x","y","z","type"))
  rcp2 <- scaleRCP(createSpat(rcp2[,c("x","y","z")]),newRadius = 0.5,oldRadius = r2)
  
  
  a <- makecluster(rcp1,rcp2,0.5,0.5,type="ppc",ppc=10,pic=1,pcp=0.06,toPlot=FALSE,showOverPts=FALSE)
  b <- makecluster(rcp1,rcp2,0.5,0.5,type="cr",cr=2,fast=FALSE,pic=1,pcp=0.06,toPlot=FALSE,showOverPts=FALSE)
  c <- makecluster(rcp1,rcp2,0.5,0.5,type="dist",d=10,pic=1,pcp=0.06,toPlot=FALSE,showOverPts=FALSE)
  d <- makecluster(rcp1,rcp2,0.5,0.5,type="cr",cr=2,fast=TRUE,pic=1,pcp=0.06,toPlot=FALSE,showOverPts=FALSE)
  
})

rm(list = ls())
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/cluster_functions.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rcp_functions.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/envelope_functions.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rapt-file.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/R/rapt-extend.R')
source('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/cluster_test.R')
library(spatstat)
library(rgl)
library(profvis)
library(parallel)

# parallel spatial testing improvement
prf1 <- profvis({
  
  # Upload RCP patterns
  r1 <- as.numeric(levels(read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/system_test_1')$V1)[2])
  r2 <- as.numeric(levels(read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/system_test_2')$V1)[2])
  rcp1 <- read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/rcp_test_1',sep=" ",col.names=c("x","y","z","type"))
  rcp1 <- scaleRCP(createSpat(rcp1[,c("x","y","z")]),newRadius = 0.5,oldRadius = r1)
  rcp2 <- read.table('C:/Users/galen/OneDrive/Documents/Research/rapt/test files/rcp_test_2',sep=" ",col.names=c("x","y","z","type"))
  rcp2 <- scaleRCP(createSpat(rcp2[,c("x","y","z")]),newRadius = 0.5,oldRadius = r2)
  rcp3 <- stitch(rcp1)
  
  # run some tests
  aa <- panomK3est(0.06,rcp1,25,correction ="iso")
  a <- panomK3est(0.06,rcp1,50,correction ="iso")
  b <- panomK3est(0.06,rcp1,75,correction ="iso")
  #c <- panomK3est(0.06,rcp1,200,correction = "iso")
  d <- panomK3est(0.06,rcp3,25,correction = "iso")
  #e <- panomK3est(0.06,rcp3,100,correction = "iso")
  #c <- panomK3est(0.06,rcp1,50,correction="bord")
  #d <- panomK3est(0.06,rcp2,50,correction="iso",toSub = a[[2]])
  #e <- panomK3est(0.06,rcp2,50,correction="trans",toSub = b[[2]])
  #f <- panomK3est(0.06,rcp2,50,correction="bord",toSub = c[[2]])
  
  
  
})

htmlwidgets::saveWidget(prf1, "panomK3est_Profile.html")
# to read a profile
browseURL("panomK3est_Profile.html")
