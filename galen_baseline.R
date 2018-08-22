# Import libraries and functions
library(spatstat)
library(parallel)

source("R/rapt-file.R")

# Read in data set and convert to pp3 (for testing purposes)
#data0 <- readATO("R45_00324-v45.ato")
#posData <- data0[,c("x","y","z")]
#data1 <- createSpat(posData)


#Function to grab a cubic subsection of the given data, given a number of points (plan on changing this to by dimensions)
subSquare <- function(n,ogPattern){
  # n = number of points you want in the box (integer)
  # og pattern = the original point pattern (pp3)

  density <- npoints(data1)/volume(as.box3(data1)); #density of OG point pattern (points/u^3)

  vol <- n/density; #required volume for new pp (u^3)
  s <- vol^(1/3) #side length of a cube

  frame <- coords(ogPattern); #extract x,y,z coordinates into a data.frame

  b <- box3(xrange=c(mean(frame$x)-s/2,mean(frame$x)+s/2),yrange=c(mean(frame$y)-s/2,mean(frame$y)+s/2),zrange=c(mean(frame$z)-s/2,mean(frame$z)+s/2));
  #create a box of dimensions that you need in order to select about the correct # of points

  newPatternFrame <- frame[inside.boxx(ogPattern,w=b),]; #select points out of original pattern and place into new data.frame
  newPattern <- createSpat(newPatternFrame);#convert new data.frame to pp3

  return(newPattern)
}

#test subset with about 1000 points
#a <- subSquare(1000,data1)

# Function to randomly relabel and select perc percent of the given point pattern
percentSelect <- function(perc,pattern){
  # perc = fraction of points you want to select (0-1)
  # pattern = the point pattern you want to draw from
  reLabel <- rlabel(pattern,labels = sample(c("A","B"),npoints(pattern),prob = c(perc,1-perc),replace = TRUE))
  inds <- which(marks(reLabel)=="A")
  newPattern <- reLabel[inds]
  return(newPattern)
}

# Function to do K3est test with random relabelings of a certain percent, a certain number of times
# returns a matrix with test results in each column, and the r values in the first column
rrK3est <- function(perc, pattern, nEvals){
  # perc = percent of original pattern to sample
  # pattern = the original pattern
  # nEvals = number of times to run K3est

  # run the test once to fill the r values
  toTest <- percentSelect(perc,pattern)
  tst <- K3est(toTest)
  tests <- tst$r
  cbind(tests,tst$iso)

  # run the tests and pull the "iso" edge condition results into the array
  for(i in 1:nEvals-1){
    toTest <- percentSelect(perc,pattern)
    tst <- K3est(toTest)
    tests <- cbind(tests,tst$iso)
  }
  return(tests)
}

# Function to plot envelope results from the rrK3est function above, based on percentiles
envPlot <- function(tests,percentiles=c(.99,.9,.5)){
  # tests = array of values returned from the rrK3est function above
  # percentiles = vector including the different percentiles you would like to see on the plot (0-1)
    # do these in descending order please

  color <- c("darksalmon","darkorchid","darkseagreen1","chocolate1","firebrick1")

  # break up data into r values and test results
  rvals <- tests[,1]
  tvals <- tests[,2:ncol(tests)]

  nTests <- ncol(tvals) # number of tests done
  prange <- percentiles*nTests # get the range of indeces for which each percentile spans

  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  percentileIndicesBig <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  percentileIndicesSmall <- round(nTests/2)-floor(prange/2) # do the same for the low end

  toPlotBigs <- NULL
  toPlotSmalls <- NULL
  # grab out the columns from the sorted test results that we will plot
  for(i in 1:length(percentiles)){
    toPlotBigs <- cbind(toPlotBigs,sortedtVals[,percentileIndicesBig[i]])
    toPlotSmalls <- cbind(toPlotSmalls,sortedtVals[,percentileIndicesSmall[i]])
  }

  # plot the envelopes from the percentile data
  plot(rvals,tvals[,1],type="n",main="Envelopes for K3est Function",xlab="r")
  for(i in 1:length(percentiles)){
    polygon(c(rvals,rev(rvals)),c(toPlotBigs[,i],rev(toPlotSmalls[,i])),col=color[i])
  }
}

# Lets get into some parallelization, shal we?
# Parellelized version of rrK3est, written above. Perfeorms the K3est function caculations in parallel
parrrK3est <- function(perc, pattern, nEvals){
  # perc = percent of original pattern to sample
  # pattern = the original pattern
  # nEvals = number of times to run K3est

  #find cores and initialize the cluster
  cores2use <- detectCores()-1
  cl <- makeCluster(cores2use)

  #initiate list to hold the different data matrices
  toTest <- vector("list",nEvals)

  # populate the list with randomly generated results
  for(i in 1:nEvals){
    toTest[[i]] <- percentSelect(perc,pattern)
  }

  # run the test once to fill the r values
  tst <- K3est(percentSelect(perc,pattern))
  tests <- tst$r
  cbind(tests,tst$iso)

  # apply K3est function to each of the pp3 patterns in parallel
  result <- parLapply(cl,toTest,K3est)

  # convert the results into the matrix tests
  for(i in 1:length(result)){
    tests <- cbind(tests,result[[i]]$iso)
  }

  # stop the cluster and revert computer to normal
  stopCluster(cl)

  return(tests)
}
