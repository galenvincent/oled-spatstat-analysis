#### localk3metrics ####
#' Extract variance of R_max over local K functions from a point pattern.
#'
#' This function takes the output of \code{\link{localK3est}}, calculates R_max
#' for the individual local K measures around each point, and returns the
#' variance of these measures.
#'
#' @param kl The output object from a \code{\link{localK3est}} function call.
#' @param start The x value index to start searching for mertics at (to avoid
#'   strange behavior near r = 0).
#' @param nsamp How many samples of R_max to take from the sample. If
#'   \code{NULL}, defaults to the number of points in the pattern (this can be
#'   very computationally expensive for large patterns).
#' @return A signle numerical value which is the variance of the R_max metric
#'   over \code{nsamp} individual points in the pattern.
#' @seealso \code{\link{localK3est}}, \code{\link[spatstat]{K3est}},
#'   \code{\link{kseries2}}

localk3metrics <- function(kl, start, nsamp = NULL){
  span <- NA
  i <- 1
  while(is.na(span)){
    r.n <- kl$r[start:nrow(kl)]
    t.n <- kl[start:nrow(kl),i+3]

    peak.info <- argmax(r.n, t.n, w = 3, span = 0.1)
    span <- (peak.info$x[1]/7)*(0.3)

    i <- i + 1
  }

  if(is.null(nsamp)){
    nsamp <- np
  }

  pks2 <- rep(0,nsamp)
  pts <- sample(1:np, nsamp, replace = FALSE)
  for(i in 1:nsamp){
    t.n <- kl[15:nrow(kl),pts[i]]
    a <- argmax(r.n, t.n, w = 3, span = span)
    pks2[i] <- a$x[1]
  }

  b <- var(pks2, na.rm = TRUE)
  return(b)
}

#### kseries ####
#' Simulate a clustered point pattern and extract the K-function metrics from
#' it.
#'
#' This function takes an RCP point pattern and a single set of parameters to
#' create a cluster and extracts the K function metrics.
#'
#' @param j The index of the RCP pattern to run.
#' @param params Vector of length 3 or 4 containing either \code{c(r, den, gbp)}
#'   or \code{c(r, den, rb, gbp)}.
#' @param maxr Maximum r to calculate K-function to.
#' @param nr Number of r values to evaluate K-function at.
#' @param toSub A vector of values to substitute for the \code{\link{anomK3est}}
#'   function.
#' @param pcp See \code{pcp} argument for \code{\link{makecluster}} function.
#' @param rcp_path String holding the file path to the directory holding the RCP
#'   'FinalConfig' and 'system' files.
#' @param s Random seed.
#'
#' @return Vector containing 5 metrics.

kseries <- function(j, params, maxr, nr, toSub, pcp = 0.06,
                    rcp_path = '~/Research/point_patterns/Final',
                    s = NULL){
  #upload
  under <- read.rcp(paste(rcp_path, '/FinalConfig', toString(j), sep=''),
                    paste(rcp_path, '/system', toString(j), sep=''),
                    scaleUp = TRUE, newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, '/FinalConfig', toString(j), sep=''),
                   paste(rcp_path, '/system', toString(j), sep=''),
                   scaleUp = TRUE, newRadius = 0.5)

  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))

  if(!is.null(s)){
    set.seed(s)
  }else{
    s <- 100
    set.seed(s)
  }

  if(length(params) == 4){
    cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",
                           pcp = pcp,
                           cr=params[1],
                           den = params[2],
                           rb = TRUE, rbp = params[1]*params[3],
                           gb = TRUE, gbp = c(0, params[4]),
                           s = s)
  }else{
    cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",
                           pcp = pcp,
                           cr=params[1],
                           den = params[2],
                           gb = TRUE, gbp = c(0, params[3]),
                           s = s)
  }
  if(is.numeric(cluster)){
    out <- c(NA, NA, NA, NA, NA)
    return(out)
  }

  result <- anomK3est(cluster[[1]],toSub,maxr,nr)
  rvals <- result$r
  tvals <- result$trans

  # get out that peak info son
  rvals.new <- rvals[15:length(rvals)]
  tvals.new <- tvals[15:length(rvals)]

  #get those metrics out
  metrics <- k3metrics(rvals.new, tvals.new, FALSE)

  out <- c(metrics[[1]], metrics[[2]], metrics[[3]], metrics[[4]], metrics[[5]])

  rm(cluster, result, rvals, tvals, rvals.new, tvals.new, over, under, over.big, under.big)
  gc()


  return(out)
}

#### kseries2 ####
#' Helper function to allow for parallelization of metric extraction from
#' cluster simulations
#'
#' This function takes a single RCP pattern and creates and tests the K function
#' metrics on different cluster properties passed into the function. Used to
#' parallelize large sets of cluster property runs.
#'
#' @param j The index of the RCP pattern to run.
#' @param p The total number of RCP patterns available.
#' @param tot A list of cluster properties to test. This should be a list of
#'   vectors containing \code{c(r, den, rb, gbp)} OR \code{c(r, den, gbp)}.
#' @param maxr Maximum r to calculate K function to.
#' @param nr Number of r values to evaluate k function at.
#' @param toSub A vector of values to substitute for the \code{\link{anomK3est}}
#'   function.
#' @param pcp See \code{pcp} argument for \code{\link{makecluster}} function.
#' @param rcp_path String holding the file path to the directory holding the RCP
#'   'FinalConfig' and 'system' files.
#' @param verbose \code{TRUE} or \code{FALSE}. Whether to output update files to
#'   a junk folder.
#' @param junk_path String holding the file path to the directory where the
#'   update files should go if \code{verbose = TRUE}.
#' @param s Random seed.
#'
#' @return Matrix containing 5 metrics for each parameter combination given in
#'   the \code{tot} list.

kseries2 <- function(j, p ,tot, maxr, nr, toSub, pcp = 0.06,
                     rcp_path = '~/Research/point_patterns/Final',
                     verbose = FALSE,
                     junk_path = '~/Research/junk/',
                     s = NULL){
  #t1 <- Sys.time()
  under.nums <- seq(2,(p+1),1)
  under.nums[length(under.nums)] <- 1
  over.nums <- seq(1,p,1)

  #upload
  under <- read.rcp(paste(rcp_path, '/FinalConfig', toString(under.nums[j]), sep=''),
                    paste(rcp_path, '/system', toString(under.nums[j]), sep=''),
                    scaleUp = TRUE,newRadius = 0.5)
  over <- read.rcp(paste(rcp_path, '/FinalConfig', toString(over.nums[j]), sep=''),
                   paste(rcp_path, '/system', toString(over.nums[j]), sep=''),
                   scaleUp = TRUE,newRadius = 0.5)

  under.big <- stitch.size(under, boxSize = c(60,60,60))
  over.big <- stitch.size(over, boxSize = c(60,60,60))

  if(!is.null(s)){
    set.seed(s)
  }

  cnt <- j*length(tot)*round(runif(length(tot), 1, 10000))

  outtemp <- matrix(NA, nrow = length(tot), ncol = 5)

  for(i in 1:length(tot)){
    #print(i)
    if(length(tot[[1]]) == 4){
      cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",
                             pcp = pcp,
                             cr=tot[[i]][1],
                             den = tot[[i]][2],
                             rb = TRUE, rbp = tot[[i]][1]*tot[[i]][3],
                             gb = TRUE, gbp = c(0, tot[[i]][4]),
                             s = cnt[i])
    }else{
      cluster <- makecluster(under.big,over.big,0.5,0.5,type="cr",speed="superfast",
                             pcp = pcp,
                             cr=tot[[i]][1],
                             den = tot[[i]][2],
                             gb = TRUE, gbp = c(0, tot[[i]][3]),
                             s = cnt[i])
    }
    if(is.numeric(cluster)){
      outtemp[i,] <- c(NA, NA, NA, NA, NA)
      a <- as.data.frame(1)
      if(verbose == TRUE){
        if(i > 1){
          file.remove(paste(junk_path, '/', toString(j), '_', toString(i-1), '.csv', sep = ''))
        }
        write.csv(a, file = paste(junk_path, '/', toString(j), '_', toString(i), '.csv', sep = ''))
      }
      next
    }
    result <- anomK3est(cluster[[1]],toSub,maxr,nr)
    rvals <- result$r
    tvals <- result$trans

    # get out that peak info son
    rvals.new <- rvals[15:length(rvals)]
    tvals.new <- tvals[15:length(rvals)]

    #get those metrics out
    metrics <- k3metrics(rvals.new, tvals.new, FALSE)

    outtemp[i,] <- c(metrics[[1]], metrics[[2]], metrics[[3]], metrics[[4]], metrics[[5]])

    rm(cluster, result, rvals, tvals, rvals.new, tvals.new)
    gc()

    cnt <- cnt + 1

    a <- as.data.frame(1)

    if(verbose == TRUE){
      if(i > 1){
        file.remove(paste(junk_path, '/', toString(j), '_', toString(i-1), '.csv', sep = ''))
      }
      write.csv(a, file = paste(junk_path, '/', toString(j), '_', toString(i), '.csv', sep = ''))
    }
  }

  rm(over, under, over.big, under.big)
  gc()

  return(outtemp)
}

