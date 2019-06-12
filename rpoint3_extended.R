#### rpoint3 ####
#' Extends \code{\link[spatstat]{rpoint}} to \code{\link[spatstat]{pp3}}.
rpoint3 <- function (n, f, im = FALSE ,fmax = 1,  win = box3(), ...,
                     giveup = 1000, verbose = FALSE, nsim = 1, drop = TRUE)
{
  if (missing(f) || (is.numeric(f) && length(f) == 1))
    return(runifpoint3(n, domain = win, nsim = nsim, drop = drop))
  if(im == TRUE){
    if(is.null(win)){
      wf <- box3(c(min(f$x$x), max(f$x$x)), c(min(f$x$y), max(f$x$y)), c(min(f$x$z), max(f$x$z)))
    }else{
      wf <- win
    }
    dx <- abs(f$eval.points[[1]][2]-f$eval.points[[1]][1])
    dy <- abs(f$eval.points[[2]][2]-f$eval.points[[2]][1])
    dz <- abs(f$eval.points[[3]][2]-f$eval.points[[3]][1])
    rxyz <- expand.grid(f$eval.points[[1]], f$eval.points[[2]], f$eval.points[[3]])
    xpix <- rxyz$Var1
    ypix <- rxyz$Var2
    zpix <- rxyz$Var3
    ppix <- vector(mode = "numeric", length = length(xpix))
    cnt <- 0
    for(i in 1:length(f$eval.points[[3]])){
      for(j in 1:length(f$eval.points[[2]])){
        for(k in 1:length(f$eval.points[[1]])){
          cnt <- cnt + 1
          ppix[cnt] <- f$estimate[k,j,i]
        }
      }
    }
    ppix[ppix < 0] <- 0
    
    result <- vector(mode = "list", length = nsim)
    for(isim in 1:nsim){
      id <- sample(length(xpix), n, replace = TRUE, prob = ppix)
      x <- xpix[id] + runif(n, min = -dx/2, max = dx/2)
      y <- ypix[id] + runif(n, min = -dy/2, max = dy/2)
      z <- zpix[id] + runif(n, min = -dz/2, max = dz/2)
      result[[isim]] <- pp3(x, y, z, window = wf)
    }
    result <- simulationresult(result, nsim, drop)
    return(result)
  }
  
  verifyclass(win, "box3")
  if (n == 0) {
    emp <- pp3(numeric(0), numeric(0), numeric(0), window = win)
    if (nsim == 1 && drop)
      return(emp)
    result <- rep(list(emp), nsim)
    names(result) <- paste("Simulation", 1:nsim)
    return(as.anylist(result))
  }
  result <- vector(mode = "list", length = nsim)
  for (isim in 1:nsim) {
    X <- pp3(numeric(0), numeric(0), numeric(0), window = win)
    pbar <- 1
    nremaining <- n
    totngen <- 0
    ntries <- 0
    repeat {
      ntries <- ntries + 1
      ngen <- nremaining/pbar + 10
      totngen <- totngen + ngen
      prop <- runifpoint3(ngen, domain = win)
      if (npoints(prop) > 0) {
        fvalues <- f(prop$data$x, prop$data$y, prop$data$z, ...)
        paccept <- fvalues/fmax
        u <- runif(npoints(prop))
        Y <- prop[u < paccept]
        if (npoints(Y) > 0) {
          X <- superimpose(X, Y, W = win)
          nX <- npoints(X)
          pbar <- nX/totngen
          nremaining <- n - nX
          if (nremaining <= 0) {
            if (verbose)
              splat("acceptance rate = ", round(100 * pbar, 2), "%")
            result[[isim]] <- if (nX == n)
              X
            else X[1:n]
            break
          }
        }
      }
      if (ntries > giveup)
        stop(paste("Gave up after", giveup * n, "trials with",
                   npoints(X), "points accepted"))
    }
  }
  if (nsim == 1 && drop)
    return(result[[1]])
  names(result) <- paste("Simulation", 1:nsim)
  return(as.anylist(result))
}
