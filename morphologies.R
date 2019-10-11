# Simulating different morphologies

library(rapt)
library(data.table)
library(parallel)

#### planes ####
#' Simulate lamellar morphology for mixed guest-host system.
#'
#' @param lambda Intensity for the background Poisson process.
#' @param frac Fraction of points to select as guest points. Between zero and
#'   one.
#' @param plane.norm Vector you wish to define normal of the planes with.
#'   Doesn't need to be normalized. Positive values only
#' @param plane.den Linear plane density.
#' @param point.den Point density within the plane volumes.
#' @param toplot \code{TRUE} or \code{FALSE}. Plot results or not.
#' @param win The simulation window.
#'
#' @return A list of: [[1]] a \code{\link[spatstat]{pp3}} object contining the
#'   guest points. [[2]] A \code{pp3} object containing the entire background
#'   underlying point pattern.

morph_lamellar <- function(lambda, 
                          frac,
                          plane.norm = c(1, 0, 0), 
                          plane.den, 
                          point.den = 1,
                          toplot = FALSE,
                          win = box3(c(0,1), c(0,1), c(0,1))){
  
  bgnd <- rpoispp3(lambda = lambda, domain = win)
  coo <- coords(bgnd)
  p.norm <- plane.norm
  p.den <- plane.den
  
  #plane density (planes per unit length)
  p.spacing <- 1/p.den
  
  #plane orientation (vector holding normal direction for planes, doesn't have to be normalized)
  p.direc <- p.norm/sqrt(sum(p.norm^2)) * p.spacing
  
  #points on planes:
  p.points <- list()
  p.points[[1]] <- c(-1, -1, -1)
  nxt <- p.points[[1]] + p.direc
  i <- 1
  while(nxt[1] <= win$xrange[2]*1.5 &
        nxt[2] <= win$yrange[2]*1.5 &
        nxt[3] <= win$zrange[2]*1.5){
    p.points[[i+1]] <- nxt
    nxt <- p.points[[i+1]] + p.direc
    i <- i + 1
  }
  
  offset <- runif(1, 0, 1)
  p.points <- lapply(p.points, function(x){x + p.norm/sqrt(sum(p.norm^2)) * offset
    print(x + p.norm/sqrt(sum(p.norm^2)) * offset)})
  
  #select points within x distance of planes to be type-A
  distmat <- matrix(FALSE, npoints(bgnd), length(p.points))
  p.norm.normed <- p.norm/sqrt(sum(p.norm^2))
  
  for(i in 1:length(p.points)){
    vecs <- t(apply(coo, 1, function(x){x - p.points[[i]]}))
    distmat[,i] <- abs(vecs %*% p.norm.normed)
  }
  
  nkeep <- round(frac*npoints(bgnd))
  mindistlist <- apply(distmat, 1, min)
  gbwhich <- order(mindistlist)[1:nkeep]
  
  if(point.den != 1){
    nin <- round(point.den*length(gbwhich))
    nout <- length(gbwhich) - nin
    points.final <- sample(gbwhich, nin)
    
    nonplane.ind <- (1:npoints(bgnd))[-gbwhich]
    noise <- sample(nonplane.ind, nout)
    total.inds <- c(points.final, noise)
  }else{
    total.inds <- gbwhich
  }
  
  p.selected <- bgnd[total.inds]
  
  if(toplot == TRUE){
    plot3d.pp3(p.selected, col = 'red', xlim = bgnd$domain$xrange, ylim = bgnd$domain$yrange, zlim = bgnd$domain$zrange)
  }
  
  return(p.selected)
}

#### rods ####
#' Simulate morphology of guest rods in a two-phase system
#'
#' @param lambda Intensity for the background Poisson process.
#' @param frac Fraction of points to select as guest points. Between zero and
#'   one.
#' @param rod.norm Either 'x', 'y', or 'z'; defines which axis the rods should
#'   be alligned with.
#' @param rod.den Rod density; 2D vector containing a linear density for both
#'   directions.
#' @param rod.spacing One of: "grid", "rcp", "hexagonal", or "random". The
#'   method used to space the rods.
#' @param rcp.path If \code{rod.spacing} is set to 'rcp', then set the file path
#'   to the folder contining the rcp FinalConfig and system files
#' @param point.den Point density within the volume of the rods.
#' @param toplot \code{TRUE} or \code{FALSE}; Plot results or not.
#' @param win The simulation window.
#'
#' @return A list of: [[1]] a \code{\link[spatstat]{pp3}} object contining the
#'   guest points. [[2]] A \code{pp3} object containing the entire background
#'   underlying point pattern.

morph_rods <- function(lambda, 
                      frac,
                      rod.norm = 'z', 
                      rod.den,
                      rod.spacing = 'grid',
                      rcp.path = 'C:/Users/galen/Documents/Research/point_patterns/2D/Final',
                      point.den = 1,
                      toplot = FALSE,
                      win = box3(c(0,1), c(0,1), c(0,1))){
  # generate points
  bgnd <- rpoispp3(lambda = lambda, domain = win)
  coo <- coords(bgnd)
  
  # grid of rods using the appropriate method and axis
  if(rod.norm == 'z') {
    comp.points <- data.frame(coo$x, coo$y)
    xdist <- diff(win$xrange)
    ydist <- diff(win$yrange)
  } 
  else if (rod.norm == 'y') {
    comp.points <- data.frame(coo$x, coo$z)
    xdist <- diff(win$xrange)
    ydist <- diff(win$zrange)
  } 
  else if (rod.norm == 'x') {
    comp.points <- data.frame(coo$y, coo$z)
    xdist <- diff(win$yrange)
    ydist <- diff(win$zrange)
  } 
  else {
    print('Please input one of \'x\', \'y\', or \'z\' for rod.norm')
  }

  # create rod spacing
  buffer <- max(comp.points)*0.1
  
  if(rod.spacing == 'grid'){
    rod.x <- seq(min(comp.points[,1]) - 1 - buffer, max(comp.points[,1]) + 1 + buffer, by = 1/rod.den[1])
    rod.y <- seq(min(comp.points[,2]) - 1 - buffer, max(comp.points[,2]) + 1 + buffer, by = 1/rod.den[2])
    offset <- c(runif(1,-1, 1), runif(1,-1,1))
    rod.x <- rod.x + offset[1]
    rod.y <- rod.y + offset[2]
    rod.x <- rod.x[rod.x > (min(comp.points[,1]) - buffer) & rod.x < (max(comp.points[,1]) + buffer)]
    rod.y <- rod.y[rod.y > (min(comp.points[,2]) - buffer) & rod.y < (max(comp.points[,2]) + buffer)]
    rod.xy <- expand.grid(rod.x, rod.y)
  } 
  else if (rod.spacing == 'rcp') {
    rcp.upload <- fread(paste(rcp.path, '/FinalConfig', sep = ''))
    rcp <- ppp(rcp.upload$V1, rcp.upload$V2)
    
    rcp.den <- rcp$n/area(rcp$window)
    scaling.factor <- sqrt(rcp.den/(rod.den[1]*rod.den[2]))
    
    rcp.upload.scaled <- coords(rcp)*scaling.factor
    rcp.scaled <- ppp(rcp.upload.scaled$x, rcp.upload.scaled$y, window = owin(c(0,max(rcp.upload.scaled$x)), c(0,max(rcp.upload.scaled$y))))
    
    xmu <- mean(rcp.scaled$window$xrange)
    ymu <- mean(rcp.scaled$window$xrange)
    
    win.select <- owin(c(xmu - xdist/2 - buffer, xmu + xdist/2 + buffer), 
                       c(ymu - ydist/2 - buffer, ymu + ydist/2 + buffer))
    rod.xy <- coords(rcp.scaled[inside.owin(rcp.scaled$x, rcp.scaled$y, win.select)])
    rod.xy$x <- rod.xy$x - xmu + xdist/2 
    rod.xy$y <- rod.xy$y - ymu + ydist/2 
  } 
  else if (rod.spacing == 'hexagonal') {
    hex <- hexgrid(owin(c(0,1), c(0,1)), 0.01)
    
    hex.den <- hex$n/area(hex$window)
    scaling.factor <- sqrt(hex.den/(rod.den[1]*rod.den[2]))
    
    hex.scaled.coo <- coords(hex)*scaling.factor
    offset <- c(runif(1, 0, 1), runif(1, 0, 1))
    
    hex.scaled <- ppp(hex.scaled.coo$x + offset[1], hex.scaled.coo$y + offset[2], window = owin(c(0,max(hex.scaled.coo$x + offset[1])), c(0,max(hex.scaled.coo$y + offset[2]))))
    
    xmu <- mean(hex.scaled$window$xrange)
    ymu <- mean(hex.scaled$window$xrange)
    
    win.select <- owin(c(xmu - xdist/2 - buffer, xmu + xdist/2 + buffer), 
                       c(ymu - ydist/2 - buffer, ymu + ydist/2 + buffer))
    rod.xy <- coords(hex.scaled[inside.owin(hex.scaled$x, hex.scaled$y, win.select)])
    rod.xy$x <- rod.xy$x - xmu + xdist/2 
    rod.xy$y <- rod.xy$y - ymu + ydist/2 
  } 
  else if (rod.spacing == 'random') {
    xrng <- win$xrange*10
    yrng <- win$yrange*10
    mat <- rMaternII((rod.den[1]*rod.den[2])*1.25, buffer*2, owin(xrng, yrng))
    
    xmu <- mean(mat$window$xrange)
    ymu <- mean(mat$window$xrange)
    
    win.select <- owin(c(xmu - xdist/2 - buffer, xmu + xdist/2 + buffer), 
                       c(ymu - ydist/2 - buffer, ymu + ydist/2 + buffer))
    rod.xy <- coords(mat[inside.owin(mat$x, mat$y, win.select)])
    rod.xy$x <- rod.xy$x - xmu + xdist/2 
    rod.xy$y <- rod.xy$y - ymu + ydist/2 
  }

  distmat <- matrix(FALSE, npoints(bgnd), nrow(rod.xy))
  for(i in 1:nrow(rod.xy)){
    distmat[,i] <- sqrt((comp.points[,1] - rod.xy[i,1])^2 + (comp.points[,2] - rod.xy[i,2])^2)
  }
  
  nkeep <- round(frac*npoints(bgnd))
  mindistlist <- apply(distmat, 1, min)
  gbwhich <- order(mindistlist)[1:nkeep]
  
  if(point.den != 1){
    nin <- round(point.den*length(gbwhich))
    nout <- length(gbwhich) - nin
    points.final <- sample(gbwhich, nin)
    
    nonplane.ind <- (1:npoints(bgnd))[-gbwhich]
    noise <- sample(nonplane.ind, nout)
    total.inds <- c(points.final, noise)
  }else{
    total.inds <- gbwhich
  }
  
  p.selected <- bgnd[total.inds]
  
  if(toplot == TRUE){
    plot3d.pp3(p.selected, col = 'red')
  }
  
  return(p.selected)
}
 
 
#### Gyroid ####
#' Simulate morphology of guest molecules on a grain boundary
#'
#' @param lambda Intensity for the background Poisson process.
#' @param frac Fraction of points to select as guest points. Between zero and
#'   one.
#' @param gyroid.scale Scale factor for the gyroid pattern. Larger = tighter
#'   spacing of the gyroid.
#' @param point.den Point density within the volume of the guest volume.
#' @param toplot \code{TRUE} or \code{FALSE}; Plot results or not.
#' @param win The simulation window.
#'
#' @return A list of: [[1]] a \code{\link[spatstat]{pp3}} object contining the
#'   guest points. [[2]] A \code{pp3} object containing the entire background
#'   underlying point pattern.

morph_gyroid <- function(lambda, 
                        frac,
                        gyroid.scale,
                        point.den = 1,
                        toplot = FALSE,
                        win = box3(c(0,1), c(0,1), c(0,1))){
  
  shift <- c(runif(1, 0, 1000), runif(1, 0, 1000), runif(1, 0, 1000))
  
  gyr <- function(x, y, z, sc){
    return(sin(8*sc*x + shift[1])*cos(8*sc*y + shift[2]) + 
             sin(8*sc*y + shift[2])*cos(8*sc*z + shift[3]) + 
             sin(8*sc*z + shift[3])*cos(8*sc*x + shift[1]))
  }
  
  bgnd <- rpoispp3(lambda = lambda, domain = win)
  coo <- coords(bgnd)
  
  coo.gyr <- gyr(coo$x, coo$y, coo$z, gyroid.scale)
  nkeep <- round(frac*npoints(bgnd))
  
  gbwhich <- tail(order(coo.gyr), n = nkeep)
  
  if(point.den != 1){
    nin <- round(point.den*length(gbwhich))
    nout <- length(gbwhich) - nin
    points.final <- sample(gbwhich, nin)
    
    nonplane.ind <- (1:npoints(bgnd))[-gbwhich]
    noise <- sample(nonplane.ind, nout)
    total.inds <- c(points.final, noise)
  }else{
    total.inds <- gbwhich
  }
  
  p.selected <- bgnd[total.inds]
  
  if(toplot == TRUE){
    plot3d.pp3(p.selected, col = 'red')
  }
  
  return(p.selected)
}

#### Grain Boundary ####
#' Simulate morphology of guest molecules on a grain boundary
#'
#' @param lambda Intensity for the background Poisson process.
#' @param frac Fraction of points to select as guest points. Between zero and
#'   one.
#' @param point.den Point density within the volume of the guest volume.
#' @param rcp.rad Radius to scale the RCP pattern to.
#' @param rcp.path File path to folder containing the 3D RCP patterns which
#'   contain both the 'FinalConfigx' and 'systemx' files.
#' @param rcp.number Which RCP file to use.
#' @param toplot \code{TRUE} or \code{FALSE}; Plot results or not.
#' @param win The simulation window.
#'
#' @return A list of: [[1]] a \code{\link[spatstat]{pp3}} object contining the
#'   guest points. [[2]] A \code{pp3} object containing the entire background
#'   underlying point pattern.

morph_gb <- function(lambda, 
                    frac,
                    point.den = 1,
                    rcp.rad,
                    rcp.path = 'C:/Users/galen/Documents/Research/point_patterns/Final',
                    rcp.number = 1,
                    toplot = FALSE,
                    win = box3(c(0,1), c(0,1), c(0,1))){
  
  if(rcp.number == 'rand'){
    rcp.number <- sample(1, 1:523)
  }
  bgnd <- rpoispp3(lambda, domain = win)
  rcp <- read.rcp(paste(rcp.path, '/FinalConfig', toString(rcp.number), sep = ''),
                  paste(rcp.path, '/system', toString(rcp.number), sep = ''),
                  scaleUp = TRUE, newRadius = rcp.rad)
  xmu <- mean(rcp$domain$xrange)
  ymu <- mean(rcp$domain$yrange)
  zmu <- mean(rcp$domain$zrange)
  
  xdist <- diff(win$xrange)
  ydist <- diff(win$yrange)
  zdist <- diff(win$zrange)
  browser()
  offset <- c(runif(1, -1, 1), runif(1, -1, 1), runif(1, -1, 1))
  
  win.select <- box3(c(xmu - xdist/2 - rcp.rad - 1 + offset[1], xmu + xdist/2 + rcp.rad + 1 + offset[1]), 
                     c(ymu - ydist/2 - rcp.rad - 1 + offset[2], ymu + ydist/2 + rcp.rad + 1 + offset[2]),
                     c(zmu - zdist/2 - rcp.rad - 1 + offset[3], zmu + zdist/2 + rcp.rad + 1 + offset[3]))
  rcp.xyz <- coords(rcp[inside.boxx(rcp, w = win.select)])
  rcp.xyz$x <- rcp.xyz$x - xmu + xdist/2 - offset[1]
  rcp.xyz$y <- rcp.xyz$y - ymu + ydist/2 - offset[2]
  rcp.xyz$z <- rcp.xyz$z - zmu + zdist/2 - offset[3]
  rcp.cut <- pp3(rcp.xyz$x, rcp.xyz$y, rcp.xyz$z, 
                 c(min(rcp.xyz$x), max(rcp.xyz$x)), 
                 c(min(rcp.xyz$y), max(rcp.xyz$y)), 
                 c(min(rcp.xyz$z), max(rcp.xyz$z)))
  
  nnc <- nncross(bgnd, rcp.cut, k = 1:2)
  
  diff <- abs(nnc$dist.1 - nnc$dist.2)
  nkeep <- round(frac*npoints(bgnd))
  gbwhich <- order(diff)[1:nkeep]
  
  if(point.den != 1){
    nin <- round(point.den*length(gbwhich))
    nout <- length(gbwhich) - nin
    points.final <- sample(gbwhich, nin)
    
    nonplane.ind <- (1:npoints(bgnd))[-gbwhich]
    noise <- sample(nonplane.ind, nout)
    total.inds <- c(points.final, noise)
  }else{
    total.inds <- gbwhich
  }
  
  p.selected <- bgnd[total.inds]
  
  if(toplot == TRUE){
    plot3d.pp3(p.selected, col = 'red')
  }
  
  return(p.selected)
}





