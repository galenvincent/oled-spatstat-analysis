# Morphology generation for Paul
library(rapt)
library(data.table)
library(dplyr)

# Random on grid
rand.grid <- function(dim = 60, guest.frac = 0.06, write.file = FALSE, filename = 'randgrid.txt'){
  # set up grid
  x <- 0:(dim-1)
  y <- 0:(dim-1)
  z <- 0:(dim-1)
  grid.xyz <- expand.grid(x, y, z)
  
  n <- nrow(grid.xyz)
  
  #random
  inds <- sample(1:n, round(n*guest.frac))
  marks <- rep(1, n)
  marks[inds] <- 2
  out.df <- data.frame(x = grid.xyz$Var3, y = grid.xyz$Var2, z = grid.xyz$Var1, type = marks)
  
  if(write.file == TRUE){
    fwrite(out.df, file = filename, row.names = FALSE, col.names = FALSE)
  }
  return(out.df)
}

# Inhibited
inhibited.grid <- function(dim = 60, type = 'regular', guest.frac = 0.06, write.file = FALSE, filename = 'randgrid.txt'){
  # set up grid
  x <- 0:(dim-1)
  y <- 0:(dim-1)
  z <- 0:(dim-1)
  grid.xyz <- expand.grid(x, y, z)
  
  n <- nrow(grid.xyz)
  nguest <- round(n*guest.frac)
  
  if(type == 'regular'){
    # Regular inhibited
    pts.per.side <- floor(nguest^(1/3))
    spacing <- ceiling(dim/pts.per.side)
    if(spacing < 2){
      print('Too high of density for this method')
      return()
    }
    
    inds.select <- seq(0, (dim-1), by = spacing)
    grid.select <- expand.grid(inds.select, inds.select, inds.select)
    grid.inds.select <- (grid.select$Var1*1 + grid.select$Var2*60 + grid.select$Var3*60^2)+1
    
    n.to.remove <- nrow(grid.select) - nguest
    inds.to.remove <- sample(1:nrow(grid.select), n.to.remove)
    
    grid.inds.select.final <- grid.inds.select[-inds.to.remove]
    
    marks <- rep(1, n)
    marks[grid.inds.select.final] <- 2
    out.df <- data.frame(x = grid.xyz$Var3, y = grid.xyz$Var2, z = grid.xyz$Var1, type = marks)
  }
  
  # More Random Inhibited
  else if(type == 'random'){
    grid.sums <- grid.xyz$Var1+grid.xyz$Var2+grid.xyz$Var3
    grid.take.tf <- grid.sums %% 2 == 0
    grid.take.inds <- which(grid.take.tf == TRUE)
    
    keeps <- grid.take.inds[sample(1:length(grid.take.inds), nguest)]
    
    marks <- rep(1, nrow(grid.xyz))
    marks[keeps] <- 2
    
    out.df <- data.frame(x = grid.xyz$Var3, y = grid.xyz$Var2, z = grid.xyz$Var1, type = marks)
  }
  
  if(write.file == TRUE){
    fwrite(out.df, file = filename, row.names = FALSE, col.names = FALSE)
  }
  return(out.df)
}

# Grid linked
linked.grid <- function(dim = 60, guest.frac = 0.06, write.file = FALSE, filename = 'randgrid.txt'){
  # set up grid
  x <- 0:(dim-1)
  y <- 0:(dim-1)
  z <- 0:(dim-1)
  grid.xyz <- expand.grid(x, y, z)
  names(grid.xyz) <- c('z', 'y', 'x')
  
  n <- nrow(grid.xyz)
  nguest <- round(n*guest.frac)
  
  f <- function(n, p, d){abs(3*d*n^2 - 2*n^3 - p)}
  a <- optimize(f, c(0,60), p = nguest, d = dim)
  ngrid <- floor(a[[1]])
  spacing <- ceiling(dim/ngrid)
  inds.select <- seq(0, dim-1, by = spacing)
  
  grid.selectxy <- grid.xyz[grid.xyz$x %in% inds.select & grid.xyz$y %in% inds.select,]
  grid.selectyz <- grid.xyz[grid.xyz$y %in% inds.select & grid.xyz$z %in% inds.select,]
  grid.selectxz <- grid.xyz[grid.xyz$x %in% inds.select & grid.xyz$z %in% inds.select,]
  
  grid.select.all <- rbind(grid.selectxy, grid.selectxz, grid.selectyz)
  grid.select.all.unique <- unique(grid.select.all)
  
  grid.inds.select <- (grid.select.all.unique$z*1 + grid.select.all.unique$y*60 + grid.select.all.unique$x*60^2)+1
  
  n.to.add <- nguest - length(grid.inds.select)
  
  marks <- rep(1, n)
  marks[grid.inds.select] <- 2
  
  out.pp3 <- pp3(grid.xyz$x, grid.xyz$y, grid.xyz$z, box3(c(0, dim), c(0, dim), c(0, dim)))
  marks(out.pp3) <- 1:npoints(out.pp3)
  
  guest.pp3 <- out.pp3[marks == 2]
  host.pp3 <- out.pp3[marks == 1]
  
  cp <- crosspairs(guest.pp3, host.pp3, rmax = 1.1, what = 'indices')
  
  inds.to.choose.from <- unique(marks(host.pp3[cp$j]))
  
  inds.to.add <- sample(inds.to.choose.from, n.to.add)
  marks[inds.to.add] <- 2
  
  out.df <- data.frame(x = grid.xyz$x, y = grid.xyz$y, z = grid.xyz$z, type = marks)
  
  
  if(write.file == TRUE){
    fwrite(out.df, file = filename, row.names = FALSE, col.names = FALSE)
  }
  return(out.df)
}

