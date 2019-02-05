helpme <- function(under,over,radius1,radius2,
                   type = "cr",
                   ppc=NULL,
                   cr=NULL,speed = "superfast",
                   d=NULL,
                   pic = 1,
                   pcp = 0.06,
                   den = 1,
                   gb = FALSE,
                   gbp = c(0,1),
                   gbmethod = 1,
                   rb = FALSE,
                   rbp,
                   rbmethod = 1,
                   s = 100,
                   toPlot=FALSE,showOverPts=FALSE){
  
  # If rb is true, we need to set a new cr to deal with spacing - Added 2/4/19
  if(rb == TRUE){
    if(rbmethod == 1){
      #scale cr to average *volume*
      cr <- (cr * (cr^2 + 3*rbp^2))^(1/3)
      #print(cr)
    }else if(rbmethod == 2){
      #scale cr to average *volume*
      cr <- (rbp[1]*rbp[2]^3 + (1-rbp[1])*rbp[3]^3)^(1/3)
    }
  }

  #real cluster percent
  set.seed(s)
  rcp <- pcp*pic
  
  under.r <- radius1
  over.r <- radius2
  over.vol <- volume(domain(over))
  under.vol <- volume(domain(under))
  
  # Calculate scaling factor for over point pattern
  # Added 2/4/19
    # There's some hairy math here that took me a while to figure out. Come see me if you need it explained, 
    # or see my reseatch notebook.
  z63 <- (((4/3)*pi*over.r^3)*npoints(over)*0.75 + ((4/3)*pi*(over.r*1.20)^3)*npoints(over)*0.25)/over.vol
  under.xdim <- domain(under)$xrange[2]
  under.ydim <- domain(under)$yrange[2]
  under.zdim <- domain(under)$zrange[2]
  
  innervol <- (under.xdim - 2 * cr)*(under.ydim - 2 * cr)*(under.zdim - 2 * cr)
  outervol <- (under.xdim + 2 * cr)*(under.ydim + 2 * cr)*(under.zdim + 2 * cr)
  middlevol <- outervol - innervol
  
  volfactor <- ((innervol/outervol) + (middlevol/outervol)*(5/12))
  
  over.rf <- cr*((under.xdim + 2*cr)*(under.ydim + 2*cr)*(under.zdim + 2*cr)*z63*volfactor/(under.vol * rcp * 1.182))^(1/3)
  
  if(rbmethod == 1){
    sdfactor <- (rbp/cr)*0.37
    over.rf <- over.rf + over.rf*sdfactor
  }

  #ncpred <- ((under.xdim + 2*cr) * (under.ydim + 2*cr) * (under.zdim + 2*cr) * z63)/((4/3)*pi*(0.75*(over.rf)^3 + 0.25*(1.2*over.rf)^3))
  #print(over.rf)
  over.sep <- over.rf*2
  over.scaled <- scaleRCP(over,newRadius = over.rf, oldRadius = over.r,win = domain(over))
  
  if(gb == TRUE){
    n <- npoints(over.scaled)
    gbval <- rgblur(n,gbp[1],gbp[2],coords = "rec", method = gbmethod)
    over.xyz <- coords(over.scaled)
    over.xyz.new <- over.xyz + gbval
    over.scaled.new <- createSpat(over.xyz.new, win = domain(over.scaled))
  }else{
    over.scaled.new <- over.scaled
  }
  
  # new addition as of 11/8/2018 - use cluster centers that can fall outside of the under pattern domain. This reduces
  # number of cluster points error significantly
  under.coo <- coords(under)
  under.new <- createSpat(under.coo + cr, win = domain(under))
  over.scaled.domain <- c(domain(under)$xrange[2]+2*cr,domain(under)$yrange[2]+2*cr,domain(under)$zrange[2]+2*cr)
  over.scaledf <- subSquare(over.scaled.new, over.scaled.domain)
  
  if(rb == TRUE){
    n <- npoints(over.scaledf)
    if(rbmethod == 1){
      crrand <- rnorm(n,mean = cr, sd = rbp)
      crrand[crrand < 0] <- 0
    } else if(rbmethod == 2){
      n1 <- round(n*rbp[1])
      n2 <- n - n1
      a <- c(rep(rbp[2], n1), rep(rbp[3], n2))
      crrand <- sample(a, replace = FALSE)
    }
  }
  
  # Deal with overlapping clusters here
  # Check for overlap, seperate to no overlap if so
  if(gb == TRUE){
    nnd <- nndist.pp3(over.scaledf)
    if(rb == TRUE){
      scrrand <- sort(crrand,decreasing = TRUE)
      comp <- scrrand[1] + scrrand[2]
    }else{
      comp <- 2*cr
    }
    
    if(any(nnd < comp)){
      check <- which(nnd < comp)
      while(!is.empty(check)) {
        nnw <- nnwhich.pp3(over.scaledf)
        direction <- (coords(over.scaledf)[nnw[check[1]],]-coords(over.scaledf)[check[1],])/nnd[check[1]]
        coords(over.scaledf)[check[1],] <- coords(over.scaledf)[check[1],]+(nnd[check[1]]-(comp+0.00001))*direction
        nnd <- nndist.pp3(over.scaledf)
        check <- which(nnd < comp)
      }
    }
  }
  
  #rb
  if(rb == TRUE){
    cluster.nnR.ind1 <- list()
    cluster.nnR.ind2 <- list()
    #t1 <- Sys.time()
    cluster.nnR.full <- crosspairs.pp3(over.scaledf,under.new,max(crrand),what="all",twice=FALSE,distinct=TRUE,neat=TRUE)
    cluster.nnR.split <- list()
    cluster.nnR.split$d <- split(cluster.nnR.full$d,cluster.nnR.full$i,drop=FALSE)
    cluster.nnR.split$j <- split(cluster.nnR.full$j,cluster.nnR.full$i,drop=FALSE)
    split.vals <- sort(unique(cluster.nnR.full$i))
    for(i in 1:length(split.vals)){
      cluster.nnR.ind2[[i]] <- cluster.nnR.split$j[[i]][cluster.nnR.split$d[[i]] < crrand[split.vals[i]]]
      cluster.nnR.ind1[[i]] <- rep(split.vals[i],length(cluster.nnR.ind2[[i]]))
    }
    #t2 <- Sys.time()
    #print(t2 - t1)
    cluster.nnR.new <- list()
    cluster.nnR.new[[1]] <- unlist(cluster.nnR.ind1)
    cluster.nnR.new[[2]] <- unlist(cluster.nnR.ind2)
  }else{
    cluster.nnR.new <- crosspairs.pp3(over.scaledf,under.new,cr,what="indices",twice=FALSE,distinct=TRUE,neat=TRUE)
  }
  
  ##### test #####
  # cluster.ind <- cluster.nnR.new[[2]]
  # more <- round(npoints(under.new)*rcp)-length(cluster.ind)
  # more.orig <- more
  # over.rf.orig <- over.rf
  # iter <- 0
  # 
  # while(abs(more) > 100){
  #   sf <- more * (-1.769e-04)
  #   print(sf)
  #   over.scaled.new <- scaleRCP(over.scaled.new, newRadius = (over.rf + sf), oldRadius = over.rf)
  #   over.scaledf <- subSquare(over.scaled.new, over.scaled.domain)
  #   over.rf <- over.rf + sf
  #   
  #   #rb
  #   if(rb == TRUE){
  #     cluster.nnR.ind1 <- list()
  #     cluster.nnR.ind2 <- list()
  #     #t1 <- Sys.time()
  #     cluster.nnR.full <- crosspairs.pp3(over.scaledf,under.new,max(crrand),what="all",twice=FALSE,distinct=TRUE,neat=TRUE)
  #     cluster.nnR.split <- list()
  #     cluster.nnR.split$d <- split(cluster.nnR.full$d,cluster.nnR.full$i,drop=FALSE)
  #     cluster.nnR.split$j <- split(cluster.nnR.full$j,cluster.nnR.full$i,drop=FALSE)
  #     split.vals <- sort(unique(cluster.nnR.full$i))
  #     for(i in 1:length(split.vals)){
  #       cluster.nnR.ind2[[i]] <- cluster.nnR.split$j[[i]][cluster.nnR.split$d[[i]] < crrand[split.vals[i]]]
  #       cluster.nnR.ind1[[i]] <- rep(split.vals[i],length(cluster.nnR.ind2[[i]]))
  #     }
  #     #t2 <- Sys.time()
  #     #print(t2 - t1)
  #     cluster.nnR.new <- list()
  #     cluster.nnR.new[[1]] <- unlist(cluster.nnR.ind1)
  #     cluster.nnR.new[[2]] <- unlist(cluster.nnR.ind2)
  #   }else{
  #     cluster.nnR.new <- crosspairs.pp3(over.scaledf,under.new,cr,what="indices",twice=FALSE,distinct=TRUE,neat=TRUE)
  #   }
  #   
  #   cluster.ind <- cluster.nnR.new[[2]]
  #   more <- round(npoints(under.new)*rcp)-length(cluster.ind)
  #   rm(cluster.nnR.ind1, cluster.nnR.ind2, cluster.nnR.full, cluster.nnR.split)
  #   iter <- iter + 1
  #   if(iter > 100){
  #     break
  #   }
  # }

  ##########
  
  if(den < 1 & den >= 0){
    cluster.ind.split <- split(cluster.nnR.new[[2]],cluster.nnR.new[[1]],drop=FALSE)
    cluster.ind.thinned <- lapply(cluster.ind.split,function(x){
      return(sample(x,round(den*length(x)),replace=FALSE))})
    cluster.ind <- unlist(cluster.ind.thinned)
    cluster.ind.split <- unlist(cluster.ind.split)
  }else {
    cluster.ind <- cluster.nnR.new[[2]]
    cluster.ind.split <- cluster.ind
  }
  
  more <- round(npoints(under.new)*pcp)-length(cluster.ind)
  if(more==0){
    
  }else if(more > 0){
    cluster.ind <- randomInsert(cluster.ind,more,npoints(under.new),s,cluster.ind.split)
  }else if(more < 0){
    cluster.ind <- randomTakeAway(cluster.ind,-1*more,npoints(under.new),s)
  }
  
  cluster.xyz <- coords(under)[cluster.ind,]
  cluster.xyz <- na.omit(cluster.xyz)
  cluster <- createSpat(cluster.xyz)
  
  over.scaledf.coo <- coords(over.scaledf)
  over.scaledf <- createSpat(over.scaledf.coo - cr,win = domain(under))
  
  cluster.nn <- nndist.pp3(over.scaledf, k = 1)
  
  if(toPlot==TRUE){
    plot3d.pp3(cluster,col="red",size=5)
    #plot3d.pp3(under,col="lightgray",add=TRUE)
    if(showOverPts==TRUE){
      plot3d.pp3(over.scaledf,size= 6,col="black",add=TRUE)
    }
  }
  
  if(rb == TRUE){
    return(list(cluster,over.scaledf,more,cluster.nn,crrand,ncpred))
  }else{
    return(list(cluster,over.scaledf,more,cluster.nn,0,ncpred))
  }
}