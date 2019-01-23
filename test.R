helpme <- function(under,over,radius1,radius2,
                   type = "ppc",
                   ppc=NULL,
                   cr=NULL,speed = "superfast",
                   d=NULL,
                   pic = 1,
                   pcp = 0.06,
                   den = 1,
                   gb = FALSE,
                   gbp = c(0,1),
                   s = 100,
                   toPlot=FALSE,showOverPts=FALSE){
rcp <- pcp*pic

under.r <- radius1
over.r <- radius2
under.vol <- volume(domain(under))

over.rf <- under.r*cr*((4*pi*npoints(under))/(3*under.vol*rcp))^(1/3)
over.sep <- over.rf*2
over.scaled <- scaleRCP(over,newRadius = over.rf, oldRadius = over.r,win = domain(over))


#what happends when clusters get into overlap land? need to figure this out
if(gb == TRUE){
  n <- npoints(over.scaled)
  gbval <- rgblur(n,gbp[1],gbp[2],coords = "rec")
  over.xyz <- coords(over.scaled)
  over.xyz.new <- over.xyz + gbval
  over.scaled.new <- createSpat(over.xyz.new, win = domain(over.scaled))
}else{
  over.scaled.new <- over.scaled
}
#browser()

#coords(under) <- coords(under) + cr
#over.scaled.domain <- c(domain(under)$xrange[2]+2*cr,domain(under)$yrange[2]+2*cr,domain(under)$zrange[2]+2*cr)
over.scaledf <- subSample(under,over.scaled.new)
# Deal with overlapping clusters here
# Check for overlap, seperate to no overlap if so
# nnd <- nndist.pp3(over.scaledf)
# if(any(nnd < 2*cr)){
#   nnw <- nnwhich.pp3(over.scaledf)
#
#   while(any(nnd < 2*cr)) {
#     check <- which(nnd < 2*cr)
#     for(i in 1:check) {
#         direction <- (coords(over.scaledf)[nnw[i],]-coords(over.scaledf)[i,])/nnd[i]
#         coords(over.scaledf)[i,] <- coords(over.scaledf)[i,]+(nnd[i]-(2*cr+0.00001))*direction
#     }
#     nnw <- nnwhich.pp3(over.scaledf)
#     nnd <- nndist.pp3(over.scaledf)
#   }
# }


cluster.nnR.new <- crosspairs.pp3(over.scaledf,under,cr,what="indices",twice=FALSE,distinct=TRUE,neat=TRUE)

if(den < 1 & den >= 0){
  set.seed(s)
  cluster.ind.split <- split(cluster.nnR.new[[2]],cluster.nnR.new[[1]],drop=FALSE)
  cluster.ind.thinned <- lapply(cluster.ind.split,function(x){
    return(sample(x,round(den*length(x)),replace=FALSE))})
  cluster.ind <- unlist(cluster.ind.thinned)
  cluster.ind.split <- unlist(cluster.ind.split)
}else {
  cluster.ind <- cluster.nnR.new[[2]]
  cluster.ind.split <- cluster.ind
}

more <- round(npoints(under)*pcp)-length(cluster.ind)
if(more==0){

}else if(more > 0){
  cluster.ind <- randomInsert(cluster.ind,more,npoints(under),s,cluster.ind.split)
}else if(more < 0){
  cluster.ind <- randomTakeAway(cluster.ind,-1*more,npoints(under),s)
}

cluster.xyz <- coords(under)[cluster.ind,]
cluster.xyz <- na.omit(cluster.xyz)
#cluster.xyz <- cluster.xyz - cr
cluster <- createSpat(cluster.xyz)

#coords(under) <- coords(under) - cr
#coords(over.scaledf) <- coords(over.scaledf) - cr

if(toPlot==TRUE){
  plot3d.pp3(cluster,col="red",size=5)
  plot3d.pp3(under,col="lightgray",add=TRUE)
  if(showOverPts==TRUE){
    plot3d.pp3(over.scaledf,size= 6,col="black",add=TRUE)
  }
}

return(list(cluster,over.scaledf,more,over.sep))}
