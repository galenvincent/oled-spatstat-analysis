# Envelope width analysis
library(data.table)
library(rapt)
library(pracma)

percentile <- .999
plot.big <- matrix(NaN,200,5)
plot.small <- matrix(NaN,200,5)
rvals <- matrix(NaN,200,5)

#### RRL Envelopes only ####
for(i in (1:5)){
  rtemp <- fread(paste('~/Research/box_size_effect/RRL',toString(i),'.csv',sep=""),select=2) # r values
  tvals <- fread(paste('~/Research/box_size_effect/RRL',toString(i),'.csv',sep=""),drop=c(1,2)) # results 
  nTests <- ncol(tvals) # number of tests done
  prange <- percentile*nTests # get the range of indeces for which each percentile spans
  
  sortedtVals <- t(apply(tvals,1,sort)) # sort the results at each r value from lowest to highest
  
  ind.big <- round(nTests/2)+floor(prange/2) # select the high end indexes based on being 1/2 of the percentile span from the middle of the tests
  ind.small <- round(nTests/2)-floor(prange/2) # do the same for the low end
  
  plot.big[,i] <- sortedtVals[,ind.big]
  plot.small[,i] <- sortedtVals[,ind.small]
  rvals[,i] <- as.numeric(rtemp$V1)
  
  rm(rtemp,tvals,nTests,prange,sortedtVals,ind.big,ind.small)
  gc()
}


rs <- seq(1.15,9.25,l=100)
widths <- matrix(NaN,length(rs),ncol(plot.big))

for(j in (1:length(rs))){
  for(i in (1:ncol(plot.big))){
    rind <- which(abs(rs[j]-rvals[,i]) == min(abs(rs[j]-rvals[,i])))
    widths[j,i] <- plot.big[rind,i] #- plot.small[rind,i]
  }
}

x1 <- c(15,20,30,40,60)
plot(x1,widths[9,1:5],pch = 1, col="blue",xlab = "Cube Side Length [nm]",ylab = "Envelope Width",main = "Cubic Volumes",ylim = c(0,max(widths[9,1:5])))

ns <- list() 
ns$val<- vector("numeric",length(rs))
ns$eb <- vector("numeric",length(rs))
models <- list()
for(i in 1:length(rs)){
  data <- widths[i,]
  rstar <- rep(rs[i],5)
  models[[i]] <- nls(data~a*x1^(-n)*rstar^n,start = list(n = 1.5, a = 12))
  #model <- lm(log(widths[i,1:5])~log(x1))
  ns$val[i] <- coefficients(models[[i]])[1]
  ns$eb[i] <- coefficients(summary(models[[i]]))[1,2]
}
plot(rs,ns$val,ylab = "Polynomial Decay Rate (n)", xlab= "Radius", pch = 19,
     main = "Envelope Width Decay Rates vs Radius",
     ylim = c(1,2), xlim = c(0,rs[length(rs)]))
arrows(rs, ns$val-ns$eb, rs, ns$val+ns$eb, length=0.05, angle=90, code=3)

x <- seq(0,60,l=100)
plot(x1,widths[100,])
lines(x,predict(models[[100]],list(x1 = x)))

sl <- c(15,20,30,40,60)
#r/l vector
rss <- rvals[,1]/sl[1]

#### skip if you want
cutoff <- 0.5
rss.fit <- rss[rss < cutoff]
n <- length(rss.fit)
rsd <- matrix(0,n,5)
for(i in 1:5){
  rsd[,i] <- rss.fit - 1/sl[i]
}

plot(rsd[,5],plot.big[rss < cutoff,5],type="n")
x <- vector("numeric",5*n)
y <- vector("numeric",5*n)
for(i in 1:5){
  lines(rsd[,i],plot.big[rss < cutoff,i] - 42*sl[i]^(-1.5),col=i)
  x[(n*(i-1)+1):(n*i)] <- rsd[,i]
  y[(n*(i-1)+1):(n*i)] <- plot.big[rss < cutoff,i] - 42*sl[i]^(-1.5)
}
xnz <- x[x >= 0]
ynz <- y[x >= 0]
plot(xnz,ynz)
mod.last <- nls(ynz~a*xnz^n, start = list(a = 12, n = 1.8))
xplot <- seq(0,0.5,l=100)
lines(xplot, predict(mod.last,list(xnz = xplot)),lwd = 2,col="red")
####

cutoff <- 0.5
rss.fit <- rss[rss < cutoff]
n <- length(rss.fit)
data <- vector("numeric",n*5)
for(i in 1:5){
  data[(n*(i-1)+1):(n*i)] <- plot.big[rss < cutoff,i]
}
rss.fit.big <- rep(rss.fit,5)
sl.fit <- c(rep(15,n),rep(20,n),rep(30,n),rep(40,n),rep(60,n))
# get rid of zeros 
data.nz <- data[data != 0]
rss.fit.big.nz <- rss.fit.big[data != 0]
sl.fit.nz <- sl.fit[data != 0]
plot(rss.fit.big.nz - 1/sl.fit.nz,data.nz,type="l")

plot(rss.fit, plot.big[rss < cutoff,1],type = "n",ylim = c(0,max(plot.big[rss < cutoff,])))
for(i in 1:5){
  lines(rss.fit - 1/rep(sl[i],n), plot.big[rss < cutoff,i],col=i)
}

m <- nls(data.nz ~ a *(rss.fit.big.nz - 1/sl.fit.nz)^n + b*(sl.fit.nz)^(-1.5), start = list(a = 12, n = 1.8, b = 42))
summary(m)

x <- seq(0,1,l=100)
plot(rss.fit,plot.big[rss < cutoff,1],ylim = c(0, max(plot.big[rss < cutoff,])),type="n")
for(i in 1:5){
  lines(rss.fit,plot.big[rss < cutoff,i],col=i)
  y <- rep(sl[i],100)
  lines(x,predict(m,list(rss.fit.big.nz = x, sl.fit.nz = y)),lwd = 2, col="red")
}

########################################
rs <- seq(1.15,6,l=100)
widths <- matrix(NaN,length(rs),ncol(plot.big))

for(j in (1:length(rs))){
  for(i in (1:ncol(plot.big))){
    rind <- which(abs(rs[j]-rvals[,i]) == min(abs(rs[j]-rvals[,i])))
    widths[j,i] <- plot.big[rind,i] #- plot.small[rind,i]
  }
}

x1 <- c(15,20,30,40,60)
plot(x1,widths[9,1:5],pch = 1, col="blue",xlab = "Cube Side Length [nm]",ylab = "Envelope Width",main = "Cubic Volumes",ylim = c(0,max(widths[9,1:5])))

ns <- list()
ns$val<- vector("numeric",length(rs))
ns$eb <- vector("numeric",length(rs))
models <- list()
for(i in 1:length(rs)){
  data <- widths[i,]
  rstar <- rep(rs[i],5)
  #x <- seq(0,60,l=100)
  #plot(x1,data)
  #lines(x, 12.74882*((rstar/x)-(1/x))^1.8 + 42.02498*x^(-1.55))
  models[[i]] <- nls(data~12.74882*((rstar/x1)-(1/x1))^n + 42.02498*x1^(-1.6), start = list(n = 1.81))
  #lines(x,predict(models[[i]],list(x1 = x)))
  ns$val[i] <- coef(summary(models[[i]]))[1]
  ns$eb[i] <- coefficients(summary(models[[i]]))[1,2]
}
plot(rs,ns$val,ylab = "Polynomial Decay Rate (n)", xlab= "Radius", pch = 19,
     main = "Envelope Width Decay Rates vs Radius",
     ylim = c(1,2), xlim = c(0,rs[length(rs)]))
arrows(rs, ns$val-ns$eb, rs, ns$val+ns$eb, length=0.05, angle=90, code=3)

data <- ns$val
mod.n <- nls(data~-a*(rs-c)^-n+d, start = list(n = 1.5,c=.5,d=1.8,a=.5))
summary(mod.n)
x <- seq(1.5,6,l=100)
lines(x, -0.75*(x-.3)^(-3)+1.8)
lines(x,predict(mod.n,list(rs = x)))


#### Aside - see what the function of shift is as a function of l ####
ris <- seq(0.05,0.25,l=100)
nis <- list() 
nis$val <- vector("numeric",100)
nis$er <-vector("numeric",100)
rsd <- matrix(0,n,5)
for(i in 1:5){
  rsd[,i] <- rss.fit - 1/sl[i]
}
for(i in 1:100){
  ri <- ris[i]
  rind <- vector("numeric",5)
  widths2 <- vector("numeric",5)
  for(j in 1:5){
    rind[j] <- which(abs(rsd[,j] - ri) == min(abs(rsd[,j] - ri)))
    widths2[j] <- plot.big[rss < cutoff,1:5][rind[j],j]
  }
  
  #plot(sl,widths2)
  #mod.hope <- nls(widths2~a*exp(-sl*n)+c,start = list(a = 1.4, n = 0.1, c = min(widths2)))
  mod.hope <- nls(widths2~a*(sl)^(-n)+c,start = list(a = 60, n = 1.5, c = min(widths2)))
  #summary(mod.hope)
  #x <- seq(0,60,l=100)
  #lines(x,predict(mod.hope,list(sl = x)))
  nis$val[i] <- coef(summary(mod.hope))[2,1]
  nis$er[i] <-coef(summary(mod.hope))[2,2]
}
abline(v = ri)
points(rep(ri,5),widths2)


plot(ris,nis$val)
arrows(ris, nis$val-nis$er, ris, nis$val+nis$er, length=0.05, angle=90, code=3)
nis.mean <- weighted.mean(x=nis$val,w=1/(nis$er))
abline(h = nis.mean)
