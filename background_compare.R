library(rapt)
library(scales)
library(data.table)

load("~/Research/oled-spatstat-analysis/bcc.res.RData")
load("~/Research/oled-spatstat-analysis/rcp.res.RData")
load("~/Research/oled-spatstat-analysis/rand.res.RData")
load("~/Research/oled-spatstat-analysis/bcc.offset.res.RData")

rvals <- bcc.res[,1]
tvals <- list()
tvals$bcc <- bcc.res[,-1]
tvals$rcp <- rcp.res[,-1]
tvals$rand <- rand.res[,-1]


percentile <- 0.999
ntests <- ncol(tvals[[1]])
prange <- percentile * ntests

tvals.sorted <- lapply(tvals, function(x){t(apply(x,1,sort))})

# Create toSub from RCP RRL
if(ntests %% 2 == 0){
  toSub <- (tvals.sorted$rcp[,(ntests/2)] + tvals.sorted$rcp[,(ntests/2 + 1)])/2
} else {
  toSub <- tvals.sorted$rcp[, ceiling(nt/2)]
}

tvals.sqrt <- lapply(tvals.sorted, function(x){sqrt(x) - sqrt((4/3)*pi*(rvals^3))})  

ind.big <- round(ntests/2) + floor(prange/2) 
ind.small <- round(ntests/2) - floor(prange/2)

plot.big <- lapply(tvals.sqrt, function(x){x[,ind.big]})
plot.small<- lapply(tvals.sqrt, function(x){x[,ind.small]})

colors <- c("black", "red", "blue")
par(mgp = c(2,1,0),mar = c(3.25,3.25,2,2))
plot(rvals, plot.big[[1]], type="n", main="99.9% AI Envelopes",
     xlab="r (arb.)", ylab=expression(sqrt('K'[3]*'(r)')*'  Anomaly'),
     ylim = c(-2.5,2.5), xlim = c(0, max(rvals)),
     cex.lab = 1, cex.main = 1, cex.axis = 1)
for(i in 2:length(plot.big)){
  a <- c(rvals, rev(rvals))
  b <- c(plot.big[[i]], rev(plot.small[[i]]))
  polygon(a, b, col=alpha(colors[i], 0.4), border="black", lwd=1)
}
lines(rvals, tvals.sqrt$bcc[,10], col = alpha(colors[1], alpha = 0.7), lwd = 1.25)
abline(h=0,lty=2,lwd=1,col="black")
legend(0, 2.75, legend = c("BCC", "RCP", "Random"),
       col = alpha(colors, alpha = 0.4), lty = 1, lwd = c(1.25,15,15), bty="n", cex = 1)




