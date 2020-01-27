# New cluster simulation code
library(rapt)

#### Testing 
# Upload RCPs for testing
under.ul <- read.rcp('~/Research/point_patterns/Final/FinalConfig1', '~/Research/point_patterns/Final/system1', scaleUp = TRUE, newRadius = 0.5)
over.ul <- read.rcp('~/Research/point_patterns/Final/FinalConfig4', '~/Research/point_patterns/Final/system4', scaleUp = TRUE, newRadius = 0.5)
under <- stitch.size(under.ul, domain(under.ul), c(60,60,60))
over <- stitch.size(over.ul, domain(over.ul), c(60,60,60))

# Clustersim test
n <- 40
cr <- runif(n, 2, 6.5)
rho1 <- runif(n, 0.2, 1)
rho2 <- runif(n, 0, 0.03)
rb <- runif(n, 0, 0.5)
pb <- runif(n, 0, 0.2)
s <- 1:n

params <- list()
for(i in 1:n){
  params[[i]] <- c(cr[i], rho1[i], rho2[i], rb[i], pb[i], s[i])
}

library(parallel)
cl <- makePSOCKcluster(8)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))
clusterExport(cl, c('params','under','over'))

a <- lapply(1:n, function(i){
  print(i)
  x <- clustersim(under, over, 0.5, pcp = 0.05114235, 
             params[[i]][1], 
             params[[i]][2], 
             params[[i]][3], 
             params[[i]][4], 
             params[[i]][5],
             tol = 0.005,
             s = params[[i]][6])
  if(is.list(x)){
    return(x[[4]])
  }else{
    return(-1)
  }
})
stopCluster(cl)
b <- unlist(a)


load('test.RData')
data <- matrix(unlist(a), ncol = 7, byrow = T)
colnames(data) <- c('pcp', 'cr', 'rho1', 'rho2', 'rb', 'pb', 's')
res <- as.data.frame(data)

hist(res$pcp)
hist(100*res$pcp[res$pcp != -1], breaks = seq(4, 6, by = 0.05), xlim = c(4.5, 5.6), main = '',xlab = 'percent cluster points')
abline(v = 5.114235, col = 'red', lwd = 2)

load('params.RData')


# Testing for parameter ranges
a <- clustersim(under, over, 0.5, pcp = 0.05114235,
                cr = 2.5,
                rho1 = 0.2, 
                rho2 = 0.025,
                rb = 0.5,
                pb = 0.2,
                s = 100,
                toplot = TRUE)




