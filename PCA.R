# practice PCA stuff
# https://www.datacamp.com/community/tutorials/pca-analysis-r

mtcars.pca <- prcomp(mtcars[,c(1:7,10,11)], center = TRUE, scale. = TRUE)
summary(mtcars.pca)

library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

ggbiplot(mtcars.pca)


# Simulated Training Data:
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/200131_params.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/200131_sim.res.RData')

# unpack data
data.unlisted <- matrix(NA, nrow = length(sim.res)*nrow(sim.res[[1]]), ncol = 10)

cnt <- 1
for(i in 1:length(sim.res)){
  for(j in 1:nrow(sim.res[[1]])){
    data.unlisted[cnt, 1:5] <- params[[i]]
    data.unlisted[cnt, 6:10] <- sim.res[[i]][j,]
    cnt <- cnt + 1
  }
}

data.unlisted <- as.data.frame(data.unlisted) 
names(data.unlisted) <- c("cr","rho1","rho2","rb","pb","Km","Rm","Rdm","Rddm","Kdm")
data.unlisted$sigma <- data.unlisted$rb * data.unlisted$cr
data.unlisted$rw <- (data.unlisted$cr^4 + 6* data.unlisted$cr^2 * data.unlisted$sigma^2 + 3*data.unlisted$sigma^4)/
  (data.unlisted$cr^3 + 3*data.unlisted$cr*data.unlisted$sigma^2)

#drop NA rows
data.unlisted <- data.unlisted[complete.cases(data.unlisted),]

data.train <- data.unlisted[,c("Km","Rm","Rdm","Rddm","Kdm")]

data.pca <- prcomp(data.train, center = TRUE, scale. = TRUE)
summary(data.pca)
ggbiplot(data.pca, labels = NA, xlim = c(-2,2))

data.pca$rotation

sdev <- data.pca$sdev
dvar <- sdev^2
prop.var <- dvar/sum(dvar)
plot(cumsum(prop.var), xlab = "prc", ylab = "percent", type = "b")

# no rb data
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/200131_params_norb.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/200131_sim.res_norb.RData')

# unpack data
data.unlisted <- matrix(NA, nrow = length(sim.res)*nrow(sim.res[[1]]), ncol = 10)
cnt <- 1
for(i in 1:length(sim.res)){
  for(j in 1:nrow(sim.res[[1]])){
    data.unlisted[cnt, 1:5] <- params[[i]]
    data.unlisted[cnt, 6:10] <- sim.res[[i]][j,]
    cnt <- cnt + 1
  }
}

data.unlisted <- as.data.frame(data.unlisted) 
names(data.unlisted) <- c("cr","rho1","rho2","rb","pb","Km","Rm","Rdm","Rddm","Kdm")
data.unlisted$rw <- data.unlisted$cr

#drop NA rows
data.unlisted <- data.unlisted[complete.cases(data.unlisted),]
data.train <- data.unlisted[,c("Km","Rm","Rdm","Rddm","Kdm")]
data.pca <- prcomp(data.train, center = TRUE, scale. = TRUE)
summary(data.pca)

# Home brewed plot
x <- data.pca$rotation[,1]
y <- data.pca$rotation[,2]

par(mfrow = c(1,1), mgp = c(3,1,0), mar = c(4.5,4.5,3,2.5))
summ <- summary(data.pca)
pc1 <- toString(round(summ$importance[2,1]*100, 1))
pc2 <- toString(round(summ$importance[2,2]*100, 1))
cols <- c("red", "blue", "forestgreen", "grey", "black")
plot(x,y, xlab = paste("Principle Component #1 - ", pc1, "% var. explained", sep = ''),
          ylab = paste("Principle Component #2 - ", pc2, "% var. explained", sep = ''),
     col = 'black', pch = 19, xlim = c(-0.5, 0.5), ylim = c(-0.9,0.9))
arrows(rep(0,5), rep(0,5), x, y , length=0.1, angle=45, code=2, col = 'black')
grid()
#legend(0.08, 0.05, legend = c("Km","Rm","Rdm","Rddm","Kdm"), col = cols, lwd = 2, lty = 1, bty = 'n')
#text(x, y, labels = c("Km","Rm","Rdm","Rddm","Kdm"))

