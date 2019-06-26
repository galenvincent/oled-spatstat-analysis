# practice PCA stuff
# https://www.datacamp.com/community/tutorials/pca-analysis-r

mtcars.pca <- prcomp(mtcars[,c(1:7,10,11)], center = TRUE, scale. = TRUE)
summary(mtcars.pca)

library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

ggbiplot(mtcars.pca)


library(data.table)

data.grid <- fread("~/Research/ML/190216_pca.csv", select = c(5:9))
data.rand <- fread("~/Research/ML/190304_pca.csv", select = c(5:9))  
data <- rbind(data.grid, data.rand)

p <- 0.75
samp <- sample(1:nrow(data),round(nrow(data)*p),replace = FALSE)

data.train <-  data[samp,]
data.test <- data[-samp,] 

colnames(data.train) <- c("Km","Rm","Rdm","Rddm","Kdm")
data.pca <- prcomp(data.train, center = TRUE, scale. = TRUE)
summary(data.pca)
ggbiplot(data.pca, labels = NA, xlim = c(-2,2))

data.pca$rotation

sdev <- data.pca$sdev
dvar <- sdev^2
prop.var <- dvar/sum(dvar)
plot(cumsum(prop.var), xlab = "prc", ylab = "percent", type = "b")

# no rb data
data <- fread("~/Research/ML/190617_norb/pca.csv",select = c(4,5,6,7,8))

p <- 0.75
samp <- sample(1:nrow(data),round(nrow(data)*p),replace = FALSE)

data.train <-  data[samp,]
data.test <- data[-samp,] 

colnames(data.train) <- c("Km","Rm","Rdm","Rddm","Kdm")
data.pca <- prcomp(data.train, center = TRUE, scale. = TRUE)


# Home brewed plot
x <- data.pca$rotation[,1]
y <- data.pca$rotation[,2]

summ <- summary(data.pca)
pc1 <- toString(round(summ$importance[2,1]*100, 1))
pc2 <- toString(round(summ$importance[2,2]*100, 1))
cols <- c("red", "blue", "forestgreen", "grey", "black")
plot(x,y, xlab = paste("Principle Component #1 - ", pc1, "% var. explained", sep = ''),
          ylab = paste("Principle Component #2 - ", pc2, "% var. explained", sep = ''),
     col = 'black', pch = 19, xlim = c(-0.5, 0.5), ylim = c(-0.8,0.8))
arrows(rep(0,5), rep(0,5), x, y , length=0.1, angle=45, code=2, col = 'black')
grid()
legend(0.08, 0.05, legend = c("Km","Rm","Rdm","Rddm","Kdm"), col = cols, lwd = 2, lty = 1, bty = 'n')
text(x, y, labels = c("Km","Rm","Rdm","Rddm","Kdm"))

