# practice PCA stuff
# https://www.datacamp.com/community/tutorials/pca-analysis-r

mtcars.pca <- prcomp(mtcars[,c(1:7,10,11)], center = TRUE, scale. = TRUE)
summary(mtcars.pca)

library(ggbiplot)

ggbiplot(mtcars.pca)


library(data.table)

data <- fread("~/Research/PCA/190216_pca.csv")
rsm <- sample(1:10000, replace = FALSE)
data.test <- data[rsm[1:1000],]
data.train <- data[rsm[1001:10000],]

data.train.pred <- data.train[,5:9]

colnames(data.train.pred) <- c("Km","Rm","Rdm","Rddm","Kdm")
data.pca <- prcomp(data.train.pred, center = TRUE, scale. = TRUE)
summary(data.pca)
ggbiplot(data.pca)

data.pca$rotation

sdev <- data.pca$sdev
dvar <- sdev^2
prop.var <- dvar/sum(dvar)
plot(cumsum(prop.var), xlab = "prc", ylab = "percent", type = "b")

train.data <- data.frame(data.train[,1], data.train[5:9])



