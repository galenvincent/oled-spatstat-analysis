# Attemp at machine learning for predictions

library(data.table)
library(caret)
library(parallel)
library(doParallel)
library(ggbiplot)


#### Grid Data ####
dep <- fread("~/Research/ML/190216_pca.csv",select = c(1:4))
names(dep) <- c("r","den","rb","gb")
met <- fread("~/Research/ML/190216_pca.csv",select = c(5,6,7,8,9))
names(met) <- c("Km","Rm","Rdm","Rddm","Kdm")

pp <- preProcess(met, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
data.pca.vals <- predict(pp, met)

# Add random data if you want
deptemp <- fread("~/Research/ML/190304_pca.csv", select = c(1:4))
dep <- as.data.frame(rbind(as.matrix(dep), as.matrix(deptemp)))

mettemp <-  fread("~/Research/ML/190304_pca.csv",select = c(5,6,7,8,9))
met <- as.data.frame(rbind(as.matrix(met), as.matrix(mettemp)))
#

pp <- preProcess(met, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
data.pca.vals <- predict(pp, met)

rm(deptemp)
rm(mettemp)
#### ####


#### Klocal data import ####
dep <- fread("~/Research/ML/190311_pca.csv", select = c(1:4))
names(dep) <- c("r","den","rb","gb")
met <- fread("~/Research/ML/190311_pca.csv",select = c(5:9))#10))
names(met) <- c("Km","Rm","Rdm","Rddm","Kdm")#,"Kl")

pp <- preProcess(met, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
data.pca.vals <- predict(pp, met)
#### ####

#### Random Data ####
dep <- fread("~/Research/ML/190304_pca.csv", select = c(1:4))
names(dep) <- c("r","den","rb","gb")
met <- fread("~/Research/ML/190304_pca.csv",select = c(5:9))
names(met) <- c("Km","Rm","Rdm","Rddm","Kdm")

pp <- preProcess(met, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
data.pca.vals <- predict(pp, met)
#### ####

#### data to predict ####
topred.dep <- fread("~/Research/ML/190216_pca.csv", select = c(1:4))
names(topred.dep) <- c("r","den","rb","gb")
topred.met <- fread("~/Research/ML/190216_pca.csv", select = c(5:9))
names(topred.met) <- c("Km","Rm","Rdm","Rddm","Kdm")

# for grid data, get rid of the edge data...
cuts <- (topred.dep$r > 2 & topred.dep$r < 6.5 &
           topred.dep$den > 0.15 & topred.dep$den < 1 & 
           topred.dep$rb > 0 & topred.dep$rb < 0.5 & 
           topred.dep$gb > 0 & topred.dep$gb< 0.6)
topred.chop.met <- topred.met[cuts,]
topred.chop.dep <- topred.dep[cuts,]
topred.chop.pca <- predict(pp, topred.chop.met)
topred.chop <- data.frame(topred.chop.met, topred.chop.dep, topred.chop.pca)
topred.chop$rw <- (topred.chop$r^4 + 6* (topred.chop$r^2) *(topred.chop$rb * topred.chop$r)^2 + 3 *(topred.chop$rb * topred.chop$r)^4)/(topred.chop$r^3 + 3*topred.chop$r*(topred.chop$rb * topred.chop$r)^2)
topred.chop$sigma <- topred.chop$rb * topred.chop$r
#


topred.pca <- predict(pp, topred.met)
topred <- data.frame(topred.met, topred.dep, topred.pca)

topred$rw <- (topred$r^4 + 6* (topred$r^2) *(topred$rb * topred$r)^2 + 3 *(topred$rb * topred$r)^4)/(topred$r^3 + 3*topred$r*(topred$rb * topred$r)^2)
topred$sigma <- topred$rb * topred$r
#### ####

## Old pca
#data.pca <- prcomp(met, center = TRUE, scale. = TRUE)
#ggbiplot(data.pca, labels = NA)
#summary(data.pca)
#data.pca.vals <- as.data.frame(data.pca$x, names = c("PC1","PC2","PC3","PC4","PC5"))#,"PC6"))

data <- data.frame(met, dep, data.pca.vals)

data$rw <- (data$r^4 + 6* (data$r^2) *(data$rb * data$r)^2 + 3 *(data$rb * data$r)^4)/(data$r^3 + 3*data$r*(data$rb * data$r)^2)
data$sigma <- data$rb * data$r


p <- 0.75
samp <- sample(1:nrow(dep),round(nrow(dep)*p),replace = FALSE)

trainSet <- data[samp,]
testSet <- data[-samp,] 


out <- 'den'
#predictors <- c("Km","Rm","Rdm","Rddm","Kdm")
predictors <- c("PC1","PC2")#,"PC3","PC4")#,"PC5")#,"PC6")

models <- list()
gc()


# Based on all
models[[1]] <- train(trainSet[,predictors],trainSet[,out],method = "glm",trControl = trainControl("cv", number = 10))

models[[2]] <- train(trainSet[,predictors],trainSet[,out],method = "knn",preProcess = c("center","scale"), tuneLength = 10)

models[[3]] <- train(trainSet[,predictors], trainSet[,out], 
               method = "glmnet", 
               trControl = trainControl("cv", number = 10),
               tuneLength = 10)
models[[4]] <- train(trainSet[,predictors],trainSet[,out],method = "brnn", trControl = trainControl("cv", number = 10))

cl <- makeCluster(detectCores())
registerDoParallel(cl)
models[[5]] <- train(trainSet[,predictors],trainSet[,out],method = "parRF", trControl = trainControl("cv", number = 10))
stopCluster(cl)
registerDoSEQ()

names(models) <- c("lm","knn","elastic-net","brnn","rf")

compall <- extractPrediction(models, testX = testSet[,predictors], testY = testSet[,out])
comp <- split(compall, compall$model)

#### predict different data ####
compall_pred <- extractPrediction(models, testX = topred[,predictors], testY = topred[,out])
comp_pred <- split(compall_pred, compall_pred$model)

# predict on the cut data, yo
compall_pred.chop <- extractPrediction(models, testX = topred.chop[,predictors], testY = topred.chop[,out])
comp_pred.chop <- split(compall_pred.chop, compall_pred.chop$model)

#Test or train "Test" or "Training"
tt <- "Test"

lapply(comp, function(x){
  xn <- x[x$dataType == tt,]
  ds <- defaultSummary(xn)
  
  plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
  abline(0,1, col = "red", lwd = 1.25)
  return(ds)
})

#plot for poster
xn <- comp[[1]][comp[[1]]$dataType == tt,]
plot(xn$obs, xn$pred, main = "Cluster Weighted Radius Model", xlab = "True Value", ylab = "Predicted Value", xlim = c(0,10),ylim = c(0,11))
abline(0,1, col = "red", lwd = 1.25)
legend(0, 11, "True Value = Predicted Value", col = "red", lty = 1, lwd = 2,bty = "n")

# Plot the predictions
lapply(comp_pred, function(x){
  xn <- x[x$dataType == "Test",]
  ds <- defaultSummary(xn)
  
  plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
  abline(0,1, col = "red", lwd = 1.25)
  return(ds)
})

# Plot the predictions - the chopped off edges of the grid data to test the 
# theory that the random model does bad predicting around the "edges" of the dataset
for(i in 1:length(models)){
  xn <- comp_pred.chop[[i]][comp_pred.chop[[i]]$dataType == "Test",]
  ds <- defaultSummary(xn)
  xo <- comp_pred[[i]][comp_pred[[i]]$dataType == "Test",]
  
  plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value", type = 'n')
  points(xo$obs, xo$pred, col = "red")
  points(xn$obs, xn$pred, col = "black")
  abline(0,1, col = "green", lwd = 1.25)
  print(toString(xn$model[1]))
  print(ds)
}



