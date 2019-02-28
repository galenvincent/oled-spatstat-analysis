# Attemp at machine learning for predictions

library(data.table)
library(caret)
library(parallel)
library(doParallel)

dep <- fread("~/Research/PCA/190216_pca.csv",select = c(1:4))
names(dep) <- c("r","den","rb","gb")
met <- fread("~/Research/PCA/190216_pca.csv",select = c(5:9))
names(met) <- c("Km","Rm","Rdm","Rddm","Kdm")

data.pca <- prcomp(met, center = TRUE, scale. = TRUE)
data.pca.vals <- as.data.frame(data.pca$x, names = c("PC1","PC2","PC3","PC4","PC5"))

data <- data.frame(met, dep, data.pca.vals)

data$rw <- (data$r^4 + 6* (data$r^2) *(data$rb * data$r)^2 + 3 *(data$rb * data$r)^4)/(data$r^3 + 3*data$r*(data$rb * data$r)^2)
data$sigma <- data$rb * data$r

p = 1
samp <- sample(1:nrow(dep),round(nrow(dep)*p),replace = FALSE)

trainSet <- data[samp,]
testSet <- data[-samp,] 


out <- 'sigma'
#predictors <- c("Km","Rm","Rdm","Rddm","Kdm")
predictors <- c("PC1","PC2","PC3","PC4","PC5")

models <- list()

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

poly <- lm(trainSet[,out] ~ poly(trainSet[,predictors[1]],
                                 trainSet[,predictors[2]],
                                 trainSet[,predictors[3]],
                                 trainSet[,predictors[4]],
                                 trainSet[,predictors[5]], 
                                 degree = 2, raw = TRUE))

names(models) <- c("lm","knn","elastic-net","brnn","rf")

compall <- extractPrediction(models)#, testX = testSet[,predictors], testY = testSet[,out])
comp <- split(compall, compall$model)
comp[[length(models) + 1]] <- data.frame(obs = trainSet[,out], pred = predict(poly), model = "polylm", dataType = "Training", object = "poly")
names(comp)[length(models) + 1] <- "polylm"

lapply(comp, function(x){
  ds <- defaultSummary(x)
  
  plot(x$obs, x$pred, main = toString(x$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
  abline(0,1, col = "red", lwd = 1.25)
  return(ds)
})


plot(varImp(models[[1]]))
plot(varImp(models[[2]]))
plot(varImp(models[[3]]))
