# Attemp at machine learning for predictions

library(data.table)
library(caret)
library(parallel)
library(doParallel)
library(ggbiplot)
library(scales)

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
## ##





#### Grid + Random no rb Data ####
dep <- fread("~/Research/ML/190617_norb/pca.csv",select = c(1:3))
names(dep) <- c("r","den","gb")
met <- fread("~/Research/ML/190617_norb/pca.csv",select = c(4,5,6,7,8))
names(met) <- c("Km","Rm","Rdm","Rddm","Kdm")

pp <- preProcess(met, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
data.pca.vals <- predict(pp, met)

## ##





#### Grid + Random norb FULL Data ####
dep.i <- fread("~/Research/ML/190617_norb/tot.csv")
dep <- dep.i
for(i in 2:101){
  dep <- rbind(dep, dep.i)
}
names(dep) <- c("r","den","gb")
rm(dep.i)

met <- fread("~/Research/ML/190617_norb/pca_full/full_1.csv")
for(i in 2:101){
  met <- rbind(met, fread(paste("~/Research/ML/190617_norb/pca_full/full_", toString(i),".csv", sep = "")))
}
names(met) <- c("Km","Rm","Rdm","Rddm","Kdm")

pp <- preProcess(met, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
data.pca.vals <- predict(pp, met)

# get rid of bad data
met <- met[!(rowSums(is.na(data.pca.vals)) > 0),]
dep <- dep[!(rowSums(is.na(data.pca.vals)) > 0),]
data.pca.vals <- data.pca.vals[!(rowSums(is.na(data.pca.vals)) > 0),]

## ##

#### Klocal data import ####
dep <- fread("~/Research/ML/190311_pca.csv", select = c(1:4))
names(dep) <- c("r","den","rb","gb")
met <- fread("~/Research/ML/190311_pca.csv",select = c(5:9))#10))
names(met) <- c("Km","Rm","Rdm","Rddm","Kdm")#,"Kl")

pp <- preProcess(met, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
data.pca.vals <- predict(pp, met)
## ##



#### Random Data ####
dep <- fread("~/Research/ML/190304_pca.csv", select = c(1:4))
names(dep) <- c("r","den","rb","gb")
met <- fread("~/Research/ML/190304_pca.csv",select = c(5:9))
names(met) <- c("Km","Rm","Rdm","Rddm","Kdm")

pp <- preProcess(met, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
data.pca.vals <- predict(pp, met)
## ##


#### Full raw data ------------------------------------------------------
load('Z:/Galen/Machine\ Learning\ Files/190710_params.RData')
load('Z:/Galen/Machine\ Learning\ Files/190710_sim.res.RData')

data.unlisted <- matrix(NA, nrow = length(sim.res)*nrow(sim.res[[1]]), ncol = 9)

cnt <- 1
for(i in 1:length(sim.res)){
  for(j in 1:nrow(sim.res[[1]])){
    data.unlisted[cnt, 1:4] <- params[[i]]
    data.unlisted[cnt, 5:9] <- sim.res[[i]][j,]
    cnt <- cnt + 1
  }
}

data.unlisted <- as.data.frame(data.unlisted)
names(data.unlisted) <- c("r","den","rb","gb","Km","Rm","Rdm","Rddm","Kdm")
data.unlisted$rw <- (data.unlisted$r^4 + 6* (data.unlisted$r^2) *(data.unlisted$rb * data.unlisted$r)^2 + 3 *(data.unlisted$rb * data.unlisted$r)^4)/(data.unlisted$r^3 + 3*data.unlisted$r*(data.unlisted$rb * data.unlisted$r)^2)
data.unlisted$sigma <- data.unlisted$rb * data.unlisted$r

test.ind <- sample(1:nrow(data.unlisted), 0.8*nrow(data.unlisted))

data.test <- data.unlisted[test.ind,]
data.train <- data.unlisted[-test.ind,]

#data.train <- data.unlisted

#drop NA rows
data.test <- data.test[!is.na(data.test$Km) & !is.na(data.test$Rdm) & !is.na(data.test$Rddm),]
data.train <- data.train[!is.na(data.train$Km) & !is.na(data.train$Rdm) & !is.na(data.train$Rddm),]

mets <- c("Km","Rm","Rdm","Rddm","Kdm")

pp.train <- preProcess(data.train[,mets], method = c("scale", "center", "pca", "BoxCox"), thresh = 1)

data.train.pca <- predict(pp.train, data.train[,mets])
data.test.pca <- predict(pp.train, data.test[,mets])

trainSet <- data.frame(data.train, data.train.pca)
testSet <- data.frame(data.test, data.test.pca)
#

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
## ##


#### Old pca ####
#data.pca <- prcomp(met, center = TRUE, scale. = TRUE)
#ggbiplot(data.pca, labels = NA)
#summary(data.pca)
#data.pca.vals <- as.data.frame(data.pca$x, names = c("PC1","PC2","PC3","PC4","PC5"))#,"PC6"))


#### ML ####
data <- data.frame(met, dep, data.pca.vals)

# only add this if there is rb data
data$rw <- (data$r^4 + 6* (data$r^2) *(data$rb * data$r)^2 + 3 *(data$rb * data$r)^4)/(data$r^3 + 3*data$r*(data$rb * data$r)^2)
data$sigma <- data$rb * data$r


p <- 0.5
samp <- sample(1:nrow(dep),round(nrow(dep)*p),replace = FALSE)

trainSet <- data[samp,]
testSet <- data[-samp,] 


out <- 'rw'
#predictors <- c("Km","Rm","Rdm","Rddm","Kdm")
predictors <- c("PC1","PC2","PC3","PC4","PC5")

models <- list()
gc()


# Based on all
models[[1]] <- train(trainSet[,predictors], trainSet[,out], method = "glm", trControl = trainControl("cv", number = 10))

models[[2]] <- train(trainSet[,predictors], trainSet[,out], method = "knn", tuneLength = 10)

models[[3]] <- train(trainSet[,predictors], trainSet[,out], 
               method = "glmnet", 
               trControl = trainControl("cv", number = 10),
               tuneLength = 10)

models[[4]] <- train(trainSet[,predictors], trainSet[,out], method = "brnn", trControl = trainControl("cv", number = 10),
                     tuneGrid = data.frame("neurons" = c(1, 3, 5, 7)))

# cl <- makeCluster(detectCores())
# registerDoParallel(cl)
# models[[5]] <- train(trainSet[,predictors],trainSet[,out],method = "parRF", trControl = trainControl("cv", number = 10))
# stopCluster(cl)
# registerDoSEQ()

names(models) <- c("lm","knn","elastic-net","brnn")#,"rf")

compall <- extractPrediction(models, testX = testSet[,predictors], testY = testSet[,out])
comp <- split(compall, compall$model)

#### predict different data ####
compall_pred <- extractPrediction(models, testX = topred[,predictors], testY = topred[,out])
comp_pred <- split(compall_pred, compall_pred$model)

#### predict on the cut data, yo ####
compall_pred.chop <- extractPrediction(models, testX = topred.chop[,predictors], testY = topred.chop[,out])
comp_pred.chop <- split(compall_pred.chop, compall_pred.chop$model)

#### See model results and comparison ####
# "Test" or "Train" - Plot the models vs the test data or on the training data
tt <- "Test"

lapply(comp, function(x){
  xn <- x[x$dataType == tt,]
  ds <- defaultSummary(xn)
  
  plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
  abline(0,1, col = "red", lwd = 1.25)
  return(ds)
})

#### plot for poster ####
xn <- comp[[1]][comp[[1]]$dataType == tt,]
plot(xn$obs, xn$pred, main = "Cluster Radius Model", xlab = "True Value", ylab = "Predicted Value", xlim = c(2,7),ylim = c(2,7))
abline(0,1, col = "red", lwd = 1.25)
legend(1.8, 7, "True Value = Predicted Value", col = "red", lty = 1, lwd = 2,bty = "n")

#### Plot the predictions ####
lapply(comp_pred, function(x){
  xn <- x[x$dataType == "Test",]
  ds <- defaultSummary(xn)
  
  plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
  abline(0,1, col = "red", lwd = 1.25)
  return(ds)
})

#### Plot the chopped predictions ####
# Look at the chopped off edges of the grid data to test the 
# theory that the random model does bad predicting around the 
# "edges" of the dataset
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





#### Paper ML work ####
# Load training data and calculate PCAs: Select either rb or norb files to load
load('Z:/Galen/Machine\ Learning\ Files/190710_params_norb.RData')
load('Z:/Galen/Machine\ Learning\ Files/190710_sim.res_norb.RData')

# RB
# data.unlisted <- matrix(NA, nrow = length(sim.res)*nrow(sim.res[[1]]), ncol = 9)
# 
# cnt <- 1
# for(i in 1:length(sim.res)){
#   for(j in 1:nrow(sim.res[[1]])){
#     data.unlisted[cnt, 1:4] <- params[[i]]
#     data.unlisted[cnt, 5:9] <- sim.res[[i]][j,]
#     cnt <- cnt + 1
#   }
# }

#NORB
data.unlisted <- matrix(NA, nrow = length(sim.res)*nrow(sim.res[[1]]), ncol = 8)

cnt <- 1
for(i in 1:length(sim.res)){
  for(j in 1:nrow(sim.res[[1]])){
    data.unlisted[cnt, 1:3] <- params[[i]]
    data.unlisted[cnt, 4:8] <- sim.res[[i]][j,]
    cnt <- cnt + 1
  }
}

data.unlisted <- as.data.frame(data.unlisted) 
#names(data.unlisted) <- c("r","den","rb","gb","Km","Rm","Rdm","Rddm","Kdm") #RB
names(data.unlisted) <- c("r","den","gb","Km","Rm","Rdm","Rddm","Kdm") #NORB

#rb only
data.unlisted$rw <- (data.unlisted$r^4 + 6* (data.unlisted$r^2) *(data.unlisted$rb * data.unlisted$r)^2 + 3 *(data.unlisted$rb * data.unlisted$r)^4)/(data.unlisted$r^3 + 3*data.unlisted$r*(data.unlisted$rb * data.unlisted$r)^2)
data.unlisted$sigma <- data.unlisted$rb * data.unlisted$r

#drop NA rows
data.unlisted <- data.unlisted[!is.na(data.unlisted$Km) & !is.na(data.unlisted$Rdm) & !is.na(data.unlisted$Rddm),]

# Select metrics here
mets <- c("Km","Rm","Rdm","Rddm","Kdm")
mets <- c("Km","Rm","Rdm","Kdm")
mets <- c("Km","Rm")
pp.train <- preProcess(data.unlisted[,mets], method = c("scale", "center", "pca", "BoxCox"), thresh = 1)


# Load in test data
load('Z:/Galen/Machine\ Learning\ Files/ml.test.results.norb.RData')
test.set.norb <- test.set
rm(test.set)
load('Z:/Galen/Machine\ Learning\ Files/ml.test.results.RData')
test.set <- as.data.frame(test.set)
test.set.norb <- as.data.frame(test.set.norb)
names(test.set) <- c("r","den","rb","gb","Km","Rm","Rdm","Rddm","Kdm")
names(test.set.norb) <- c("r","den","gb","Km","Rm","Rdm","Rddm","Kdm")
test.set <- test.set[!is.na(test.set$Km) & !is.na(test.set$Rdm) & !is.na(test.set$Rddm),]
test.set.norb <- test.set.norb[!is.na(test.set.norb$Km) & !is.na(test.set.norb$Rdm) & !is.na(test.set.norb$Rddm),]
test.set$rw <- (test.set$r^4 + 6* (test.set$r^2) *(test.set$rb * test.set$r)^2 + 3 *(test.set$rb * test.set$r)^4)/(test.set$r^3 + 3*test.set$r*(test.set$rb * test.set$r)^2)
test.set$sigma <- test.set$rb * test.set$r

# Make sure the pca is evaled up above in "full raw data"
test.pca <- predict(pp.train, test.set[,mets])
test.pca.norb <- predict(pp.train, test.set.norb[,mets])


##Load in a model
load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_rw.RData')
load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_4met_rw.RData')
load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_2met_rw.RData')

load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_den.RData')
load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_4met_den.RData')
load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_2met_den.RData')

load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_norb_r.RData')
load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_norb_4met_r.RData')
load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_norb_2met_r.RData')

load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_norb_den.RData')
load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_norb_4met_den.RData')
load('Z:/Galen/Machine\ Learning\ Files/Models/ml.models_norb_2met_den.RData')

predictors <- sapply(1:length(mets),function(x){paste('PC',toString(x),sep = '')})
out <- 'r'

test <- data.frame(test.set.norb, test.pca.norb)

##Predict test data with model
toRound <- 6
par(mar = c(4,4,3,3))
model.res <- lapply(models, function(x){
  a <- predict(x, newdata = test[,predictors])
  RMSE <- mean((a-test[,out])^2, na.rm = TRUE)
  plot(test[,out], a, main = x$method, xlab = "Input Simulation Value", ylab = "Estimated Value")
  abline(0,1, col = "red", lwd = 1.25)
  print(x$method)
  print(paste('RMSE: ', toString(round(RMSE,toRound)), sep = ''))
})

#Plot for paper
par(mgp = c(2, 1, 0), mar = c(3.5, 3.5, 2, 2))
i <- 3
a <- predict(models[[i]], newdata = test[,predictors])
plot(test[,out], a, main = models[[i]]$method, xlab = "Input Simulation Value", ylab = "Estimated Value")
abline(0,1, col = "red", lwd = 1.25)
legend(1.8, 10, "True Value = Predicted Value", col = "red", lty = 1, lwd = 2,bty = "n")
#


