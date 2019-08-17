# ML analysis for paper - is it better to have more random points, or more
# realizations at each point? Also, what is the difference between training with
# the mean and training with the raw values

library(caret)
library(data.table)
library(rapt)
library(parallel)
library(doParallel)
library(RColorBrewer)

#
# Generate Neccesary Data -------------------------------------------------
# 100 params with 100 each
# 500 params with 20 each
# 1000 params with 10 each
# 2500 params with 4 each (do you have this? - yes but lets re-do it in a more structured way.)
# 5000 params with 2 each
# 10000 params with 1 each

cluster.size <- c(60,60,60)

# Number of random simulations desired
nrand <- 10000 # Number of random parameters to generate
nparam <- c(1, 100, 500, 1000, 2500, 5000, 10000)
nclust <- c(100, 20, 10, 4, 2, 1) 

#nrand <- 10
#nparam <- c(0, 1, 3, 10)
#nclust <- c(3, 2, 1)

s <- 113 # random seed (change between runs)

#add in the random values
set.seed(s)
params <- list()
r <- runif(nrand, min = 2, max = 6.5)
den <- runif(nrand, min = 0.15, max = 1) 
rbp <- runif(nrand, min = 0, max = 0.5)
gbp <- runif(nrand, min = 0, max = 0.6)
for(i in 1:nrand){
  params[[i]] <- c(r[i], den[i], rbp[i], gbp[i])
}

#Set up the final data structure
sim.res <- list()
for(x in 1:(length(nparam)-1)){
  for(i in (nparam[x]+1):nparam[x+1]){
    sim.res[[i]] <- matrix(NA, nrow = nclust[x], ncol = 5)
  }
}

save(params, file = 'params.RData')
print("Parameter set initialized...")
print(paste("Total parameters to test: ", toString(nrand), sep = ''))

# Upload cube RRL files
toSub <- fread('~/Research/K_cluster_series/cubetoSub_r35.csv', drop=1)
env.r <- fread('~/Research/K_cluster_series/cube_r35.csv', select=2)

# HPC
#toSub <- fread('cubetoSub_r35.csv', drop=1)
#env.r <- fread('cube_r35.csv', select=2)

# Max r value and number of r values to go to in the k tests
maxr <- max(env.r)
nr <- nrow(env.r)

cores2use <- detectCores() 

print(paste("Using ", toString(cores2use), " cores.", sep = ''))
print("starting simulations...")

cl <- makePSOCKcluster(cores2use)
clusterEvalQ(cl, library(rapt))
clusterEvalQ(cl, library(zoo))

for(k in 1:length(nclust)){
  t1 <- Sys.time()
  out <- parLapply(cl, 1:nclust[k], kseries2, nclust[k], params[(nparam[k]+1):nparam[k+1]], maxr, nr, toSub, verbose = FALSE)
  t2 <- Sys.time()
  print(t2-t1)
  
  for(i in 1:nclust[k]){
    for(j in 1:nrow(out[[1]])){
      sim.res[[nparam[k]+j]][i,] <- out[[i]][j,] 
    }
  }
}

stopCluster(cl)
save(sim.res, file = "sim.res.RData")




# What matters more, number of random parameter combos, or number of repeats at each? -----------------------------------------
rm(list = ls())
gc()

load('Z:/Galen/Machine\ Learning\ Files/Model\ Fitting\ Sweep/params.RData')
load('Z:/Galen/Machine\ Learning\ Files/Model\ Fitting\ Sweep/sim.res.RData')
s.r <- sim.res
p <- params
rm(sim.res, params)

load('Z:/Galen/Machine\ Learning\ Files/190710_sim.res.RData')
load('Z:/Galen/Machine\ Learning\ Files/190710_params.RData')
s.r.t <- sim.res
p.t <- params
rm(sim.res, params)

#### Train with means, predict means and all ####
nparam <- c(100, 500, 1000, 2500, 5000, 10000)
nclust <- c(100, 20, 10, 4, 2, 1) 

temp.mean <- list()
for(i in 1:length(nparam)){
  if(nclust[i] != 1){
    a <- lapply(s.r[1:nparam[i]], function(x){
      b <- x[1:nclust[i],]
      return(apply(b, 2, mean, na.rm = TRUE))})
  } else {
    a <- lapply(s.r[1:nparam[i]], function(x){
      b <- x[1:nclust[i],]
      return(b)})
  }
  temp.mean[[i]] <- data.frame(matrix(unlist(a), nrow = nparam[i], ncol = 5, byrow = TRUE))
  names(temp.mean[[i]]) <- c('Km','Rm','Rdm','Rddm','Kdm')
}

train <- list(params = data.frame(matrix(unlist(p), nrow = 10000, ncol = 4, byrow = TRUE)), mets = temp.mean)
names(train$params) <- c('r','den','rb','gb')
train$params$rw <- (train$params$r^4 + 6* (train$params$r^2) *(train$params$rb * train$params$r)^2 + 3 *(train$params$rb * train$params$r)^4)/(train$params$r^3 + 3*train$params$r*(train$params$rb * train$params$r)^2)
train$params$sigma <- train$params$rb * train$params$r


temp.mean <- lapply(s.r.t, function(x){apply(x, 2, mean, na.rm = TRUE)})
s.r.t$mean <- data.frame(matrix(unlist(temp.mean), nrow = 2500, ncol = 5, byrow = TRUE))
names(s.r.t$mean) <- c('Km','Rm','Rdm','Rddm','Kdm')

test <- list(params = data.frame(matrix(unlist(p.t), nrow = 2500, ncol = 4, byrow = TRUE)), mets = s.r.t$mean)
names(test$params) <- c('r','den','rb','gb')
test$params$rw <- (test$params$r^4 + 6* (test$params$r^2) *(test$params$rb * test$params$r)^2 + 3 *(test$params$rb * test$params$r)^4)/(test$params$r^3 + 3*test$params$r*(test$params$rb * test$params$r)^2)
test$params$sigma <- test$params$rb * test$params$r

s.r.t$full <- matrix(NA, nrow = 2500*50, ncol = 5)
cnt <- 1
for(i in 1:2500){
  for(j in 1:50){
    s.r.t$full[cnt, ] <- s.r.t[[i]][j,]
    cnt <- cnt + 1
  }
}
s.r.t$full <- as.data.frame(s.r.t$full)
names(s.r.t$full) <- c('Km','Rm','Rdm','Rddm','Kdm')

params <- matrix(NA, nrow = 2500*50, ncol = 4)
cnt <- 1
for(i in 1:2500){
  for(j in 1:50){
    params[cnt,] <- p.t[[i]]
    cnt <- cnt + 1
  }
}
params <- as.data.frame(params)
names(params) <- c('r','den','rb','gb')

test.all <- list(params = params, mets = s.r.t$full)
test.all$params$rw <- (test.all$params$r^4 + 6* (test.all$params$r^2) *(test.all$params$rb * test.all$params$r)^2 + 3 *(test.all$params$rb * test.all$params$r)^4)/(test.all$params$r^3 + 3*test.all$params$r*(test.all$params$rb * test.all$params$r)^2)
test.all$params$sigma <- test.all$params$rb * test.all$params$r

mets <- c('Km','Rm','Rdm','Rddm','Kdm')

predictors <- c("PC1","PC2")

out <- c('rw','den')

res.mean <- list()
res.all <- list()

for(k in 1:length(out)){
  print(k)
  res.mat.mean <- matrix(NA, nrow = length(nclust), ncol = 4)
  res.mat.all <- matrix(NA, nrow = length(nclust), ncol = 4)
  for(i in 1:length(nclust)){
    
    print(paste('STARTING # ', toString(i), '/6',sep = ''))
    train$mets.cut <- train$mets[[i]][,mets]
    train$params.cut <- train$params[1:nparam[i],]
    
    train$params.cut <- train$params.cut[complete.cases(train$mets.cut),]
    train$mets.cut <- train$mets.cut[complete.cases(train$mets.cut),]
    
    pp <- preProcess(train$mets.cut, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
    pcas <- predict(pp, train$mets.cut)
    pcas.test <- predict(pp, test$mets[,mets])
    pcas.test.all <- predict(pp, test.all$mets[,mets])
    
    train$full <- data.frame(train$params.cut, train$mets.cut, pcas)
    test$full <- data.frame(test$params, test$mets[,mets], pcas.test)
    test.all$full <- data.frame(test.all$params, test.all$mets[,mets], pcas.test.all)
    
    # get rid of na rows
    test.all$full <- test.all$full[complete.cases(test.all$full),]
    
    models <- list()
    gc()
    
    # Based on all
    
    models[[1]] <- train(train$full[,predictors], train$full[,out[k]], method = "knn", tuneLength = 10)
    
    models[[2]] <- train(train$full[,predictors], train$full[,out[k]], 
                         method = "glmnet", 
                         trControl = trainControl("cv", number = 10),
                         tuneLength = 10)
    
    models[[3]] <- train(train$full[,predictors], train$full[,out[k]], method = "brnn", trControl = trainControl("cv", number = 10), verbose = FALSE)
    
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    models[[4]] <- train(train$full[,predictors],train$full[,out[k]],method = "parRF", trControl = trainControl("cv", number = 10))
    stopCluster(cl)
    registerDoSEQ()
    
    names(models) <- c("knn","elastic-net","brnn","rf")
    
    # Test on means
    compall <- extractPrediction(models, testX = test$full[,predictors], testY = test$full[,out[k]])
    comp <- split(compall, compall$model)
    
    res <- lapply(comp, function(x){
      xn <- x[x$dataType == 'Test',]
      ds <- defaultSummary(xn)
      
      #plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
      #abline(0,1, col = "red", lwd = 1.25)
      return(ds)
    })
    
    for(j in 1:4){
      res.mat.mean[i,j] <- res[[j]][1] 
    }
    
    # Test on all
    compall <- extractPrediction(models, testX = test.all$full[,predictors], testY = test.all$full[,out[k]])
    comp <- split(compall, compall$model)
    
    res <- lapply(comp, function(x){
      xn <- x[x$dataType == 'Test',]
      ds <- defaultSummary(xn)
      
      #plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
      #abline(0,1, col = "red", lwd = 1.25)
      return(ds)
    })
    
    for(j in 1:4){
      res.mat.all[i,j] <- res[[j]][1] 
    }
    
  }
  
  res.df.mean <- as.data.frame(res.mat.mean)
  res.df.all <- as.data.frame(res.mat.all)
  names(res.df.mean) <- c('brnn', 'glmnet', 'knn', 'rf')
  names(res.df.all) <- c('brnn', 'glmnet', 'knn', 'rf')
  
  res.mean[[k]] <- res.df.mean
  res.all[[k]] <- res.df.all
  
  rm(res.df.all, res.df.mean, res.mat.all, res.mat.mean)
}
names(res.mean) <- c('rw', 'den')
names(res.all) <- c('rw', 'den')

save(res.mean, file = 'MEAN.res.mean.RData')
save(res.all, file = 'MEAN.res.all.RData')

pchs <- c(15, 16, 17, 18)
cols <- c("red","black","blue","magenta")
par(mar = c(3.5, 3.5, 2, 2), mgp = c(2, 1, 0))

plot(res.mean[[1]]$brnn, type= 'n', ylim = c(0,0.45), 
     xlab = 'Training Data', ylab = "RMSE", 
     main = paste(out[1], ": Mean -> Mean", sep = ''))
ylines <- seq(0,0.5,by = 0.025)
for(i in 1:length(ylines)){
  abline(h = ylines[i], col = 'grey', lwd = 0.5)
}
for(i in 1:4){
  points(res.mean[[1]][,i], col = cols[i], pch = pchs[i], cex = 1.25)
}
legend(1, 0.15, legend = c('brnn','glm','knn','rf'), col = cols, pch = pchs, cex = 1.25)


plot(res.mean[[2]]$brnn, type= 'n', ylim = c(0,0.06), 
     xlab = 'Training Data', ylab = "RMSE", 
     main = paste(out[2], ": Mean -> Mean", sep = ''))
ylines <- seq(0,0.08,by = 0.0025)
for(i in 1:length(ylines)){
  abline(h = ylines[i], col = 'grey', lwd = 0.5)
}
for(i in 1:4){
  points(res.mean[[2]][,i], col = cols[i], pch = pchs[i], cex = 1.25)
}
legend(1, 0.02, legend = c('brnn','glm','knn','rf'), col = cols, pch = pchs, cex = 1.25)


plot(res.all[[1]]$brnn, type= 'n', ylim = c(0.3,max(res.all[[1]])), 
     xlab = 'Training Data', ylab = "RMSE", 
     main = paste(out[1], ": Mean -> All", sep = ''))
ylines <- seq(0,0.6,by = (0.025/2))
for(i in 1:length(ylines)){
  abline(h = ylines[i], col = 'grey', lwd = 0.5)
}
for(i in 1:4){
  points(res.all[[1]][,i], col = cols[i], pch = pchs[i], cex = 1.25)
}
legend(1, 0.4, legend = c('brnn','glm','knn','rf'), col = cols, pch = pchs, cex = 1.25)


plot(res.all[[2]]$brnn, type= 'n', ylim = c(0,max(res.all[[2]])), 
     xlab = 'Training Data', ylab = "RMSE", 
     main = paste(out[2], ": Mean -> All", sep = ''))
ylines <- seq(0,0.09,by = 0.005)
for(i in 1:length(ylines)){
  abline(h = ylines[i], col = 'grey', lwd = 0.5)
}
for(i in 1:4){
  points(res.all[[2]][,i], col = cols[i], pch = pchs[i], cex = 1.25)
}
legend(1, 0.04, legend = c('brnn','glm','knn','rf'), col = cols, pch = pchs, cex = 1.25)


#### Train with all, predict means and all ####
nparam <- c(100, 500, 1000, 2500, 5000, 10000)
nclust <- c(100, 20, 10, 4, 2, 1)

temp.all <- list()
temp.params <- list()
for(i in 1:length(nparam)){
  a <- lapply(s.r[1:nparam[i]], function(x){x[1:nclust[i],]})
  
  temp.all[[i]] <- matrix(NA, nrow = nparam[i]*nclust[i], ncol = 5)
  cnt <- 1
  if(nclust[i] != 1){
    for(j in 1:nparam[i]){
      for(k in 1:nclust[i]){
        temp.all[[i]][cnt,] <- a[[j]][k,]
        cnt <- cnt + 1
      }
    }} else {
      for(j in 1:nparam[i]){
        for(k in 1:nclust[i]){
          temp.all[[i]][cnt,] <- a[[j]]
          cnt <- cnt + 1
        }
      }
    }
  temp.all[[i]] <- as.data.frame(temp.all[[i]])
  names(temp.all[[i]]) <- c('Km','Rm','Rdm','Rddm','Kdm')
  
  temp.params[[i]] <- matrix(NA, nrow = nparam[i]*nclust[i], ncol = 4)
  cnt <- 1
  for(j in 1:nparam[i]){
    for(k in 1:nclust[i]){
      temp.params[[i]][cnt,] <- p[[j]]
      cnt <- cnt + 1
    }
  }
  temp.params[[i]] <- as.data.frame(temp.params[[i]])
  names(temp.params[[i]]) <- c('r','den','rb','gb')
  temp.params[[i]]$rw <- (temp.params[[i]]$r^4 + 6* (temp.params[[i]]$r^2) *(temp.params[[i]]$rb * temp.params[[i]]$r)^2 + 3 *(temp.params[[i]]$rb * temp.params[[i]]$r)^4)/(temp.params[[i]]$r^3 + 3*temp.params[[i]]$r*(temp.params[[i]]$rb * temp.params[[i]]$r)^2)
  temp.params[[i]]$sigma <- temp.params[[i]]$rb * temp.params[[i]]$r
}

train <- list(params = temp.params, mets = temp.all)


temp.mean <- lapply(s.r.t, function(x){apply(x, 2, mean, na.rm = TRUE)})
s.r.t$mean <- data.frame(matrix(unlist(temp.mean), nrow = 2500, ncol = 5, byrow = TRUE))
names(s.r.t$mean) <- c('Km','Rm','Rdm','Rddm','Kdm')

test <- list(params = data.frame(matrix(unlist(p.t), nrow = 2500, ncol = 4, byrow = TRUE)), mets = s.r.t$mean)
names(test$params) <- c('r','den','rb','gb')
test$params$rw <- (test$params$r^4 + 6* (test$params$r^2) *(test$params$rb * test$params$r)^2 + 3 *(test$params$rb * test$params$r)^4)/(test$params$r^3 + 3*test$params$r*(test$params$rb * test$params$r)^2)
test$params$sigma <- test$params$rb * test$params$r

s.r.t$full <- matrix(NA, nrow = 2500*50, ncol = 5)
cnt <- 1
for(i in 1:2500){
  for(j in 1:50){
    s.r.t$full[cnt, ] <- s.r.t[[i]][j,]
    cnt <- cnt + 1
  }
}
s.r.t$full <- as.data.frame(s.r.t$full)
names(s.r.t$full) <- c('Km','Rm','Rdm','Rddm','Kdm')

params <- matrix(NA, nrow = 2500*50, ncol = 4)
cnt <- 1
for(i in 1:2500){
  for(j in 1:50){
    params[cnt,] <- p.t[[i]]
    cnt <- cnt + 1
  }
}
params <- as.data.frame(params)
names(params) <- c('r','den','rb','gb')

test.all <- list(params = params, mets = s.r.t$full)
test.all$params$rw <- (test.all$params$r^4 + 6* (test.all$params$r^2) *(test.all$params$rb * test.all$params$r)^2 + 3 *(test.all$params$rb * test.all$params$r)^4)/(test.all$params$r^3 + 3*test.all$params$r*(test.all$params$rb * test.all$params$r)^2)
test.all$params$sigma <- test.all$params$rb * test.all$params$r

mets <- c('Km','Rm','Rdm','Rddm','Kdm')

predictors <- c("PC1","PC2")

out <- c('rw','den')

res.mean <- list()
res.all <- list()

for(k in 1:length(out)){
  print(k)
  res.mat.mean <- matrix(NA, nrow = length(nclust), ncol = 4)
  res.mat.all <- matrix(NA, nrow = length(nclust), ncol = 4)
  for(i in 1:length(nclust)){
    
    print(paste('STARTING # ', toString(i), '/6',sep = ''))
    train$mets.cut <- train$mets[[i]][,mets]
    train$params.cut <- train$params[[i]]
    
    train$params.cut <- train$params.cut[complete.cases(train$mets.cut),]
    train$mets.cut <- train$mets.cut[complete.cases(train$mets.cut),]
    
    pp <- preProcess(train$mets.cut, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
    pcas <- predict(pp, train$mets.cut)
    pcas.test <- predict(pp, test$mets[,mets])
    pcas.test.all <- predict(pp, test.all$mets[,mets])
    
    train$full <- data.frame(train$params.cut, train$mets.cut, pcas)
    test$full <- data.frame(test$params, test$mets[,mets], pcas.test)
    test.all$full <- data.frame(test.all$params, test.all$mets[,mets], pcas.test.all)
    
    # get rid of na rows
    test.all$full <- test.all$full[complete.cases(test.all$full),]
    
    models <- list()
    gc()
    
    # Based on all
    
    models[[1]] <- train(train$full[,predictors], train$full[,out[k]], method = "knn", tuneLength = 10)
    
    models[[2]] <- train(train$full[,predictors], train$full[,out[k]], 
                         method = "glmnet", 
                         trControl = trainControl("cv", number = 10),
                         tuneLength = 10)
    
    models[[3]] <- train(train$full[,predictors], train$full[,out[k]], method = "brnn", trControl = trainControl("cv", number = 10), 
                         tuneGrid = data.frame('neurons' = c(3, 5, 7, 9)))
    
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    models[[4]] <- train(train$full[,predictors],train$full[,out[k]],method = "parRF", trControl = trainControl("cv", number = 10))
    stopCluster(cl)
    registerDoSEQ()
    
    names(models) <- c("knn","elastic-net","brnn","rf")
    
    # Test on means
    compall <- extractPrediction(models, testX = test$full[,predictors], testY = test$full[,out[k]])
    comp <- split(compall, compall$model)
    
    res <- lapply(comp, function(x){
      xn <- x[x$dataType == 'Test',]
      ds <- defaultSummary(xn)
      
      #plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
      #abline(0,1, col = "red", lwd = 1.25)
      return(ds)
    })
    
    for(j in 1:4){
      res.mat.mean[i,j] <- res[[j]][1] 
    }
    
    # Test on all
    compall <- extractPrediction(models, testX = test.all$full[,predictors], testY = test.all$full[,out[k]])
    comp <- split(compall, compall$model)
    
    res <- lapply(comp, function(x){
      xn <- x[x$dataType == 'Test',]
      ds <- defaultSummary(xn)
      
      plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
      abline(0,1, col = "red", lwd = 1.25)
      return(ds)
    })
    
    for(j in 1:4){
      res.mat.all[i,j] <- res[[j]][1] 
    }
    
  }
  
  res.df.mean <- as.data.frame(res.mat.mean)
  res.df.all <- as.data.frame(res.mat.all)
  names(res.df.mean) <- c('brnn', 'glmnet', 'knn', 'rf')
  names(res.df.all) <- c('brnn', 'glmnet', 'knn', 'rf')
  
  res.mean[[k]] <- res.df.mean
  res.all[[k]] <- res.df.all
  
  rm(res.df.mean, res.df.all, res.mat.all, res.mat.mean)
}
names(res.mean) <- c('rw', 'den')
names(res.all) <- c('rw', 'den')
# Run the above on the HPC

load('~/Research/ML/ALL.res.all.RData')
load('~/Research/ML/ALL.res.mean.RData')
out <- c('rw', 'den')


pchs <- c(15, 16, 17, 18)
cols <- c("red","black","blue","magenta")
par(mar = c(3.5, 3.5, 2, 2), mgp = c(2, 1, 0))


plot(res.mean[[1]]$brnn, type= 'n', ylim = c(0,max(res.mean[[1]])), 
     xlab = 'Training Data', ylab = "RMSE", 
     main = paste(out[1], ": All -> Mean", sep = ''))
ylines <- seq(0,0.5,by = 0.025)
for(i in 1:length(ylines)){
  abline(h = ylines[i], col = 'grey', lwd = 0.5)
}
for(i in 1:4){
  points(res.mean[[1]][,i], col = cols[i], pch = pchs[i], cex = 1.25)
}
legend(1, 0.15, legend = c('brnn','glm','knn','rf'), col = cols, pch = pchs, cex = 1.25)


plot(res.mean[[2]]$brnn, type= 'n', ylim = c(0,max(res.mean[[2]])), 
     xlab = 'Training Data', ylab = "RMSE", 
     main = paste(out[2], ": All -> Mean", sep = ''))
ylines <- seq(0,0.08,by = 0.0025)
for(i in 1:length(ylines)){
  abline(h = ylines[i], col = 'grey', lwd = 0.5)
}
for(i in 1:4){
  points(res.mean[[2]][,i], col = cols[i], pch = pchs[i], cex = 1.25)
}
legend(1, 0.03, legend = c('brnn','glm','knn','rf'), col = cols, pch = pchs, cex = 1.25)


plot(res.all[[1]]$brnn, type= 'n', ylim = c(0.3,max(res.all[[1]])), 
     xlab = 'Training Data', ylab = "RMSE", 
     main = paste(out[1], ": All -> All", sep = ''))
ylines <- seq(0,0.6,by = (0.025/2))
for(i in 1:length(ylines)){
  abline(h = ylines[i], col = 'grey', lwd = 0.5)
}
for(i in 1:4){
  points(res.all[[1]][,i], col = cols[i], pch = pchs[i], cex = 1.25)
}
legend(1, 0.4, legend = c('brnn','glm','knn','rf'), col = cols, pch = pchs, cex = 1.25)


plot(res.all[[2]]$brnn, type= 'n', ylim = c(0,max(res.all[[2]])), 
     xlab = 'Training Data', ylab = "RMSE", 
     main = paste(out[2], ": All -> All", sep = ''))
ylines <- seq(0,0.09,by = 0.005)
for(i in 1:length(ylines)){
  abline(h = ylines[i], col = 'grey', lwd = 0.5)
}
for(i in 1:4){
  points(res.all[[2]][,i], col = cols[i], pch = pchs[i], cex = 1.25)
}
legend(1, 0.04, legend = c('brnn','glm','knn','rf'), col = cols, pch = pchs, cex = 1.25)


#

# Sweep through how much multiple realizations at a value impact performance --------
load('Z:/Galen/Machine\ Learning\ Files/2500x50\ sims/190710_sim.res.RData')
load('Z:/Galen/Machine\ Learning\ Files/2500x50\ sims/190710_params.RData')

test <- list(params = data.frame(matrix(unlist(params[1:500]), nrow = 500, byrow = TRUE)), mets = sim.res[1:500])
names(test$params) <- c('r','den','rb','gb')
test$params$rw <- (test$params$r^4 + 6* (test$params$r^2) *(test$params$rb * test$params$r)^2 + 3 *(test$params$rb * test$params$r)^4)/(test$params$r^3 + 3*test$params$r*(test$params$rb * test$params$r)^2)
test$params$sigma <- test$params$rb * test$params$r

train <- list(params = data.frame(matrix(unlist(params[501:2500]), nrow = 2000, byrow = TRUE)), mets = sim.res[501:2500])
names(train$params) <- c('r','den','rb','gb')
train$params$rw <- (train$params$r^4 + 6* (train$params$r^2) *(train$params$rb * train$params$r)^2 + 3 *(train$params$rb * train$params$r)^4)/(train$params$r^3 + 3*train$params$r*(train$params$rb * train$params$r)^2)
train$params$sigma <- train$params$rb * train$params$r

rm(params, sim.res)
mets <- c('Km','Rm','Rdm','Rddm','Kdm')

nclust <- c(50,25,10,5,2)

res.mat <- matrix(NA, nrow = length(nclust), ncol = 2)

#### Train with means, predict means ####
for(i in 1:length(nclust)){
  print(paste('STARTING # ', toString(i), '/5',sep = ''))
  train$mets.cut <- list()
  test$mets.cut <- list()
  
  for(j in 1:length(train$mets)){
    train$mets.cut[[j]] <- train$mets[[j]][1:nclust[i], 1:5]
  }
  for(j in 1:length(test$mets)){
    test$mets.cut[[j]] <- test$mets[[j]][1:nclust[i], 1:5]
  }
  
  train$mets.mean <- lapply(train$mets.cut, apply, 2, mean, na.rm = TRUE)
  test$mets.mean <- lapply(test$mets.cut, apply, 2, mean, na.rm = TRUE)
  
  temp <- data.frame(matrix(unlist(train$mets.mean), ncol = 5, byrow = TRUE))
  temp.test <- data.frame(matrix(unlist(test$mets.mean), ncol = 5, byrow = TRUE))
  
  names(temp) <- mets
  names(temp.test) <- mets
  
  pp <- preProcess(temp, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
  pcas <- predict(pp, temp)
  pcas.test <- predict(pp, temp.test)
  
  train$full <- data.frame(train$params, temp, pcas)
  test$full <- data.frame(test$params, temp.test, pcas.test)
  
  predictors <- c("PC1","PC2")
  out <- 'rw'
  
  models <- list()
  gc()
  
  # Based on all
  models[[1]] <- train(train$full[,predictors], train$full[,out], method = "glm", trControl = trainControl("cv", number = 10))
  
  models[[2]] <- train(train$full[,predictors], train$full[,out], method = "knn", tuneLength = 10)
  
  models[[3]] <- train(train$full[,predictors], train$full[,out], 
                       method = "glmnet", 
                       trControl = trainControl("cv", number = 10),
                       tuneLength = 10)
  
  models[[4]] <- train(train$full[,predictors], train$full[,out], method = "brnn", trControl = trainControl("cv", number = 10))
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  models[[5]] <- train(train$full[,predictors],train$full[,out],method = "parRF", trControl = trainControl("cv", number = 10))
  stopCluster(cl)
  registerDoSEQ()
  
  names(models) <- c("lm","knn","elastic-net","brnn","rf")
  
  compall <- extractPrediction(models, testX = test$full[,predictors], testY = test$full[,out])
  comp <- split(compall, compall$model)
  
  res <- lapply(comp, function(x){
    xn <- x[x$dataType == 'Test',]
    ds <- defaultSummary(xn)
    
    #plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
    #abline(0,1, col = "red", lwd = 1.25)
    return(ds)
  })
  
  for(j in 1:5){
    res.mat[i,j] <- res[[j]][1] 
  } 
}

save(res.mat, file = 'START_HERE.RData')

#### Train with means, predict all data ####
#### Train with all data, predict means ####
#### Train with all data, predict all data ####
for(i in 1:length(nclust)){
  print(paste('STARTING # ', toString(i), '/5',sep = ''))
  train$mets.cut <- list()
  test$mets.cut <- list()
  
  for(j in 1:length(train$mets)){
    train$mets.cut[[j]] <- train$mets[[j]][1:nclust[i], 1:5]
  }
  for(j in 1:length(test$mets)){
    test$mets.cut[[j]] <- test$mets[[j]][1:nclust[i], 1:5]
  }
  
  temp <- as.data.frame(do.call(rbind, train$mets.cut))
  temp.test <- as.data.frame(do.call(rbind, test$mets.cut))
  
  names(temp) <- mets
  names(temp.test) <- mets
  
  pp <- preProcess(temp, method = c("scale", "center", "pca", "BoxCox"), thresh = 1)
  pcas <- predict(pp, temp)
  pcas.test <- predict(pp, temp.test)
  
  temp.params <- train$params[rep(seq(nrow(train$params)), each = nclust[i]),]
  temp.params.test <- test$params[rep(seq(nrow(test$params)), each = nclust[i]),]
  
  train$full <- data.frame(temp.params, temp, pcas)
  test$full <- data.frame(temp.params.test, temp.test, pcas.test)
  
  train$full <- train$full[complete.cases(train$full),]
  test$full <- test$full[complete.cases(test$full),]
  
  predictors <- c("PC1","PC2")
  out <- 'rw'
  
  models <- list()
  gc()
  
  # Based on all
  models[[1]] <- train(train$full[,predictors], train$full[,out], 
                       method = "glmnet", 
                       trControl = trainControl("cv", number = 10),
                       tuneLength = 10)
  
  models[[2]] <- train(train$full[,predictors], train$full[,out], method = "brnn", trControl = trainControl("cv", number = 10),
                       tuneGrid = data.frame('neurons' = c(5)))

  
  names(models) <- c('glmnet',"brnn")
  
  compall <- extractPrediction(models, testX = test$full[,predictors], testY = test$full[,out])
  comp <- split(compall, compall$model)
  
  res <- lapply(comp, function(x){
    xn <- x[x$dataType == 'Test',]
    ds <- defaultSummary(xn)
    
    #plot(xn$obs, xn$pred, main = toString(xn$model[1]), xlab = "Actual Value", ylab = "Predicted Value")
    #abline(0,1, col = "red", lwd = 1.25)
    return(ds)
  })
  
  for(j in 1:2){
    res.mat[i,j] <- res[[j]][1] 
  } 
}





# Actual paper analysis and plots -----------------------------------------
#### NO RB ####
# Load training data and calculate PCAs: Select either rb or norb files to load
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/190725_params_norb.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/190725_sim.res_norb.RData')

# unpack data
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
names(data.unlisted) <- c("r","den","gb","Km","Rm","Rdm","Rddm","Kdm")
data.unlisted$rw <- data.unlisted$r

#drop NA rows
data.unlisted <- data.unlisted[complete.cases(data.unlisted),]

# Select metrics here
mets <- c("Km","Rm","Rdm","Rddm","Kdm")
mets <- c("Km","Rm","Rdm","Kdm")
mets <- c("Km","Rm")
pp.train <- preProcess(data.unlisted[,mets], method = c("scale", "center", "pca", "BoxCox"), thresh = 1)

# Upload test data
load('Z:/Galen/Machine\ Learning\ Files/Test\ Data/ml.test.results.norb.RData')
test.set <- as.data.frame(test.set)
names(test.set) <- c("r","den","gb","Km","Rm","Rdm","Rddm","Kdm")
test.set <- test.set[complete.cases(test.set),]
test.set$rw <- test.set$r

test.pca <- predict(pp.train, test.set[,mets])
test <- data.frame(test.set, test.pca)

# Upload model of choice
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_norb_r.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_norb_4met_r.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_norb_2met_r.RData')

load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_norb_den.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_norb_4met_den.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_norb_2met_den.RData')

predictors <- sapply(1:length(mets),function(x){paste('PC',toString(x),sep = '')})
out <- 'den'

# Test models
toRound <- 6
model.res <- lapply(models, function(x){
  a <- predict(x, newdata = test[,predictors])
  RMSE <- sqrt(mean((a-test[,out])^2))
  plot(test[,out], a, main = x$method, xlab = "Input Simulation Value", ylab = "Estimated Value")
  abline(0,1, col = "red", lwd = 1.25)
  print(x$method)
  print(paste('RMSE: ', toString(round(RMSE,toRound)), sep = ''))
})

# upload the next model
rm(models)
gc()


#### RB ####
# Load training data and calculate PCAs: Select either rb or norb files to load
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/190725_params.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/190725_sim.res.RData')

# unpack data
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


#drop NA rows
data.unlisted <- data.unlisted[complete.cases(data.unlisted),]

# Select metrics here
mets <- c("Km","Rm","Rdm","Rddm","Kdm")
mets <- c("Km","Rm","Rdm","Kdm")
mets <- c("Km","Rm")
pp.train <- preProcess(data.unlisted[,mets], method = c("scale", "center", "pca", "BoxCox"), thresh = 1)

# Upload test data
load('Z:/Galen/Machine\ Learning\ Files/Test\ Data/ml.test.results.RData')
test.set <- as.data.frame(test.set)
names(test.set) <- c("r","den","rb","gb","Km","Rm","Rdm","Rddm","Kdm")
test.set <- test.set[complete.cases(test.set),]
test.set$rw <- (test.set$r^4 + 6* (test.set$r^2) *(test.set$rb * test.set$r)^2 + 3 *(test.set$rb * test.set$r)^4)/(test.set$r^3 + 3*test.set$r*(test.set$rb * test.set$r)^2)
test.set$sigma <- test.set$rb * test.set$r

test.pca <- predict(pp.train, test.set[,mets])
test <- data.frame(test.set, test.pca)

# Upload model of choice
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_rw.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_4met_rw.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_2met_rw.RData')

load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_den.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_4met_den.RData')
load('Z:/Galen/Machine\ Learning\ Files/10000x10\ sims\ and\ models/ml.models_2met_den.RData')


predictors <- sapply(1:length(mets),function(x){paste('PC',toString(x),sep = '')})
out <- 'rw'

# Test models
toRound <- 6
model.res <- lapply(models, function(x){
  a <- predict(x, newdata = test[,predictors])
  RMSE <- sqrt(mean((a-test[,out])^2))
  plot(test[,out], a, main = x$method, xlab = "Input Simulation Value", ylab = "Estimated Value")
  abline(0,1, col = "red", lwd = 1.25)
  print(x$method)
  print(paste('RMSE: ', toString(round(RMSE,toRound)), sep = ''))
})

# upload the next model
rm(models)
gc()

rm(pp.train, test)

# Residual analysis and plots ---------------------------------------------
pred.full <- extractPrediction(models[3], testX = test[,predictors], testY = test[,out])
pred <- pred.full[pred.full$dataType == 'Test',]
defaultSummary(pred)
diff <- pred$obs - pred$pred
diff.perc <- diff*100/pred$obs

hist(diff, breaks = 100, xlim = c(-2,2))
hist(diff.perc, breaks = 100, xlab = 'Percent Error', ylab = 'Frequency')
qqnorm(diff)
qqnorm(diff.perc)
qqline(diff)

# Percent error error plot
percs <- c(50, 75, 90)
cols <- brewer.pal(3, "Set1")
downsample.ratio <- 0.5
par(mar = c(4, 4, 2, 2), mgp = c(2.3, 1, 0))

n <- length(diff.perc)
diff.perc.sorted <- sort(abs(diff.perc))
xs <- (1:n)*100/n

plot(xs, diff.perc.sorted, 
     xlab = 'Sorted Data Percentile', ylab = 'Percent Error',
     type = 'n')
vals <- list('all' = data.frame('perc' = xs, 'error' = diff.perc.sorted))
ind <- rep(0, length(percs)+1)
for(i in 1:length(percs)){
  ind[i+1] <- round(n*percs[i]/100)
  vals[[i+1]] <- vals$all[(ind[i]+1):ind[i+1],]
  n.i <- nrow(vals[[i+1]])
  sample.inds <- sample(1:n.i, round(n.i*downsample.ratio))
  points(vals[[i+1]]$perc[sample.inds], vals[[i+1]]$error[sample.inds], pch = 16, col = cols[i], cex = 1)
  segments(-10, tail(vals[[i+1]]$error, n = 1), tail(vals[[i+1]]$perc, n = 1), tail(vals[[i+1]]$error, n = 1),
           col = 'red', lwd = 1.5)
  #text(0, tail(vals[[i+1]]$error, n = 1) + 1.25, paste(toString(percs[i]), '% of data'), pos = 4, cex = 0.75, col = 'red')
}
q <- (tail(ind, n = 1 )+1):n
ind.last <- sample(q, length(q)*downsample.ratio)
points(vals$all$perc[ind.last], vals$all$error[ind.last], col = 'black', pch = 16, cex = 1)

# 1:1 points envelope plot
sorted <- data.frame('obs' = pred$obs[order(abs(diff.perc))], 'pred' = pred$pred[order(abs(diff.perc))])
n <- nrow(sorted)

plot(sorted$obs, sorted$pred, col = 'black', pch = 16, type = 'n',
     xlab = 'Input Simiulation Value', ylab = 'Estimated Value')

vals.subed <- list('all' = sorted)
for(i in 1:length(percs)){
  vals.subed[[i+1]] <- sorted[(ind[i]+1):ind[i+1],]
  n.i <- nrow(vals.subed[[i+1]])
  sample.inds <- sample(1:n.i, round(n.i*downsample.ratio))
  points(vals.subed[[i+1]]$obs[sample.inds], vals.subed[[i+1]]$pred[sample.inds], pch = 16, col = cols[i], cex = 0.35)
}
q <- (tail(ind, n =1 )+1):n
ind.last <- sample(q, length(q)*downsample.ratio)
points(vals.subed$all$obs[ind.last], vals.subed$all$pred[ind.last], col = 'black', pch = 16, cex = 0.35)

leg <- rep('', length(percs)+1)
for(i in 1:(length(leg)-1)){
  leg[i] <- paste(toString(percs[i]), '%', sep = '')
}
leg[length(leg)] <- '100%'

legend(.2, 1, legend = leg, col = c(cols, 'black'), pch = 16, cex = 1, title = 'Sorted Data Percentile:', bty = 'n')
abline(0, 1, col = 'black', lty = 2, lwd = 1.5)



