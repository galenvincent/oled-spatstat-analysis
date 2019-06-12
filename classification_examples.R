# Practice with multinomial regression
library(foreign)
library(nnet)
library(ggplot2)
library(reshape2)

ml <- read.dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")
head(ml)

table(ml$ses, ml$prog)
do.call(rbind, tapply(ml$write, ml$prog, function(x){c(M = mean(x), SD = sd(x))}))

ml$prog2 <- relevel(ml$prog, ref = "academic")
test <- multinom(prog2 ~ ses + write, data = ml)
summary(test)


# example 2
rm(list = ls())
gc()
cmcData <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/cmc/cmc.data", stringsAsFactors=FALSE, header=F)
colnames(cmcData) <- c("wife_age", "wife_edu", "hus_edu", "num_child", "wife_rel", "wife_work", "hus_occu", "sil", "media_exp", "cmc")
head(cmcData)

cmcData$wife_edu <- factor(cmcData$wife_edu, levels=sort(unique(cmcData$wife_edu)))
cmcData$hus_edu <- factor(cmcData$hus_edu, levels=sort(unique(cmcData$hus_edu)))
cmcData$wife_rel <- factor(cmcData$wife_rel, levels=sort(unique(cmcData$wife_rel)))
cmcData$wife_work <- factor(cmcData$wife_work, levels=sort(unique(cmcData$wife_work)))
cmcData$hus_occu <- factor(cmcData$hus_occu, levels=sort(unique(cmcData$hus_occu)))
cmcData$sil <- factor(cmcData$sil, levels=sort(unique(cmcData$sil)))
cmcData$media_exp <- factor(cmcData$media_exp, levels=sort(unique(cmcData$media_exp)))
cmcData$cmc <- factor(cmcData$cmc, levels=sort(unique(cmcData$cmc)))


trainingRows <- sample(1:nrow(cmcData), 0.7*nrow(cmcData))
train <- cmcData[trainingRows,]
test <- cmcData[-trainingRows,]

model <- multinom(cmc ~ ., data = train)
summary(model)

pred <- predict(model, test, "probs")
pred_class <- predict(model, test)
head(pred)
head(pred_class)

table(pred_class, test$cmc)
