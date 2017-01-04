# Name: Q1_diversity_3m.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 04/01/2017
# Desc: Shannon diversity measured at 3 months is associated with eczema status at 3 months.


## data import
dfImport = read.csv('Data_external/merged_dataset_UN_Manual.csv', header=T)
dfImport = dfImport[dfImport$age == '3m',]
dfImport = dfImport[,-c(1,2, 3)] # drop the subject id, sample id, age
dim(dfImport)

# which columns have the highest NAs
i = sapply(1:ncol(dfImport), function(x) sum(is.na(dfImport[,x])))
# remove columns with higher than ~ 35 NAs
f = i > 35
table(f)
# remove 35 columns
dfImport = dfImport[,!f]

#### perform a quick random forest calculation to look at important covariates
library(randomForest)
dim(na.omit(dfImport))
set.seed(123)
fit.rf = randomForest(shannonmean ~ ., data=na.omit(dfImport))
# get variables importance
varImpPlot(fit.rf)
dfRF = data.frame(importance(fit.rf))
head(dfRF)
ivScore = dfRF$IncNodePurity
names(ivScore) = rownames(dfRF)
ivScore = sort(ivScore, decreasing = T)
# which covariates have the abundances
i = grep('__', names(ivScore))
# remove these
ivScore = ivScore[-i]
# remove selected covariates which are not of interest
data.frame(names(ivScore))
i = c(1:13, 17, 19,24,28,29,30,35,41,44,47,55,56,57,60,76,79,80,81,82,83)
data.frame(names(ivScore)[i])
ivScore = ivScore[-i]

### select only these variables
dfData = data.frame(dfImport[,names(ivScore)])

## find correlated variables
m = NULL;

for (i in 1:ncol(dfData)){
  m = cbind(m, dfData[,i])
}
colnames(m) = colnames(dfData)
mCor = cor(m, use="na.or.complete")
library(caret)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.7, names=T)
data.frame(n)
n = n[-9]

i = which(names(ivScore) %in% n)
ivScore = ivScore[-i]
data.frame(ivScore)
### perform another round of random forest
dfData = data.frame(Shannon=dfImport$shannonmean, dfImport[,names(ivScore)])
set.seed(123)
fit.rf = randomForest(Shannon ~ ., data=na.omit(dfData))
# get variables importance
varImpPlot(fit.rf)
dfRF = data.frame(importance(fit.rf))
head(dfRF)
ivScore = dfRF$IncNodePurity
names(ivScore) = rownames(dfRF)
ivScore = sort(ivScore, decreasing = T)
data.frame(names(ivScore))
i = c(20,31,38)

## check correlations of these 3 different measures of eczema status
dfData = data.frame(dfImport[,names(ivScore)[i]])

## find correlated variables
m = NULL;

for (i in 1:ncol(dfData)){
  m = cbind(m, dfData[,i])
}
colnames(m) = colnames(dfData)
cor(m, use="na.or.complete")

## keep the eczema status only
data.frame(names(ivScore))
i = c(20,31)
data.frame(names(ivScore)[i])
ivScore = ivScore[-i]

### check variable combinations 
## download the ccrossvalidation class to select variables of importance
# this section of code causes errors due to liner dependencies - skip it
# if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')
# 
# url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
# download(url, 'CCrossValidation.R')
# 
# # load the required packages
# source('CCrossValidation.R')
# # delete the file after source
# unlink('CCrossValidation.R')
# 
# dfData = data.frame(Shannon=dfImport$shannonmean, dfImport[,names(ivScore)])
# set.seed(123)
# dfData.sub = na.omit(dfData)
# dim(dfData.sub)
# oVar.s = CVariableSelection.ReduceModel(dfData.sub[,-1], dfData.sub$Shannon , boot.num = 50, cvMethod = 'forward')
# plot.var.selection(oVar.s)
# dfData.sub = dfData.sub[,-1]
# sapply(seq_along(1:ncol(dfData.sub)), function(x) CVariableSelection.ReduceModel.getMinModel(oVar.s, x, cvMethod = 'forward'))

dfData = data.frame(Shannon=dfImport$shannonmean, dfImport[,names(ivScore)])
fm01 = lm(Shannon ~ ., data=dfData)
summary(fm01)
library(car)
Anova(fm01)
anova(fm01, test='Chisq')

## keep only the interesting covariates
fm01.a = lm(Shannon ~ exam3mecz*frozen + DaysOldAtSampling+ parecz3mg+Siblings+MouldAtHome+nonhigherpatsch_q3mgen+
              anyabxupto3m+nonhighermatsch_q3mgen, data=dfData)

summary(fm01.a)
Anova(fm01.a)

