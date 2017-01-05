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
# remove exam3mcz from this list as we want this variable in the model
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
# the 3 eczema variables
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

### univariate report
i = which(names(ivScore) %in% c('frozen', 'exam3mecz'))
dfData = data.frame(Shannon=dfImport$shannonmean, Frozen=dfImport$frozen, Eczema3m=dfImport$exam3mecz, dfImport[,names(ivScore)[-i]])
str(dfData)
cn = names(ivScore[-i])

lFm = lapply(cn, function(x) {
  df = dfData[,c('Shannon', 'Frozen', 'Eczema3m', x)]
  fm = lm(Shannon ~ Eczema3m*Frozen + . , data=df)
})

names(lFm) = cn
library(car)
f_get.coef.pvalue = function(x){
  s = Anova(lFm[[x]])
  round(s$`Pr(>F)`[3], 3)
}

m = sapply(cn, function(x) {
  tryCatch(expr = f_get.coef.pvalue(x), error=function(e) NULL)
})

dfUnivariate = data.frame(m)
rownames(dfUnivariate) = names(m)
Sig = ifelse(m < 0.05, 'Yes', 'No')
dfUnivariate$Significant = Sig
colnames(dfUnivariate) = c('P-Value', 'Significant')

dir.create('Temp')
write.csv(dfUnivariate, file='Temp/uni.csv')

lapply(lFm[Sig == 'Yes'], summary)

dfData = data.frame(Shannon=dfImport$shannonmean, dfImport[,names(ivScore)])
fm01 = lm(Shannon ~ exam3mecz*frozen + ., data=dfData)
summary(fm01)

Anova(fm01)
anova(fm01, test='Chisq')

## keep only the interesting covariates based on univariate and multivariate results
fm01.a = lm(Shannon ~ exam3mecz*frozen + DaysOldAtSampling+ parecz3mg+MouldAtHome+nonhigherpatsch_q3mgen+
              nonhighermatsch_q3mgen+numrespinf3m+abxmnth3m, data=dfData)

summary(fm01.a)
Anova(fm01.a)

## subset model further 
fm01.b = lm(Shannon ~ exam3mecz*frozen + DaysOldAtSampling+ MouldAtHome+numrespinf3m+abxmnth3m, data=dfData)

summary(fm01.b)
Anova(fm01.b)

# convert num respiratory infections variable to numeric
x = dfData$numrespinf3m
class(x)
str(x)
x2 = rep(NA, length=length(x))
levels(x)
x2[x == '4'] = 4
x2[x == '5'] = 5
x2[x == "At least three episodes"] = 3
x2[x == "One episode"] = 1
x2[x == "Two episodes"] = 2
x2[x == "None"] = 0

dfData$numrespinf3m = x2
# refit after converting it to numeric
fm01.b = lm(Shannon ~ exam3mecz*frozen + DaysOldAtSampling+ MouldAtHome+numrespinf3m+abxmnth3m, data=dfData)

summary(fm01.b)
Anova(fm01.b)

## subset the model again
fm01.c = lm(Shannon ~ exam3mecz*frozen + MouldAtHome+abxmnth3m, data=dfData)
summary(fm01.c)
Anova(fm01.c)

library(lattice)

xyplot(ifelse(dfData$exam3mecz == "Yes", 1, 0) ~ Shannon | frozen, data = dfData,
       type = c("p", "smooth"),
       ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity given Frozen')

xyplot(ifelse(dfData$exam3mecz == "Yes", 1, 0) ~ Shannon, data = dfData,
       type = c("p"),
       ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity')

xyplot(ifelse(dfData$exam3mecz == "Yes", 1, 0) ~ Shannon | frozen, data = dfData,
       type = c("p"),
       ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity given Frozen')

xyplot(Shannon ~ abxmnth3m, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=T,
       type=c('p', 'g', 'r'), groups=exam3mecz)

xyplot(Shannon ~ abxmnth3m | ifelse(MouldAtHome == 'Yes', 'Mould Yes', 'Mould No'), data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=T,
       type=c('p', 'g', 'r'), groups=ifelse(exam3mecz == 'Yes', 'Eczema Yes', 'Eczema No'),
       main='Diversity vs Antibiotic treatment Given Presense of Mould at home and eczema', xlab='Antibiotics 3 months')

xyplot(ifelse(dfData$MouldAtHome == "Yes", 1, 0) ~ Shannon , data = dfData,
       type = c("p"),
       ylab = "Mould at home", xlab = "Shannon Diversity", main='Mould at vs Diversity')

densityplot(~ Shannon , data=dfData, auto.key=T, groups=ifelse(MouldAtHome == 'Yes', 'Mould Yes', 'Mould No'),
            main='Shannon Diversity 3m given Mould at home')


xyplot(shannonmean ~ age | id, index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',
       data=dfData, type=c('g', 'p', 'r'), par.strip.text=list(cex=0.4))#, layout=c(12,5))