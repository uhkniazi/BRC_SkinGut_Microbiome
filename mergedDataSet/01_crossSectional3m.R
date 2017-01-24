# Name: 01_crossSectional3m.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 23/01/2017
# Desc: relationships of diversity with covariates


p.old = par()

## data import
dfImport = read.csv('mergedDataSet/Data_external/3month_data_jan_2017.csv', header=T)
rownames(dfImport) = dfImport$id
i = sapply(c('id', 'mean', 'anyteeth1isyes'), grep, colnames(dfImport))
colnames(dfImport)[unlist(i)]
# remove these columns apart from shannonmean
# [1] "id"              "shannonmean"     "chaomean"        "pdwholetreemean" "simpsonmean"    
# [6] "obsotumean"      "anyteeth1isyes" 
i = unlist(i)[-2]
dfImport = dfImport[,-i]
dfImport = dfImport[,-2]
dfImport = dfImport[,-1]
dim(dfImport)
str(dfImport)
summary(dfImport)
# which columns have the highest NAs
i = sapply(1:ncol(dfImport), function(x) sum(is.na(dfImport[,x])))
table(i)

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
ivScore

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
sapply(n, function(x) {
  (abs(mCor[,x]) >= 0.7)
})

n = c('eczemaSeverity', 'scorad.at.sample', 'Ige3m0.10', 'eczema.at.sample', 'Ige12m0.35', 'skinPrick36m', 'scorad.at.12m', 
      'ecz.severity.at.12m', 'furryAnimal', 'ruralCategorical')

i = which(names(ivScore) %in% n)
ivScore = ivScore[-i]
data.frame(ivScore)

## check correlations again
dfData = data.frame(dfImport[,names(ivScore)])
m = NULL;

for (i in 1:ncol(dfData)){
  m = cbind(m, dfData[,i])
}
colnames(m) = colnames(dfData)
mCor = cor(m, use="na.or.complete")
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.7, names=T)
data.frame(n)
sapply(n, function(x) {
  (abs(mCor[,x]) >= 0.7)
})
n = c('ecz.severity.at.3m', 'weight12m')

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

## check variable combinations
# download the ccrossvalidation class to select variables of importance

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

dfData = data.frame(Shannon=dfImport$shannonmean, dfImport[,names(ivScore)])
dfData$abxMonth = factor(ifelse(dfData$abxMonth == 1, 'Yes', 'No'))
set.seed(123)
dfData.sub = na.omit(dfData)
dim(dfData.sub)
oVar.s = CVariableSelection.ReduceModel(dfData.sub[,-1], dfData.sub$Shannon , boot.num = 50)
plot.var.selection(oVar.s)
sapply(seq_along(1:29), function(x) CVariableSelection.ReduceModel.getMinModel(oVar.s, x))

### univariate report with and without interaction with frozen
dfData = data.frame(Shannon=dfImport$shannonmean, dfImport[,names(ivScore)])
dfData$abxMonth = factor(ifelse(dfData$abxMonth == 1, 'Yes', 'No'))
str(dfData)
cn = names(ivScore)

lFm = lapply(cn, function(x) {
  df = dfData[,c('Shannon', x)]
  fm = lm(Shannon ~ . , data=df)
})

names(lFm) = cn
library(car)
f_get.coef.pvalue = function(x){
  s = Anova(lFm[[x]])
  round(s$`Pr(>F)`[1], 3)
}

m = sapply(cn, function(x) {
  tryCatch(expr = f_get.coef.pvalue(x), error=function(e) NULL)
})

dfUnivariate = data.frame(m)
rownames(dfUnivariate) = names(m)
Sig = ifelse(m < 0.05, 'Yes', 'No')
dfUnivariate$Significant = Sig
colnames(dfUnivariate) = c('P-Value', 'Significant')

dir.create('mergedDataSet/Temp')
write.csv(dfUnivariate, file='mergedDataSet/Temp/uni.csv')

lapply(lFm[Sig == 'Yes'], summary)

### repeat but with frozen included in the model
str(dfData)
cn = names(ivScore)
i = which(cn == 'frozen')
cn = cn[-i]

lFm = lapply(cn, function(x) {
  df = dfData[,c('Shannon', 'frozen', x)]
  fm = lm(Shannon ~ frozen + . , data=df)
})

names(lFm) = cn

f_get.coef.pvalue = function(x){
  s = Anova(lFm[[x]])
  round(s$`Pr(>F)`[2], 3)
}

m = sapply(cn, function(x) {
  tryCatch(expr = f_get.coef.pvalue(x), error=function(e) NULL)
})

dfUnivariate = data.frame(m)
rownames(dfUnivariate) = names(m)
Sig = ifelse(m < 0.05, 'Yes', 'No')
dfUnivariate$Significant = Sig
colnames(dfUnivariate) = c('P-Value', 'Significant')

write.csv(dfUnivariate, file='mergedDataSet/Temp/bivariate.csv')

lapply(lFm[Sig == 'Yes'], summary)

### repeat but with frozen included and an interaction with frozen in the model
str(dfData)
cn = names(ivScore)
i = which(cn == 'frozen')
cn = cn[-i]

lFm = lapply(cn, function(x) {
  df = dfData[,c('Shannon', 'frozen', x)]
  fm = lm(Shannon ~ frozen*. , data=df)
})

names(lFm) = cn

f_get.coef.pvalue = function(x){
  s = Anova(lFm[[x]])
  round(s$`Pr(>F)`[1:3], 3)
}

m = sapply(cn, function(x) {
  tryCatch(expr = f_get.coef.pvalue(x), error=function(e) NULL)
})

dfUnivariate = data.frame(t(m))
colnames(dfUnivariate) = c('Frozen', 'Covariate-pvalue', 'Interaction-pvalue')
Sig = ifelse(m[3,] < 0.05, 'Yes', 'No')
dfUnivariate$Significant = Sig

write.csv(dfUnivariate, file='mergedDataSet/Temp/bivariate_interaction.csv')

lapply(lFm[Sig == 'Yes'], summary)

## the reason eczema needs an interaction with frozen as a lot of samples with eczema were frozen
xtabs(~ eczema3m + frozen, data=dfData)

## interaction only with eczema status at 3m
dfData = data.frame(Shannon=dfImport$shannonmean, dfImport[,names(ivScore)])
dfData$abxMonth = factor(ifelse(dfData$abxMonth == 1, 'Yes', 'No'))
fm01 = lm(Shannon ~ eczema3m*frozen + ., data=dfData)
summary(fm01)

Anova(fm01)
anova(fm01, test='Chisq')
s = Anova(fm01)
df = data.frame(s)
dfMultivariate = data.frame(PValue=round(df$Pr..F., 3))
rownames(dfMultivariate) = rownames(df)
write.csv(dfMultivariate, file='mergedDataSet/Temp/multivariate.csv')

## keep only the interesting covariates based on univariate and multivariate results
fm01.a = lm(Shannon ~ eczema3m*frozen + daysOldAtSample + dogOrCat + parentEczema + siblings +
              diarrhoeaCategorical3m + MaternalSchooling + abxMonth , data=dfData)

summary(fm01.a)
Anova(fm01.a)

s = Anova(fm01.a)
df = data.frame(s)
dfMultivariate = data.frame(PValue=round(df$Pr..F., 3))
rownames(dfMultivariate) = rownames(df)
write.csv(dfMultivariate, file='mergedDataSet/Temp/multivariate_submodel.csv')
## subset model further 
fm01.b = lm(Shannon ~ eczema3m*frozen + diarrhoeaCategorical3m + abxMonth , data=dfData)

summary(fm01.b)
Anova(fm01.b)
print(summary(fm01.b), digits=3)

## make some plots of interest

library(lattice)

xyplot(ifelse(dfData$eczema3m == "Yes", 1, 0) ~ Shannon | frozen, data = dfData,
       type = c("p", "smooth"),
       ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity given Frozen')

xyplot(ifelse(dfData$eczema3m == "Yes", 1, 0) ~ Shannon, data = dfData,
       type = c("p"),
       ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity')

xyplot(ifelse(dfData$eczema3m == "Yes", 1, 0) ~ Shannon | frozen, data = dfData,
       type = c("p"),
       ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity given Frozen')

xyplot(Shannon ~ abxMonth, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=T,
       type=c('p', 'g', 'r'), groups=eczema3m)

xyplot(Shannon ~ eczema3m+ifelse(abxMonth == 'Yes', 'AB Yes', 'AB No') | ifelse(frozen == 'Yes', 'Frozen', 'Not Frozen'), data=dfData, xlab='Antibiotic Yes/No, Eczema No Yes',
       main='Affect on Shannon Diversity of Eczema status and Antibiotics')

