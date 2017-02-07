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
# remove these columns apart from chaomean
# [1] "id"              "shannonmean"     "chaomean"        "pdwholetreemean" "simpsonmean"    
# [6] "obsotumean"      "anyteeth1isyes" 
i = unlist(i)[-3]
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
fit.rf = randomForest(chaomean ~ ., data=na.omit(dfImport))
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
dfData = data.frame(Shannon=dfImport$chaomean, dfImport[,names(ivScore)])
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

dfData = data.frame(Chao=log(dfImport$chaomean), dfImport[,names(ivScore)])
dfData$abxMonth = factor(ifelse(dfData$abxMonth == 1, 'Yes', 'No'))
set.seed(123)
dfData.sub = na.omit(dfData)
dim(dfData.sub)
oVar.s = CVariableSelection.ReduceModel(dfData.sub[,-1], dfData.sub$Chao , boot.num = 50)
plot.var.selection(oVar.s)
sapply(seq_along(1:29), function(x) CVariableSelection.ReduceModel.getMinModel(oVar.s, x))

### univariate report with and without interaction with frozen
dfData = data.frame(Chao=dfImport$chaomean, dfImport[,names(ivScore)])
dfData$abxMonth = factor(ifelse(dfData$abxMonth == 1, 'Yes', 'No'))
str(dfData)
cn = names(ivScore)

lFm = lapply(cn, function(x) {
  df = dfData[,c('Chao', x)]
  fm = glm(Chao ~ . , data=df, family = Gamma(link='log'))
})

names(lFm) = cn
library(car)
f_get.coef.pvalue = function(x){
  s = Anova(lFm[[x]])
  round(s$`Pr(>Chisq)`[1], 3)
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
  df = dfData[,c('Chao', 'frozen', x)]
  fm = glm(Chao ~ frozen + . , data=df, family=Gamma(link='log'))
})

names(lFm) = cn

f_get.coef.pvalue = function(x){
  s = Anova(lFm[[x]])
  round(s$`Pr(>Chisq)`[2], 3)
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
  df = dfData[,c('Chao', 'frozen', x)]
  fm = glm(Chao ~ frozen*. , data=df, family=Gamma(link='log'))
})

names(lFm) = cn

f_get.coef.pvalue = function(x){
  s = Anova(lFm[[x]])
  round(s$`Pr(>Chisq)`[1:3], 3)
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
dfData = data.frame(Chao=dfImport$chaomean, dfImport[,names(ivScore)])
dfData$abxMonth = factor(ifelse(dfData$abxMonth == 1, 'Yes', 'No'))
fm01 = glm(Chao ~ eczema3m*frozen + ., data=dfData, family=Gamma(link='log'))
summary(fm01)

Anova(fm01)
anova(fm01, test='Chisq')
s = Anova(fm01)
df = data.frame(s)
dfMultivariate = data.frame(PValue=round(df$Pr..Chisq., 3))
rownames(dfMultivariate) = rownames(df)
write.csv(dfMultivariate, file='mergedDataSet/Temp/multivariate.csv')

## keep only the interesting covariates based on univariate and multivariate results
fm01.a = glm(Chao ~ eczema3m*frozen + daysOldAtSample + dogOrCat + parentEczema + siblings +
              diarrhoeaCategorical3m + MaternalSchooling + abxMonth , data=dfData, family=Gamma(link='log'))

summary(fm01.a)
Anova(fm01.a)

s = Anova(fm01.a)
df = data.frame(s)
dfMultivariate = data.frame(PValue=round(df$Pr..Chisq., 3))
rownames(dfMultivariate) = rownames(df)
write.csv(dfMultivariate, file='mergedDataSet/Temp/multivariate_submodel.csv')
## subset model further 
fm01.b = glm(Chao ~ eczema3m*frozen + diarrhoeaCategorical3m + abxMonth , data=dfData, family=Gamma(link='log'))

summary(fm01.b)
Anova(fm01.b)
print(summary(fm01.b), digits=3)

## make some plots of interest

library(lattice)

xyplot(ifelse(dfData$eczema3m == "Yes", 1, 0) ~ Chao | frozen, data = dfData,
       type = c("p", "smooth"),
       ylab = "Eczema at 3 months", xlab = "Chao Diversity", main='Eczema at 3m vs Diversity given Frozen')

xyplot(ifelse(dfData$eczema3m == "Yes", 1, 0) ~ Chao, data = dfData,
       type = c("p"),
       ylab = "Eczema at 3 months", xlab = "Chao Diversity", main='Eczema at 3m vs Diversity')

xyplot(ifelse(dfData$eczema3m == "Yes", 1, 0) ~ Chao | frozen, data = dfData,
       type = c("p"),
       ylab = "Eczema at 3 months", xlab = "Chao Diversity", main='Eczema at 3m vs Diversity given Frozen')

xyplot(Chao ~ abxMonth, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=T,
       type=c('p', 'g', 'r'), groups=eczema3m)

xyplot(Chao ~ eczema3m+ifelse(abxMonth == 'Yes', 'AB Yes', 'AB No') | ifelse(frozen == 'Yes', 'Frozen', 'Not Frozen'), data=dfData, xlab='Antibiotic Yes/No, Eczema No Yes',
       main='Affect on Chao Diversity of Eczema status and Antibiotics')

#############################################################
## MCMC RStan

getStanSD = function(obj){
  return(apply(extract(obj)$betas, 2, sd))
}
getStanMean = function(obj){
  return(apply(extract(obj)$betas, 2, mean))
}
getStanPValue = function(obj){
  pnorm(-abs(getStanMean(obj)/getStanSD(obj)))*2
}

library(rstan)
stanDso = rstan::stan_model(file='mergedDataSet/robustRegressionT.stan')

## prepare data for input
mModMatrix = model.matrix(weight ~ height, data=dfData)
lStanData = list(Ntotal=nrow(mModMatrix), Ncol=ncol(mModMatrix), X=mModMatrix, y=dfData$weight, meanY=mean(dfData$weight), 
                 sdY=sd(dfData$weight))

fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=4, pars=c('betas', 'nu', 'sigma'))
print(fit.stan)

getStanMean(fit.stan)
getStanSD(fit.stan)


#############################################################
## try a bayesian approach for this 
mylogpost = function(theta, data){
  nu = exp(theta['nu']) ## normality parameter for t distribution
  sigma = exp(theta['sigma']) # scale parameter for t distribution
  betas = theta[3:4] # vector of betas i.e. regression coefficients
  dfData = data$dfData # predictors and response 
  # function to use to use scale parameter
  ## see here https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
  dt_ls = function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
  ## likelihood function
  lf = function(dat, pred){
    return(log(dt_ls(dat, nu, pred, sigma)))
  }
  # the predicted value
  mModMatrix = model.matrix(y ~ ., data=dfData)
  mu = mModMatrix %*% betas ## link is identity
  # likelihood function
  val = sum(lf(dfData$y, mu))
  val = val + dexp(nu, 1/29, log = T) + dunif(sigma, 1, 1e+4, log = T)
  return(val)
}

mylogpost.inv = function(theta, data){return(1/mylogpost(theta, data))}
start = c('nu'=1, 'sigma'=1, 'betas'=c(1, 1))
nlm(mylogpost.inv, start, lData)

## load some data and libraries
library(LearnBayes)
dfData = read.csv(file.choose(), header=T)

## prepare input data
lData = list('dfData'=data.frame(y=dfData$weight, height=dfData$height))
# starting values
f = lm(y ~ height, data=lData$dfData)
start = c('nu'=0, 'sigma'=log(sd(lData$dfData$y)), 'betas'=coef(f))
summary(f)
mylogpost(start, lData)
fit = laplace(mylogpost, start, lData)
fit
fit$mode
se = sqrt(diag(fit$var))
round(pnorm(-abs(fit$mode/se))*2,3)

# mylogpost = function(theta, data){
#   #nu = exp(theta['nu']) ## normality parameter for t distribution
#   sigma = exp(lData$sigma) # scale parameter for t distribution
#   betas = theta[1:2] # vector of betas i.e. regression coefficients
#   dfData = data$dfData # predictors and response 
#   # function to use to use scale parameter
#   ## see here https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
#   dt_ls = function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
#   ## likelihood function
#   lf = function(dat, pred){
#     return(log(dnorm(dat, pred, sigma)))
#   }
#   # the predicted value
#   mModMatrix = model.matrix(y ~ ., data=dfData)
#   mu = mModMatrix %*% betas ## link is identity
#   # likelihood function
#   val = sum(lf(dfData$y, mu))
#   val = val 
#   return(val)
# }
# 
# start = c('betas'=c('beta0'=-100, 'beta1'=0))
# lData$sigma = log(sd(lData$dfData$y))
# 
# fit = laplace(mylogpost, start, lData)
# se = sqrt(diag(fit$var))


