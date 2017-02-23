# Name: 01_crossSectional3m_P2.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 23/02/2017
# Desc: relationships of diversity with covariates


p.old = par()

## data import
dfImport = read.csv(file.choose(), header=T)
## keep only the 3 m data
dfImport = dfImport[dfImport$age == '3m',]
# remove some duplicate and other columns
i = sapply(c('id$', 'mean', 'anyteeth1isyes'), grep, colnames(dfImport))
colnames(dfImport)[unlist(i)]
# remove these columns apart from chaomean
# [1] "id"              "shannonmean"     "chaomean"        "pdwholetreemean"
# [5] "simpsonmean"     "obsotumean"   
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
# which column has 70 NA's remove those
dfImport = dfImport[,-which(i > 60)]
dim(dfImport)
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
# remove lowest scores 
length(ivScore)
ivScore = ivScore[1:51]
tail(ivScore)
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

cvKeep = c('dogOrCat3m', 'eczema3m', 'foodAllergy', 'respInfection3m', 'countryside')
n = n[!(n%in% cvKeep)]
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
n = c('eczemaSev3m', 'skinPrick36m', 'catAtSample', 'urti3m', 'urban')

i = which(names(ivScore) %in% n)
ivScore = ivScore[-i]
data.frame(ivScore)

### perform another round of random forest
dfData = data.frame(Chao=dfImport$chaomean, dfImport[,names(ivScore)])
set.seed(123)
fit.rf = randomForest(Chao ~ ., data=na.omit(dfData))
# get variables importance
varImpPlot(fit.rf)
dfRF = data.frame(importance(fit.rf))
head(dfRF)
ivScore = dfRF$IncNodePurity
names(ivScore) = rownames(dfRF)
ivScore = sort(ivScore, decreasing = T)
data.frame(names(ivScore))


## keep only the interesting covariates based on previous analysis
dfData = data.frame(Chao=dfImport$chaomean, dfImport[,names(ivScore)])
dfData$abxMonth = factor(ifelse(dfData$abxMonth == 1, 'Yes', 'No'))
dfData$abx3Months = factor(ifelse(dfData$abx3Months == 1, 'Yes', 'No'))
library(heavy)
fm01.a = heavyLm(Chao ~ eczema3m*frozen + daysOldAtSample + weight3m + dogOrCat3m +
              parentEczema  + caesarean+ numSiblings3m +
              diarrhoea3m + MaternalSchooling + abxMonth , data=dfData, family=Student(df=3))

summary(fm01.a)
## subset model further 
fm01.b = heavyLm(Chao ~ eczema3m*frozen +  dogOrCat3m +
              caesarean+  
              diarrhoea3m  + abxMonth , data=dfData, family=Student(df=3))

summary(fm01.b)

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
dfData.full = dfData
dfData = dfData.full[,c('Chao', 'eczema3m', 'frozen', 'dogOrCat3m',
                          'caesarean', 'diarrhoea3m', 'abxMonth')]
dfData = na.omit(dfData)
mModMatrix = model.matrix(Chao ~ eczema3m*frozen +  dogOrCat3m +
                            caesarean+  
                            diarrhoea3m  + abxMonth , data=dfData)
lStanData = list(Ntotal=nrow(mModMatrix), Ncol=ncol(mModMatrix), X=mModMatrix, y=dfData$Chao, meanY=mean(dfData$Chao), 
                 sdY=sd(dfData$Chao))

fit.stan = sampling(stanDso, data=lStanData, iter=10000, chains=4, pars=c('betas', 'nu', 'sigma'))
print(fit.stan)

getStanMean(fit.stan)
getStanSD(fit.stan)
round(getStanPValue(fit.stan), 3)

#############################################################
## try another bayesian approach for this 
library(LearnBayes)
mylogpost = function(theta, data){
  nu = exp(theta['nu']) ## normality parameter for t distribution
  sigma = exp(theta['sigma']) # scale parameter for t distribution
  betas = theta[-c(1:2)] # vector of betas i.e. regression coefficients
  dfData = data$dfData # predictors and response 
  # function to use to use scale parameter
  ## see here https://grollchristian.wordpress.com/2013/04/30/students-t-location-scale/
  dt_ls = function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
  ## likelihood function
  lf = function(dat, pred){
    return(log(dt_ls(dat, nu, pred, sigma)))
  }
  # the predicted value
  mModMatrix = model.matrix(Chao ~ eczema3m*frozen +  dogOrCat3m +
                                           caesarean+  
                                           diarrhoea3m  + abxMonth, data=dfData)
  mu = mModMatrix %*% betas ## link is identity
  # likelihood function
  val = sum(lf(dfData$Chao, mu))
  val = val + dexp(nu, 1/29, log = T) + dunif(sigma, 1, 1e+4, log = T)
  return(val)
}

#mylogpost.inv = function(theta, data){return(1/mylogpost(theta, data))}
fitdistr(dfData$Chao, densfun = 't')
ivBetas = coef(lm(Chao ~ eczema3m*frozen +  dogOrCat3m +
                    caesarean+  
                    diarrhoea3m  + abxMonth , data=dfData))
start = c('nu'=log(3), 'sigma'=log(119), 'betas'=ivBetas)
lData = list('dfData'=dfData)
#nlm(mylogpost.inv, start, lData)
mylogpost(start, lData)
## get starting values
op = optim(start, mylogpost, control = list(fnscale = -1, maxit=10000), data=lData)
op
start = op$par
fit = laplace(mylogpost, start, lData)
fit
## redefine laplace to take iterations
mylaplace = function (logpost, mode, ...) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL, ..., hessian = TRUE, 
              control = list(fnscale = -1, maxit=10000))
  options(warn = 0)
  mode = fit$par
  h = -solve(fit$hessian)
  p = length(mode)
  int = p/2 * log(2 * pi) + 0.5 * log(det(h)) + logpost(mode, 
                                                        ...)
  stuff = list(mode = mode, var = h, int = int, converge = fit$convergence == 
                 0)
  return(stuff)
}

fit = mylaplace(mylogpost, start, lData)
fit

fit$mode
se = sqrt(diag(fit$var))
round(pnorm(-abs(fit$mode/se))*2,3)
exp(fit$mode['nu'])
exp(fit$mode['sigma'])

# 
# 
# 
# ## load some data and libraries
# library(LearnBayes)
# dfData = read.csv(file.choose(), header=T)
# 
# ## prepare input data
# lData = list('dfData'=data.frame(y=dfData$weight, height=dfData$height))
# # starting values
# f = lm(y ~ height, data=lData$dfData)
# start = c('nu'=0, 'sigma'=log(sd(lData$dfData$y)), 'betas'=coef(f))
# summary(f)
# mylogpost(start, lData)
# fit = laplace(mylogpost, start, lData)
# fit
# fit$mode
# se = sqrt(diag(fit$var))
# round(pnorm(-abs(fit$mode/se))*2,3)

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


