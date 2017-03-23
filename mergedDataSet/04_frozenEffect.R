# Name: 04_frozenEffect.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 20/03/2017
# Desc: the effect of freezing on diversity


p.old = par()

## data import
dfImport = read.csv(file.choose(), header=T)
## keep only the 3 m data
dfImport = dfImport[dfImport$age == '3m',]
# remove some duplicate and other columns
i = sapply(c('id$', 'mean', 'anyteeth1isyes'), grep, colnames(dfImport))
colnames(dfImport)[unlist(i)]
# remove these columns apart from shannonmean
# [1] "id"              "shannonmean"     "chaomean"        "pdwholetreemean"
# [5] "simpsonmean"     "obsotumean"   
i = unlist(i)[-2]
dfImport = dfImport[,-i]
dfImport = dfImport[,-2]
dfImport = dfImport[,-1]
dim(dfImport)
str(dfImport)
summary(dfImport)

### covariates of interest for analysis 
dfData = data.frame(resp=dfImport$shannonmean, frozen=dfImport$frozen, 
                    abxMonth=ifelse(dfImport$abxMonth == 0, 'No', 'Yes'),
                    dogOrCat3m=dfImport$dogOrCat3m,
                    eczema=dfImport$eczema3m,
                    diarrhoea=dfImport$diarrhoea3m,
                    siblings=dfImport$sibling3m,
                    caes = dfImport$caesarean)

dfData = na.omit(dfData)


### write the functions for the analysis
library(LearnBayes)
library(MASS)
library(numDeriv)
mylogpost = function(theta, data){
  sigma = exp(theta['sigma']) # scale parameter for normal distribution
  # hyperparameter for the hierarchichal standard deviation parameter of regression coefficients
  betaSigma = exp(theta['betaSigma'])
  betas = theta[-c(1:2)] # vector of betas i.e. regression coefficients
  dfData = data$dfData # predictors and response 
  # hyper-hyperparameters for the hierarchichal standard deviation parameter
  ivBetaParam = data$betaParam
  
  ## likelihood function
  lf = function(dat, pred){
    return(log(dnorm(dat, pred, sigma)))
  }
  # the predicted value
  mModMatrix = model.matrix(resp ~ ., data=dfData)
  mu = mModMatrix %*% betas ## link is identity
  # likelihood function
  val = sum(lf(dfData$resp, mu))
  # add the priors
  val = val  + dunif(sigma, 0.1, 1e+4, log = T) +
    sum(dcauchy(betas, 0, betaSigma, log = T)) + 
    dgamma(betaSigma, shape = ivBetaParam['shape'], rate = ivBetaParam['rate'], log=T)
    # dcauchy(betas[1], 0, 1e2, log=T) + 
    # sum(dcauchy(betas[2:length(betas)], 0, betaSigma, log=T)) + 
    # dgamma(betaSigma, shape = ivBetaParam['shape'], rate = ivBetaParam['rate'], log=T)
  return(val)
}

mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=10000), data=data)
  # calculate hessian
  fit$hessian = (hessian(logpost, fit$par, data=data))
  colnames(fit$hessian) = names(mode)
  rownames(fit$hessian) = names(mode)
  options(warn = 0)
  mode = fit$par
  h = -solve(fit$hessian)
  stuff = list(mode = mode, var = h, converge = fit$convergence == 
                 0)
  return(stuff)
}

getOptimizedSummary = function(obj){
  se = sqrt(abs(diag(obj$var)))
  m = obj$mode
  es = m/se
  p = pnorm(-abs(m/se))*2
  ret  = cbind('Coef' = round(m, 3), 'SE'=round(se, 3), 'ES'=es, 'P-Value'=signif(p, 3))
  colnames(ret) = c('Coef', 'SE', 'ES', 'P-Value')
  return(ret)
}


fit.model = function(dfData){
  #### initial values for the optimizer
  ivBetas = coef(lm(resp ~ . , data=dfData))
  start = c('sigma'=log(sd(dfData$resp)), 'betaSigma'=log(sd(dfData$resp)), 'betas'=ivBetas)
  # parameters and data
  lData = list('dfData'=dfData)
  # jeffery's prior for gamma distribution, hyper-prior distribution for scale of coefficients
  lData$betaParam = c('shape'=0.5, 'rate'=0.0001)
  fit = laplace(mylogpost, start, lData)
  return(getOptimizedSummary(fit))
}

########### subset the data based on frozen and fresh
dfFresh = dfData[dfData$frozen == 'No',]
dfFrozen = dfData[dfData$frozen == 'Yes',]
dim(dfFresh); dim(dfFrozen)

#### check each covariate with fresh and frozen samples
colnames(dfData)

fit.model.split = function(dfData){
  mFresh = fit.model(dfData[dfData$frozen == 'No',c(1,3)])
  mFrozen = fit.model(dfData[dfData$frozen == 'Yes',c(1,3)])
  return(list('fresh' = mFresh, 'frozen'=mFrozen))
}

labx = fit.model.split(dfData[,c(1, 2, 3)])
lpet = fit.model.split(dfData[,c(1, 2, 4)])
leczema = fit.model.split(dfData[,c(1, 2, 5)])
ldia = fit.model.split(dfData[,c(1, 2, 6)])

summary(lm(resp ~ diarrhoea, data=dfData[dfData$frozen == 'Yes',]))

#########################################################################################
####### redefine logposterior function and use a constant hyperparameter
mylogpost = function(theta, data){
  sigma = exp(theta['sigma']) # scale parameter for normal distribution
  # hyperparameter for the hierarchichal standard deviation parameter of regression coefficients
  betaSigma = 10#exp(theta['betaSigma'])
  betas = theta[-c(1)] # vector of betas i.e. regression coefficients
  dfData = data$dfData # predictors and response 
  # hyper-hyperparameters for the hierarchichal standard deviation parameter
  #ivBetaParam = data$betaParam
  
  ## likelihood function
  lf = function(dat, pred){
    return(log(dnorm(dat, pred, sigma)))
  }
  # the predicted value
  mModMatrix = model.matrix(resp ~ ., data=dfData)
  mu = mModMatrix %*% betas ## link is identity
  # likelihood function
  val = sum(lf(dfData$resp, mu))
  # add the priors
  val = val  + dunif(sigma, 0.1, 1e+4, log = T) +
    sum(dcauchy(betas, 0, betaSigma, log = T))  
    #dgamma(betaSigma, shape = ivBetaParam['shape'], rate = ivBetaParam['rate'], log=T)
  # dcauchy(betas[1], 0, 1e2, log=T) + 
  # sum(dcauchy(betas[2:length(betas)], 0, betaSigma, log=T)) + 
  # dgamma(betaSigma, shape = ivBetaParam['shape'], rate = ivBetaParam['rate'], log=T)
  return(val)
}

fit.model = function(dfData){
  #### initial values for the optimizer
  ivBetas = coef(lm(resp ~ . , data=dfData))
  start = c('sigma'=log(sd(dfData$resp)), 'betas'=ivBetas)
  # parameters and data
  lData = list('dfData'=dfData)
  fit = mylaplace(mylogpost, start, lData)
  return(getOptimizedSummary(fit))
}

fit.model.split = function(dfData){
  mFresh = fit.model(dfData[dfData$frozen == 'No',c(1,3)])
  print(summary(lm(resp ~ ., data=dfData[dfData$frozen == 'No',c(1,3)])))
  ## select equal sizes for frozen
  df = dfData[dfData$frozen == 'Yes',c(1,3)]
  df = df[sample(1:182, size = 81, replace = F),]
  print(summary(lm(resp ~ ., data=df)))
  mFrozen = fit.model(df)
  return(list('fresh' = mFresh, 'frozen'=mFrozen))
}



labx = fit.model.split(dfData[,c(1, 2, 3)])
lpet = fit.model.split(dfData[,c(1, 2, 4)])
leczema = fit.model.split(dfData[,c(1, 2, 5)])
ldia = fit.model.split(dfData[,c(1, 2, 6)])
lsib = fit.model.split(dfData[,c(1, 2, 7)])
lcae = fit.model.split(dfData[,c(1, 2, 8)])

dfFresh = data.frame(rbind(labx$fresh[3,], lpet$fresh[3,], leczema$fresh[3,], ldia$fresh[3,], lsib$fresh[3,], lcae$fresh[3,]))
dfFrozen = data.frame(rbind(labx$frozen[3,], lpet$frozen[3,], leczema$frozen[3,], ldia$frozen[3,], lsib$frozen[3,], lcae$frozen[3,]))
rownames(dfFresh) = c('antibiotics', 'pets', 'eczema', 'diahorrea', 'siblings', 'caesarean')
rownames(dfFrozen) = c('antibiotics', 'pets', 'eczema', 'diahorrea', 'siblings', 'caesarean')

dfFresh[order(dfFresh$ES),]
dfFrozen[order(dfFrozen$ES),]

