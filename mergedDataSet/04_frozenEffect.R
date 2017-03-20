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
dfData = data.frame(resp=dfImport$shannonmean, #frozen=dfImport$frozen, 
                    abxMonth=ifelse(dfImport$abxMonth == 0, 'No', 'Yes'),
                    dogOrCat3m=dfImport$dogOrCat3m,
                    eczema=dfImport$eczema3m,
                    diarrhoea=dfImport$diarrhoea3m)

dfData = na.omit(dfData)


### write the functions for the analysis
library(LearnBayes)
library(MASS)
library(numDeriv)
mylogpost = function(theta, data){
  sigma = exp(theta['sigma']) # scale parameter for normal distribution
  # hyperparameter for the hierarchichal standard deviation parameter
  betaSigma = 1; #exp(theta['betaSigma'])
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


#### initial values
ivBetas = coef(lm(resp ~ . , data=dfData))

#start = c('sigma'=log(sd(dfData$resp)), 'betaSigma'=log(sd(dfData$resp)), 'betas'=ivBetas)
start = c('sigma'=log(sd(dfData$resp)), 'betas'=ivBetas)
lData = list('dfData'=dfData)
#lData$betaParam = c('shape'=0.5, 'rate'=0.0001)

#fit = laplace(mylogpost, start, lData)
fit2 = mylaplace(mylogpost, start, lData)
fit3 = (lm(resp ~ . , data=dfData))

getOptimizedSummary = function(obj){
  se = sqrt(abs(diag(obj$var)))
  m = obj$mode
  es = m/se
  p = pnorm(-abs(m/se))*2
  ret  = cbind('Coef' = round(m, 3), 'SE'=round(se, 3), 'ES'=es, 'P-Value'=signif(p, 3))
  colnames(ret) = c('Coef', 'SE', 'ES', 'P-Value')
  return(ret)
}

mFrozenNo = round(getOptimizedSummary(fit2), 3)
summary(fit3)
#####################################################################
## with frozen
dfData = data.frame(resp=dfImport$shannonmean, frozen=ifelse(dfImport$frozen=='Yes', 1, 0), 
                    abxMonth=ifelse(dfImport$abxMonth == 0, 'No', 'Yes'),
                    dogOrCat3m=dfImport$dogOrCat3m,
                    eczema=dfImport$eczema3m,
                    diarrhoea=dfImport$diarrhoea3m)

dfData = na.omit(dfData)
## initial values
ivBetas = coef(lm(resp ~ . , data=dfData))
start = c('sigma'=log(sd(dfData$resp)), 'betas'=ivBetas)
lData = list('dfData'=dfData)

fit2.frozen = mylaplace(mylogpost, start, lData)
mFrozenYes = round(getOptimizedSummary(fit2.frozen), 3)

# 
# ## make some plots of interest
# 
# library(lattice)
# 
# xyplot(ifelse(dfData$eczema3m == "Yes", 1, 0) ~ Shannon | frozen, data = dfData,
#        type = c("p", "smooth"),
#        ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity given Frozen')
# 
# xyplot(ifelse(dfData$eczema3m == "Yes", 1, 0) ~ Shannon, data = dfData,
#        type = c("p"),
#        ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity')
# 
# xyplot(ifelse(dfData$eczema3m == "Yes", 1, 0) ~ Shannon | frozen, data = dfData,
#        type = c("p"),
#        ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity given Frozen')
# 
# xyplot(Shannon ~ abxMonth, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=T,
#        type=c('p', 'g', 'r'), groups=eczema3m)
# 
# xyplot(Shannon ~ eczema3m+ifelse(abxMonth == 'Yes', 'AB Yes', 'AB No') | ifelse(frozen == 'Yes', 'Frozen', 'Not Frozen'), data=dfData, xlab='Antibiotic Yes/No, Eczema No Yes',
#        main='Affect on Shannon Diversity of Eczema status and Antibiotics')

