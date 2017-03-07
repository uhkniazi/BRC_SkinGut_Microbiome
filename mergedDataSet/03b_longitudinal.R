# Name: 03b_longitudinal.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 02/03/2017
# Desc: longitudinal analysis of EAT dataset


p.old = par()
dfImport = read.csv(file.choose(), header=T)

dim(dfImport)

# setup data
dfData = dfImport
dfData$age = as.numeric(gsub('m', '', dfData$age))
rownames(dfData) = dfImport$samples
dfData$samples = factor(dfData$samples)

str(dfData)

## extract repeated samples
ids = which(duplicated(as.character(dfData$id)))
ids = unique(as.character(dfData$id[ids]))
dfData = dfData[dfData$id %in% ids,]
dfData = droplevels.data.frame(dfData)
# sanity checks
table(dfData$id)


library(lattice)
library(MASS)
library(car)
library(lme4)
library(lmerTest)

# check data distribution
hist(dfData$chaomean)
fitdistr(dfData$chaomean, 'gamma')
qqPlot(dfData$chaomean, 'gamma', shape=6, rate=0.012,  ylab='Chao1 diversity')

### trying model fits with covariates
# model with just random intercept
fm01 = glmer(chaomean ~ 1 + (1 | id), data=dfData, family=Gamma(link='identity'))
summary(fm01)
fm01.cor = update(fm01, chaomean ~ 1 + (1 + age | id))
summary(fm01.cor)
# compare the 2 models
anova(fm01, fm01.cor)

# try uncorrelated random effects
fm01.uncor = update(fm01, chaomean ~ 1 + (1 | id) + (0 + age | id))
summary(fm01.uncor)
# compare 3 models
anova(fm01.uncor, fm01.cor, fm01 )

## choose the model with correlated intercept and slope
fm01 = fm01.cor

## update this model with one covariate after adding time
fm02 = update(fm01, chaomean ~ 1 + age + (1 + age | id))
summary(fm02)
anova(fm01, fm02)

# test the four models in a sequence
fm00 = update(fm01, chaomean ~ 1 + (1 + age | id))
fm01 = update(fm01, chaomean ~ 1 + age + (1 + age | id))
fm02 = update(fm01, chaomean ~ 1 + age + intervention + (1 + age | id))
fm03 = update(fm01, chaomean ~ 1 + age*intervention + (1 + age | id), family=Gamma(link='log'))
anova(fm00, fm01, fm02, fm03)

# fit the get p-values for coefficients
fm = glmer(chaomean ~ 1 + age + intervention + (1 + age | id), data=dfData, family=Gamma(link='log'))
summary(fm)

# clean up data first before proceeding
# which columns have the highest NAs
i = sapply(1:ncol(dfData), function(x) sum(is.na(dfData[,x])))
table(i)
# which column has 40 NA's remove those
dfData = dfData[,-which(i > 40)]
dim(dfData)

cvCols = colnames(dfData)
# drop some of the columns that are redundant
cvCols = cvCols[!(cvCols %in% c('samples', "shannonmean", "pdwholetreemean",
                                "simpsonmean", "obsotumean"))]
dfData = dfData[,cvCols]

### finding correlated covariates
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
cvKeep = c('foodAllergy', 'dogOrCat3m', 'dogOrCatAtSample', 'diet6FoodsAtSample', 'dietDiversityAtSample',
             'eczema3m', 'eczema12m')
n = n[!(n%in% cvKeep)]

i = which(cvCols %in% n)
cvCols = cvCols[-i]
data.frame(cvCols)

## check correlations again
dfData = dfData[,cvCols]
dim(dfData)
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
n = c('eczemaSev3m', 'scorad3m', 'skinPrickAero12m', 'skinPrickAero36m', 'DogsAndCatsAtSampleImputed5m',
      'furryAnimal3m', 'catsCount3m', 'catAtSample', 'catAtSample', 'eczemaSev12m', 
      'scorad12m')

i = which(cvCols %in% n)
cvCols = cvCols[-i]
data.frame(cvCols)

dfData = dfData[,cvCols]
dim(dfData)


# ## fit a model with each covariate 
# f_get.coef.pvalue = function(fit){
#   s = summary(fit)
#   return(s$coefficients[,'Pr(>|t|)'])
# }
# 
# f_get.coef.estimate = function(fit){
#   s = summary(fit)
#   return(s$coefficients[,'Estimate'])
# }
# 
# cvCov = colnames(dfData)[4:34]
# lFits = vector('list', length=length(cvCov))
# 
# for(i in 1:length(lFits)){
#   x = cvCov[i]
#   print(x)
#   df = dfData[,c('id', 'age', 'chaomean', x)]
#   cvFormula = paste('chaomean ~ 1 + ', x,' + (1 + age | id)', sep='')
#   ## get stats
#   getstats = function() {s = summary(lmerTest::lmer(cvFormula, data=df, REML=T));
#   return(c(f_get.coef.pvalue(s), f_get.coef.estimate(s)))}
#   
#   lFits[[i]] = tryCatch(expr = getstats() , error=function(e) NULL) 
#   #return(lmerTest::lmer(cvFormula, data=df, REML=F))
# }
# 
# names(lFits) = cvCov
# ## which models did not fit
# bNotFit = sapply(lFits, is.null)
# which(bNotFit)
# #lFits = lFits[!bNotFit]

cvCov = colnames(dfData)[4:34]
lFits = vector('list', length=length(cvCov))

for(i in 1:length(lFits)){
  x = cvCov[i]
  print(x)
  df = dfData[,c('id', 'age', 'chaomean', x)]
  cvFormula = paste('chaomean ~ 1 + ', x,' + (1 + age | id)', sep='')
  lFits[[i]] = Anova(glmer(cvFormula, data=df, family=Gamma(link='identity')))
}

names(lFits) = cvCov

f_get.coef.pvalue = function(s){
  round(s$`Pr(>Chisq)`[1], 3)
}

m = sapply(lFits, function(x) {
  tryCatch(expr = f_get.coef.pvalue(x), error=function(e) NULL)
})

dfUnivariate = data.frame(m)
rownames(dfUnivariate) = names(m)
Sig = ifelse(m < 0.05, 'Yes', 'No')
dfUnivariate$Significant = Sig
colnames(dfUnivariate) = c('P-Value', 'Significant')

dir.create('mergedDataSet/Temp')
write.csv(dfUnivariate, file='mergedDataSet/Temp/uni_identity.csv')

cvCov = rownames(dfUnivariate[dfUnivariate$Significant == 'Yes',])

for(i in seq_along(cvCov)){
  x = cvCov[i]
  cat('-------------------------------------------\n')
  df = dfData[,c('id', 'age', 'chaomean', x)]
  cvFormula = paste('chaomean ~ 1 + ', x,' + (1 + age | id)', sep='')
  print(summary(glmer(cvFormula, data=df, family=Gamma(link='identity'))))
  cat('------------------------------------------\n\n')
}

### multivariate with all covariates
which(colnames(dfData) %in% 'abxWeek')
dfData = dfData[,-28]
dfData$abxMonth = factor(ifelse(dfData$abxMonth == 1, 'Yes', 'No'))
dfData$abx3Months = factor(ifelse(dfData$abx3Months == 1, 'Yes', 'No'))

fm01 = glmer(chaomean ~ 1 + age + intervention +
                    (1 + age | id), data=dfData, family=Gamma(link='log'))
summary(fm01)

fm02 = glmer(chaomean ~ 1 + age + intervention + 
                        daysRoomTemp + 
                        (1 + age | id), data=dfData, family=Gamma(link='log'))

summary(fm02)

fm03 = glmer(chaomean ~ 1 + age + intervention + 
                        daysRoomTemp + diet6FoodsAtSample +
                        (1 + age | id), data=dfData, family=Gamma(link='log'))

summary(fm03)

fm04 = glmer(chaomean ~ 1 + age + intervention + 
                        daysRoomTemp + diet6FoodsAtSample + dietDiversityAtSample +
                        (1 + age | id), data=dfData, family=Gamma(link='log'))
summary(fm04)
## diet variables are correlated with age and including them makes 
## coefficients unstable
summary(glmer(chaomean ~ 1 + age + intervention + 
                        daysRoomTemp + caesarean + sibling3m + frozen +abxMonth + diarrDaysPreSample + sterilise +
                        (1 + age | id), data=dfData, family=Gamma(link='log')))

fm.mv = glmer(chaomean ~ 1 + age + intervention + 
                daysRoomTemp + caesarean + sibling3m + frozen +abxMonth + diarrDaysPreSample + sterilise +
                (1 + age | id), data=dfData, family=Gamma(link='log'))

Anova(fm.mv)
anova(fm.mv, test='Chisq')



###########################################################################################################################
########################## custom regression functions - just playing around
## custom optimiser
library(LearnBayes)
library(numDeriv)
library(MASS)

mylaplace = function (logpost, mode, data) 
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL,  
              control = list(fnscale = -1, maxit=1000), method='CG', data=data)
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

######### add a custom random effects model with gamma response and log link
# myloglike = function(theta, data){
#   ## parameters to track/estimate
#   sigmaPop = exp(theta['sigmaPop'])
#   sigmaRan = exp(theta['sigmaRan'])
#   betas = theta[3:4] # vector of betas i.e. regression coefficients for population
#   #iGroupsJitter = theta[-c(1:4)]
#   ## data
#   dfData = data$dfData # predictors and response
#   ## likelihood function
#   lf = function(dat, fitted){
#     alpha = (fitted^2)/sigmaPop
#     beta = fitted/sigmaPop
#     return(log(dgamma(dat, shape = alpha, rate = beta)))
#   }
#   ## random effect jitter for the population intercept
#   ## how many groups
#   iGroupsCount = unique(dfData$groupIndex)
#   # each group contributes a jitter centered on 0
#   iGroupsJitter = rnorm(iGroupsCount, 0, sigmaRan)
#   # population slope + random jitter
#   ivBetaRand = betas[1] + iGroupsJitter
#   # create model matrix
#   mModMatrix = model.matrix(chaomean ~ 1+age, data=dfData)
#   # create a matrix of betas
#   ivIntercept = ivBetaRand[dfData$groupIndex]
#   mBetas = matrix(rep(betas, times=length(ivIntercept)), ncol = length(betas), byrow = T)
#   mBetas[,1] = ivIntercept
#   mBetas = t(mBetas)
#   iFitted = diag(mModMatrix %*% mBetas)
#   # using log link
#   # write the likelihood function
#   val = sum(lf(dfData$resp, exp(iFitted)))
#   val = val +  dunif(sigmaPop, 0, 1e10, log=T) + sum(dnorm(iGroupsJitter, 0, sigmaRan, log=T)) +
#     dunif(sigmaRan, 0, 1e3, log=T)
#   return(val)
# }

dfData.re = droplevels.data.frame(na.omit(dfData[,c('id', 'age', 'chaomean', 'intervention', 'daysRoomTemp', 'caesarean', 'sibling3m')]))
dfData.re$groupIndex = as.numeric(dfData.re$id)
dfData.re$resp = dfData.re$chaomean
# # set starting values
# 
# start = c('sigmaPop'=log(5), 'sigmaRan'=log(2), rep(0, times=2))#, rep(0, times=length(unique(dfData.re$groupIndex)))) 
# 
# dfData.re$resp = dfData.re$chaomean
# lData = list('dfData'=dfData.re)
# 
# myloglike(start, lData)
# 
# op = optim(start, myloglike, control = list(fnscale = -1, maxit=10000), data=lData, method='SANN', hessian=T)
# start = as.numeric(op$par)
# names(start)[1:2] = c('sigmaPop', 'sigmaRan')
# h = hessian(myloglike, start, data=lData)
# # 
# # start = as.numeric(op$par)
# # names(start)[1:2] = c('sigmaPop', 'sigmaRan')
# # start
# # op = optim(start, myloglike, control = list(fnscale = -1, maxit=20000), data=lData, hessian=T)
# 
# h = -solve(h)
# op$var = h
# fit = op
# # 
# # start[3:4] = op$par[3:4]
# # names(start) = NULL
# # names(start)[1:2] = c('sigmaPop', 'sigmaRan')
# # fit = mylaplace(myloglike, start, lData)
# 
# se = sqrt(abs(diag(fit$var)))
# m = fit$par

fit.lm = glmer(resp ~ 1+ age + intervention + daysRoomTemp+caesarean+sibling3m+ (1 | id), data=dfData.re, family=Gamma(link='log'))
summary(fit.lm)
# m[3:4]
# se[3:4]

myloglike = function(theta, data){
  ## parameters to track/estimate
  sigmaPop = exp(theta['sigmaPop'])
  sigmaRan = exp(theta['sigmaRan'])
  betas = theta[3:7] # vector of betas i.e. regression coefficients for population
  iGroupsJitter = theta[-c(1:7)]
  ## data
  dfData = data$dfData # predictors and response
  ## likelihood function
  lf = function(dat, fitted){
    alpha = (fitted^2)/sigmaPop
    beta = fitted/sigmaPop
    return(log(dgamma(dat, shape = alpha, rate = beta)))
  }
  ## random effect jitter for the population intercept
  ## how many groups
  #iGroupsCount = unique(dfData$groupIndex)
  # each group contributes a jitter centered on 0
  #iGroupsJitter = rnorm(iGroupsCount, 0, sigmaRan)
  # population slope + random jitter
  ivBetaRand = betas[1] + iGroupsJitter
  # create model matrix
  mModMatrix = model.matrix(resp ~ 1+ age + daysRoomTemp+caesarean+sibling3m, data=dfData)
  # create a matrix of betas
  ivIntercept = ivBetaRand[dfData$groupIndex]
  mBetas = matrix(rep(betas, times=length(ivIntercept)), ncol = length(betas), byrow = T)
  mBetas[,1] = ivIntercept
  mBetas = t(mBetas)
  iFitted = diag(mModMatrix %*% mBetas)
  # using log link
  # write the likelihood function
  val = sum(lf(dfData$resp, exp(iFitted)))
  val = val +  dunif(sigmaPop, 0, 1e10, log=T) + sum(dnorm(iGroupsJitter, 0, sigmaRan, log=T)) + 
    sum(dcauchy(betas, 0, 10, log=T)) +
    dunif(sigmaRan, 0, 1e3, log=T) 
  return(val)
}

# set starting values
start = c('sigmaPop'=log(5), 'sigmaRan'=log(2), rep(0, times=5), rep(0, times=length(unique(dfData.re$groupIndex)))) 

lData = list('dfData'=dfData.re)

myloglike(start, lData)

# op = optim(start, myloglike, control = list(fnscale = -1, maxit=10000), data=lData, method='SANN', hessian=T)
# start = as.numeric(op$par)
# names(start)[1:2] = c('sigmaPop', 'sigmaRan')
# h = hessian(myloglike, start, data=lData)
# # 
# # start = as.numeric(op$par)
# # names(start)[1:2] = c('sigmaPop', 'sigmaRan')
# # start
# # op = optim(start, myloglike, control = list(fnscale = -1, maxit=20000), data=lData, hessian=T)
# 
# h = -solve(h)
# op$var = h
# fit = op
# # 
# # start[3:4] = op$par[3:4]
# # names(start) = NULL
# # names(start)[1:2] = c('sigmaPop', 'sigmaRan')
# # fit = mylaplace(myloglike, start, lData)
# 
# se = sqrt(abs(diag(fit$var)))
# m = fit$par

# fit.lm = glmer(chaomean ~ 1+age + (1 | id), data=dfData, family=Gamma(link='log'))
# summary(fit.lm)
# m[3:4]
# se[3:4]
# m[3:4]/se[3:4]

#start = c('sigmaPop'=log(5), 'sigmaRan'=log(2), rep(0, times=2), rep(0, times=length(unique(dfData.re$groupIndex)))) 
fit = mylaplace(myloglike, start, lData)
fit$mode[3:7]
se = sqrt(abs(diag(fit$var)))
se[3:7]

############################## try without tracking RE intercepts
myloglike = function(theta, data){
  ## parameters to track/estimate
  sigmaPop = exp(theta['sigmaPop'])
  sigmaRan = exp(theta['sigmaRan'])
  betas = theta[3:7] # vector of betas i.e. regression coefficients for population
  #iGroupsJitter = theta[-c(1:7)]
  ## data
  dfData = data$dfData # predictors and response
  ## likelihood function
  lf = function(dat, fitted){
    alpha = (fitted^2)/sigmaPop
    beta = fitted/sigmaPop
    return(log(dgamma(dat, shape = alpha, rate = beta)))
  }
  ## random effect jitter for the population intercept
  ## how many groups
  iGroupsCount = unique(dfData$groupIndex)
  # each group contributes a jitter centered on 0
  iGroupsJitter = rnorm(iGroupsCount, 0, sigmaRan)
  # population slope + random jitter
  ivBetaRand = betas[1] + iGroupsJitter
  # create model matrix
  mModMatrix = model.matrix(resp ~ 1+ age + daysRoomTemp+caesarean+sibling3m, data=dfData)
  # create a matrix of betas
  ivIntercept = ivBetaRand[dfData$groupIndex]
  mBetas = matrix(rep(betas, times=length(ivIntercept)), ncol = length(betas), byrow = T)
  mBetas[,1] = ivIntercept
  mBetas = t(mBetas)
  iFitted = diag(mModMatrix %*% mBetas)
  # using log link
  # write the likelihood function
  val = sum(lf(dfData$resp, exp(iFitted)))
  val = val +  dunif(sigmaPop, 0, 1e10, log=T) + #sum(dnorm(iGroupsJitter, 0, sigmaRan, log=T)) +
    dunif(sigmaRan, 0, 1e3, log=T)
  return(val)
}

# set starting values
start = c('sigmaPop'=log(5), 'sigmaRan'=log(2), rep(0, times=5))
fit.old = fit
myloglike(start, lData)

op = optim(start, myloglike, control = list(fnscale = -1, maxit=20000), data=lData, method='SANN', hessian=T)
start = as.numeric(op$par)
names(start)[1:2] = c('sigmaPop', 'sigmaRan')
h = hessian(myloglike, start, data=lData)
# 
# start = as.numeric(op$par)
# names(start)[1:2] = c('sigmaPop', 'sigmaRan')
# start
# op = optim(start, myloglike, control = list(fnscale = -1, maxit=20000), data=lData, hessian=T)

h = -solve(h)
op$var = h

fit = mylaplace(myloglike, start, lData)
fit$mode[3:7]
se = sqrt(abs(diag(fit$var)))
se[3:7]

se = sqrt(abs(diag(fit$var)))
m = fit$par

############################# test with stan
library(rstan)
stanDso = rstan::stan_model(file='mergedDataSet/glmerRegression.stan')

mModMatrix = model.matrix(resp ~ 1+ age + daysRoomTemp+caesarean+sibling3m, data=dfData.re)
lStanData = list(Ntotal=nrow(mModMatrix), Nclusters=length(unique(dfData.re$id)), NgroupMap=dfData.re$groupIndex,
                 Ncol=ncol(mModMatrix), X=mModMatrix, y=dfData.re$resp)

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('betas', 'sigmaPop', 'sigmaRan'))
print(fit.stan)

extract()

######################################################################
# #### section with lattice plots
# dfData.bk = dfData
# # sort on age for plotting
# dfData = dfData[order(dfData$age),]
# 
# # xyplot(Shannon ~ age | intervention+ifelse(rural == 'Yes', 'Rural=Y', 'Rural=N')+ifelse(sterlise == 'Yes', 'Ster=Y', 'Ster=N'),
# #        data=dfData, type=c('g', 'b', 'p'), groups=id, layout=c(4,2))
# # xyplot(Shannon ~ age | sterlise+rural+intervention , data=dfData, type=c('g', 'b', 'p'), groups=id, layout=c(4,2))
# # xyplot(Shannon ~ age | caesarean , data=dfData, type=c('g', 'b', 'p'), groups=id)
# # xyplot(Shannon ~ age | rural , data=dfData, type=c('g', 'b', 'p'), groups=id)
# # xyplot(Shannon ~ age | intervention , data=dfData, type=c('g', 'b', 'p'), groups=id)
# # xyplot(Shannon ~ age | intervention , data=dfData, type=c('g', 'r', 'p'), groups=id)
# # xyplot(Shannon ~ age | pet , data=dfData, type=c('g', 'r'), groups=id)
# # xyplot(Shannon ~ age | ethnicity , data=dfData, type=c('g', 'r', 'p'), groups=id)
# # xyplot(Shannon ~ age | siblings , data=dfData, type=c('g', 'r', 'p'), groups=id)
# 
# xyplot(chaomean ~ age | intervention+ifelse(rural_q3mgen == 'Yes', 'Rural=Y', 'Rural=N')+ifelse(sterlise == 'Yes', 'Ster=Y', 'Ster=N')+
#          ifelse(Caesarean == 'Yes', 'caes=Y', 'caes=N'),
#        data=dfData, type=c('g', 'smooth', 'r'), index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy', layout=c(8,2), 
#        par.strip.text=list(cex=0.7))
# 
# xyplot(chaomean ~ age | interventiongroup,
#        data=dfData, type=c('g', 'r'), index.cond = function(x,y) coef(lm(y ~ x))[1],
#        par.strip.text=list(cex=0.7), groups=id)
# 
# 
# xyplot(chaomean ~ age | id, index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',
#        data=dfData, type=c('g', 'p', 'r'), par.strip.text=list(cex=0.4))#, layout=c(12,5))
# 
# xyplot(chaomean ~ age | id, index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',
#        data=dfData, type=c('g', 'p', 'r'), par.strip.text=list(cex=0.4), groups=interventiongroup, auto.key = list(columns=2))#, layout=c(12,5))
# 
# xyplot(chaomean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], 
#        type=c('p', 'g', 'r'), groups=interventiongroup)
# 
# df = aggregate(dfData$chaomean, by=list(age=dfData$age, intervention=dfData$interventiongroup), mean, na.rm = T)
# xyplot(x ~ age, data=df, auto.key=list(columns=2), xlab='Shannon Diversity', main='Average Shannon Diversity Given Intervention',
#        type=c('o', 'g'), groups=intervention)
# 
# ## make some plots with covariates
# xyplot(chaomean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
#        ylab='Shannon Diversity', main='Shannon Diversity Given Rural',
#        type=c('p', 'g', 'r'), groups=rural_q3mgen)
# 
# xyplot(chaomean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
#        ylab='Shannon Diversity', main='Shannon Diversity Given Cat or Dog',
#        type=c('p', 'g', 'r'), groups=catOrDog3m)
# 
# xyplot(chaomean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
#        ylab='Shannon Diversity', main='Shannon Diversity Given Caesarean',
#        type=c('p', 'g', 'r'), groups=Caesarean)
# 
# xyplot(chaomean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
#        ylab='Shannon Diversity', main='Shannon Diversity Given No Siblings 3m',
#        type=c('p', 'g', 'r'), groups=NumSiblings3m)
# 
# xyplot(chaomean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
#        ylab='Shannon Diversity', main='Shannon Diversity Given Childminder12m',
#        type=c('p', 'g', 'r'), groups=Childminder12m)
# 
# xyplot(chaomean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
#        ylab='Shannon Diversity', main='Shannon Diversity Given Oral AB past month',
#        type=c('p', 'g', 'r'), groups=OralABmonth)
# 
# xyplot(chaomean ~ FirstTooth | age, data=dfData, index.cond = list(as.table=TRUE), auto.key=list(columns=2),
#        ylab='Shannon Diversity', main='Shannon Diversity Given Age',
#        type=c('p', 'g', 'r'))
# 
# xyplot(chaomean ~ FirstTooth, data=dfData, auto.key=list(columns=2),
#        ylab='Shannon Diversity', main='Shannon Diversity Given First Tooth',
#        type=c('p', 'g', 'r'))
# 
# xyplot(chaomean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
#        ylab='Shannon Diversity', main='Shannon Diversity Given sterlise',
#        type=c('p', 'g', 'r'), groups=sterlise)


# densityplot(~ chaomean | interventiongroup, data=dfData, main='Intervention')
# densityplot(~ chaomean | Caesarean, data=dfData, main='Caesarean')
# densityplot(~ chaomean | sterlise, data=dfData, main='Sterlise')
# densityplot(~ chaomean | rural_q3mgen, data=dfData, main='Rural')
# densityplot(~ chaomean | catOrDog3m, data=dfData, main='Pet')
# densityplot(~ chaomean | NumSiblings3m, data=dfData, main='Siblings')
# densityplot(~ chaomean | OralABmonth, data=dfData, main='Oral AB last month')

# xyplot(chaomean ~ age | interventiongroup+sterlise+Caesarean, index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',
#        data=dfData, type=c('g', 'p', 'r'), layout=c(4,2))#, par.strip.text=list(cex=0.4))#, layout=c(12,5))




































