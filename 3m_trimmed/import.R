# Name: import.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 3/10/2016
# Desc: import the 3 month old data with covariates and save data object after cleanup


### functions
getalphabeta.poisson = function(lambda){
  m = mean(lambda)
  v = var(lambda)
  alpha = (m^2)/v
  beta = alpha/m
  return(c(alpha=alpha, beta=beta))
}

getalphabeta.beta = function(m, v){
  al.be = (m * (1-m) / v) - 1
  al = al.be * m
  be = al.be * (1-m)
  return(c(alpha=al, beta=be))
}

logit = function(p) log(p/(1-p))
logit.inv = function(p) {exp(p)/(exp(p)+1) }
###

dfImport = read.csv(file.choose(), header=T)

str(dfImport)
head(dfImport)

# data frame for diversity indices
dfDiversity = data.frame(cbind(dfImport$shannonmean, dfImport$simpsonmean, dfImport$chaomean))#, dfImport$pdwholetreemean, dfImport$obsotumean))
colnames(dfDiversity) = c('Shannon', 'Simpson', 'Chao')#, 'TreeMean', 'ObsotuMean')
str(dfDiversity)

# check the corrleations
pairs(dfDiversity)

# distributions for the data
library(car)
library(MASS)

f_normqq = function(x, t){
  hist(x, main=paste(t))
  d = fitdistr(x, 'normal')
  print(d)
  qqPlot(x, 'norm', main = paste(t)) 
}

f_gammaqq = function(x, t){
  d = getalphabeta.poisson(x)
  print(d)
  qqPlot(x, 'gamma', shape=d['alpha'], rate=d['beta'], main = paste(t)) 
}

f_betaqq = function(x, t){
  d = getalphabeta.beta(mean(x), var(x))
  print(d)
  if (d['alpha'] < 0) return(0)
  qqPlot(x, 'beta', shape1=d['alpha'], shape2=d['beta'], main = paste(t)) 
}

cn = colnames(dfDiversity)
temp = sapply(cn, function(x){
  y = dfDiversity[,x]
  f_normqq(y, x)
})

temp = sapply(cn, function(x){
  y = dfDiversity[,x]
  f_gammaqq(y, x)
})

temp = sapply(cn, function(x){
  y = dfDiversity[,x]
  f_betaqq(y, x)
})

## transform some variables and try again
dfDiversity.tr = dfDiversity
dfDiversity.tr$Simpson = logit(dfDiversity.tr$Simpson)
dfDiversity.tr$Chao = log(dfDiversity.tr$Chao)
#dfDiversity.tr$TreeMean = sqrt(dfDiversity.tr$TreeMean)
#dfDiversity.tr$ObsotuMean = sqrt(dfDiversity.tr$ObsotuMean)
pairs(dfDiversity.tr)

cn = colnames(dfDiversity.tr)
temp = sapply(cn, function(x){
  y = dfDiversity.tr[,x]
  f_normqq(y, x)
})

temp = sapply(cn, function(x){
  y = dfDiversity.tr[,x]
  f_gammaqq(y, x)
})

temp = sapply(cn, function(x){
  y = dfDiversity.tr[,x]
  f_betaqq(y, x)
})

## original variables i.e. 
## shannon = normal
## chao = gamma
## simpson = beta

## transformed scale
## shannon = no transformation
## chao = log = normal
## simpson = logit = normal

library(lattice)
dfData = data.frame(Shannon=dfImport$shannonmean, Frozen=factor(ifelse(dfImport$frozen1yes == 1, 'Yes', 'No')), Male=dfImport$male,
                                                                Caesarean=dfImport$Caesarean.born, Siblings=dfImport$anysibs,
                                                                Pets=factor(dfImport$catordog3m), Eczema3m=dfImport$ecz3m,
                                                                Eczema3m.sev=factor(dfImport$eczsev3m))

densityplot(~ Shannon | Eczema3m.sev, data=dfData, groups=Eczema3m, auto.key=T)


