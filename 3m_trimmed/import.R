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
                    Caesarean=dfImport$Caesarean.born, Siblings=dfImport$anysibs, Pets=factor(dfImport$catordog3m),
                    Eczema3m=dfImport$ecz3m, Eczema3m.sev=factor(dfImport$eczsev3m),
                    Eczema12m=dfImport$exam12mecz,
                    Eczema12m.sev=dfImport$exam12meczsev,
                    FoodAllergy=dfImport$foodallergy,
                    AntiBiotics = factor(dfImport$anyrouteantbxlastmonth),
                    ige009.3m = dfImport$ige009at3m,
                    ige035.12m = dfImport$ige035at12m,
                    anyfoodsensitised_cv12m = dfImport$anyfoodsensitised_cv12m,
                    ige035.36m = dfImport$ige035at36m
                    )

densityplot(~ Shannon | Eczema3m.sev, data=dfData, groups=Eczema3m, auto.key=T)

## how does the density compare with eczema at 3 months and covariates
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=Frozen, main='Eczema at 3 months given Frozen')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=Male, main='Eczema at 3 months given Male')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=Caesarean, main='Eczema at 3 months given Caesarean')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=Siblings, main='Eczema at 3 months given Siblings')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=Pets, main='Eczema at 3 months given Pets')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=Eczema3m.sev, main='Eczema at 3 months given Eczema 3m Severity')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=Eczema12m, main='Eczema at 3 months given Eczema 12m')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=FoodAllergy, main='Eczema at 3 months given FoodAllergy')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=AntiBiotics, main='Eczema at 3 months given AntiBiotics')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=ige009.3m, main='Eczema at 3 months given ige009.3m')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=ige035.12m, main='Eczema at 3 months given ige035.12m')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=anyfoodsensitised_cv12m, 
            main='Eczema at 3 months given anyfoodsensitised_cv12m')
densityplot(~ Shannon | Eczema3m, data=dfData, auto.key=T, groups=ige035.36m, main='Eczema at 3 months given ige035.36m')

## organism abundance scores
cn = colnames(dfImport)
i = grep('__', cn)
cn = cn[i]

sapply(cn, function(x){
  org = logit(dfImport[,x])
  p = xyplot(dfData$Shannon ~ org | dfData$Eczema3m, main=x)
  print(p)
})

sapply(cn, function(x){
  org = logit(dfImport[,x])
  p = densityplot( ~ org, main=x, groups = dfData$Eczema3m, auto.key=T, xlab = 'Diversity')
  print(p)
})

dfStack = dfImport[,cn]
str(dfStack)
df = stack(dfStack)
str(df)
df$Eczema = dfData$Eczema3m
densityplot( ~ logit(values) | ind, data=df, groups=Eczema, type='n', auto.key=T)

# xyplot(ifelse(Eczema3m == "Yes", 1, 0) ~ Shannon | some, data = dfData,
#         groups = FoodAllergy, type = c("g", "smooth"),
#         auto.key = list(space = "top", points = FALSE,
#                         lines = TRUE, columns = 4),
#         ylab = "Eczema at 3 months", xlab = "Shannon Diversity")

xyplot(ifelse(Eczema3m == "Yes", 1, 0) ~ Shannon, data = dfData,
       type = c("p"), groups=FoodAllergy,
       auto.key = T,
       ylab = "Eczema at 3 months", xlab = "Shannon Diversity")




# split the Samples into 2 groups based on the last quantile of Shannon or Chao diversity distribution
# use the shannon index
iCut.pt = qnorm(0.95, mean(dfImport$shannonmean), sd(dfImport$shannonmean))
i = which(dfImport$shannonmean > iCut.pt)
fSubGroups = rep('Low', times=length(dfImport$shannonmean))
fSubGroups[i] = 'High'
fSubGroups = factor(fSubGroups, levels=c('Low', 'High'))

dfData$fSubGroups = fSubGroups

densityplot(~ Shannon | Eczema3m, data=dfData, groups=fSubGroups, auto.key=T)

densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=Frozen, main='Eczema at 3 months given Frozen')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=Male, main='Eczema at 3 months given Male')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=Caesarean, main='Eczema at 3 months given Caesarean')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=Siblings, main='Eczema at 3 months given Siblings')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=Pets, main='Eczema at 3 months given Pets')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=Eczema3m.sev, main='Eczema at 3 months given Eczema 3m Severity')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=Eczema12m, main='Eczema at 3 months given Eczema 12m')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=FoodAllergy, main='Eczema at 3 months given FoodAllergy')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=AntiBiotics, main='Eczema at 3 months given AntiBiotics')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=ige009.3m, main='Eczema at 3 months given ige009.3m')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=ige035.12m, main='Eczema at 3 months given ige035.12m')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=anyfoodsensitised_cv12m, 
            main='Eczema at 3 months given anyfoodsensitised_cv12m')
densityplot(~ Shannon | Eczema3m+fSubGroups, data=dfData, auto.key=T, groups=ige035.36m, main='Eczema at 3 months given ige035.36m')

## organism abundance scores
cn = colnames(dfImport)
i = grep('__', cn)
cn = cn[i]

sapply(cn, function(x){
  org = logit(dfImport[,x])
  p = xyplot(dfData$Shannon ~ org | dfData$Eczema3m, groups = dfData$fSubGroups, main=x, auto.key=T)
  print(p)
})

dfStack = dfImport[,cn]
str(dfStack)
df = stack(dfStack)
str(df)
df$Ezcema = dfData$Eczema3m
df$Group = dfData$fSubGroups
densityplot( ~ logit(values) | ind, data=df, groups=Group, type='n', auto.key=T)


## calculate prior and posterior probabilities
## P(Ec = Y | Group=High), P(Ec = N | Group=High)
t = xtabs( ~ Eczema3m + fSubGroups, data=dfData)
# confusion matrix
mConf = t
mConf = mConf/sum(rowSums(mConf))
pec = rowSums(mConf)
pgr = colSums(mConf)
mProbEc.Given.Group = sweep(mConf, 2, pgr, '/')
t(outer(pgr, pec))
mConf
mProbEc.Given.Group
# prior
pec
# post
mProbEc.Given.Group[,'High']
mProbEc.Given.Group[,'Low']

