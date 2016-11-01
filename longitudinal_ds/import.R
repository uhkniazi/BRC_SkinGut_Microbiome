# Name: import.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 31/10/2016
# Desc: import the longitudinal dataset and cleanup


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
dfData = data.frame(Shannon=dfImport$shannonmean, Frozen=factor(ifelse(dfImport$frozen1yes == 1, 'Frozen.Yes', 'Frozen.No')), 
                    Male=factor(ifelse(dfImport$male == 'Yes', 'Male', 'Female')),
                    Caesarean=factor(ifelse(dfImport$Caesarean.born == 'Yes', 'Caesarean.Yes', 'Caesarean.No')),
                    Siblings=factor(ifelse(dfImport$anysibs == 'Yes', 'Siblings.Yes', 'Siblings.No')),
                    Pets=factor(ifelse(dfImport$catordog3m == 1, 'Pet.Yes', 'Pet.No')),
                    Eczema3m=dfImport$ecz3m, Eczema3m.sev=factor(dfImport$eczsev3m),
                    Eczema12m=factor(ifelse(dfImport$exam12mecz == 'Yes', 'Eczema12m.Yes', 'Eczema12m.No')),
                    Eczema12m.sev=dfImport$exam12meczsev,
                    FoodAllergy=factor(ifelse(dfImport$foodallergy == 'Yes', 'FoodAllergy.Yes', 'FoodAllergy.No')),
                    AntiBiotics = factor(ifelse(dfImport$anyrouteantbxlastmonth == 1, 'Antibiotics.Yes', 'Antibiotics.No')),
                    ige009.3m = factor(ifelse(dfImport$ige009at3m == 'Yes', 'ige009at3m.Yes', 'ige009at3m.No')),
                    ige035.12m = factor(ifelse(dfImport$ige035at12m == 'Yes', 'ige035at12m.Yes', 'ige035at12m.No')),
                    anyfoodsensitised_cv12m = factor(ifelse(dfImport$anyfoodsensitised_cv12m == 'Yes', 'anyfoodsensitised_cv12m.Yes',
                                                            'anyfoodsensitised_cv12m.No')),
                    ige035.36m = factor(ifelse(dfImport$ige035at36m == 'Yes', 'ige035at36m.Yes', 'ige035at36m.No'))
                    )

densityplot(~ Shannon | Eczema3m.sev, data=dfData, groups=Eczema3m, auto.key=T)

## how does the density compare with eczema at 3 months and covariates
densityplot(~ Shannon | Frozen, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Frozen')
densityplot(~ Shannon | Male, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Male')
densityplot(~ Shannon | Caesarean, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Caesarean')
densityplot(~ Shannon | Siblings, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Siblings')
densityplot(~ Shannon | Pets, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Pets')
densityplot(~ Shannon | Eczema3m.sev, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Eczema 3m Severity')
densityplot(~ Shannon | Eczema12m, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Eczema 12m')
densityplot(~ Shannon | FoodAllergy, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given FoodAllergy')
densityplot(~ Shannon | AntiBiotics, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given AntiBiotics')
densityplot(~ Shannon | ige009.3m, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given ige009.3m')
densityplot(~ Shannon | ige035.12m, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given ige035.12m')
densityplot(~ Shannon | anyfoodsensitised_cv12m, data=dfData, auto.key=T, groups=Eczema3m, 
            main='Eczema at 3 months given anyfoodsensitised_cv12m')
densityplot(~ Shannon | ige035.36m, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given ige035.36m')

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
dfData$Eczema3m.Sub = factor(paste(dfData$Eczema3m, fSubGroups, sep='.'))

densityplot(~ Shannon, data=dfData, groups=Eczema3m.Sub, auto.key=list(columns=4), 
            main='Distribution of Diversity given Eczema3m (No/Yes) and Group (High/Low)')
dfData.bk = dfData
dfData = subset(dfData, dfData$Eczema3m.Sub != 'Yes.High')

## how does the density compare with eczema at 3 months and covariates
densityplot(~ Shannon | Frozen, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given Frozen')
densityplot(~ Shannon | Male, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given Male')
densityplot(~ Shannon | Caesarean, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given Caesarean')
densityplot(~ Shannon | Siblings, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given Siblings')
densityplot(~ Shannon | Pets, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given Pets')
densityplot(~ Shannon | Eczema3m.sev, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given Eczema 3m Severity')
densityplot(~ Shannon | Eczema12m, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given Eczema 12m')
densityplot(~ Shannon | FoodAllergy, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given FoodAllergy')
densityplot(~ Shannon | AntiBiotics, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given AntiBiotics')
densityplot(~ Shannon | ige009.3m, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given ige009.3m')
densityplot(~ Shannon | ige035.12m, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given ige035.12m')
densityplot(~ Shannon | anyfoodsensitised_cv12m, data=dfData, auto.key=T, groups=Eczema3m.Sub, 
            main='Eczema at 3 months given anyfoodsensitised_cv12m')
densityplot(~ Shannon | ige035.36m, data=dfData, auto.key=T, groups=Eczema3m.Sub, main='Eczema at 3 months given ige035.36m')

## organism abundance scores
cn = colnames(dfImport)
i = grep('__', cn)
cn = cn[i]

sapply(cn, function(x){
  org = logit(dfImport[,x])
  p = xyplot(dfData$Shannon ~ org | dfData$Eczema3m, groups = dfData$fSubGroups, main=x, auto.key=T)
  print(p)
})

dfData = dfData.bk
dfStack = dfImport[,cn]
str(dfStack)
df = stack(dfStack)
str(df)
df$Ezcema = dfData$Eczema3m
df$Group = dfData$Eczema3m.Sub
df = subset(df, df$Group != 'Yes.High')
densityplot( ~ logit(values) | ind, data=df, groups=Group, type='n', auto.key=list(columns=4), xlab='logit Diversity')

## pairs plots
col.p = c('black', 'red')
col = col.p[as.numeric(dfData$fSubGroups)]
pairs(dfDiversity, pch=20, col=col)



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

