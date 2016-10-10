# Name: univariate_report.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 10/10/2016
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

## original variables i.e. 
## shannon = normal
## chao = gamma
## simpson = beta

## transformed scale
## shannon = no transformation
## chao = log = normal
## simpson = logit = normal

library(lattice)
library(car)
dfData = data.frame(Shannon=dfImport$shannonmean, Frozen=factor(ifelse(dfImport$frozen1yes == 1, 'Frozen.Yes', 'Frozen.No')),
                    DaysOldWhenSampled= dfImport$daysoldwhensamplepassed,
                    Male=factor(ifelse(dfImport$male == 'Yes', 'Male', 'Female')),
                    Caesarean=factor(ifelse(dfImport$Caesarean.born == 'Yes', 'Caesarean.Yes', 'Caesarean.No')),
                    Siblings=factor(ifelse(dfImport$anysibs == 'Yes', 'Siblings.Yes', 'Siblings.No')),
                    Pets=factor(ifelse(dfImport$catordog3m == 1, 'Pet.Yes', 'Pet.No')),
                    Eczema3m=dfImport$ecz3m, Eczema3m.sev=factor(dfImport$eczsev3m),
                    Scorad_cv3m = dfImport$scorad_cv3m, 
                    Eczema12m=factor(ifelse(dfImport$exam12mecz == 'Yes', 'Eczema12m.Yes', 'Eczema12m.No')),
                    Eczema12m.sev=dfImport$exam12meczsev,
                    Scorad_cv12m = dfImport$scorad_cv12m,
                    FoodAllergy=factor(ifelse(dfImport$foodallergy == 'Yes', 'FoodAllergy.Yes', 'FoodAllergy.No')),
                    AntiBiotics.any = factor(ifelse(dfImport$anyrouteantbxlastmonth == 1, 'Antibiotics.Yes', 'Antibiotics.No')),
                    AntiBiotics.oral = factor(ifelse(dfImport$oralantibioticslastweek == 1, 'Antibiotics.Yes', 'Antibiotics.No')),
                    ige009.3m = factor(ifelse(dfImport$ige009at3m == 'Yes', 'ige009at3m.Yes', 'ige009at3m.No')),
                    ige035.12m = factor(ifelse(dfImport$ige035at12m == 'Yes', 'ige035at12m.Yes', 'ige035at12m.No')),
                    anyfoodsensitised_cv12m = factor(ifelse(dfImport$anyfoodsensitised_cv12m == 'Yes', 'anyfoodsensitised_cv12m.Yes',
                                                            'anyfoodsensitised_cv12m.No')),
                    ige035.36m = factor(ifelse(dfImport$ige035at36m == 'Yes', 'ige035at36m.Yes', 'ige035at36m.No'))
                    )


cn = colnames(dfData)
cn = cn[2:length(cn)]

lFm = lapply(cn, function(x) {
  df = dfData[,c('Shannon', x)]
  fm = lm(Shannon ~ ., data=df)
})

names(lFm) = cn

m = sapply(lFm, function(x){
  s = Anova(x)
  s$`Pr(>F)`
})

dfShannon = round((m[1,]), 3)
dfShannon = data.frame(dfShannon)
Sig = ifelse(dfShannon[,1] < 0.05, 'Yes', 'No')
dfShannon$Significant = Sig
colnames(dfShannon) = c('P-Value', 'Significant')

## repeat for each diversity proportion of organisms
## organism abundance scores
cn2 = colnames(dfImport)
i = grep('__', cn2)
cn2 = cn2[i]

lFm02 = lapply(cn2, function(x){
  print(x)
  org = logit(abs(jitter(dfImport[,x])))
  cn = colnames(dfData)
  cn = cn[2:length(cn)]
  
  lFm = lapply(cn, function(x) {
    df = data.frame(org=org, x=dfData[,x])
    fm = lm(org ~ ., data=df)
  })
  
  names(lFm) = cn
  
  m = sapply(lFm, function(x){
    s = Anova(x)
    s$`Pr(>F)`
  })
  
  dfShannon = round((m[1,]), 3)
  dfShannon = data.frame(dfShannon)
  Sig = ifelse(dfShannon[,1] < 0.05, 'Yes', 'No')
  dfShannon$Significant = Sig
  colnames(dfShannon) = c('P-Value', 'Significant')
  return(dfShannon)
})

names(lFm02) = cn2


## chao diversity
dfData$Shannon = dfImport$chaomean

cn = colnames(dfData)
cn = cn[2:length(cn)]
lFm.g = lapply(cn, function(x) {
  df = dfData[,c('Shannon', x)]
  fm = glm(Shannon ~ ., data=df, family = Gamma(link='inverse'))
})

names(lFm.g) = cn

m = sapply(lFm.g, function(x){
  s = Anova(x)
  s$`Pr(>Chisq)`
})

dfChao = round(m, 3)
dfChao = data.frame(dfChao)
Sig = ifelse(dfChao[,1] < 0.05, 'Yes', 'No')
dfChao$Significant = Sig
colnames(dfChao) = c('P-Value', 'Significant')

### box plots and density plots for diversity
dfData = data.frame(Shannon=dfImport$shannonmean, Chao=dfImport$chaomean, Eczema3m=dfImport$ecz3m, Eczema3m.sev=factor(dfImport$eczsev3m))
str(dfData)

boxplot(Shannon ~ Eczema3m, data=dfData)

t = dfData$Shannon[dfData$Eczema3m == 'Yes']
# calculate the mid points for histogram/discrete distribution
h = hist(t, plot=F)
dn = dnorm(h$mids, mean(t), sd(t))
# which distribution can approximate the frequency
hist(t, prob=T, xlab='Shannon Diversity', ylab='', ylim=c(0, max(dn, h$density)), main='Distribution of Shannon Diversity in Disease',
     xlim=c(1.5, 5.5))
# parameterized on the means
lines(h$mids, dn, col='black', type='b')
points(qnorm(0.95, mean(t), sd(t)), 0, pch=20, col='red', cex=2)
legend('topright', legend =c('Normal Density Curve', 'P-Value 0.05 Cutoff'), fill = c('black', 'red'))

t = dfData$Shannon[dfData$Eczema3m == 'No']
# calculate the mid points for histogram/discrete distribution
h = hist(t, plot=F)
dn = dnorm(h$mids, mean(t), sd(t))
# which distribution can approximate the frequency
hist(t, prob=T, xlab='Shannon Diversity', ylab='', ylim=c(0, max(dn, h$density)), main='Distribution of Shannon Diversity in Healthy')
# parameterized on the means
lines(h$mids, dn, col='black', type='b')
points(qnorm(0.95, mean(t), sd(t)), 0, pch=20, col='red', cex=2)
legend('topright', legend =c('Normal Density Curve', 'P-Value 0.05 Cutoff'), fill = c('black', 'red'))

boxplot(Shannon ~ Eczema3m, data=dfData)
abline(h = qnorm(0.95, mean(t), sd(t)), col='red', lwd=2)

iCut.pt = qnorm(0.95, mean(t), sd(t))
i = which(dfData$Shannon > iCut.pt)
fSubGroups = rep('Low', times=length(dfData$Shannon))
fSubGroups[i] = 'High'
fSubGroups = factor(fSubGroups, levels=c('Low', 'High'))
dfData$fSubGroups = fSubGroups
# drop one value i.e outlier in the data
# i = which(dfData$Eczema3m == 'Yes' & dfData$fSubGroups == 'High')
# dfData = dfData[-i,]
f = paste(dfData$Eczema3m, dfData$fSubGroups, sep=' ')
f = gsub('Yes (Low|High)', 'Disease', f)
f = gsub('No (\\w+)', 'Healthy \\1', f)
dfData$Eczema3m.Diversity = factor(f)

boxplot(Shannon ~ Eczema3m.Diversity, data=dfData, ylab='Shannon Diversity', main='High diversity has a protective effect')
abline(h = qnorm(0.95, mean(t), sd(t)), col='red', lwd=2)




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

