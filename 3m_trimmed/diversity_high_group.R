# Name: diversity_high_group.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 01/11/2016
# Desc: multinomial regression for 3 month diversity split into groups


### functions
# getalphabeta.poisson = function(lambda){
#   m = mean(lambda)
#   v = var(lambda)
#   alpha = (m^2)/v
#   beta = alpha/m
#   return(c(alpha=alpha, beta=beta))
# }
# 
# getalphabeta.beta = function(m, v){
#   al.be = (m * (1-m) / v) - 1
#   al = al.be * m
#   be = al.be * (1-m)
#   return(c(alpha=al, beta=be))
# }
# 
# logit = function(p) log(p/(1-p))
# logit.inv = function(p) {exp(p)/(exp(p)+1) }
###

dfImport = read.csv(file.choose(), header=T)

str(dfImport)
head(dfImport)

library(lattice)
library(car)



### density plots for shannon diversity in eczema3m yes and no groups
dfData = data.frame(Shannon=dfImport$shannonmean, Frozen=factor(ifelse(dfImport$frozen1yes == 1, 'Yes', 'No')),
                    DaysOldWhenSampled= dfImport$daysoldwhensamplepassed,
                    Male=factor(ifelse(dfImport$male == 'Yes', 'Male', 'Female')),
                    Caesarean=factor(ifelse(dfImport$Caesarean.born == 'Yes', 'Yes', 'No')),
                    Siblings=factor(ifelse(dfImport$anysibs == 'Yes', 'Yes', 'No')),
                    Pets=factor(ifelse(dfImport$catordog3m == 1, 'Yes', 'No')),
                    Eczema3m=dfImport$ecz3m, Eczema3m.sev=factor(dfImport$eczsev3m),
                    Scorad_cv3m = dfImport$scorad_cv3m, 
                    Eczema12m=factor(ifelse(dfImport$exam12mecz == 'Yes', 'Yes', 'No')),
                    Eczema12m.sev=dfImport$exam12meczsev,
                    Scorad_cv12m = dfImport$scorad_cv12m,
                    FoodAllergy=factor(ifelse(dfImport$foodallergy == 'Yes', 'Yes', 'No')),
                    AntiBiotics.any = factor(ifelse(dfImport$anyrouteantbxlastmonth == 1, 'Yes', 'No')),
                    AntiBiotics.oral = factor(ifelse(dfImport$oralantibioticslastweek == 1, 'Yes', 'No')),
                    ige009.3m = factor(ifelse(dfImport$ige009at3m == 'Yes', 'Yes', 'No')),
                    ige035.12m = factor(ifelse(dfImport$ige035at12m == 'Yes', 'Yes', 'No')),
                    anyfoodsensitised_cv12m = factor(ifelse(dfImport$anyfoodsensitised_cv12m == 'Yes', 'Yes',
                                                            'No')),
                    ige035.36m = factor(ifelse(dfImport$ige035at36m == 'Yes', 'Yes', 'No'))
)


# dfData = data.frame(Shannon=dfImport$shannonmean, Eczema3m=dfImport$ecz3m, 
#                     Eczema12m=dfImport$exam12mecz)
str(dfData)

t = dfData$Shannon[dfData$Eczema3m == 'Yes']
# calculate the mid points for histogram/discrete distribution
h = hist(t, plot=F)
dn = dnorm(h$mids, mean(t), sd(t))
# which distribution can approximate the frequency
hist(t, prob=T, xlab='Shannon Diversity', ylab='', ylim=c(0, max(dn, h$density)), main='Distribution of Shannon Diversity in Disease',
     xlim=c(1.5, 5.5))
# parameterized on the means
lines(h$mids, dn, col='black', type='b')
#points(qnorm(0.95, mean(t), sd(t)), 0, pch=20, col='red', cex=2)
legend('topright', legend =c('Normal Density Curve'), fill = c('black'))


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

## diversity higher than this cutoff is Diversity High group
iCut.pt = qnorm(0.95, mean(t), sd(t))
i = which(dfData$Shannon > iCut.pt)
fSubGroups = rep('Low', times=length(dfData$Shannon))
fSubGroups[i] = 'High'
fSubGroups = factor(fSubGroups, levels=c('Low', 'High'))
dfData$DiversityClass = fSubGroups
# drop one value i.e outlier in the data
# i = which(dfData$Eczema3m == 'Yes' & dfData$fSubGroups == 'High')
# dfData = dfData[-i,]
f = paste(dfData$Eczema3m, dfData$DiversityClass, sep=' ')
# f = gsub('Yes (Low|High)', 'Disease', f)
# f = gsub('No (\\w+)', 'Healthy \\1', f)
# dfData$Condition = factor(f, levels=c('Disease', 'Healthy Low', 'Healthy High'))
table(f)
dfData$Condition = factor(f, levels=c('No High', 'No Low', 'Yes High', 'Yes Low'))

str(dfData)

boxplot(Shannon ~ Condition, data=dfData, ylab='Shannon Diversity 3M', main='High diversity has a protective effect')
abline(h = qnorm(0.95, mean(t), sd(t)), col='red', lwd=2)

### perform logistic regression where response is DiversityClass
### find the associated variables
## download the ccrossvalidation class to select variables of importance
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

# check what is predictive of eczema at DiversityClass
as.data.frame(colnames(dfData))
dfData.sub = dfData
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
dim(dfData.sub)
fGroups = dfData.sub$DiversityClass
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-21]
str(dfData.sub)

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups, boot.num = 1000)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
dfRF
rn = rownames(dfRF)
data.frame(rn)
## look at the top variables and perform subset selection
i = which(colnames(dfData.sub) %in% rn[c(5:19)])
dfData.sub = dfData.sub[,i]
str(dfData.sub)

oVar.s = CVariableSelection.ReduceModel(dfData.sub, fGroups, boot.num = 100)
plot.var.selection(oVar.s)

sapply(seq_along(1:ncol(dfData.sub)), function(x) CVariableSelection.ReduceModel.getMinModel(oVar.s, x))

cn = CVariableSelection.ReduceModel.getMinModel(oVar.s, 4)

## perform logistic regression analysis on these variables
# create formula
cvFormula = paste0('DiversityClass ~ 1 + ', paste(cn, collapse = ' + '))
cvFormula
fm01 = glm(cvFormula, family=binomial(), data=dfData, control = list(maxit=100))
library(car)
Anova(fm01)
anova(fm01, test='Chisq')
summary(fm01)



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

