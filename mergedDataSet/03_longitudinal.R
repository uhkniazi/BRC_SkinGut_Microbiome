# Name: 03_longitudinal.R
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
library(lmerTest)

# check data distribution
hist(dfData$shannonmean)
fitdistr(dfData$shannonmean, 'normal')
qqPlot(dfData$shannonmean, 'norm',ylab='Shannon diversity')

### trying model fits with covariates
# model with just random intercept
fm01 = lmer(shannonmean ~ 1 + (1 | id), data=dfData, REML=F)
summary(fm01)
fm01.cor = update(fm01, shannonmean ~ 1 + (1 + age | id))
summary(fm01.cor)
# compare the 2 models
anova(fm01, fm01.cor)

# try uncorrelated random effects
fm01.uncor = update(fm01, shannonmean ~ 1 + (1 | id) + (0 + age | id))
summary(fm01.uncor)
# compare 3 models
anova(fm01.uncor, fm01.cor, fm01 )

## choose the model with correlated intercept and slope
fm01 = fm01.cor

## update this model with one covariate after adding time
fm02 = update(fm01, shannonmean ~ 1 + age + (1 + age | id))
summary(fm02)
anova(fm01, fm02)

# test the four models in a sequence
fm00 = update(fm01, shannonmean ~ 1 + (1 + age | id))
fm01 = update(fm01, shannonmean ~ 1 + age + (1 + age | id))
fm02 = update(fm01, shannonmean ~ 1 + age + intervention + (1 + age | id))
fm03 = update(fm01, shannonmean ~ 1 + age*intervention + (1 + age | id))
anova(fm00, fm01, fm02, fm03)

# fit the model with lmer test to get p-values for coefficients
fm = lmerTest::lmer(shannonmean ~ 1 + age + intervention + (1 + age | id), data=dfData)
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
cvCols = cvCols[!(cvCols %in% c('samples', "chaomean", "pdwholetreemean",
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
n = c('eczemaSev3m', 'scorad3m', 'skinPrickAero12m', 'skinPrickAero36m',
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
#   df = dfData[,c('id', 'age', 'shannonmean', x)]
#   cvFormula = paste('shannonmean ~ 1 + ', x,' + (1 + age | id)', sep='')
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
  df = dfData[,c('id', 'age', 'shannonmean', x)]
  cvFormula = paste('shannonmean ~ 1 + ', x,' + (1 + age | id)', sep='')
  lFits[[i]] = Anova(lmerTest::lmer(cvFormula, data=df, REML=T))
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
write.csv(dfUnivariate, file='mergedDataSet/Temp/uni.csv')


## create a file with the results
oFile = file('longitudinal_ds/Results/univariate_summary.txt', 'wt')
p1 = c('||', 'Covariate || Coefficient || P-Value', ' ||', sep='')
writeLines(paste(p1, collapse = ''), oFile)

for(i in 1:length(lPVals)){
  pv = lPVals[[i]][-1]
  co = lCoef[[i]][-1]
  # check length and print as some may have multiple levels
  for (l in 1:length(pv)){
    p2 = paste('||', names(lPVals)[i], ' ', names(pv)[l], '||', round(co[l], 3), '||', round(pv[l], 3), '||', sep='')
    writeLines(p2, oFile)
  }
}

close(oFile)


#### section with lattice plots
dfData.bk = dfData
# sort on age for plotting
dfData = dfData[order(dfData$age),]

# xyplot(Shannon ~ age | intervention+ifelse(rural == 'Yes', 'Rural=Y', 'Rural=N')+ifelse(sterlise == 'Yes', 'Ster=Y', 'Ster=N'),
#        data=dfData, type=c('g', 'b', 'p'), groups=id, layout=c(4,2))
# xyplot(Shannon ~ age | sterlise+rural+intervention , data=dfData, type=c('g', 'b', 'p'), groups=id, layout=c(4,2))
# xyplot(Shannon ~ age | caesarean , data=dfData, type=c('g', 'b', 'p'), groups=id)
# xyplot(Shannon ~ age | rural , data=dfData, type=c('g', 'b', 'p'), groups=id)
# xyplot(Shannon ~ age | intervention , data=dfData, type=c('g', 'b', 'p'), groups=id)
# xyplot(Shannon ~ age | intervention , data=dfData, type=c('g', 'r', 'p'), groups=id)
# xyplot(Shannon ~ age | pet , data=dfData, type=c('g', 'r'), groups=id)
# xyplot(Shannon ~ age | ethnicity , data=dfData, type=c('g', 'r', 'p'), groups=id)
# xyplot(Shannon ~ age | siblings , data=dfData, type=c('g', 'r', 'p'), groups=id)

xyplot(shannonmean ~ age | interventiongroup+ifelse(rural_q3mgen == 'Yes', 'Rural=Y', 'Rural=N')+ifelse(sterlise == 'Yes', 'Ster=Y', 'Ster=N')+
         ifelse(Caesarean == 'Yes', 'caes=Y', 'caes=N'),
       data=dfData, type=c('g', 'smooth', 'r'), index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy', layout=c(8,2), 
       par.strip.text=list(cex=0.7))

xyplot(shannonmean ~ age | interventiongroup,
       data=dfData, type=c('g', 'r'), index.cond = function(x,y) coef(lm(y ~ x))[1],
       par.strip.text=list(cex=0.7), groups=id)


xyplot(shannonmean ~ age | id, index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',
       data=dfData, type=c('g', 'p', 'r'), par.strip.text=list(cex=0.4))#, layout=c(12,5))

xyplot(shannonmean ~ age | id, index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',
       data=dfData, type=c('g', 'p', 'r'), par.strip.text=list(cex=0.4), groups=interventiongroup, auto.key = list(columns=2))#, layout=c(12,5))

xyplot(shannonmean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], 
       type=c('p', 'g', 'r'), groups=interventiongroup)

df = aggregate(dfData$shannonmean, by=list(age=dfData$age, intervention=dfData$interventiongroup), mean, na.rm = T)
xyplot(x ~ age, data=df, auto.key=list(columns=2), xlab='Shannon Diversity', main='Average Shannon Diversity Given Intervention',
       type=c('o', 'g'), groups=intervention)

## make some plots with covariates
xyplot(shannonmean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
       ylab='Shannon Diversity', main='Shannon Diversity Given Rural',
       type=c('p', 'g', 'r'), groups=rural_q3mgen)

xyplot(shannonmean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
       ylab='Shannon Diversity', main='Shannon Diversity Given Cat or Dog',
       type=c('p', 'g', 'r'), groups=catOrDog3m)

xyplot(shannonmean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
       ylab='Shannon Diversity', main='Shannon Diversity Given Caesarean',
       type=c('p', 'g', 'r'), groups=Caesarean)

xyplot(shannonmean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
       ylab='Shannon Diversity', main='Shannon Diversity Given No Siblings 3m',
       type=c('p', 'g', 'r'), groups=NumSiblings3m)

xyplot(shannonmean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
       ylab='Shannon Diversity', main='Shannon Diversity Given Childminder12m',
       type=c('p', 'g', 'r'), groups=Childminder12m)

xyplot(shannonmean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
       ylab='Shannon Diversity', main='Shannon Diversity Given Oral AB past month',
       type=c('p', 'g', 'r'), groups=OralABmonth)

xyplot(shannonmean ~ FirstTooth | age, data=dfData, index.cond = list(as.table=TRUE), auto.key=list(columns=2),
       ylab='Shannon Diversity', main='Shannon Diversity Given Age',
       type=c('p', 'g', 'r'))

xyplot(shannonmean ~ FirstTooth, data=dfData, auto.key=list(columns=2),
       ylab='Shannon Diversity', main='Shannon Diversity Given First Tooth',
       type=c('p', 'g', 'r'))

xyplot(shannonmean ~ age, data=dfData, index.cond = function(x,y) coef(lm(y ~ x))[1], auto.key=list(columns=2),
       ylab='Shannon Diversity', main='Shannon Diversity Given sterlise',
       type=c('p', 'g', 'r'), groups=sterlise)


densityplot(~ shannonmean | interventiongroup, data=dfData, main='Intervention')
densityplot(~ shannonmean | Caesarean, data=dfData, main='Caesarean')
densityplot(~ shannonmean | sterlise, data=dfData, main='Sterlise')
densityplot(~ shannonmean | rural_q3mgen, data=dfData, main='Rural')
densityplot(~ shannonmean | catOrDog3m, data=dfData, main='Pet')
densityplot(~ shannonmean | NumSiblings3m, data=dfData, main='Siblings')
densityplot(~ shannonmean | OralABmonth, data=dfData, main='Oral AB last month')

# xyplot(shannonmean ~ age | interventiongroup+sterlise+Caesarean, index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',
#        data=dfData, type=c('g', 'p', 'r'), layout=c(4,2))#, par.strip.text=list(cex=0.4))#, layout=c(12,5))

