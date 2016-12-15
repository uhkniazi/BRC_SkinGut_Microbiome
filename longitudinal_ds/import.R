# Name: import.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 31/10/2016
# Desc: import the merged data set and extract longitudinal data to make some plots and summaries


dfImport = read.csv(file.choose(), header=T)

str(dfImport)
head(dfImport)
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

oFile = file('longitudinal_ds/Results/data_summary.txt', 'wt')
p1 = c('||', paste(colnames(dfData), ' ||', sep=''))
writeLines(paste(p1, collapse = ''), oFile)

f_paste_factor = function(f){
  s = summary(f)
  p = paste(names(s), '=', s, collapse = ' ')
  return(p)
}

f_paste_other = function(f){
  return(format(mean(f, na.rm = T), digits = 3))
}

## print the information for each time point
f_print_summary = function(df){
  cn = colnames(df)
  cn = cn[4:length(cn)]
  p2 = lapply(cn, function(x) ifelse(is.factor(df[,x]), f_paste_factor(df[,x]), f_paste_other(df[,x])))
  p2 = c(format(length(df$id)), format(unique(df$age)), format(length(df$samples)), do.call(c, p2)) 
  p2 = paste(c('||', paste(p2, '||', sep='', collapse = '')), collapse = '')
  return(p2)
}

writeLines(f_print_summary(dfData[dfData$age == 3,]), oFile)
writeLines(f_print_summary(dfData[dfData$age == 5,]), oFile)
writeLines(f_print_summary(dfData[dfData$age == 12,]), oFile)
close(oFile)


library(lattice)
library(MASS)
library(car)
#library(lme4)
library(lmerTest)
hist(dfData$shannonmean)
fitdistr(dfData$shannonmean, 'normal')
qqp(dfData$shannonmean, 'norm', mean=3.4, sd=0.94, ylab='Shannon diversity')
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
fm00 = update(fm01, shannonmean ~ 1 + (1 + age | id))
fm01 = update(fm01, shannonmean ~ 1 + age + (1 + age | id))
fm02 = update(fm01, shannonmean ~ 1 + age + interventiongroup + (1 + age | id))
fm03 = update(fm01, shannonmean ~ 1 + age*interventiongroup + (1 + age | id))
anova(fm00, fm01, fm02, fm03)

# fit the model with lmer test to get p-values for coefficients
fm = lmerTest::lmer(shannonmean ~ 1 + age*interventiongroup + (1 + age | id), data=dfData)
summary(fm)

## fit a model with each covariate 
f_get.coef.pvalue = function(fit){
  s = summary(fit)
  return(s$coefficients[,'Pr(>|t|)'])
}

f_get.coef.estimate = function(fit){
  s = summary(fit)
  return(s$coefficients[,'Estimate'])
}

cvCov = colnames(dfData)[4:214]

lFits = sapply(cvCov, function(x){
  print(x)
  df = dfData[,c('id', 'age', 'shannonmean', x)]
  cvFormula = paste('shannonmean ~ 1 + ', x,' + (1 + age | id)', sep='')
  tryCatch(expr = lmerTest::lmer(cvFormula, data=df, REML=T), error=function(e) NULL) 
  #return(lmerTest::lmer(cvFormula, data=df, REML=F))
})

## which models did not fit
bNotFit = sapply(lFits, is.null)
which(bNotFit)
lFits = lFits[!bNotFit]

lPVals = lapply(lFits, function(x) tryCatch(expr = f_get.coef.pvalue(x), error=function(e) f_get.coef.estimate(x)))
lCoef = lapply(lFits, function(x) tryCatch(expr = f_get.coef.estimate(x), error=function(e) NULL))

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

