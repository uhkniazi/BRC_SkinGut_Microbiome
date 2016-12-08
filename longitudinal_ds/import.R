# Name: import.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 31/10/2016
# Desc: import the longitudinal dataset and cleanup


dfImport = read.csv(file.choose(), header=T)

str(dfImport)
head(dfImport)

# setup data
dfData = dfImport
dfData$age = as.numeric(gsub('m', '', dfData$age))
rownames(dfData) = dfImport$samples
dfData$samples = factor(dfData$samples)
dfData$OralABweek = factor(dfData$OralABweek)
dfData$OralABmonth = factor(dfData$OralABmonth)

str(dfData)

# merge the 3 columns into one
df = data.frame(m3=dfImport$currentlysterilise3months, 
                m5=dfImport$currentlysterilise5months, 
                m12=dfImport$currentlysterilise12months, 
                age=dfData$age, 
                id=dfData$id)

fm = function(df){
  id = as.character(df$id)
  v = rep(NA, length.out=length(id))
  names(v) = id
  for (i in 1:length(id)){
    a = df$age[i]
    am = paste0('m',a)
    v[i] = ifelse(df[i,am] == '', NA, as.character(df[i,am]))
  }
  return(v)
}

dfData$sterlise = fm(df)
dfData$sterlise = factor(dfData$sterlise)

# drop the sterlize columns as they have been merged
as.data.frame(colnames(dfData))
dfData = dfData[,-c(16:18)]
str(dfData)
summary(dfData)

## import the 3 month new dataset and merge
dfData.3m = read.csv(file.choose(), header=T)
## match the data ids
i = match(dfData$id, dfData.3m$id)

rt = dfData.3m$dayssampleatroomtemperature[i]
dfData$SampleAtRoomTemp[dfData$age==3] = rt[dfData$age==3]

## other covariates
dfData$eczema3m = dfData.3m$exam3mecz[i]
dfData$eczema12m = dfData.3m$exam12mecz[i]
dfData$foodallergy = dfData.3m$foodallergy[i]
dfData$peanutallergy = dfData.3m$peanutallergy[i]
dfData$eggallergy = dfData.3m$eggallergy[i]
dfData$ige009at3m = dfData.3m$ige009at3m[i]
dfData$ige035at12m = dfData.3m$ige035at12m[i]
dfData$ige035at36m = dfData.3m$ige035at36m[i]
dfData$Frozen3m = factor(ifelse(dfData.3m$frozen1yes[i] == 1, 'Yes', 'No'))
wt = dfData.3m$weight_cv3m[i]
dfData$Weight = NA
dfData$Weight[dfData$age == 3] = wt[dfData$age == 3]
wt = dfData.3m$weight_cv12m[i]
dfData$Weight[dfData$age == 12] = wt[dfData$age == 12]
dfData$Weight36m = dfData.3m$weight_cv36m[i]
dfData$BMI36m = dfData.3m$bmi_cv36m[i]
dfData$Sex = dfData.3m$male_cv3m[i]
dfData$homeMould3m = dfData.3m$homemould_q3mgen[i]
dfData$diarrhea3m = dfData.3m$anydiarrdays3m[i]
dfData$respInfections3m = dfData.3m$numrespinf3m[i]


## import the second 3month dataset with additional covariates
dfData.3m = read.csv(file.choose(), header=T, sep='\t')
head(dfData.3m)
str(dfData.3m)
## match the data ids
i = match(dfData$id, dfData.3m$id)

dfData$MotherAge = dfData.3m$MotherAge[i]
dfData$diarrhea12m = dfData.3m$diarrhea12m[i]
dfData$ParentSmoke = dfData.3m$ParentSmoke[i]
dfData$ParentEczema = dfData.3m$ParentEczema[i]
dfData$MotherEducation = dfData.3m$MotherEducation[i]
dfData$FatherEducation = dfData.3m$FatherEducation[i]
# create a longitudinal version of diarrhea variable
dia = rep(NA, length=nrow(dfData))
dia[dfData$age == 3] = ifelse(dfData$diarrhea3m[dfData$age == 3] == 'Yes', 'Yes', 'No')
dia[dfData$age == 12] = ifelse(dfData$diarrhea12m[dfData$age == 12] == 'Yes', 'Yes', 'No')
dfData$diarrhea = factor(dia)

## save this data 
write.csv(dfData, file='longitudinal_ds/Data_external/longitudinal_dataset_manual_UN_2.csv')

oFile = file('longitudinal_ds/Results/data_summary.txt', 'wt')
p1 = c('||', paste(colnames(dfData), ' ||', sep=''))
writeLines(paste(p1, collapse = ''), oFile)

f_paste_factor = function(f){
  s = summary(f)
  p = paste(names(s), '=', s, collapse = ' ')
  return(p)
}

## print the information for each time point
f_print_summary = function(df){
  p2 = c(format(length(df$id)), format(unique(df$age)), format(length(df$samples)), format(mean(df$shannonmean, na.rm = T), digits = 3), 
         format(mean(df$chaomean, na.rm = T), digits = 3), f_paste_factor(df$rural_q3mgen), f_paste_factor(df$catOrDog3m),
         format(mean(df$catsOwned3m), digits = 3), format(mean(df$dogsOwned3m, na.rm = T), digits = 3), f_paste_factor(df$interventiongroup),
         f_paste_factor(df$Caesarean), f_paste_factor(df$NonCaucasian), format(mean(df$NumSiblings3m, na.rm = T), digits = 3), 
         f_paste_factor(df$Childminder12m), format(mean(df$SampleAtRoomTemp, na.rm = T), digits = 3), 
         f_paste_factor(df$OralABweek), f_paste_factor(df$OralABmonth), format(mean(df$FirstTooth, na.rm = T), digits = 3),
         f_paste_factor(df$sterlise),
         f_paste_factor(df$eczema3m),
         f_paste_factor(df$eczema12m),
         f_paste_factor(df$foodallergy),
         f_paste_factor(df$peanutallergy),
         f_paste_factor(df$eggallergy),
         f_paste_factor(df$ige009at3m),
         f_paste_factor(df$ige035at12m),
         f_paste_factor(df$ige035at36m),
         f_paste_factor(df$Frozen3m),
         format(mean(df$Weight, na.rm = T), digits=3),
         format(mean(df$Weight36m, na.rm = T), digits=3),
         format(mean(df$BMI36m, na.rm = T), digits=3),
         f_paste_factor(df$Sex),
         f_paste_factor(df$homeMould3m),
         f_paste_factor(df$diarrhea3m),
         format(mean(df$respInfections3m, na.rm = T), digits=3),
         format(mean(df$MotherAge, na.rm = T), digits=3),
         f_paste_factor(df$diarrhea12m),
         f_paste_factor(df$ParentSmoke),
         f_paste_factor(df$ParentEczema),
         f_paste_factor(df$MotherEducation),
         f_paste_factor(df$FatherEducation),
         f_paste_factor(df$diarrhea)
         ) 
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
library(lme4)
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

cvCov = colnames(dfData)[6:42]

lFits = sapply(cvCov, function(x){
  print(x)
  df = dfData[,c('id', 'age', 'shannonmean', x)]
  cvFormula = paste('shannonmean ~ 1 + ', x,' + (1 + age | id)', sep='')
  tryCatch(expr = lmerTest::lmer(cvFormula, data=df, REML=F), error=function(e) NULL) 
  #return(lmerTest::lmer(cvFormula, data=df, REML=F))
})

## which models did not fit
bNotFit = sapply(lFits, is.null)
which(bNotFit)
lFits = lFits[!bNotFit]

lPVals = lapply(lFits, f_get.coef.pvalue)
lCoef = lapply(lFits, f_get.coef.estimate)

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

