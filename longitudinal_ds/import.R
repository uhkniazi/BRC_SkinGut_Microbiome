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
         f_paste_factor(df$sterlise)
         ) 
  p2 = paste(c('||', paste(p2, '||', sep='', collapse = '')), collapse = '')
  return(p2)
}

writeLines(f_print_summary(dfData[dfData$age == 3,]), oFile)
writeLines(f_print_summary(dfData[dfData$age == 5,]), oFile)
writeLines(f_print_summary(dfData[dfData$age == 12,]), oFile)
close(oFile)

write.csv(dfData, file='longitudinal_ds/Data_external/longitudinal_dataset_manual_UN_2.csv')


library(lattice)
library(MASS)
library(car)
hist(dfData$shannonmean)
fitdistr(dfData$shannonmean, 'normal')
qqp(dfData$shannonmean, 'norm', mean=3.4, sd=0.94, ylab='Shannon diversity')
qqPlot(dfData$shannonmean, 'norm',ylab='Shannon diversity')

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

