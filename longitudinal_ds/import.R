# Name: import.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 31/10/2016
# Desc: import the longitudinal dataset and cleanup


dfImport = read.csv(file.choose(), header=T)

str(dfImport)
head(dfImport)

# setup data
dfData = data.frame(id=dfImport$id)
dfData$age = as.numeric(gsub('m', '', dfImport$age))
dfData$Shannon = dfImport$shannonmean
rownames(dfData) = dfImport$samples
dfData$rural = dfImport$rural_q3mgen
dfData$pet = dfImport$petownfur
dfData$intervention = dfImport$interventiongroup
dfData$caesarean = dfImport$mode
dfData$ethnicity = dfImport$eth
dfData$siblings = dfImport$numsibs

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

# ## antibiotic data needs to be coded correctly
# df = data.frame(m3=dfImport$abxmnth3m,
#                 m5=dfImport$abxmnth5m,
#                 m12=dfImport$abxmnth12m,
#                 age=dfData$age,
#                 id=dfData$id)
# 
# #
# # fm2 = function(df){
# #   id = as.character(df$id)
# #   v = rep(NA, length.out=length(id))
# #   names(v) = id
# #   for (i in 1:length(id)){
# #     a = df$age[i]
# #     am = paste0('m',a)
# #     v[i] = ifelse(is.na(df[i,am]), NA, as.character(df[i,am]))
# #   }
# #   return(v)
# # }

str(dfData)

library(lattice)
library(MASS)
library(car)
hist(dfData$Shannon)
fitdistr(dfData$Shannon, 'normal')
qqp(dfData$Shannon, 'norm', mean=3.4, sd=0.94, ylab='Shannon diversity')
qqPlot(dfData$Shannon, 'norm',ylab='Shannon diversity')

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

xyplot(Shannon ~ age | intervention+ifelse(rural == 'Yes', 'Rural=Y', 'Rural=N')+ifelse(sterlise == 'Yes', 'Ster=Y', 'Ster=N')+
         ifelse(caesarean == 'Yes', 'caes=Y', 'caes=N'),
       data=dfData, type=c('g', 'smooth', 'r'), index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy', layout=c(8,2), 
       par.strip.text=list(cex=0.7))

xyplot(Shannon ~ age | intervention,
       data=dfData, type=c('g', 'r'), index.cond = function(x,y) coef(lm(y ~ x))[1],
       par.strip.text=list(cex=0.7), groups=id)


xyplot(Shannon ~ age | id, index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',
       data=dfData, type=c('g', 'p', 'r'), par.strip.text=list(cex=0.4))#, layout=c(12,5))

xyplot(Shannon ~ age | id, index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',
       data=dfData, type=c('g', 'p', 'r'), par.strip.text=list(cex=0.4), groups=intervention, auto.key = list(columns=2))#, layout=c(12,5))

densityplot(~ Shannon | intervention, data=dfData, main='Intervention')
densityplot(~ Shannon | caesarean, data=dfData, main='Caesarean')
densityplot(~ Shannon | sterlise, data=dfData, main='Sterlise')
densityplot(~ Shannon | rural, data=dfData, main='Rural')
densityplot(~ Shannon | pet, data=dfData, main='Pet')
densityplot(~ Shannon | siblings, data=dfData, main='Siblings')
densityplot(~ Shannon | ethnicity, data=dfData, main='Ethnicity')

xyplot(Shannon ~ age | intervention+sterlise+caesarean, index.cond = function(x,y) coef(lm(y ~ x))[1], aspect='xy',
       data=dfData, type=c('g', 'p', 'r'), layout=c(4,2))#, par.strip.text=list(cex=0.4))#, layout=c(12,5))

