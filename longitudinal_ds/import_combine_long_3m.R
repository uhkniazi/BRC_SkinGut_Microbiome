# Name: import_combine_long_3m.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 31/10/2016
# Desc: import the longitudinal dataset and 3m datasets and merge


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

dfData.long = dfData

## 3 month dataset
dfImport = read.csv(file.choose(), header=T)

str(dfImport)
head(dfImport)


dfData = data.frame(id = dfImport$id, Shannon=dfImport$shannonmean, Frozen=factor(ifelse(dfImport$frozen1yes == 1, 'Yes', 'No')),
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

dfData.3m = dfData

i = match(as.character(dfData.3m$id), as.character(dfData.long$id))
dfData = dfData.3m
dfData$intervention = dfData.long$intervention[i]  

i1 = dfData[,c('id', 'intervention')]
i1 = na.omit(i1)
i2 = dfData.long[,c('id', 'intervention')]
i2 = i2[!duplicated(i2$id),]
table(i1$id %in% i2$id)

str(dfData)

# the response variable to model is eczema status at 12 months
# remove NA values from this 
i = is.na(dfData$Eczema12m)
table(i)

# 26 values missing
dfData = dfData[!i,]


## download the ccrossvalidation class to select variables of importance
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')
library(corrplot)
# url = 'http://www.sthda.com/upload/rquery_cormat.r'
# download(url, 'rquery_cormat.r')
source('rquery_cormat.r')
# drop the id
dfData = dfData[,-c(1)]

# check what is predictive of eczema at 12 months
# make a correlation plot for the variables to see their association
as.data.frame(colnames(dfData))

m = NULL;

for (i in 1:ncol(dfData)){
  m = cbind(m, dfData[,i])
}
colnames(m) = colnames(dfData)

rquery.cormat(m)

str(dfData)
dfData.sub = dfData
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
dim(dfData.sub)
fGroups = dfData.sub$Eczema12m
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-11]
str(dfData.sub)

# oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups, boot.num = 3000)
# plot.var.selection(oVar.r)
# dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
rf = randomForest(Eczema12m ~., data=na.omit(dfData), importance = TRUE, ntree = 500)
varImpPlot(rf)
dfRF = importance(rf)
dfRF = dfRF[order(dfRF[,3],decreasing = T ),]
dfRF
rn = rownames(dfRF)
as.data.frame(rn)
## choose the 3 + intervention variable we used earlier in eczema12m_diversity.R script for 3m data analysis
rn = rn[c(4, 5, 14)] 
cn = c('Eczema12m', rn)
str(dfData[,cn])
dfData = dfData[,cn]
summary(dfData)

# perform a multivariate and univariate logitic regression
fm03 = glm(Eczema12m ~ 1 + ., family=binomial(), data=dfData)
library(car)
Anova(fm03)
anova(fm03, test='Chisq')
summary(fm03)
exp(coef(fm03))

df = na.omit(dfData)
str(df)
oCV = CCrossValidation.LDA(test.dat = (df[,-1]), train.dat = (df[,-1]), test.groups = df$Eczema12m,
                           train.groups = df$Eczema12m, level.predict = 'Yes', boot.num = 100, k.fold = 10)

plot.cv.performance(oCV)

## univariate analysis
df = dfData
str(df)

fm03.s = glm(Eczema12m ~ 1 + intervention, data=df, family=binomial)
summary(fm03.s)
exp(coef(fm03.s))

df = na.omit(dfData[,c('Eczema12m', 'intervention')])
str(df)

oCV = CCrossValidation.LDA(test.dat = data.frame(g=df[,-1]), train.dat = data.frame(g=df[,-1]), test.groups = df$Eczema12m,
                           train.groups = df$Eczema12m, level.predict = 'Yes', boot.num = 100, k.fold = 10)

plot.cv.performance(oCV)


