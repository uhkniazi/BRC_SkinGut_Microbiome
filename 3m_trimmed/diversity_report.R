# Name: diversity_report.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 02/11/2016
# Desc: multiple regression analysis for diversity at 3 months



dfImport = read.csv(file.choose(), header=T)

str(dfImport)
head(dfImport)

library(lattice)
library(car)



### organize the data
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


str(dfData)

### find the associated variables
## download the ccrossvalidation class to select variables of importance
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

# check what is predictive of diversity 
as.data.frame(colnames(dfData))
fit.rf = randomForest(Shannon ~., data=na.omit(dfData), importance = TRUE, ntree = 500)
# get variables importance
varImpPlot(fit.rf)
dfRF = importance(fit.rf)
dfRF = dfRF[order(dfRF[,1],decreasing = T ),]
dfRF

dfData.sub = dfData
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
dim(dfData.sub)
fGroups = dfData.sub$Shannon
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-1]
str(dfData.sub)

rn = rownames(dfRF)
data.frame(rn)
## look at the covariate variables and perform subset selection
i = which(colnames(dfData.sub) %in% rn[c(1, 2, 6, 7, 9, 10, 11, 14, 15, 17, 18)])
str(dfData.sub[,i])
dfData.sub = dfData.sub[,i]
str(dfData.sub)

oVar.s = CVariableSelection.ReduceModel(dfData.sub, fGroups, boot.num = 50)
plot.var.selection(oVar.s)

sapply(seq_along(1:ncol(dfData.sub)), function(x) CVariableSelection.ReduceModel.getMinModel(oVar.s, x))

cn = CVariableSelection.ReduceModel.getMinModel(oVar.s, 11)
cn
## perform regression analysis on all covariates variables
# create formula
cvFormula = paste0('Shannon ~ 1 + ', paste(cn, collapse = ' + '))
cvFormula
fm01 = lm(cvFormula, data=dfData)
Anova(fm01)
anova(fm01, test='Chisq')
summary(fm01)

## the 5 variable model from subset selection and after regression
cn = CVariableSelection.ReduceModel.getMinModel(oVar.s, 3)
cn


library(lattice)
str(dfData)
## how does the density compare with eczema at 3 months and covariates
densityplot(~ Shannon, data=dfData, auto.key=T, groups=Frozen, main='Shannon Diversity 3m given Frozen')
densityplot(~ Shannon, data=dfData, auto.key=T, groups=Pets, main='Shannon Diversity 3m given Pets')
densityplot(~ Shannon, data=dfData, auto.key=T, groups=Male, main='Shannon Diversity 3m given Male')
densityplot(~ Shannon, data=dfData, auto.key=T, groups=Siblings, main='Shannon Diversity 3m given Siblings')

xyplot(Shannon ~ DaysOldWhenSampled | Frozen, data=dfData, groups=Eczema3m, auto.key = T)

densityplot(~ Shannon, data=dfData, auto.key=T, groups=Eczema3m, main='Shannon Diversity 3m Eczema Status 3m')

## choose a 4 variable model 
cn = c('Frozen', 'Pets', 'Siblings', 'Eczema3m')
## perform regression analysis on all covariates variables
# create formula
cvFormula = paste0('Shannon ~ ', paste(cn, collapse = ' + '))
cvFormula
fm02 = lm(cvFormula, data=dfData)
Anova(fm02)
anova(fm02, test='Chisq')
summary(fm02)

fm03 = lm(Shannon ~ 1 + Eczema3m, data=dfData)
summary(fm03)

xyplot(ifelse(dfData$Eczema3m == "Yes", 1, 0) ~ Shannon | Frozen, data = dfData,
       type = c("p", "smooth"),
       ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity given Frozen')


xyplot(ifelse(dfData$Eczema3m == "Yes", 1, 0) ~ Shannon, data = dfData,
       type = c("p"),
       ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity')

xyplot(ifelse(dfData$Eczema3m == "Yes", 1, 0) ~ Shannon | Frozen, data = dfData,
       type = c("p"),
       ylab = "Eczema at 3 months", xlab = "Shannon Diversity", main='Eczema at 3m vs Diversity given Frozen')

## add an interaction term after observing this pattern
fm03 = update(fm02, Shannon ~ 1 + Frozen + Pets + Siblings + Eczema3m + Eczema3m:Frozen)
summary(fm03)
anova(fm02, fm03, test='Chisq')



# 
# 
# densityplot(~ Shannon | Frozen, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Frozen')
# densityplot(~ Shannon | Male, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Male')
# densityplot(~ Shannon | Caesarean, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Caesarean')
# densityplot(~ Shannon | Siblings, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Siblings')
# densityplot(~ Shannon | Pets, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Pets')
# densityplot(~ Shannon | Eczema3m.sev, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Eczema 3m Severity')
# densityplot(~ Shannon | Eczema12m, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given Eczema 12m')
# densityplot(~ Shannon | FoodAllergy, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given FoodAllergy')
# densityplot(~ Shannon | AntiBiotics, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given AntiBiotics')
# densityplot(~ Shannon | ige009.3m, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given ige009.3m')
# densityplot(~ Shannon | ige035.12m, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given ige035.12m')
# densityplot(~ Shannon | anyfoodsensitised_cv12m, data=dfData, auto.key=T, groups=Eczema3m, 
#             main='Eczema at 3 months given anyfoodsensitised_cv12m')
# densityplot(~ Shannon | ige035.36m, data=dfData, auto.key=T, groups=Eczema3m, main='Eczema at 3 months given ige035.36m')









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

