# Name: eczema12m_diversity.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 13/10/2016
# Desc: explore the relationship between eczema12m and 3 months and diversity


dfImport = read.csv(file.choose(), header=T)

str(dfImport)
head(dfImport)

## original variables i.e. 
## shannon = normal
## chao = gamma
## simpson = beta

## transformed scale
## shannon = no transformation
## chao = log = normal
## simpson = logit = normal

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

# the response variable to model is eczema status at 12 months
# remove NA values from this 
i = is.na(dfData$Eczema12m)
table(i)

# 26 values missing
dfData = dfData[!i,]

t = xtabs(~ Eczema12m + Eczema3m, dfData)
t
## perform McNemar's chi-squared test to see if the change between conditions is significant
## difference between changes in opposite directions
mcnemar.test(t)
mcnemar.test(t, correct=F)

## proportion of change
(36+38)/sum(rowSums(t))
se = sqrt(abs(0.28 - (1-0.28))/261)
0.28+2*se
0.28-2*se


## download the ccrossvalidation class to select variables of importance
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

# check what is predictive of eczema at 12 months
as.data.frame(colnames(dfData))
dfData.sub = dfData
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
dim(dfData.sub)
fGroups = dfData.sub$Eczema12m
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-11]
str(dfData.sub)

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups, boot.num = 100)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
dfRF

## remove the top 2 predictors i.e. score and eczema severity at 12 months
rn = rownames(dfRF)[1:2]
i = which(colnames(dfData) %in% rn)
dfData = dfData[,-i]
str(dfData)

### repeat with reduced variables
# select a subset of the data, to check what is predictive of eczema at 12 months
as.data.frame(colnames(dfData))
dfData.sub = dfData
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
fGroups = dfData.sub$Eczema12m
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-11]
str(dfData.sub)

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups, boot.num = 100)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
dfRF

# create formula
cvFormula = paste0('Eczema12m ~ 1 + ', paste(rownames(dfRF), collapse = ' + '))
cvFormula
fm01 = glm(cvFormula, family=binomial(), data=dfData)
library(car)
Anova(fm01)
anova(fm01, test='Chisq')
summary(fm01)

## remove eczema severity from the data
i = which(colnames(dfData) %in% 'Eczema3m.sev')
dfData = dfData[,-i]
str(dfData)

### repeat with reduced variables
as.data.frame(colnames(dfData))
dfData.sub = dfData
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
fGroups = dfData.sub$Eczema12m
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-10]
str(dfData.sub)

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups, boot.num = 100)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
dfRF

# create formula
cvFormula = paste0('Eczema12m ~ 1 + ', paste(rownames(dfRF), collapse = ' + '))
cvFormula
fm02 = glm(cvFormula, family=binomial(), data=dfData)
Anova(fm02)
anova(fm02, test='Chisq')
summary(fm02)

# select a subset of the important variables
a = Anova(fm02)
i = which(a$`Pr(>Chisq)` < 0.1)
coef(fm02)[i+1]
rn = rownames(data.frame(a))[i]

#### further reduce the model size
####
rn
i = which(colnames(dfData) %in% c(rn, 'Eczema12m'))
str(dfData[,i])
dfData = dfData[,i]

### repeat with reduced variables
# select a subset of the data, to check what is predictive of eczema at 12 months
as.data.frame(colnames(dfData))
dfData.sub = dfData
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
fGroups = dfData.sub$Eczema12m
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-3]
str(dfData.sub)

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups, boot.num = 100)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
dfRF

# create formula
cvFormula = paste0('Eczema12m ~ 1 + ', paste(rownames(dfRF), collapse = ' + '))
cvFormula
fm03 = glm(cvFormula, family=binomial(), data=dfData)
Anova(fm03)
anova(fm03, test='Chisq')
summary(fm03)
# remove food allergy
rn = rownames(dfRF)
rn = rn[c(1, 2, 3)]
cn = c('Eczema12m', rn)
str(dfData[,cn])
dfData = dfData[,cn]
summary(dfData)

fm03 = glm(Eczema12m ~ 1 + ., family=binomial(), data=dfData)
Anova(fm03)
anova(fm03, test='Chisq')
summary(fm03)

df = na.omit(dfData)

oCV = CCrossValidation.LDA(test.dat = (df[,-1]), train.dat = (df[,-1]), test.groups = df$Eczema12m,
                           train.groups = df$Eczema12m, level.predict = 'Yes', boot.num = 100, k.fold = 10)

plot.cv.performance(oCV)

s = summary(fm03)$coefficients
c = format(s[,1], digi=3)
s.e = format(s[,2], digi=3)
p = format(s[,4], digi=3)
o = format(exp(s[,1]), digi=3)

oFile = file('3m_trimmed/Results/eczema12m_covariates.xls', 'wt')

writeLines('Prediction of Eczema status at 12 months using logistic regression', oFile)
writeLines('Covariate,Coefficient,Std_err,P-Value,Odds', oFile)
p1 = paste(rownames(s), c, s.e, p, o, sep=',')
writeLines(p1, oFile)
writeLines('\n\n', oFile)
close(oFile)


## try multinomial regression before binomial
############################################################################
# what about an interaction between eczema12m and eczema3m, which people
# had eczema at 3 months and don't have it at 12 and vice versa
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

# the response variable to model is eczema status at 12 months
# remove NA values from this 
i = is.na(dfData$Eczema12m)
table(i)

# 26 values missing
dfData = dfData[!i,]

table(dfData$Eczema3m:dfData$Eczema12m)
dfData$EczemaChange = factor(dfData$Eczema3m:dfData$Eczema12m)
dfData$EczemaChange = relevel(dfData$EczemaChange, 'No:No')

# which variables are the most predictive
df = na.omit(dfData)
rf = randomForest(EczemaChange ~ ., data=df, importance=T)
dfImp = importance(rf)
dfImp = dfImp[order(dfImp[,'MeanDecreaseAccuracy'], decreasing = T),]

## remove the highest scoring variables
dfImp
rn = rownames(dfImp)[-c(1:6)]
str(dfData[,c('EczemaChange', rn)])
dfData = dfData[,c('EczemaChange', rn)]
# remove na 
dfData = na.omit(dfData)
## fit the full model
library(nnet)
fm04 = multinom(EczemaChange ~ ., data=(dfData))
summary(fm04)
# wald test p-values
z = summary(fm04)$coefficients/summary(fm04)$standard.errors
# calculate p-values from z scores
fm04.pvalues = (1 - pnorm(abs(z), 0, 1))*2
fm04.pvalues

# put these in a dataframe and examine
c = format(t(summary(fm04)$coefficients), digi=3)
s = format(t(summary(fm04)$standard.errors), digi=3)
p = format(t(fm04.pvalues), digi=3)
o = format(exp(t(summary(fm04)$coefficients)), digi=3)

oFile = file('3m_trimmed/Results/change_in_eczema_status_3_to_12_months_covariates_2.xls', 'wt')

writeLines('No:Yes - developed eczema from 3 to 12 months', oFile)
writeLines('Covariate,Coefficient,Std_err,P-Value,Odds', oFile)
p1 = paste(rownames(c), c[,1], s[,1], p[,1], o[,1], sep=',')
writeLines(p1, oFile)
writeLines('\n\n', oFile)

writeLines('Yes:No - lost eczema from 3 to 12 months', oFile)
writeLines('Covariate,Coefficient,Std_err,P-Value,Odds', oFile)
p1 = paste(rownames(c), c[,2], s[,2], p[,2], o[,2], sep=',')
writeLines(p1, oFile)
writeLines('\n\n', oFile)

writeLines('Yes:Yes - stayed eczema from 3 to 12 months', oFile)
writeLines('Covariate,Coefficient,Std_err,P-Value,Odds', oFile)
p1 = paste(rownames(c), c[,3], s[,3], p[,3], o[,3], sep=',')
writeLines(p1, oFile)
writeLines('\n\n', oFile)

close(oFile)

## repeat after using only significant variables
rn = rownames(dfImp)
as.data.frame(rn)
rn = c('ige009.3m', 'Male', 'DaysOldWhenSampled')#rn[c(7,15)]

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

# the response variable to model is eczema status at 12 months
# remove NA values from this 
i = is.na(dfData$Eczema12m)
table(i)

# 26 values missing
dfData = dfData[!i,]

table(dfData$Eczema3m:dfData$Eczema12m)
dfData$EczemaChange = factor(dfData$Eczema3m:dfData$Eczema12m)
dfData$EczemaChange = relevel(dfData$EczemaChange, 'No:No')

str(dfData[,c('EczemaChange', rn)])
dfData = dfData[,c('EczemaChange', rn)]
summary(dfData)
dim(na.omit(dfData))
dim(dfData)
dfData = na.omit(dfData)

fm04 = multinom(EczemaChange ~ ., data=(dfData))
summary(fm04)
# wald test p-values
z = summary(fm04)$coefficients/summary(fm04)$standard.errors
# calculate p-values from z scores
fm04.pvalues = (1 - pnorm(abs(z), 0, 1))*2
fm04.pvalues

# put these in a dataframe and examine
c = format(t(summary(fm04)$coefficients), digi=3)
s = format(t(summary(fm04)$standard.errors), digi=3)
p = format(t(fm04.pvalues), digi=3)
o = format(exp(t(summary(fm04)$coefficients)), digi=3)

oFile = file('3m_trimmed/Results/change_in_eczema_status_3_to_12_months_covariates_3.xls', 'wt')

writeLines('No:Yes - developed eczema from 3 to 12 months', oFile)
writeLines('Covariate,Coefficient,Std_err,P-Value,Odds', oFile)
p1 = paste(rownames(c), c[,1], s[,1], p[,1], o[,1], sep=',')
writeLines(p1, oFile)
writeLines('\n\n', oFile)

writeLines('Yes:No - lost eczema from 3 to 12 months', oFile)
writeLines('Covariate,Coefficient,Std_err,P-Value,Odds', oFile)
p1 = paste(rownames(c), c[,2], s[,2], p[,2], o[,2], sep=',')
writeLines(p1, oFile)
writeLines('\n\n', oFile)

writeLines('Yes:Yes - stayed eczema from 3 to 12 months', oFile)
writeLines('Covariate,Coefficient,Std_err,P-Value,Odds', oFile)
p1 = paste(rownames(c), c[,3], s[,3], p[,3], o[,3], sep=',')
writeLines(p1, oFile)
writeLines('\n\n', oFile)

close(oFile)


############################################################################
# what about an interaction between eczema12m and eczema3m, which people
# had eczema at 3 months and don't have it at 12 and vice versa
## use 2 group comparisons as a binomial problem
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

# the response variable to model is eczema status at 12 months
# remove NA values from this 
i = is.na(dfData$Eczema12m)
table(i)

# 26 values missing
dfData = dfData[!i,]

table(dfData$Eczema3m:dfData$Eczema12m)

fEczemaNoYes = rep('No:Yes', length=nrow(dfData))
f = dfData$Eczema3m:dfData$Eczema12m
i = which(f != 'No:Yes')
fEczemaNoYes[i] = 'Other'
fEczemaNoYes = factor(fEczemaNoYes, levels=c('Other', 'No:Yes'))
table(fEczemaNoYes)
levels(fEczemaNoYes)

dfData$EczemaNoYes = fEczemaNoYes
# select a subset of the data 
as.data.frame(colnames(dfData))
dfData.sub = dfData#[,-c(8, 9, 10, 11, 12, 13)]
rm(fEczemaNoYes)
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
fGroups = dfData.sub$EczemaNoYes
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-21]
str(dfData.sub)

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
dfRF
# drop the top terms i.e. eczema status at 12m 3m and scores
rn = rownames(dfRF)
rn = rn[-c(1:6)]
# create formula
cvFormula = paste0('EczemaNoYes ~ 1 + ', paste(rn, collapse = ' + '))
cvFormula
fm04 = glm(cvFormula, family=binomial(), data=dfData)
Anova(fm04)
anova(fm04, test='Chisq')
summary(fm04)

## drop these terms from the data.frame
cn = c('EczemaNoYes', rn)
str(dfData[,cn])
dfData = dfData[,cn]

df = na.omit(data.frame(days=dfData$DaysOldWhenSampled, EczemaNoYes=dfData$EczemaNoYes))
df2 = data.frame(days=df$days)
oCV = CCrossValidation.LDA(test.dat = df2, train.dat = df2, test.groups = df$EczemaNoYes,
                           train.groups = df$EczemaNoYes, level.predict = 'No:Yes', boot.num = 100, k.fold = 10)
plot.cv.performance(oCV)



# refit the model with only the 1 variable
fm05 = update(fm04, EczemaNoYes ~ DaysOldWhenSampled)
Anova(fm05)
summary(fm05)
plogis(coef(fm05)[2])
table(dfData$EczemaNoYes)
library(lattice)
## what is the relationship of days old when sampled
xyplot(ifelse(dfData$EczemaNoYes == "No:Yes", 1, 0) ~ DaysOldWhenSampled, data = dfData,
       type = c("p", "smooth"),
       ylab = "Developed Eczema from 3 to 12 months", xlab = "DaysOldWhenSampled")


####### repeat the analysis in other direction
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

# the response variable to model is eczema status at 12 months
# remove NA values from this 
i = is.na(dfData$Eczema12m)
table(i)

# 26 values missing
dfData = dfData[!i,]

table(dfData$Eczema3m:dfData$Eczema12m)

fEczemaYesNo = rep('Yes:No', length=nrow(dfData))
f = dfData$Eczema3m:dfData$Eczema12m
i = which(f != 'Yes:No')
fEczemaYesNo[i] = 'Other'
fEczemaYesNo = factor(fEczemaYesNo, levels=c('Other', 'Yes:No'))
table(fEczemaYesNo)
levels(fEczemaYesNo)
dfData$EczemaYesNo = fEczemaYesNo
rm(fEczemaYesNo)
str(dfData)

as.data.frame(colnames(dfData))
dfData.sub = dfData
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
fGroups = dfData.sub$EczemaYesNo
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-21]
str(dfData.sub)

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
dfRF

# drop the top terms i.e. eczema status at 12m 3m and scores
rn = rownames(dfRF)
rn = rn[-c(1:6)]
# create formula
cvFormula = paste0('EczemaYesNo ~ 1 + ', paste(rn, collapse = ' + '))
cvFormula
fm06 = glm(cvFormula, family=binomial(), data=dfData)
Anova(fm06)
anova(fm06, test='Chisq')
summary(fm06)

## drop these terms from the data.frame
cn = c('EczemaYesNo', rn)
str(dfData[,cn])
dfData = dfData[,cn]


# refit the model with only the 1 variable
fm07 = update(fm06, EczemaYesNo ~ AntiBiotics.any)
Anova(fm07)
summary(fm07)
table(dfData$EczemaYesNo)

