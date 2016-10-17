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

# # try a logistic regression with diversity as covariate
# library(lattice)
# xyplot(ifelse(Eczema3m == "Yes", 1, 0) ~ Shannon | Eczema12m , data = dfData,
#         type = c("p", "smooth"),
#         auto.key = list(space = "top", points = FALSE,
#                         lines = TRUE, columns = 4),
#         ylab = "Eczema at 3 months", xlab = "Shannon Diversity")
# 
# xyplot(ifelse(Eczema3m == "Yes", 1, 0) ~ Shannon , data = dfData,
#        type = c("p", "smooth"),
#        auto.key = list(space = "top", points = FALSE,
#                        lines = TRUE, columns = 4),
#        ylab = "Eczema at 3 months", xlab = "Shannon Diversity")
# 
# 
# xyplot(ifelse(Eczema12m == "Yes", 1, 0) ~ Shannon | Eczema3m , data = dfData,
#        type = c("p", "smooth"),
#        auto.key = list(space = "top", points = FALSE,
#                        lines = TRUE, columns = 4),
#        ylab = "Eczema at 12 months", xlab = "Shannon Diversity")
# 
# xyplot(ifelse(Eczema12m == "Yes", 1, 0) ~ Shannon , data = dfData,
#        type = c("g", 'smooth'),
#        auto.key = T,
#        groups = Eczema3m,
#        ylab = "Eczema at 12 months", xlab = "Shannon Diversity at 3 months")
# 
# xyplot(ifelse(Eczema12m == "Yes", 1, 0) ~ Scorad_cv3m , data = dfData,
#        type = c("p"),
#        auto.key = T,
#        groups = NULL,
#        ylab = "Eczema at 12 months", xlab = "Scorad 3 months")
# 
# xyplot(ifelse(Eczema12m == "Yes", 1, 0) ~ Scorad_cv12m , data = dfData,
#        type = c("p"),
#        auto.key = T,
#        groups = NULL,
#        ylab = "Eczema at 12 months", xlab = "Scorad 12 months")
# 
# xyplot(ifelse(Eczema3m == "Yes", 1, 0) ~ Scorad_cv3m , data = dfData, jitter=T,
#        type = c("p"),
#        auto.key = T,
#        groups = NULL,
#        ylab = "Eczema at 3 months", xlab = "Scorad 3 months")


## download the ccrossvalidation class to select variables of importance
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

# # select a subset of the data 
# dfData.sub = dfData[,-c(8,9,10,12, 13)]
# dfData.sub = na.omit(dfData.sub)
# dfData.sub = droplevels.data.frame(dfData.sub)
# fGroups = dfData.sub$Eczema12m
# dfData.sub = dfData.sub[,-8]
# str(dfData.sub)
# 
# oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups)
# plot.var.selection(oVar.r)
# dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# 
# oVar.r.12m = oVar.r

# select a subset of the data, to check what is predictive of eczema at 12 months
as.data.frame(colnames(dfData))
dfData.sub = dfData#[,-c(8,9,10,12, 13)]
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

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups, boot.num = 1000)
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

### quality check
plot(profile(fm03))

library(lattice)
xyplot(ifelse(Eczema12m == "Yes", 1, 0) ~ (Scorad_cv3m) | FoodAllergy , data = dfData, 
       type = c("p"),
       auto.key = T,
       groups = ige009.3m, pch=20,
       ylab = "Eczema at 12 months", xlab = "Scorad 3 months")

df = na.omit(dfData)

oCV = CCrossValidation.LDA(test.dat = (df[,-3]), train.dat = (df[,-3]), test.groups = df$Eczema12m,
                           train.groups = df$Eczema12m, level.predict = 'Yes', boot.num = 100, k.fold = 10)

plot.cv.performance(oCV)

p = predict(fm03, type='response', newdata=df)
bEczema12m.p = rep('No', length.out=length(p))
bEczema12m.p[p > 0.5] = 'Yes'
table(bEczema12m.p)
table(orig=df$Eczema12m, pred=bEczema12m.p)
mean(df$Eczema12m != bEczema12m.p)
mean(df$Eczema12m == bEczema12m.p)

## drop everything apart from ige
str(dfData)
dfData = dfData[,c(3, 5)]
fm03 = glm(Eczema12m ~ ., family=binomial(), data=dfData)
Anova(fm03)
anova(fm03, test='Chisq')
summary(fm03)

### quality check
plot(profile(fm03))

df = na.omit(dfData)
str(df)
oCV = CCrossValidation.LDA(test.dat = data.frame(ige=df$ige009.3m), train.dat = data.frame(ige=df$ige009.3m), test.groups = df$Eczema12m,
                           train.groups = df$Eczema12m, level.predict = 'Yes', boot.num = 100, k.fold = 10)

plot.cv.performance(oCV)

p = predict(fm03, type='response', newdata=df)
bEczema12m.p = rep('No', length.out=length(p))
bEczema12m.p[p > 0.5] = 'Yes'
table(bEczema12m.p)
table(orig=df$Eczema12m, pred=bEczema12m.p)
mean(df$Eczema12m != bEczema12m.p)
mean(df$Eczema12m == bEczema12m.p)
# 
# 
# # select a subset of the data 
# as.data.frame(colnames(dfData))
# dfData.sub = dfData[,-c(11,9,10,12, 13)]
# dfData.sub = na.omit(dfData.sub)
# dfData.sub = droplevels.data.frame(dfData.sub)
# fGroups = dfData.sub$Eczema3m
# as.data.frame(colnames(dfData.sub))
# dfData.sub = dfData.sub[,-8]
# str(dfData.sub)
# 
# oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups)
# plot.var.selection(oVar.r)
# dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# dfRF
# 
# # create formula
# cvFormula = paste0('Eczema3m ~ 1 + ', paste(rownames(dfRF), collapse = ' + '))
# 
# fm02 = glm(cvFormula, family=binomial(), data=dfData)
# library(car)
# Anova(fm02)
# anova(fm02, test='Chisq')
# summary(fm02)

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
## fit the full model
library(nnet)
fm04 = multinom(EczemaChange ~ ., data=(dfData))
summary(fm04)
# wald test p-values
z = summary(fm04)$coefficients/summary(fm04)$standard.errors
# calculate p-values from z scores
fm04.pvalues = (1 - pnorm(abs(z), 0, 1))*2
fm04.pvalues

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
rn = rn[-c(1:3)]
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

df = na.omit(dfData)

oCV = CCrossValidation.LDA(test.dat = (df[,-1]), train.dat = (df[,-1]), test.groups = df$EczemaNoYes,
                           train.groups = df$EczemaNoYes, level.predict = 'No:Yes', boot.num = 100, k.fold = 10)
#### doesnt fit, reduce variable count
##
as.data.frame(colnames(dfData))
dfData.sub = dfData
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
fGroups = dfData.sub$EczemaNoYes
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-1]
str(dfData.sub)

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
dfRF
# drop the top terms i.e. scorad 3m and eczema severity 3m
rn = rownames(dfRF)
rn = rn[-c(1:2)]
# create formula
cvFormula = paste0('EczemaNoYes ~ 1 + ', paste(rn, collapse = ' + '))
cvFormula
fm04 = glm(cvFormula, family=binomial(), data=dfData)
Anova(fm04)
anova(fm04, test='Chisq')
summary(fm04)

p = predict(fm04, type='response', newdata=na.omit(dfData))
bEczema12m.p = rep('Other', length.out=length(p))
bEczema12m.p[p > 0.5] = 'No:Yes'
table(bEczema12m.p)
table(orig=na.omit(dfData)$EczemaNoYes, pred=bEczema12m.p)
mean(na.omit(dfData)$EczemaNoYes != bEczema12m.p)
mean(na.omit(dfData)$EczemaNoYes == bEczema12m.p)

## this has no predictive power, drop some variables
rn = rownames(dfRF)
rn = rn[-c(1:3)]
cn = c('EczemaNoYes', rn)
str(dfData[,cn])
dfData = dfData[,cn]
as.data.frame(colnames(dfData))
dfData.sub = dfData
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
fGroups = dfData.sub$EczemaNoYes
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-1]
str(dfData.sub)

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
dfRF
rn = rownames(dfRF)
# create formula
cvFormula = paste0('EczemaNoYes ~ 1 + ', paste(rn, collapse = ' + '))
cvFormula
fm04 = glm(cvFormula, family=binomial(), data=dfData)
Anova(fm04)
anova(fm04, test='Chisq')
summary(fm04)

## only food allergy and days old are significant
p = predict(fm04, type='response', newdata=na.omit(dfData))
bEczema12m.p = rep('Other', length.out=length(p))
bEczema12m.p[p > 0.5] = 'No:Yes'
table(bEczema12m.p)
table(orig=na.omit(dfData)$EczemaNoYes, pred=bEczema12m.p)
mean(na.omit(dfData)$EczemaNoYes != bEczema12m.p)
mean(na.omit(dfData)$EczemaNoYes == bEczema12m.p)

# refit the model with only the 2 variables
fm05 = update(fm04, EczemaNoYes ~ DaysOldWhenSampled + FoodAllergy)
Anova(fm05)
summary(fm05)
plogis(coef(fm05)[2])


## what is the relationship of days old when sampled
xyplot(ifelse(dfData$EczemaNoYes == "No:Yes", 1, 0) ~ DaysOldWhenSampled | FoodAllergy, data = dfData,
       type = c("p", "smooth"),
       ylab = "Developed Eczema from 3 to 12 months", xlab = "DaysOldWhenSampled")


## try in other direction, i.e. how many got better from 3 to 12 months
table(dfData$Eczema3m:dfData$Eczema12m)

fEczemaYesNo = rep('Yes:No', length=nrow(dfData))
f = dfData$Eczema3m:dfData$Eczema12m
i = which(f != 'Yes:No')
fEczemaYesNo[i] = 'Other'
fEczemaYesNo = factor(fEczemaYesNo, levels=c('Other', 'Yes:No'))
table(fEczemaYesNo)
levels(fEczemaYesNo)

# select a subset of the data 
as.data.frame(colnames(dfData))
dfData.sub = dfData[,-c(8, 9, 10, 11, 12, 13)]
dfData.sub$fEczemaYesNo = fEczemaYesNo
dfData.sub = na.omit(dfData.sub)
dfData.sub = droplevels.data.frame(dfData.sub)
fGroups = dfData.sub$fEczemaYesNo
as.data.frame(colnames(dfData.sub))
dfData.sub = dfData.sub[,-15]
str(dfData.sub)

oVar.r = CVariableSelection.RandomForest(data = dfData.sub, fGroups)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
dfRF
dfData.sub = dfData
dfData.sub$fEczemaYesNo = fEczemaYesNo
rm(fEczemaYesNo)

# create formula
cvFormula = paste0('fEczemaYesNo ~ 1 + ', paste(rownames(dfRF), collapse = ' + '))
cvFormula
fm04 = glm(cvFormula, family=binomial(), data=dfData.sub)
library(car)
Anova(fm04)
anova(fm04, test='Chisq')
summary(fm04)

## what is the relationship of days old when sampled
xyplot(ifelse(fEczemaYesNo == "Yes:No", 1, 0) ~ DaysOldWhenSampled, data = dfData.sub,
       type = c("p", "smooth"),
       auto.key = list(space = "top", points = FALSE,
                       lines = TRUE, columns = 4),
       ylab = "Eczema healed from 3 to 12 months", xlab = "DaysOldWhenSampled")



dfData.sub$change = dfData.sub$Eczema3m:dfData.sub$Eczema12m
xyplot(ifelse(Eczema3m == "Yes", 1, 0) ~ DaysOldWhenSampled | change, data = dfData.sub,
       type = c("p"),
       auto.key = list(space = "top", points = FALSE,
                       lines = TRUE, columns = 4),
       ylab = "Eczema healed from 3 to 12 months", xlab = "DaysOldWhenSampled")





















fm02 = update(fm01, Eczema12m ~ Shannon + Eczema3m)
summary(fm02)

p = predict(fm02, data=dfData, type='response')
plot(dfData$Shannon, p)

fm02 = update(fm01, Eczema12m ~ Eczema3m)
summary(fm02)

fm03 = glm(Eczema3m ~ Shannon, family=binomial(), data=dfData)
summary(fm03)
fm04 = update(fm03, Eczema3m ~ Shannon)
summary(fm04)

p = predict(fm01, data=dfData, type='response')
plot(dfData$Shannon, p, pch=20)

str(dfData[,-1])

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

dfReport = data.frame(Shannon = dfShannon$`P-Value`)
rownames(dfReport) = rownames(dfShannon)

df = data.frame(lFm02[[1]][1])
for(i in 2:length(lFm02)){
  df = cbind(df, data.frame(lFm02[[i]][,1]))
}
colnames(df) = names(lFm02)

dfReport = cbind(dfReport, df)
write.csv(dfReport, file='3m_trimmed/Results/Covariate_p_values.xls')
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
dfData = data.frame(Shannon=dfImport$shannonmean, Chao=dfImport$chaomean, Eczema3m=dfImport$ecz3m, Eczema3m.sev=factor(dfImport$eczsev3m),
                    Eczema12m=dfImport$exam12mecz)
str(dfData)

boxplot(Shannon ~ Eczema3m, data=dfData, ylab='3M Shannon Diversity', main='3M Shannon diversity and Eczema status at 3 months')
boxplot(Shannon ~ Eczema12m, data=dfData, ylab='3M Shannon Diversity', main='3M Shannon diversity and Eczema status at 12 months')

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

boxplot(Shannon ~ Eczema3m.Diversity, data=dfData, ylab='Shannon Diversity 3M', main='High diversity has a protective effect')
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

