# scratch.R

dfData = data.frame(dfImport[,n])

## find correlated variables
m = NULL;

for (i in 1:ncol(dfData)){
  m = cbind(m, dfData[,i])
}
colnames(m) = colnames(dfData)
mCor = cor(m, use="na.or.complete")
library(caret)

i = findCorrelation(abs(mCor), cutoff = 0.7)
n = colnames(mCor)[i]