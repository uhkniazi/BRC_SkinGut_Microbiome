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
## some coda tests


library(coda)
m = as.matrix(fit.stan)
mc = mcmc(m)
head(mc)
ncol(fit.stan)
plot(fit.stan)
pairs(fit.stan)
m = extract(fit.stan)
length(m)
names(m)
sapply(m, length)
head(mc)
library(lattice)
xyplot(mc)
autocorr.plot(mc)
plot(cumsum(mc[,1])/1:20000)
plot(cumsum(mc[,2])/1:20000)
plot(cumsum(mc[,3])/1:20000)
plot(cumsum(mc[,4])/1:20000)
plot(cumsum(mc[,5])/1:20000)
# cumplot(mc) # not run slow
acfplot(mc)
acf(mc)
heidel.diag(mc)
geweke.diag(mc)
geweke.plot(mc)

mc2 = extract(fit.stan, permuted=F)
dim(mc2)
# second dimension is the chains
# get samples split into chains
mc2 = mcmc.list( lapply( 1:ncol(fit.stan) , 
                         function(x) { mcmc(as.array(fit.stan)[,x,]) } ) )
gelman.diag(mc2)
gelman.plot(mc2)
heidel.diag(mc2)
geweke.diag(mc2)
geweke.plot(mc2)
autocorr.diag(mc)
autocorr.diag(mc2)
effectiveSize(mc)
effectiveSize(mc2)