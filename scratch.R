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

#################################################
t = dfData$Chao
h = hist(t, plot=F)
dn = dt_ls(h$mids, df = 2.5, mu = 403, a = 118)
# which distribution can approximate the frequency
hist(t, prob=T, xlab='Shannon Diversity', ylab='', main='Distribution of Shannon Diversity in Disease')
# parameterized on the means
lines(h$mids, dn, col='black', type='b')
#points(qnorm(0.95, mean(t), sd(t)), 0, pch=20, col='red', cex=2)
legend('topright', legend =c('Normal Density Curve'), fill = c('black'))

library(heavy)
hfit = heavyLm(Chao ~ eczema3m*frozen + diarrhoeaCategorical3m + abxMonth , data=dfData, family=Student(df=2.95))
summary(hfit)

