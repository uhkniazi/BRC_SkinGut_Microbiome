data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Nclusters; // number of groups for random intercepts
    int<lower=1> NgroupMap[Ntotal]; // mapping variable to map each observation to a group 
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    real y[Ntotal]; // response variable gamma distributed
}
// transformed data {
// }
parameters {
  // parameters to estimate in the model
    row_vector[Ncol] betas; // regression parameters
    //real betaRandom; // random intercept
    real<lower=0> sigmaPop; // population standard deviation
    real<lower=0> sigmaRan; // random effect standard deviation
    vector[Nclusters] rGroupsJitter; // number of random jitters for each cluster member
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  matrix[Ntotal, Ntotal] muMat; 
  vector[Ntotal] alpha;// shape parameter
  vector[Ntotal] beta;// rate parameter
  vector[Ntotal] rNewIntercept;
  matrix[Ntotal, Ncol] mBetas;
  // add this random jitter to the population intercept
  vector[Nclusters] rBetaRand;
  rBetaRand = rGroupsJitter + betas[1];
  // new intercept variable equal to number of observations
  rNewIntercept = rBetaRand[NgroupMap];
  // fill this matrix
  for (r in 1:Ntotal) mBetas[r] = betas;
  // replace first column i.e. original intercept with group
  // intercept
  mBetas[,1] = rNewIntercept;
  muMat = X * (mBetas)'; 
  # get the diagonal of this matrix
  for(i in 1:Ntotal) mu[i] = muMat[i,i];
  mu = exp(mu);
  // convert this fitted value to the gamma parameters
  alpha = mu .* mu / sigmaPop; 
  beta = mu/sigmaPop; 
}
model {
  // using non-informative priors to start with
  //sigmaPop ~ uniform(0, 1e10);
  //sigmaRan ~ uniform(0, 1e3);
  // random effects sample
  rGroupsJitter ~ normal(0, sigmaRan);
  
  // likelihood function
  y ~ gamma(alpha, beta);
}