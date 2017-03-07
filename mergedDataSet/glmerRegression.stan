data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Nclusters; // number of groups for random intercepts
    int<lower=1> NgroupMap[Ntotal]; // mapping variable to map each observation to a group 
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    real y[Ntotal]; // response variable gamma distributed
}
// transformed data {
//   real unifLo;
//   real unifHi; // the lower and upper bounds for the t-distrubtion scale
//   real expLambda; // rate parameter for the exponential distribution for nu
//   unifLo = sdY/1000;
//   unifHi = sdY*1000;
//   expLambda = 1/29.0 ;
// }
parameters {
  // parameters to estimate in the model
    row_vector[Ncol] betas; // regression parameters
    //real betaRandom; // random intercept
    real<lower=0> sigmaPop; // population standard deviation
    real<lower=0> sigmaRan; // random effect standard deviation
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  vector[Nclusters] rGroupsJitter; // number of random jitters for each cluster member
  matrix[Ntotal, Ntotal] muMat; 
  vector[Ntotal] alpha;// shape parameter
  vector[Ntotal] beta;// rate parameter
  // add this random jitter to the population intercept
  vector[Nclusters] rBetaRand = rGroupsJitter + betas[1];
  # convert these betas to a matrix form
  vector[Ntotal] rNewIntercept = rBetaRand[NgroupMap];
  matrix[Ntotal, Ncol] mBetas;
  // fill this matrix
  for (i in 1:Ntotal) mBetas[i,] = betas;
  mBetas[,1] = rNewIntercept;
  mBetas = (mBetas)'; // transpose the matrix
  muMat = X * mBetas; 
  # get the diagonal of this matrix
  for(i in 1:Ntotal) mu[i] = muMat[i,i];
  mu = exp(mu);
  // convert this fitted value to the gamma parameters
  alpha = mu .* mu / sigmaPop; 
  beta = mu/sigmaPop; 
}
model {
  // using non-informative priors to start with
  sigmaPop ~ uniform(0, 1e10);
  sigmaRan ~ uniform(0, 1e3);
  // random effects sample
  rGroupsJitter ~ normal(0, sigmaRan);
  
  // likelihood function
  y ~ gamma(alpha, beta);
}