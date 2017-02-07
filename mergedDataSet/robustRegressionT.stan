data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    real y[Ntotal]; // response variable - t distributed
  }
parameters {
    vector[Ncol] betas; // regression parameters
    real<lower=0> nu; // normality parameter for t distribution or degree of freedom 
    real<lower=0> sigma; // scale parameter  
  }
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = X*betas; // fitted value using identity link
}
model {
  y ~ student_t( nu , mu , sigma);
}
