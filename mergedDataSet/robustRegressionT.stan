data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    real y[Ntotal]; // response variable - t distributed
    real meanY; // the mean value for response
    real sdY; // the scale for response
  }
transformed data {
  real unifLo;
  real unifHi; // the lower and upper bounds for the t-distrubtion scale
  real expLambda; // rate parameter for the exponential distribution for nu
  unifLo = sdY/1000;
  unifHi = sdY*1000;
  expLambda = 1/29.0 ;
}
parameters {
    vector[Ncol] betas; // regression parameters
    real<lower=0.1> nu; // normality parameter for t distribution or degree of freedom 
    real<lower=unifLo, upper=unifHi> sigma; // scale parameter  
  }
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = X*betas; // fitted value using identity link
}
model {
  nu ~ exponential(expLambda);
  y ~ student_t( nu , mu , sigma);
}
