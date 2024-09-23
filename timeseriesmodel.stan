data {
  int T;
  array[T] int y;
  vector[T] law;
  vector[T] kms;
  vector[T] time;
}
transformed data{
  vector[T] log_y;
  for(i in 1:T){
    log_y[i] = log(y[i]);
  }
}
parameters {
  // intercept
  real alpha;
  
  // AR parameters
  real beta_lag1;
  real beta_lag2;
  real beta_lag12;
  
  // trend
  real beta_time;
  
  // covariates
  real beta_law;
  real beta_kms;
}
model {
  // priors
  alpha ~ normal(5.7, 10);
  beta_lag1 ~ normal(0, 0.5);
  beta_lag2 ~ normal(0, 0.5);
  beta_lag12 ~ normal(0, 0.5);
  beta_law~ normal(0, 0.5);
  
  // likelihood
  for(t in 13:T){
    y[t] ~ poisson_log(alpha+beta_law*law[t]+beta_kms*kms[t]+beta_lag1*log_y[t-1]+beta_lag2*log_y[t-2]+beta_lag12*log_y[t-12] + beta_time*time);
  }
}
