data {
  int T;
  array[T] int y;
  vector[T] law;
  vector[T] dist;
  vector[T] time;
}
transformed data{
  vector[T] log_y;
  vector[T] log_time;
  for(i in 1:T){
    log_y[i] = log(y[i]+1);
    log_time[i] = log(time[i]);
  }
}
parameters {
  // intercept
  real alpha;
  
  // AR parameters
  real gamma_lag1;
  real gamma_lag2;
  real gamma_lag12;
  
  // covariates
  real beta_law;
  real beta_dist;
  
  // trend component
  real delta_trend;
  
  // negbin dispersion parameter
  real<lower=0> reciprocal_phi;
}
transformed parameters{
  real<lower=0> phi;
  phi = 1/reciprocal_phi;
  
}
model {
  // priors
  alpha ~ normal(5.11, 10);
  beta_law~ normal(0, 0.5);
  beta_dist~ normal(0, 0.5);
  gamma_lag1 ~ normal(0,0.5);
  gamma_lag2 ~ normal(0,0.5);
  gamma_lag12 ~ normal(0,0.5);
  delta_trend ~ normal(0,10);
  
  // likelihood
  for(t in 13:T){
    y[t] ~ neg_binomial_2_log(alpha+beta_law*law[t]+beta_dist*dist[t]+gamma_lag1*log_y[t-1]+gamma_lag2*log_y[t-2]+gamma_lag12*log_y[t-12]+delta_trend*log_time[t], phi);
  }
}
generated quantities {
 vector[T-12] loglik;
  for(t in 13:T){
    loglik[t-12] = neg_binomial_2_log_lpmf(y[t]|alpha+beta_law*law[t]+beta_dist*dist[t]+gamma_lag1*log_y[t-1]+gamma_lag2*log_y[t-2]+gamma_lag12*log_y[t-12]+delta_trend*log_time[t], phi); 
  }
}