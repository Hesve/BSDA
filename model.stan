data {
  int N;
  array[N] int y;
  vector[N] law;
  vector[N] kms;
  vector [N] Februari;
  vector [N] Mars;
  vector [N] April;
  vector [N] May;
  vector [N] Juli;
  vector [N] June;
  vector [N] August;
  vector [N] September;
  vector [N] October;
  vector [N] November;
  vector [N] December;
  vector[N] time;
}
parameters {
  real alpha;
  real beta_law;
  real beta_kms;
  real beta_Februari;
  real beta_Mars;
  real beta_April;
  real beta_Juli;
  real beta_June;
  real beta_August;
  real beta_September;
  real beta_October;
  real beta_November;
  real beta_December;
  real trend;
}
transformed parameters{
  vector[N] lambda = exp(alpha + beta_law*law + beta_kms*kms + beta_Februari* Februari + beta_Mars * Mars
  + beta_April * April + beta_Juli * Juli + beta_June* June + beta_August * August + beta_September * September  + beta_October *October+ beta_November * November + beta_December * December + trend*time);
}
model {
  y ~ poisson(lambda);
}
