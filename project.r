## importing data
data("Seatbelts")
library(posterior)
library(xtable)
library(tidyverse)
library(forecast)
forecast::ggtsdisplay(Seatbelts[,"DriversKilled"])
# using all cores for stan
options(mc.cores = parallel::detectCores())

# standardizing kms for easier model fit
kms_std <- (Seatbelts[,"kms"] / 1000)

# saving data in stan-format
stan_data <- list(
  y = Seatbelts[,"DriversKilled"], # number of driver deaths each month
  T = nrow(Seatbelts), # number of observations
  law = Seatbelts[,"law"], # binary law variable of interest
  dist = kms_std, # standardized kms
  time = seq(1:192) # monthly index for trend component
)

# fitting poisson model
poisson_model <- rstan::stan_model(file = "./poisson_model.stan")

poisson_fit <- rstan::sampling(poisson_model, data = stan_data,
                         iter = 4000, seed = 1402)

poisson_res <- summarise_draws(poisson_fit)
poisson_res[1:7,c(1,2, 6,7,8)] %>% 
  xtable(caption = "Posterioir statistics for the poisson model", label="poisson_results") %>% 
  print(caption.placement="top")

# fitting negative binomial model
negbin_model <- rstan::stan_model(file = "./negbin_model.stan")
negbin_fit <- rstan::sampling(negbin_model, data = stan_data,
                             iter = 4000, seed = 1402)
negbin_res <- summarise_draws(negbin_fit)
negbin_res[1:9,c(1,2, 6,7,8)] %>% 
  xtable(caption = "Posterioir statistics for the negative binomial model", 
         label="negbin_results") %>% 
  print(caption.placement="top")



digits_func <- function(value)  {
  #function to calculate the how many digits with 0 until the first 
  #non-zero digit using regexp.
  #used for MSCE rounding.
  if (value > 1){
    return(0)
  }
  else if (value < 0.001){
    value <- as.character(format(value,scientific=FALSE))
  }
  n_digits <- attr(regexpr("(?<=\\.)0+|$", value, perl = TRUE), "match.length")
  return(n_digits)
}

mean_and_quantile_func <- function(model,variable, probs = 0.9){
  #function to calculate and present the mean and quantiles for a parameter 
  #given a long with the mcse for a given model. 
  data <- extract_variable_matrix(model, variable=variable)
  #this extracts the result for the variable with all chains as columns
  lower <- (1-probs)/2
  upper <- (1-lower)
  mcse_mean <- mcse_mean(data)
  mean_digits <- digits_func(mcse_mean)
  mean <- round(mean(data), digits = mean_digits)
  
  q_5_mcse <- mcse_quantile(data, probs=lower)
  q_5_digits <- digits_func(q_5_mcse)
  q_5 <- round(quantile(data, probs=lower), digits = q_5_digits)
  
  
  q_95_mcse <- mcse_quantile(data, probs=upper)
  q_95_digits <- digits_func(q_95_mcse)
  q_95 <- round(quantile(data, probs=upper), digits = q_95_digits)
  rhat <- round(posterior::rhat(data), digits=3)
  
  results <- matrix(cbind(mean, q_5, q_95, 
                          round(mcse_mean, digits= mean_digits + 2), 
                          round(q_5_mcse, digits = q_5_digits+2),
                          round(q_95_mcse, digits= q_95_digits + 2), rhat), 
                    ncol=7)
  colnames(results) <- c("Mean", paste("q", lower, sep=""), 
                         paste("q", upper, sep=""),"MCSE mean", "MCSE q5",
                         "MCSE q95", "Rhat")
  row_name<- variable
  rownames(results) <- row_name
  return(results)
}

vars = list("alpha", "gamma_lag1", "gamma_lag2", "gamma_lag12", "delta_trend", "beta_law", "beta_dist")

mapply(FUN= mean_and_quantile_func, variable=vars, MoreArgs = list(model = poisson_fit), SIMPLIFY =FALSE)

vars_neg <- list("alpha", "gamma_lag1", "gamma_lag2", "gamma_lag12", "delta_trend","beta_law", "beta_dist", "phi")

mapply(FUN= mean_and_quantile_func, variable=vars_neg, MoreArgs = list(model = negbin_fit), SIMPLIFY =FALSE)

###obs när man sammanställer resultaten i en matris / data frame fyller den ut med nollor på slutet
mean_and_quantile_table <- function(model, variable_list, name_of_model){
  results <- t(mapply(FUN= mean_and_quantile_func, variable=variable_list, MoreArgs = list(model = model), SIMPLIFY = TRUE)) 
  rownames(results) <- unlist(variable_list)
  colnames(results) <- c("Mean", "q5", 
                         "q95","MCSE mean", "MCSE q5",
                         "MCSE q95", "Rhat")
 results %>% 
   xtable(caption = paste("Results for", name_of_model), label= paste(name_of_model, "table"), sep="_") %>% 
   print(caption.placement="top")
}
mean_and_quantile_table(model=poisson_fit, variable_list=vars, name_of_model = "Poisson model")

mean_and_quantile_func(poisson_fit, variable="alpha")
mean_and_quantile_mcse <- function(model, variable, probs=0.9){
 res <- mean_and_quantile_func(model, variable, probs=probs)[,c(4,5,6)]
 return(res)
}
#######################
  mcse_func <- function(model,variable, probs = 0.9){
    #function to calculate and present the mean and quantiles for a parameter 
    #given a long with the mcse for a given model. 
    data <- extract_variable_matrix(model, variable=variable)
    #this extracts the result for the variable with all chains as columns
    lower <- (1-probs)/2
    upper <- (1-lower)
    mcse_mean <- mcse_mean(data)
    mean_digits <- digits_func(mcse_mean)
    mean <- round(mean(data), digits = mean_digits)
    
    q_5_mcse <- mcse_quantile(data, probs=lower)
    q_5_digits <- digits_func(q_5_mcse)
    q_5 <- round(quantile(data, probs=lower), digits = q_5_digits)
    
    
    q_95_mcse <- mcse_quantile(data, probs=upper)
    q_95_digits <- digits_func(q_95_mcse)
    q_95 <- round(quantile(data, probs=upper), digits = q_95_digits)
    rhat <- round(posterior::rhat(data), digits=3)
    
    results <- matrix(cbind(mean, q_5, q_95, 
                            round(mcse_mean, digits= mean_digits + 2), 
                            round(q_5_mcse, digits = q_5_digits+2),
                            round(q_95_mcse, digits= q_95_digits + 2), rhat), 
                      ncol=7)
    colnames(results) <- c("Mean", paste("q", lower, sep=""), 
                           paste("q", upper, sep=""),"MCSE mean", "MCSE q5",
                           "MCSE q95", "Rhat")
    row_name<- variable
    rownames(results) <- row_name
    return(results[,-c(1,2,3,7)])
  }
  
mcse_func(poisson_fit, variable="alpha")  
###############
MCSE_table <- function(model, variable_list, name_of_model){
  results <- t(mapply( FUN= mcse_func, variable = variable_list, MoreArgs = list(model=model), SIMPLIFY = TRUE))
  rownames(results) <- unlist(variable_list)
  colnames(results) <- c("MSCE mean", "MSCE q5", "MSCE q95 ")
  results %>% 
    xtable(caption = paste("Results for", name_of_model), label= paste(name_of_model, "table"), sep="_", digits = c(5,5,5,5)) %>% 
    print(caption.placement="top")
}

MCSE_table(model= poisson_fit, variable_list = vars, name_of_model = "Poisson model MSCE")

MCSE_table(model= negbin_fit, variable_list = vars_neg, name_of_model = "Negative binomial MSCE")

######
######
######

poisson_loglik = extract_log_lik(poisson_fit, parameter_name = "loglik")

poisson_loglik = poisson_loglik[,-c(1:12)]

r_eff_poisson <- relative_eff(exp(poisson_loglik), chain_id = rep(1:4, each= 2000))

poisson_loo = loo(poisson_loglik, r_eff = r_eff_poisson)
plot(poisson_loo, main="")

negbin_loglik = extract_log_lik(negbin_fit, parameter_name = "loglik")

negbin_loglik = poisson_loglik[,-c(1:12)]

r_eff_negbin <- relative_eff(exp(negbin_loglik), chain_id = rep(1:4, each= 2000))

negbin_loo = loo(negbin_loglik, r_eff = r_eff_negbin)

plot(negbin_loo, main="")
