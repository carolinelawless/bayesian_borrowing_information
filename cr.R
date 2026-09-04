remove(list = ls())
setwd("~/Documents/travail/dmd_when_next_clinical_trial?/bayesian_borrowing_information")
#setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("functions.R")




M <- 100
B <- 2000

coverage_interval <- 0.95
thres1 <- 0.1

a_theta <- b_theta <- 0.5
mean_theta <- sd_theta <- sd_obs <- 0.5


binomial_params <- list(
  a_theta = a_theta,
  b_theta = b_theta
)

gaussian_params <- list(
  mean_theta = mean_theta,
  sd_theta = sd_theta,
  sd_obs = sd_obs
)

thres1 <- 0.1


model <- "gaussian"
model_params <- gaussian_params
# model <- "binomial"
# model_params <- binomial_params


#scenario <- "stable"
#scenario <- "gradual"
scenario <- "abrupt"


lambda <- 100
n_versions <- 10
lambdas <- rep(lambda - 1, n_versions)

if(scenario == "stable"){
  params <- rep(0.7, n_versions)
}else if(scenario == "gradual"){
  params <- seq(0.6, 0.9, length = n_versions)
}else if(scenario == "abrupt"){
  params <- c(rep(0.6, n_versions - 1), 0.9)
}


CSD <- 1
epsilon_const <- 0
res <- posterior_sim(params, M, B, lambdas, epsilon_const , CSD, thres1, model, model_params, coverage_interval)
coverage_TEA_0 <- sapply(res$coverage, mean)

epsilon_const <- 0.5
res <- posterior_sim(params, M, B, lambdas, epsilon_const , CSD, thres1, model, model_params, coverage_interval)
coverage_TEA_0.5 <- sapply(res$coverage, mean)

epsilon_const <- 1
res <- posterior_sim(params, M, B, lambdas, epsilon_const , CSD, thres1, model, model_params, coverage_interval)
coverage_TEA_1 <- sapply(res$coverage, mean)

CSD <- 0.2
epsilon_const <- 0.5
res <- posterior_sim(params, M, B, lambdas, epsilon_const , CSD, thres1, model, model_params, coverage_interval)
coverage_TEA_EB_0.5 <- sapply(res$coverage, mean)

epsilon_const <- 1
res <- posterior_sim(params, M, B, lambdas, epsilon_const , CSD, thres1, model, model_params, coverage_interval)
coverage_TEA_EB_1 <- sapply(res$coverage, mean)

plot_coverage(n_versions,
scenario, 
coverage_TEA_0,
coverage_TEA_0.5,
coverage_TEA_1,
coverage_TEA_EB_0.5,
coverage_TEA_EB_1) 






