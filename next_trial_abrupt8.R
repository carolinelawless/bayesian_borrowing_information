remove(list = ls())
setwd("~/Documents/travail/dmd_when_next_clinical_trial?/bayesian_borrowing_information")
source("functions8.R")

t1 <- Sys.time()

M <- 100
B <- 500

n_versions <- 10
lambda <- 20
lambdas <- rep(lambda - 1, n_versions)
epsilon_const <- 0.5
a_theta <- b_theta <- 0.5
thres <- 0.2 #0.1

params <- rep(0.7, n_versions)
#params <- seq(0.6, 0.9, length = n_versions)
#params <- c(rep(0.6, n_versions - 1), 0.9)

J <- length(params)
true_diff <- params[J] - params[1]



proba_TEA_0 <- vector()
proba_TEA_0.5 <- vector()
proba_TEA_1 <- vector()
proba_TEA_EB_0.5 <- vector()
proba_TEA_EB_1 <- vector()


lambda_vector <- (1:25)*2
for(lambda in lambda_vector){
  lambdas <- rep(lambda - 1, n_versions)
  
  CSD <- 1 #TEA without EB
  epsilon_const <- 0
  res<- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_const, CSD)
  diffs <- res$diffs
  proba_new_trial_TEA_0 <- sum(diffs > thres)/B

  epsilon_const <- 0.5
  res<- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_const, CSD)
  diffs <- res$diffs
  proba_new_trial_TEA_0.5 <- sum(diffs > thres)/B

  epsilon_const <- 1
  res<- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_const, CSD)
  diffs <- res$diffs
  proba_new_trial_TEA_1 <- sum(diffs > thres)/B
  
  
  CSD <- 0.2 #TEA with EB
  epsilon_const <- 0
  res<- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_const, CSD)
  diffs <- res$diffs
  proba_new_trial_TEA_EB_0 <- sum(diffs > thres)/B
  
  epsilon_const <- 0.5
  res<- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_const, CSD)
  diffs <- res$diffs
  proba_new_trial_TEA_EB_0.5 <- sum(diffs > thres)/B
  
  epsilon_const <- 1
  res<- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_const, CSD)
  diffs <- res$diffs
  proba_new_trial_TEA_EB_1 <- sum(diffs > thres)/B
  
  proba_TEA_0 <- c(proba_TEA_0, proba_new_trial_TEA_0)
  proba_TEA_0.5 <- c(proba_TEA_0.5, proba_new_trial_TEA_0.5)
  proba_TEA_1 <- c(proba_TEA_1, proba_new_trial_TEA_1)
  proba_TEA_EB_0.5 <- c(proba_TEA_EB_0.5, proba_new_trial_TEA_EB_0.5)
  proba_TEA_EB_1 <- c(proba_TEA_EB_1, proba_new_trial_TEA_EB_1)
  
}

t2 <- sys.time()



