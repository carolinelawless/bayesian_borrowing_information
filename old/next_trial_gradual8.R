remove(list = ls())
setwd("~/Documents/travail/dmd_when_next_clinical_trial?/bayesian_borrowing_information")
source("functions8.R")

n_versions <- 200

M <- 100
B <- 2000
#lambda <- 9
lambda <- 20
lambdas <- rep(lambda - 1, n_versions)
#CSD <- 0.2
CSD <- 1
epsilon_const <- 1
a_theta <- b_theta <- 0.5
thres <- 0.1

#params <- rep(0.7, n_versions)
#params <- c(rep(0.6, n_versions - 1), 0.9)
params <- seq(0.6, 0.9, length = n_versions)


proba_new_trial <- 0
res_EB<- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_const, CSD)



plot_trajectories(res_EB$thetas, res_EB$epsilons, params, epsilon_type = "TEA")

J <- length(params)
true_diff <- params[J] - params[1]
bias <- mean(res_EB$diffs) - true_diff
variance <- var(res_EB$diffs)
round(bias, 3)
round(variance, 3)



j <- 1
while(proba_new_trial < 0.8 & j < J){
  j <- j + 1



  diff_mean <- res_EB$thetas[[j]] - res_EB$thetas[[1]]
  
  proba_new_trial <- sum(diff_mean > thres)/B

}

print(paste0("bias =", round(bias, 3)))
print(paste0("variance =", round(variance, 3)))
print(paste0("delta stop =", round(params[j] - params[1], 3)))


###

lambda <- 9
lambdas <- rep(lambda - 1, n_versions)

res_EB<- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_const, CSD)
plot_trajectories(res_EB$thetas, res_EB$epsilons, params, epsilon_type = "TEA")
diff_mean <- res_EB$thetas[[j]] - res_EB$thetas[[1]]

proba_new_trial <- sum(diff_mean > thres)/B
proba_new_trial
###
