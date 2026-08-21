remove(list = ls())
setwd("~/Documents/travail/dmd_when_next_clinical_trial?/bayesian_borrowing_information")
#setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("functions9.R")




M <- 100
B <- 500

n_versions <- 5

a_theta <- b_theta <- 0.5


#plot(density(rbeta(1e5, a_theta, b_theta)))

thres1 <- 0.1
thres2 <- 0.8

#scenario <- "stable"
#scenario <- "gradual"
scenario <- "abrupt"


if(scenario == "stable"){
  params <- rep(0.7, n_versions)
}else if(scenario == "gradual"){
  params <- seq(0.6, 0.9, length = n_versions)
}else if(scenario == "abrupt"){
  params <- c(rep(0.6, n_versions - 1), 0.9)
}

lambda <- 100
lambdas <- rep(lambda - 1, n_versions)
CSD <- 0.2 #TEA with EB
epsilon_const <- 0.5
res<- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_const, CSD, thres1)
plot_estimates(res$thetas, res$epsilons, params, lambda)


