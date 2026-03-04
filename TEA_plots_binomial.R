remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_functions_binomial.R")


M <- 100
B <- 1000
lambda <- 9

params <- c(0.5, 0.55, 0.6, 0.65, 0.7)   # true means per version

a_theta <- 5
b_theta <- 5
p_eps <- 0.5


res <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, p_eps)

plot_binomial_trajectories(res$thetas, res$epsilons, true_theta = params)
