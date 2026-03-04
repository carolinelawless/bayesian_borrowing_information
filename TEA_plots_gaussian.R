remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_functions_gaussian.R")


M <- 100
B <- 1000
lambda <- 9

params <- rep(0.7, 20)
params <- c(seq(0.6, 0.9, length = 10), seq(0.9, 0.6, length = 10))
params <- seq(0.6, 0.9, length = 20)   # true means per version
params <- c(rep(0.6, 19), 0.9)

mean_theta <- 0.5                          # prior mean
sd_theta   <- 0.5                          # prior sd
sigma      <- 0.5                          # measurement noise
p_eps <- 0.5

# Suppose you ran posterior_sim_gaussian()
res <- posterior_sim_gaussian(params, M, B, lambda, mean_theta, sd_theta, sigma, p_eps)

plot_gaussian_trajectories(res$thetas, res$epsilons, true_theta = params)

