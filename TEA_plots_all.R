remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_functions.R")


M <- 100
B <- 1000
lambda <- 9

params <- rep(0.7, 20)
#params <- c(seq(0.6, 0.9, length = 10), seq(0.9, 0.6, length = 10))
#params <- seq(0.6, 0.9, length = 20)   # true means per version
#params <- c(rep(0.6, 19), 0.9)

a_theta <- 1
b_theta <- 1
mean_theta <- 0.5                          # prior mean
sd_theta   <- 0.5                          # prior sd
sigma      <- 0.5                          # measurement noise
p_eps <- 0.5

###1) Trajectories

# Suppose you ran posterior_sim_gaussian()
res <- posterior_sim_gaussian(params, M, B, lambda, mean_theta, sd_theta, sigma, p_eps)
plot_trajectories(res$thetas, res$epsilons, params)

res <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, p_eps)
plot_trajectories(res$thetas, res$epsilons, params, 0)


###2) Power curves

#results from 'TEA_power_curves.R'

lambdas <- 1:50 * 4

params_gradual <- seq(0.6, 0.9, length = 20)
params_jump    <- c(rep(0.6, 19), 0.9)
params_drift <- c(seq(0.6, 0.9, length = 10), seq(0.9, 0.6, length = 10))
params_stable <- rep(0.7, 20)

params_list <- list(params_gradual, params_jump, params_drift, params_stable)

p_binomial_gradual <- c( 0.34, 0.54, 0.54, 0.46, 0.52, 0.62, 0.6, 0.64, 0.58, 0.54, 0.54, 0.6, 0.66, 0.68, 0.56, 0.54, 0.56, 0.5, 0.56, 0.56, 0.54, 0.5, 0.78, 0.5, 0.62, 0.44, 0.62, 0.52, 0.6, 0.56, 0.58, 0.62, 0.62, 0.64, 0.7, 0.6, 0.62, 0.68, 0.52, 0.8, 0.62, 0.68, 0.56, 0.66, 0.74, 0.64, 0.66, 0.84, 0.68, 0.8 )
p_binomial_jump <- c( 0.32, 0.34, 0.42, 0.44, 0.32, 0.46, 0.56, 0.52, 0.52, 0.58, 0.68, 0.5, 0.56, 0.58, 0.7, 0.66, 0.58, 0.62, 0.48, 0.58, 0.72, 0.56, 0.8, 0.48, 0.76, 0.76, 0.54, 0.6, 0.68, 0.62, 0.64, 0.68, 0.66, 0.72, 0.68, 0.68, 0.62, 0.7, 0.64, 0.68, 0.66, 0.66, 0.6, 0.78, 0.72, 0.66, 0.82, 0.76, 0.76, 0.72 )
p_binomial_drift <- c( 0.2, 0.16, 0.06, 0.04, 0.02, 0.02, 0, 0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
p_binomial_stable <- c( 0.18, 0.14, 0.14, 0.1, 0, 0.04, 0.02, 0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
plot_power_curves(lambdas, cbind(p_binomial_gradual, p_binomial_jump), c("gradual", "abrupt"), 1)
plot_power_curves(lambdas, cbind(p_binomial_drift, p_binomial_stable), c("drift", "stable"), 0)


p_gaussian_gradual <- c( 0.5, 0.8, 0.9, 1, 1, 1, 1, 0.9, 1, 1, 1, 1, 1, 1, 0.9, 0.9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 )
p_gaussian_jump <- c( 0.8, 0.9, 0.8, 1, 0.8, 0.8, 0.9, 0.9, 0.9, 0.9, 0.7, 0.9, 1, 1, 0.9, 1, 0.8, 0.7, 0.9, 1, 1, 0.8, 1, 0.9, 1, 1, 0.8, 1, 1, 1, 1, 0.9, 0.8, 0.9, 0.8, 1, 0.9, 0.9, 1, 1, 0.8, 1, 1, 0.7, 0.9, 1, 0.9, 1, 1, 0.8 )
p_gaussian_drift <- c( 0.8, 0.7, 0.5, 0.5, 0.9, 0.4, 0.5, 0.6, 0.4, 0.4, 0, 0.5, 0.6, 0.5, 0.7, 0.1, 0.2, 0.4, 0.3, 0.1, 0.1, 0.4, 0.3, 0.4, 0.1, 0.5, 0, 0.2, 0, 0.2, 0.3, 0.4, 0.3, 0.2, 0.4, 0, 0, 0.1, 0.1, 0.3, 0.2, 0.3, 0.1, 0.3, 0.2, 0, 0.3, 0.5, 0.1, 0 )
p_gaussian_stable <- c( 0.4, 0.7, 0.8, 0.6, 0.5, 0.8, 0.4, 0.4, 0.3, 0.3, 0.5, 0.4, 0.4, 0.3, 0.3, 0.5, 0, 0.1, 0.3, 0.5, 0.2, 0.2, 0.4, 0.3, 0.2, 0.2, 0.2, 0.1, 0.3, 0.1, 0, 0, 0.1, 0.2, 0.3, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.2, 0, 0.3, 0, 0.1, 0.1, 0.3, 0.2 )
 
p_gaussian_gradual <- c( 0.4, 0.68, 0.42, 0.4, 0.56, 0.4, 0.54, 0.34, 0.52, 0.52, 0.42, 0.54, 0.34, 0.46, 0.38, 0.44, 0.44, 0.38, 0.44, 0.42, 0.4, 0.46, 0.42, 0.36, 0.36, 0.34, 0.5, 0.5, 0.24, 0.38, 0.48, 0.22, 0.48, 0.36, 0.46, 0.36, 0.42, 0.46, 0.38, 0.44, 0.3, 0.32, 0.3, 0.34, 0.36, 0.26, 0.48, 0.34, 0.32, 0.22 )
p_gaussian_jump <- c( 0.3, 0.36, 0.28, 0.32, 0.32, 0.26, 0.4, 0.34, 0.44, 0.52, 0.36, 0.52, 0.36, 0.48, 0.46, 0.44, 0.38, 0.54, 0.58, 0.48, 0.42, 0.44, 0.52, 0.5, 0.44, 0.46, 0.46, 0.48, 0.58, 0.48, 0.46, 0.46, 0.42, 0.4, 0.46, 0.44, 0.46, 0.46, 0.42, 0.6, 0.56, 0.38, 0.42, 0.34, 0.6, 0.44, 0.58, 0.5, 0.58, 0.58 )
p_gaussian_drift <- c( 0.22, 0.12, 0.08, 0.02, 0.06, 0.02, 0.02, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
p_gaussian_stable <- c( 0.26, 0.14, 0.12, 0.06, 0.02, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
plot_power_curves(lambdas, cbind(p_gaussian_gradual, p_gaussian_jump), c("gradual", "abrupt"), 1)
plot_power_curves(lambdas, cbind(p_gaussian_drift, p_gaussian_stable), c("drift", "stable"), 0)


###3) Variance, bias, MSE

