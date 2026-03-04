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
plot_trajectories(res$thetas, res$epsilons, params)


###2) Power curves

res <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, p_eps)
diffs <- res$diffs
threshold <- quantile(diffs, 0.9) #gaussian threshold = 0.324, binomial threshold = 0.264
threshold


 ###

params_gradual <- seq(0.6, 0.9, length = 20)
params_jump    <- c(rep(0.6, 19), 0.9)

lambdas <- 1:50 * 4


threshold <- 0.264
out <- plot_power_curves_binomial(
  params_gradual,
  params_jump,
  label1 = "Gradual drift",
  label2 = "Abrupt drift (last version)",
  lambdas,
  threshold,
  M, B,
  a_theta, b_theta,
  p_eps
)



threshold <- 0.324
out <- plot_power_curves_gaussian(
  params_gradual,
  params_jump,
  label1 = "Gradual drift",
  label2 = "Abrupt drift",
  lambdas,
  threshold,
  M, B,
  mean_theta, sd_theta,
  sigma,
  p_eps
)

###3) Variance, bias, MSE

