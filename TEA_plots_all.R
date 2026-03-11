remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_functions.R")


M <- 100
B <- 1000
lambda <- 9


params_gradual <- seq(0.6, 0.9, length = 20)   # true means per version
params_abrupt <- c(rep(0.6, 19), 0.9)
params_drift <- c(seq(0.6, 0.9, length = 10), seq(0.9, 0.6, length = 10))
params_stable <- rep(0.7, 20)

a_theta <- 1
b_theta <- 1
mean_theta <- 0.5                          # prior mean
sd_theta   <- 0.5                          # prior sd
sigma      <- 0.5                          # measurement noise
p_eps <- 0.5

###1) Trajectories

par <- params_drift
res <- posterior_sim_binomial(par, M, B, lambda, a_theta, b_theta, p_eps)
plot_trajectories(res$thetas, res$epsilons, par)

res <- posterior_sim_gaussian(par, M, B, lambda, mean_theta, sd_theta, sigma, p_eps)
plot_trajectories(res$thetas, res$epsilons, par)


par(mfrow = c(1,1))
par <- params_gradual
res <- posterior_sim_binomial(par, M, B, lambda, a_theta, b_theta, p_eps)
plot_trajectories(res$thetas, res$epsilons, par)
par <- params_abrupt
res <- posterior_sim_binomial(par, M, B, lambda, a_theta, b_theta, p_eps)
plot_trajectories(res$thetas, res$epsilons, par)
par <- params_drift
res <- posterior_sim_binomial(par, M, B, lambda, a_theta, b_theta, p_eps)
plot_trajectories(res$thetas, res$epsilons, par)
par <- params_stable
res <- posterior_sim_binomial(par, M, B, lambda, a_theta, b_theta, p_eps)
plot_trajectories(res$thetas, res$epsilons, par)

###2) Power curves

#results from 'TEA_power_curves.R'
source("TEA_functions.R")
lambdas <- 1:50 * 4

# params_gradual <- seq(0.6, 0.9, length = 20)
# params_jump    <- c(rep(0.6, 19), 0.9)
# params_drift <- c(seq(0.6, 0.9, length = 10), seq(0.9, 0.6, length = 10))
# params_stable <- rep(0.7, 20)

# lambda control = 19 (for calculating threshold, type-1 error at 0.1)
# "time elapsed =50.3986024498939"
# "M =100"
# "B = 1000"
# "a_theta =1"
# "b_theta =1"
# "mean_theta =0.5"
# "sd_theta =0.5"
# "sigma =0.5"
# "p_eps =0.5"
# "threshold_binomial <-0.200181073877029"
p_binomial_gradual <- c( 0.598, 0.673, 0.685, 0.733, 0.767, 0.766, 0.792, 0.804, 0.823, 0.848, 0.847, 0.881, 0.861, 0.884, 0.892, 0.883, 0.904, 0.921, 0.904, 0.918, 0.927, 0.924, 0.94, 0.942, 0.95, 0.948, 0.943, 0.961, 0.963, 0.96, 0.963, 0.945, 0.969, 0.967, 0.964, 0.967, 0.975, 0.979, 0.971, 0.974, 0.983, 0.985, 0.977, 0.984, 0.979, 0.983, 0.99, 0.983, 0.99, 0.986 )
p_binomial_jump <- c( 0.476, 0.529, 0.552, 0.607, 0.653, 0.674, 0.706, 0.757, 0.78, 0.816, 0.824, 0.853, 0.849, 0.867, 0.874, 0.895, 0.894, 0.91, 0.929, 0.931, 0.932, 0.931, 0.93, 0.936, 0.932, 0.945, 0.936, 0.954, 0.958, 0.965, 0.966, 0.964, 0.972, 0.968, 0.976, 0.984, 0.967, 0.965, 0.977, 0.98, 0.982, 0.98, 0.982, 0.986, 0.99, 0.984, 0.989, 0.983, 0.99, 0.99 )
p_binomial_drift <- c( 0.357, 0.256, 0.196, 0.16, 0.116, 0.09, 0.08, 0.058, 0.073, 0.043, 0.045, 0.033, 0.015, 0.017, 0.013, 0.01, 0.015, 0.008, 0.009, 0.004, 0.003, 0.005, 0.004, 0.004, 0.001, 0.006, 0.002, 0, 0.001, 0, 0, 0, 0, 0.001, 0, 0, 0, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
p_binomial_stable <- c( 0.284, 0.231, 0.158, 0.115, 0.092, 0.068, 0.06, 0.037, 0.029, 0.017, 0.012, 0.011, 0.017, 0.01, 0.005, 0.006, 0.005, 0.002, 0.002, 0.001, 0.002, 0.001, 0, 0.001, 0.001, 0, 0, 0.001, 0, 0, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
# "threshold_gaussian <-0.211111187737812"
p_gaussian_gradual <- c( 0.639, 0.649, 0.691, 0.697, 0.715, 0.723, 0.725, 0.756, 0.777, 0.793, 0.776, 0.796, 0.797, 0.82, 0.817, 0.828, 0.837, 0.838, 0.855, 0.835, 0.852, 0.877, 0.862, 0.88, 0.881, 0.89, 0.874, 0.9, 0.888, 0.896, 0.899, 0.888, 0.894, 0.922, 0.899, 0.916, 0.914, 0.914, 0.937, 0.921, 0.927, 0.927, 0.932, 0.948, 0.937, 0.956, 0.937, 0.952, 0.951, 0.949 )
p_gaussian_jump <- c( 0.52, 0.512, 0.549, 0.551, 0.599, 0.617, 0.623, 0.668, 0.652, 0.682, 0.717, 0.715, 0.709, 0.76, 0.768, 0.786, 0.804, 0.824, 0.847, 0.838, 0.855, 0.878, 0.844, 0.859, 0.866, 0.872, 0.872, 0.884, 0.892, 0.884, 0.895, 0.904, 0.898, 0.899, 0.93, 0.915, 0.922, 0.931, 0.919, 0.918, 0.923, 0.929, 0.927, 0.936, 0.938, 0.935, 0.933, 0.937, 0.924, 0.943 )
p_gaussian_drift <- c( 0.393, 0.31, 0.201, 0.169, 0.133, 0.094, 0.08, 0.068, 0.049, 0.04, 0.03, 0.021, 0.023, 0.009, 0.009, 0.015, 0.006, 0.005, 0.003, 0.005, 0.002, 0.003, 0.001, 0.002, 0.003, 0.002, 0, 0.003, 0.001, 0, 0.002, 0.001, 0, 0, 0, 0.001, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
p_gaussian_stable <- c( 0.407, 0.273, 0.205, 0.157, 0.111, 0.068, 0.066, 0.046, 0.033, 0.029, 0.018, 0.009, 0.019, 0.007, 0.008, 0.005, 0.008, 0.006, 0.003, 0.003, 0, 0.001, 0.002, 0, 0.001, 0.001, 0, 0.002, 0, 0, 0, 0, 0, 0, 0, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )



 

plot_power_curves(lambdas, cbind(p_binomial_gradual, p_binomial_jump), c("gradual", "abrupt"), 1)
plot_power_curves(lambdas, cbind(p_binomial_drift, p_binomial_stable), c("drift", "stable"), 0)
plot_power_curves(lambdas, cbind(p_gaussian_gradual, p_gaussian_jump), c("gradual", "abrupt"), 1)
plot_power_curves(lambdas, cbind(p_gaussian_drift, p_gaussian_stable), c("drift", "stable"), 0)


###3) Variance, bias, MSE


remove(list = ls())

source("TEA_functions.R")

params <- seq(0.6, 0.9, length = 20)
M <- 500
B <- 5000
lambdas <- 1:10*10
a_theta <- 1
b_theta <- 1
p_eps <- 0.5


stats <- TEA_eval_binomial(params, M, B, lambdas, a_theta, b_theta, p_eps)

tea_values <- stats$variance_TEA
naive_values <- stats$variance_naive
stat_name <- "variance"
plot_model_comparison(lambdas, tea_values, naive_values, stat_name)


tea_values <- stats$bias_TEA
naive_values <- stats$bias_naive
stat_name <- "bias"
plot_model_comparison(lambdas, tea_values, naive_values, stat_name)

tea_values <- stats$MSE_TEA
naive_values <- stats$MSE_naive
stat_name <- "MSE"
plot_model_comparison(lambdas, tea_values, naive_values, stat_name)

###4)Bayesian decision

