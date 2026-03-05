remove(list = ls())
#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("TEA_functions.R")

start_time <- Sys.time()

M <- 500
B <- 1000
lambda <- 49

params0 <- rep(0.7, 20)

#params <- rep(0.7, 20)
#params <- c(seq(0.6, 0.9, length = 10), seq(0.9, 0.6, length = 10))
#params <- seq(0.6, 0.9, length = 20)   # true means per version
#params <- c(rep(0.6, 19), 0.9)

a_theta <- 1
b_theta <- 1
mean_theta <- 0.5                          # prior mean
sd_theta   <- 0.5                          # prior sd
sigma      <- 0.5                          # measurement noise
p_eps <- 0.5

res <- posterior_sim_binomial(params0, M, B, lambda, a_theta, b_theta, p_eps)
diffs <- res$diffs
threshold_binomial <- quantile(diffs, 0.9) #gaussian threshold = 0.324, binomial threshold = 0.264


res <- posterior_sim_gaussian(params0, M, B, lambda, mean_theta, sd_theta, sigma, p_eps)
diffs <- res$diffs
threshold_gaussian <- quantile(diffs, 0.9)


threshold_binomial
threshold_gaussian
###

params_gradual <- seq(0.6, 0.9, length = 20)
params_jump    <- c(rep(0.6, 19), 0.9)
params_drift <- c(seq(0.6, 0.9, length = 10), seq(0.9, 0.6, length = 10))
params_stable <- rep(0.7, 20)

params_list <- list(params_gradual, params_jump, params_drift, params_stable)

lambdas <- 1:50 * 4

out_binomial <- compute_power_binomial(
  params_list,
  lambdas,
  threshold_binomial,
  M, B,
  a_theta, b_theta,
  p_eps
)


p_binomial_gradual <- out_binomial[[1]]
p_binomial_jump <- out_binomial[[2]]
p_binomial_drift <- out_binomial[[3]]
p_binomial_stable <- out_binomial[[4]]




out_gaussian <- compute_power_gaussian(
  params_list,
  lambdas,
  threshold_gaussian,
  M, B,
  mean_theta, sd_theta,
  sigma,
  p_eps
)


p_gaussian_gradual <- out_gaussian[[1]]
p_gaussian_jump <- out_gaussian[[2]]
p_gaussian_drift <- out_gaussian[[3]]
p_gaussian_stable <- out_gaussian[[4]]

end_time <- Sys.time()

paste0("time elapsed =", end_time - start_time)
paste0("M =", M)
paste0("B = ", B)
paste0("a_theta =", a_theta)
paste0("b_theta =", b_theta)
paste0("mean_theta =", mean_theta)
paste0("sd_theta =", sd_theta)
paste0("sigma =", sigma)
paste0("p_eps =", p_eps)


paste0("threshold_binomial <-", threshold_binomial)
cat("p_binomial_gradual <- c(", paste(p_binomial_gradual, collapse = ", "), ")\n")
cat("p_binomial_jump <- c(", paste(p_binomial_jump, collapse = ", "), ")\n")
cat("p_binomial_drift <- c(", paste(p_binomial_drift, collapse = ", "), ")\n")
cat("p_binomial_stable <- c(", paste(p_binomial_stable, collapse = ", "), ")\n")

paste0("threshold_gaussian <-", threshold_gaussian)
cat("p_gaussian_gradual <- c(", paste(p_gaussian_gradual, collapse = ", "), ")\n")
cat("p_gaussian_jump <- c(", paste(p_gaussian_jump, collapse = ", "), ")\n")
cat("p_gaussian_drift <- c(", paste(p_gaussian_drift, collapse = ", "), ")\n")
cat("p_gaussian_stable <- c(", paste(p_gaussian_stable, collapse = ", "), ")\n")
