remove(list = ls())
#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("TEA_functions.R")


M <- 50
B <- 10
lambda <- 9

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

lambdas <- 1:50 * 4



out_binomial <- plot_power_curves_binomial(
  params_gradual,
  params_jump,
  label1 = "Gradual drift",
  label2 = "Abrupt drift (last version)",
  lambdas,
  threshold_binomial,
  M, B,
  a_theta, b_theta,
  p_eps
)

p1_binomial <- out_binomial$power1
p2_binomial <- out_binomial$power2


out_gaussian <- plot_power_curves_gaussian(
  params_gradual,
  params_jump,
  label1 = "Gradual drift",
  label2 = "Abrupt drift",
  lambdas,
  threshold_gaussian,
  M, B,
  mean_theta, sd_theta,
  sigma,
  p_eps
)


p1_gaussian <- out_gaussian$power1
p2_gaussian <- out_gaussian$power2


paste0("threshold_binomial <-", threshold_binomial)
cat("p1_binomial <- c(", paste(p1_binomial, collapse = ", "), ")\n")
cat("p2_binomial <- c(", paste(p2_binomial, collapse = ", "), ")\n")

paste0("threshold_gaussian <-", threshold_gaussian)
cat("p1_gaussian <- c(", paste(p1_gaussian, collapse = ", "), ")\n")
cat("p2_gaussian <- c(", paste(p2_gaussian, collapse = ", "), ")\n")
