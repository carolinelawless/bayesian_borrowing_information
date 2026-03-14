remove(list = ls())
setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("functions.R")

params <- seq(0.6, 0.9, length = 20)
#params <- rep(0.7, 20)
M <- 500
B <- 1000
lambdas <- 1:50
a_theta <- 1
b_theta <- 1
mean_theta <- 0.5
sd_theta <- 0.5
sigma <- 0.5
p_eps <- 0.5


start_time <- Sys.time()


stats_binomial <- TEA_eval_binomial(params, M, B, lambdas, a_theta, b_theta, p_eps)
stats_gaussian <- TEA_eval_gaussian(params, M, B, lambdas, mean_theta, sd_theta, sigma, p_eps)


tea_variance_binomial <- stats_binomial$variance_TEA
naive_variance_binomial <- stats_binomial$variance_naive
tea_bias_binomial <- stats_binomial$bias_TEA
naive_bias_binomial <- stats_binomial$bias_naive
tea_mse_binomial <- stats_binomial$MSE_TEA
naive_mse_binomial <- stats_binomial$MSE_naive

tea_variance_gaussian <- stats_gaussian$variance_TEA
naive_variance_gaussian <- stats_gaussian$variance_naive
tea_bias_gaussian <- stats_gaussian$bias_TEA
naive_bias_gaussian <- stats_gaussian$bias_naive
tea_mse_gaussian <- stats_gaussian$MSE_TEA
naive_mse_gaussian <- stats_gaussian$MSE_naive


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



cat("params <- c(", paste(params, collapse = ", "), ")\n")


model <- "binomial"
stat_name <- "bias"
cat(paste0("tea_",stat_name,"_",model)," <- c(", paste(get(paste0("tea_",stat_name,"_",model)), collapse = ", "), ")\n")
stat_name <- "variance"
cat(paste0("tea_",stat_name,"_",model)," <- c(", paste(get(paste0("tea_",stat_name,"_",model)), collapse = ", "), ")\n")
stat_name <- "mse"
cat(paste0("tea_",stat_name,"_",model)," <- c(", paste(get(paste0("tea_",stat_name,"_",model)), collapse = ", "), ")\n")
stat_name <- "bias"
cat(paste0("naive_",stat_name,"_",model)," <- c(", paste(get(paste0("naive_",stat_name,"_",model)), collapse = ", "), ")\n")
stat_name <- "variance"
cat(paste0("naive_",stat_name,"_",model)," <- c(", paste(get(paste0("naive_",stat_name,"_",model)), collapse = ", "), ")\n")
stat_name <- "mse"
cat(paste0("naive_",stat_name,"_",model)," <- c(", paste(get(paste0("naive_",stat_name,"_",model)), collapse = ", "), ")\n")




model <- "gaussian"
stat_name <- "bias"
cat(paste0("tea_",stat_name,"_",model)," <- c(", paste(get(paste0("tea_",stat_name,"_",model)), collapse = ", "), ")\n")
stat_name <- "variance"
cat(paste0("tea_",stat_name,"_",model)," <- c(", paste(get(paste0("tea_",stat_name,"_",model)), collapse = ", "), ")\n")
stat_name <- "mse"
cat(paste0("tea_",stat_name,"_",model)," <- c(", paste(get(paste0("tea_",stat_name,"_",model)), collapse = ", "), ")\n")
stat_name <- "bias"
cat(paste0("naive_",stat_name,"_",model)," <- c(", paste(get(paste0("naive_",stat_name,"_",model)), collapse = ", "), ")\n")
stat_name <- "variance"
cat(paste0("naive_",stat_name,"_",model)," <- c(", paste(get(paste0("naive_",stat_name,"_",model)), collapse = ", "), ")\n")
stat_name <- "mse"
cat(paste0("naive_",stat_name,"_",model)," <- c(", paste(get(paste0("naive_",stat_name,"_",model)), collapse = ", "), ")\n")






