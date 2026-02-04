remove(list = ls())


setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
source("TEA_functions.R")
source("TEA_scenarios.R")



M <- 1e2

a_theta <- 1
b_theta <- 1
p_eps <- 0.5

B <- 1e3 #number of estimates

n_versions <- 5
params_grad <- seq(0.1, 0.5, length.out = n_versions)
params_sudden <- c(rep(0.1, n_versions -1), 0.5)
params_constant <- rep(0.1, n_versions)
params <- params_constant
params


x <- 1:50

K <- length(params)

bias_naive_vec <- vector()
bias_tea_vec <- vector()
var_naive_vec <- vector()
var_tea_vec <- vector()
mse_naive_vec <- vector()
mse_tea_vec <- vector()

lambda <- 4
tea <- tea_diff_mean(params, lambda, M, B)
naive <- naive_diff(params, lambda, M, B)


for(lambda in x){
  print(lambda)
  tea <- tea_diff_mean(params, lambda, M, B)
  naive <- naive_diff(params, lambda, M, B)
  
  bias_naive <- abs(mean(naive) - (params[K] - params[1]))
  bias_tea <- abs(mean(tea) - (params[K] - params[1]))
  var_naive <- var(naive)
  var_tea <- var(tea)
  mse_naive <- bias_naive*bias_naive + var_naive
  mse_tea <- bias_tea*bias_tea + var_tea
 
  bias_naive_vec <- c(bias_naive_vec, bias_naive)
  bias_tea_vec <- c(bias_tea_vec, bias_tea)
  var_naive_vec <- c(var_naive_vec, var_naive)
  var_tea_vec <- c(var_tea_vec, var_tea)
  mse_naive_vec <- c(mse_naive_vec, mse_naive)
  mse_tea_vec <- c(mse_tea_vec, mse_tea)
}


cols <- c("#4E79A7", "#59A14F", "#9C755F", "#B07AA1", "#76B7B2")


ylim_var <- range(c(var_naive_vec, var_tea_vec), na.rm = TRUE)

plot(x, var_naive_vec, type = "l", col = cols[1], lwd = 2,
     xlab = "lambda", ylab = "variance", ylim = ylim_var, main = "")
lines(x, var_tea_vec, col = cols[2], lwd = 2)

legend("topright",
       legend = c("Naive", "TEA"),
       col = c(cols[1], cols[2]),
       lwd = 2)


ylim_bias <- range(c(bias_naive_vec, bias_tea_vec), na.rm = TRUE)

plot(x, bias_naive_vec, type = "l", col = cols[1], lwd = 2,
     xlab = "lambda", ylab = "bias", ylim = ylim_bias, main = "")
lines(x, bias_tea_vec, col = cols[2], lwd = 2)

legend("topright",
       legend = c("Naive", "TEA"),
       col = c(cols[1], cols[2]),
       lwd = 2)


ylim_mse <- range(c(mse_naive_vec, mse_tea_vec), na.rm = TRUE)

plot(x, mse_naive_vec, type = "l", col = cols[1], lwd = 2,
     xlab = "lambda", ylab = "MSE", ylim = ylim_mse, main = "")
lines(x, mse_tea_vec, col = cols[2], lwd = 2)

legend("topright",
       legend = c("Naive", "TEA"),
       col = c(cols[1], cols[2]),
       lwd = 2)




