remove(list = ls())
source("functions7.R")



M <- 500
B <- 100
lambda <- 9


params <- seq(0.6, 0.9, length = 10)
params <- c(rep(0.6, 9), 0.9)
params <- rep(0.7, 10)
CSD <- 0.2
epsilon_const <- 1
a_theta <- b_theta <- 0.5
#plot(density(rbeta(1e4, a_theta, b_theta)))
res_EB <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_const, CSD)

B <- 40

diff_mean <- vector()
params <- rep(0.7, 10)
for(b in 1:B){
  res_EB<- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_const, CSD)
  diff_mean[b]<- mean(res_EB$thetas[[10]]) - mean(res_EB$thetas[[1]])
  print(b)
}

alpha <- 0.05
thres_upper <- quantile(diff_mean, 1 - alpha/2)
thres_lower <- quantile(diff_mean, alpha/2)
thres <- quantile(diff_mean, 1 - alpha)
thres
length(which(diff_mean > thres))/length(diff_mean)
