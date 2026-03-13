remove(list = ls())
#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("functions.R")

start_time <- Sys.time()

M <- 500
B <- 5000
a_theta <- 0.5
b_theta <- 0.5
p_eps <- 0.5
mean_theta <- 0.5
sd_theta <- 0.5
sigma <- 0.5
model <- "binomial"
#model <- "gaussian"


thres1 <- 0.1
thres2 <- 0.8


params <- seq(0.6, 0.9, length = 20)
#params <- rep(0.7, 20)

K <- length(params)


tea_stops <- vector()
naive_stops <- vector()


lambdas <- 1:100/2


stat_tea <- vector(length = length(lambdas))
stat_naive <- vector(length = length(lambdas))



for(i in 1:length(lambdas)){
  lambda <- lambdas[i]
  print(lambda)
  if(model == "binomial"){
    res_tea <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, p_eps)
    res_naive <- posterior_sim_naive_binomial(params, M, B, lambda, a_theta, b_theta)
  }else if(model == "gaussian"){
    res_tea <- posterior_sim_gaussian(params, M, B, lambda, mean_theta, sd_theta, sigma, p_eps)
    res_naive <- posterior_sim_naive_gaussian(params, M, B, lambda, mean_theta, sd_theta, sigma)
  }

  stat_tea[i] <- sum(res_tea$thetas[[K]] - res_tea$thetas[[1]] > thres1)/B
  stat_naive[i] <- sum(res_naive$thetas[[K]] - res_naive$thetas[[1]] > thres1)/B
}


tea_stops <- vector()
naive_stops <- vector()
for(lambda in lambdas){
  print(lambda)
  if(model == "binomial"){
    res_tea <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, p_eps)
    res_naive <- posterior_sim_naive_binomial(params, M, B, lambda, a_theta, b_theta)
  }else if(model == "gaussian"){
    res_tea <- posterior_sim_gaussian(params, M, B, lambda, mean_theta, sd_theta, sigma, p_eps)
    res_naive <- posterior_sim_naive_gaussian(params, M, B, lambda, mean_theta, sd_theta, sigma)
  }
  k_tea <- 0
  stop_tea <- 0
  while(stop_tea == 0 & k_tea < K){
    k_tea <- k_tea + 1
    stat <- sum(res_tea$thetas[[k_tea]] - res_tea$thetas[[1]] > thres1)/B
    if(stat > thres2){
      stop_tea <- 1
    }
  }
  k_naive <- 0
  stop_naive <- 0
  while(stop_naive == 0 & k_naive < K){
    k_naive <- k_naive + 1
    stat <- sum(res_naive$thetas[[k_naive]] - res_naive$thetas[[1]] > thres1)/B
    if(stat > thres2){
      stop_naive <- 1
    }
  }
  
  
  tea_stops <- c(tea_stops, k_tea)
  naive_stops <- c(naive_stops, k_naive)
}



end_time <- Sys.time()

paste0("time elapsed =", end_time - start_time)
paste0("M =", M)
paste0("B = ", B)
paste0("a_theta =", a_theta)
paste0("b_theta =", b_theta)
paste0("p_eps =", p_eps)
paste0("threshold1 = ", thres1)
paste0("threshold2 = ", thres2)

paste0("model = ", model)
cat("params <- c(", paste(params, collapse = ", "), ")\n")

cat("lambdas <- c(", paste(lambdas, collapse = ", "), ")\n")
cat("stat_tea <- c(", paste(stat_tea, collapse = ", "), ")\n")
cat("stat_naive <- c(", paste(stat_naive, collapse = ", "), ")\n")

cat("lambdas <- c(", paste(lambdas, collapse = ", "), ")\n")
cat("naive_stops <- c(", paste(naive_stops, collapse = ", "), ")\n")
cat("tea_stops <- c(", paste(tea_stops, collapse = ", "), ")\n")



