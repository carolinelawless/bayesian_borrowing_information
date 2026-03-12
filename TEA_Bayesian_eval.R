remove(list = ls())
#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("TEA_functions.R")

start_time <- Sys.time()

M <- 500
B <- 1000
a_theta <- 0.5
b_theta <- 0.5
p_eps <- 0.5

source("TEA_functions.R")

thres1 <- 0.1
thres2 <- 0.8

lambda <- 20
params <- seq(0.6, 0.9, length = 20)


K <- length(params)


tea_stops <- vector()
naive_stops <- vector()


lambdas <- 1:50

for(lambda in lambdas){
  print(lambda)
  res_tea <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, p_eps)
  res_naive <- posterior_sim_naive_binomial(params, M, B, lambda, a_theta, b_theta)
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


cat("params <- c(", paste(params, collapse = ", "), ")\n")
cat("lambdas <- c(", paste(lambdas, collapse = ", "), ")\n")
cat("naive_stops <- c(", paste(naive_stops, collapse = ", "), ")\n")
cat("tea_stops <- c(", paste(tea_stops, collapse = ", "), ")\n")




