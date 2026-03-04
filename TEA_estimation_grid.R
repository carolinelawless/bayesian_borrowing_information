##

remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_functions.R")
source("TEA_scenarios.R")


M <- 1e2
B <- 1e3 #number of estimates

a_theta <- 1
b_theta <- 1
p_eps <- 0.5


#params <- scenarios[[8]] 

params <- rep(0.5, 2)
params <- rep(0.5, 5)
params <- scenarios[[25]]


K <- length(params)

n <- rep(10, K)
#n[2] <- 100



res <- posterior_sim(params, M, B)
thetas <- res[[1]]
epsilons <- res[[2]]
