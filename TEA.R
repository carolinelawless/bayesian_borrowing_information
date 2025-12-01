remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_functions.R")

n <- 1e2
params <- rep(0.001, 5)
V <- sapply(params, function(p) rbinom(n, 1, p))
V <- t(V)   # now V is K x n (row = dataset)
M <- 1e3




# Simulate under H0 (same params for all datasets)
S_H0_1_all <- simulate_S(rep(0.01, 5), n)
S_H0_1 <- S_H0_1_all[[1]]
epsilons_H0_1 <- S_H0_1_all[[2]]
threshold1 <- quantile(S_H0_1, 0.05)  # 5% quantile for 5% false positive rate

S_H0_2_all <- simulate_S(rep(0.5, 5), n)
S_H0_2 <- S_H0_2_all[[1]]
epsilons_H0_2 <- S_H0_2_all[[2]]
threshold2 <- quantile(S_H0_2, 0.05)  # 5% quantile for 5% false positive rate

epsilons_H0_1
epsilons_H0_2

threshold <- (threshold1 + threshold2)/2
threshold <- threshold2 #the most conservative threshold


# Simulate under H1 (all params different)
params_H1_1 <- c(0.5, 0.55, 0.6, 0.65, 0.7)  # all params different
S_H1_1_all <- simulate_S(params_H1_1, n)
S_H1_1 <- S_H1_1_all[[1]]
epsilons_H1_1 <- S_H1_1_all[[2]]
epsilons_H1_1 #0.75, 0.68, 0.65, 0.64
# Check how often S_H1 < threshold to estimate power
length(which(S_H1_1 < threshold))/length(S_H1_1) 

# Simulate under H1 (last param different)
params_H1_2 <- c(0.5, 0.5, 0.5, 0.5, 0.7)  # last param different
S_H1_2_all <- simulate_S(params_H1_2, n)
S_H1_2 <- S_H1_2_all[[1]]
epsilons_H1_2 <- S_H1_2_all[[2]]
epsilons_H1_2 #0.77, 0.78, 0.79, 0.23
# Check how often S_H1 < threshold to estimate power
length(which(S_H1_1 < threshold))/length(S_H1_2) 

length(which(S_H1_2 < threshold))/length(S_H0_2) 


mean(S_H0_1)
mean(S_H0_2)
mean(S_H1_1)
mean(S_H1_2)



##Posterior predictive
n <- 1e3
params <- rep(0.001, 5)
V <- sapply(params, function(p) rbinom(n, 1, p))
V <- t(V)   # now V is K x n (row = dataset)

res1 <- smc_sampler(V)
res2 <- posterior_predictive_diff(res1, n)
quantile(res2, probs = c(0.025, 0.975)) # -0.487, 0.46 

  


