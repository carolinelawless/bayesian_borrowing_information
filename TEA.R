remove(list = ls())

print("test1")

#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")

source("TEA_functions.R")

print("test2")

n <- 1e3
nsim <- 1e3


S_H0_all <- simulate_S(rep(0.5, 5), n, n_sim = nsim)
S_H0 <- S_H0_all[[1]]
epsilons_H0 <- S_H0_all[[2]]
threshold <- quantile(S_H0, 0.05)  # 5% quantile for 5% false positive rate


print(threshold)
print("test3")



# Simulate under H1 (all params different)
params_H1_1 <- c(0.5, 0.55, 0.6, 0.65, 0.7)  # all params different
S_H1_1_all <- simulate_S(params_H1_1, n, n_sim = nsim)
S_H1_1 <- S_H1_1_all[[1]]
epsilons_H1_1 <- S_H1_1_all[[2]]
# epsilons_H1_1 #0.75, 0.68, 0.65, 0.64
# Check how often S_H1 < threshold to estimate power
power1 <- length(which(S_H1_1 < threshold))/length(S_H1_1) #0.491 #0.461 0.461

print("test4")

# Simulate under H1 (last param different)
params_H1_2 <- c(0.5, 0.5, 0.5, 0.5, 0.7)  # last param different
S_H1_2_all <- simulate_S(params_H1_2, n, n_sim = nsim)
S_H1_2 <- S_H1_2_all[[1]]
epsilons_H1_2 <- S_H1_2_all[[2]]
# epsilons_H1_2 #0.77, 0.78, 0.79, 0.23
# Check how often S_H1 < threshold to estimate power
power2 <- length(which(S_H1_2 < threshold))/length(S_H1_2) #0.76 0.723 0.71

print("test5")

print(power1)
print(power2)





# ##Posterior predictive
# n <- 1e3
# params <- rep(0.001, 5)
# V <- sapply(params, function(p) rbinom(n, 1, p))
# V <- t(V)   # now V is K x n (row = dataset)
# 
# res1 <- smc_sampler(V)
# res2 <- posterior_predictive_diff(res1, n)
# quantile(res2, probs = c(0.025, 0.975)) # -0.487, 0.46 

  


