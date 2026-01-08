remove(list = ls())

#test
#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")


source("TEA_functions.R")


lambda <- 1
M <- 1000 #number of particles
a <- 1#beta prior
b <- 1 #beta prior
J <- 5 #number of versions
n <- rep(2e1, 5) #number of observations per version
nsim <- 1e3 #number of simulations for Monte Carlo S estimate


params_H0 <- rep(0.5, 5)
S_H0_all <- simulate_S(params_H0, n, M, nsim, lambda)
S_H0 <- S_H0_all[[1]]
epsilons_H0 <- S_H0_all[[2]]
threshold <- quantile(S_H0, 0.10)  # 5% quantile for 5% false positive rate


# Type 1 error
S_H0_all <- simulate_S(params_H0, n, M, nsim, lambda)
S_H0 <- S_H0_all[[1]]
type1_error <- length(which(S_H0 < threshold))/length(S_H0) 


# Simulate under H1 (all params different)
params_H1_1 <- c(0.5, 0.6, 0.7, 0.6, 0.5)  # all params different
S_H1_1_all <- simulate_S(params_H1_1, n, M, nsim, lambda)
S_H1_1 <- S_H1_1_all[[1]]
epsilons_H1_1 <- S_H1_1_all[[2]]
# epsilons_H1_1 #0.75, 0.68, 0.65, 0.64
# Check how often S_H1 < threshold to estimate power
power1 <- length(which(S_H1_1 < threshold))/length(S_H1_1) #0.491 #0.461 0.461



# Simulate under H1 (last param different)
params_H1_2 <- c(0.5, 0.5, 0.9, 0.5, 0.5)  # last param different
S_H1_2_all <- simulate_S(params_H1_2, n, M, nsim, lambda)
S_H1_2 <- S_H1_2_all[[1]]
epsilons_H1_2 <- S_H1_2_all[[2]]
# epsilons_H1_2 #0.77, 0.78, 0.79, 0.23
# Check how often S_H1 < threshold to estimate power
power2 <- length(which(S_H1_2 < threshold))/length(S_H1_2) #0.76 0.723 0.71


# Simulate under H1 (last param different)
params_H1_3 <- c(0.5, 0.5, 0.7, 0.5, 0.5)  
S_H1_3_all <- simulate_S(params_H1_3, n, M, nsim, lambda)
S_H1_3 <- S_H1_3_all[[1]]
epsilons_H1_3 <- S_H1_3_all[[2]]
# epsilons_H1_2 #0.77, 0.78, 0.79, 0.23
# Check how often S_H1 < threshold to estimate power
power3 <- length(which(S_H1_3 < threshold))/length(S_H1_3) #0.76 0.723 0.71

print(paste0("n = ", paste(n, collapse = ", ")))
print(paste0("M =", M))
print(paste0("lambda =", lambda))
print(paste0("prior = Beta(",a, ",", b, ")"))

print(paste0("params H0 = ", paste(params_H0, collapse = ", ")))
print(paste0("params H1_1 = ", paste(params_H1_1, collapse = ", ")))
print(paste0("params H1_2 = ", paste(params_H1_2, collapse = ", ")))
print(paste0("params H1_3 = ", paste(params_H1_3, collapse = ", ")))

print(paste0("threshold =",threshold))
print(paste0("type1 error =",type1_error))
print(paste0("power1 =",power1))
print(paste0("power2 =",power2))
print(paste0("power3 =",power3))


print(paste0("mean_S =", round(mean(S_H0), 2),",", round(mean(S_H1_1), 2), ",", round(mean(S_H1_2), 2))  ) 
print(paste0("sd_S =", round(sqrt(var(S_H0)), 2),",", round(sqrt(var(S_H1_1)), 2), ",", round(sqrt(var(S_H1_2)), 2)  ) )

# ##Posterior predictive
# n <- 1e3
# params <- rep(0.001, 5)
# V <- sapply(params, function(p) rbinom(n, 1, p))
# V <- t(V)   # now V is K x n (row = dataset)
# 
# res1 <- smc_sampler(V, M, lambda)
# res2 <- posterior_predictive_diff(res1, n)
# quantile(res2, probs = c(0.025, 0.975)) # -0.487, 0.46 

  


