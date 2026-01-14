remove(list = ls())

setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
source("TEA_functions.R")

a <- 1
b <- 1
M <- 1000
lambda <- 1
J <- 5
n <- rep(1e2, J)
params0 <- rep(0.5, J)
params1 <- c(0.5, 0.6, 0.7, 0.8, 0.9)
params2 <- c(0.5, 0.5, 0.5, 0.5, 0.9)
npred <- 5e3
delta <- 0.2
V <- list()

for(j in 1:J){
  V[[j]] <- rbinom(n[j], 1, params0[j])
}

smc_out <- smc_sampler(V, M, lambda)

diffs <- posterior_predictive_diff(smc_out, npred)

gamma <- sum(abs(diffs) > delta)/length(diffs)
round(gamma, 2)

