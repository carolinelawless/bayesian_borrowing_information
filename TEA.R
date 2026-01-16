remove(list = ls())


setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
source("TEA_functions.R")

M <- 1e5
J <- 20

a_theta <- 1
b_theta <- 1

a_eps <- 1
b_eps <- 1
p_eps <- 0.5


n <- rep(1e1, J)

B <- 1e2


params_H0 <- rep(0.1, J)
params_H1 <- c(rep(0.1, (J-1)), 0.4)



V_H0 <- list()
for(j in 1:J){
  V_H0[[j]] <- rbinom(n[j], 1, params_H0[j])
}
smc_out_H0 <- smc_sampler(V_H0, M)

V_H1 <- list()
for(j in 1:J){
  V_H1[[j]] <- rbinom(n[j], 1, params_H1[j])
}
smc_out_H1 <- smc_sampler(V_H1, M)

mean(smc_out_H1$epsilon_all[[J]])
mean(smc_out_H0$epsilon_all[[J]])


res_H0 <- posterior_predictive_diff(smc_out_H0, B)
res_H1 <- posterior_predictive_diff(smc_out_H1, B)

mean(res_H0)
mean(res_H1)




