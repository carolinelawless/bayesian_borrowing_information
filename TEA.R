remove(list = ls())


setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
source("TEA_functions.R")





M <- 1e3
J <- 5

# prior_mean <- 0.5
# prior_ess <- 50
# a_theta <- prior_mean*prior_ess
# b_theta <- (1-prior_mean)*prior_ess

a_theta <- 1
b_theta <- 1

a_eps <- 1
b_eps <- 1
p_eps <- 0.5


n <- rep(5e1, J)

B <- 1e2


params <- c(rep(0.5, 4), 0.9)
J <- length(params)
V <- list()
for(j in 1:J){
  V[[j]] <- rbinom(50, 1, params[j])
}


smc_out <- smc_sampler(V, M)
thetas <- smc_out$theta_all
epsilons <- smc_out$epsilon_all



params_H0 <- rep(0.1, 5)
params_H1 <- c(0.1, 0.1, 0.1, 0.1, 0.2)



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


res_H0 <- posterior_predictive_diff(smc_out_H0, B)
res_H1 <- posterior_predictive_diff(smc_out_H1, B)

mean(res_H0)
mean(res_H1)




