remove(list = ls())
setwd("~/Documents/travail/dmd_when_next_clinical_trial?/bayesian_borrowing_information")
source("functions7.R")


M <- 100
B <- 100
lambdas <- c(rep(4,9),49)


params <- seq(0.6, 0.9, length = 10)
#params <- c(rep(0.6, 9), 0.9)
#params <- rep(0.7, 10)
CSD <- 0.3
epsilon_const <- 0
a_theta <- b_theta <- 0.5
#plot(density(rbeta(1e4, a_theta, b_theta)))
res_EB <- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_const, CSD)
plot_trajectories(res_EB$thetas, res_EB$epsilons, params, epsilon_type = "TEA")

a_theta <- b_theta <- 0.5
res_EB0.05 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_const, CSD)
plot_trajectories(res_EB0.05$thetas, res_EB0.05$epsilons, params, epsilon_type = "TEA")
a_theta <- b_theta <- 10
res_EB100 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_const, CSD)
plot_trajectories(res_EB100$thetas, res_EB100$epsilons, params, epsilon_type = "TEA")


a_theta <- b_theta <- 10
plot(density(rbeta(1e5, a_theta, b_theta)))
res_TEA0.5 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0.5, CSD)
plot_trajectories(res_TEA0.5$thetas, res_TEA0.5$epsilons, params, epsilon_type = "TEA")

a_eps <- b_eps <- 0.01
res_TEA0.5 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0.5, CSD)
plot_trajectories(res_TEA0.5$thetas, res_TEA0.5$epsilons, params, epsilon_type = "TEA")


res_TEA0.5 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0.5, CSD)
plot_trajectories(res_TEA0.5$thetas, res_TEA0.5$epsilons, params, epsilon_type = "TEA")
res_TEA0 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0, CSD)
plot_trajectories(res_TEA0$thetas, res_TEA0$epsilons, params, epsilon_type = "TEA")
res_TEA1 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 1, CSD)
plot_trajectories(res_TEA1$thetas, res_TEA1$epsilons, params, epsilon_type = "TEA")

lambda <- 19
res_TEA0 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0, CSD)
res_TEA0.5 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0.5, CSD)
res_TEA1 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 1, CSD)
res_EB <- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "EB", a_eps, b_eps, epsilon_const = 0.5, CSD)
plot_trajectories(res_EB$thetas, res_EB$epsilons, params, epsilon_type = "TEA")
res_ATEA <- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "ATEA", a_eps, b_eps, epsilon_const, CSD)


plot_trajectories(res_TEA0$thetas, res_TEA0$epsilons, params, epsilon_type = "TEA")
plot_trajectories(res_TEA0.5$thetas, res_TEA0.5$epsilons, params, epsilon_type = "TEA")
plot_trajectories(res_TEA1$thetas, res_TEA1$epsilons, params, epsilon_type = "TEA")
plot_trajectories(res_EB$thetas, res_EB$epsilons, params, epsilon_type = "EB")
plot_trajectories(res_ATEA$thetas, res_ATEA$epsilons, params, epsilon_type = "ATEA")

