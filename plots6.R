remove(list = ls())
source("functions7.R")

scenario <- "Gradual"
scenario <- "Abrupt"
#scenario <- "Drift"
scenario <- "Stable"

if(scenario == "Gradual"){
  params <- seq(0.6, 0.9, length = 10)
}else if(scenario == "Abrupt"){
  params <- c(rep(0.6, 9), 0.9)
}else if(scenario == "Drift"){
  params <- c(seq(0.6, 0.9, length = 5), seq(0.9, 0.6, length = 5))
}else if(scenario == "Stable"){
  params <- rep(0.7, 10)}

a_theta <- 0.5
b_theta <- 0.5

M <- 500
B <- 100


CSD <- 0.2
epsilon_const <- 0.5
#thres <- 0.8
#lambdas <- 1:50
lambda <- 19

a_theta <- b_theta <- 0.005
res_EB0.05 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_const, CSD)
plot_trajectories(res_EB0.05$thetas, res_EB0.05$epsilons, params, epsilon_type = "TEA")
a_theta <- b_theta <- 100
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

