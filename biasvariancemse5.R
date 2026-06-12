remove(list = ls())
source("functions5.R")

scenario <- "Gradual"
scenario <- "Abrupt"
scenario <- "Drift"
scenario <- "Stable"

if(scenario == "Gradual"){
  params <- seq(0.6, 0.9, length = 10)
}else if(scenario == "Abrupt"){
  params <- c(rep(0.6, 9), 0.9)
}else if(scenario == "Drift"){
  params <- c(seq(0.6, 0.9, length = 5), seq(0.9, 0.6, length = 5))
}else if(scenario == "Stable"){
  params <- rep(0.7, 10)}



M <- 100
B <- 2000
a_theta <- 1
b_theta <- 1
a_eps <- 1
b_eps <- 1
CSD <- 0.1
lambdas <- 1:50





res_TEA0 <- TEA_eval_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0, CSD)
res_TEA0.5 <- TEA_eval_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0.5, CSD)
res_TEA1 <- TEA_eval_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 1, CSD)
res_EB <- TEA_eval_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "EB", a_eps, b_eps, epsilon_const, CSD)
res_ATEA <- TEA_eval_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "ATEA", a_eps, b_eps, epsilon_const, CSD)


plot(
  lambdas,
  res_TEA0$bias,
  type = "l",
  lwd = 2,
  col = "black",
  xlab = expression(lambda),
  ylab = "Bias",
  ylim = c(
    0,
    max(
      res_TEA0$bias,
      res_TEA0.5$bias,
      res_TEA1$bias,
      res_EB$bias,
      res_ATEA$bias
    )
  ),
  main = scenario
)

lines(lambdas, res_TEA0.5$bias,
      col = "purple",
      lwd = 2)

lines(lambdas, res_TEA1$bias,
      lty = 2,
      lwd = 2)

lines(lambdas, res_EB$bias,
      col = "blue",
      lwd = 2)

lines(lambdas, res_ATEA$bias,
      col = "red",
      lwd = 2)


legend(
  "topright",
  legend = c(
    expression("TEA ("*epsilon*" = 0)"),
    expression("TEA ("*epsilon*" = 0.5)"),
    expression("TEA ("*epsilon*" = 1)"),
    "Empirical Bayes",
    "A-TEA"
  ),
  col = c("black", "purple", "black", "blue", "red"),
  lty = c(1, 1, 2, 1, 1),
  lwd = 2,
  bty = "n"
)


plot(
  lambdas,
  res_TEA0$variance,
  type = "l",
  lwd = 2,
  col = "black",
  xlab = expression(lambda),
  ylab = "Variance",
  ylim = c(
    0,
    max(
      res_TEA0$variance,
      res_TEA0.5$variance,
      res_TEA1$variance,
      res_EB$variance,
      res_ATEA$variance
    )
  ),
  main = scenario
)

lines(lambdas, res_TEA0.5$variance,
      col = "purple",
      lwd = 2)

lines(lambdas, res_TEA1$variance,
      lty = 2,
      lwd = 2)

lines(lambdas, res_EB$variance,
      col = "blue",
      lwd = 2)

lines(lambdas, res_ATEA$variance,
      col = "red",
      lwd = 2)


legend(
  "topright",
  legend = c(
    expression("TEA ("*epsilon*" = 0)"),
    expression("TEA ("*epsilon*" = 0.5)"),
    expression("TEA ("*epsilon*" = 1)"),
    "Empirical Bayes",
    "A-TEA"
  ),
  col = c("black", "purple", "black", "blue", "red"),
  lty = c(1, 1, 2, 1, 1),
  lwd = 2,
  bty = "n"
)



plot(
  lambdas,
  res_TEA0$MSE,
  type = "l",
  lwd = 2,
  col = "black",
  xlab = expression(lambda),
  ylab = "MSE",
  ylim = c(
    0,
    max(
      res_TEA0$MSE,
      res_TEA0.5$MSE,
      res_TEA1$MSE,
      res_EB$MSE,
      res_ATEA$MSE
    )
  ),
  main = scenario
)

lines(lambdas, res_TEA0.5$MSE,
      col = "purple",
      lwd = 2)

lines(lambdas, res_TEA1$MSE,
      lty = 2,
      lwd = 2)

lines(lambdas, res_EB$MSE,
      col = "blue",
      lwd = 2)

lines(lambdas, res_ATEA$MSE,
      col = "red",
      lwd = 2)


legend(
  "topright",
  legend = c(
    expression("TEA ("*epsilon*" = 0)"),
    expression("TEA ("*epsilon*" = 0.5)"),
    expression("TEA ("*epsilon*" = 1)"),
    "Empirical Bayes",
    "A-TEA"
  ),
  col = c("black", "purple", "black", "blue", "red"),
  lty = c(1, 1, 2, 1, 1),
  lwd = 2,
  bty = "n"
)

lambda <- 19
res_TEA0 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0, CSD)
res_TEA0.5 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0.5, CSD)
res_TEA1 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 1, CSD)
res_EB <- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "EB", a_eps, b_eps, epsilon_const, CSD)
res_ATEA <- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "ATEA", a_eps, b_eps, epsilon_const, CSD)


plot_trajectories(res_TEA0$thetas, res_TEA0$epsilons, params, epsilon_type = "TEA")
plot_trajectories(res_TEA0.5$thetas, res_TEA0.5$epsilons, params, epsilon_type = "TEA")
plot_trajectories(res_TEA1$thetas, res_TEA1$epsilons, params, epsilon_type = "TEA")
plot_trajectories(res_EB$thetas, res_EB$epsilons, params, epsilon_type = "EB")
plot_trajectories(res_ATEA$thetas, res_ATEA$epsilons, params, epsilon_type = "ATEA")
