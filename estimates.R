remove(list = ls())
setwd("~/Documents/travail/dmd_when_next_clinical_trial?/bayesian_borrowing_information")
#setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("functions_EB.R")




M <- 100
B <- 500

n_versions <- 10

a_theta <- b_theta <- 0.5
mean_theta <- sd_theta <- sd_obs <- 0.5


binomial_params <- list(
  a_theta = a_theta,
  b_theta = b_theta
)

gaussian_params <- list(
  mean_theta = mean_theta,
  sd_theta = sd_theta,
  sd_obs = sd_obs
)


#plot(density(rbeta(1e5, a_theta, b_theta)))

thres <- 0.1
thres2 <- 0.8
coverage_interval <- 0.95

scenario <- "stable"
scenario <- "gradual"
scenario <- "abrupt"


if(scenario == "stable"){
  params <- rep(0.7, n_versions)
}else if(scenario == "gradual"){
  params <- seq(0.6, 0.9, length = n_versions)
}else if(scenario == "abrupt"){
  params <- c(rep(0.6, n_versions - 1), 0.9)
}

lambda <- 100
#lambda <- 100
lambdas <- rep(lambda - 1, n_versions)
CSD <- 0.2 #TEA with EB
epsilon_const <- 0


# 
# res <- posterior_sim(
#   params = params,
#   M = M,
#   B = B,
#   lambdas = lambdas,
#   epsilon_const = epsilon_const,
#   CSD = CSD,
#   thres = thres,
#   model = "binomial",
#   model_params = binomial_params
# )




res <- posterior_sim(
  params = params,
  M = M,
  B = B,
  lambdas = lambdas,
  epsilon_const = epsilon_const,
  CSD = CSD,
  thres = thres,
  model = "gaussian",
  model_params = gaussian_params, 
  coverage_interval = coverage_interval
)


plot_estimates(res$thetas, res$epsilons, params, lambda)


x <- 1:length(params)

plot(1:length(params), params,
     ylim = c(0, 1),
     pch = 16,
     col = rgb(0, 0, 1, 0.4),
     cex = 2,
     ylab = "θ",
     xlab = "Version",
     xaxt = "n")

axis(1, at = x, labels = x)

