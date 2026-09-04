remove(list = ls())
#setwd("~/Documents/travail/dmd_when_next_clinical_trial?/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("functions_EB.R")

t1 <- Sys.time()


M <- 100
B <- 20000

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

thres1 <- 0.1
thres2 <- 0.8
coverage_interval <- 0.75

model <- "gaussianl"
model_params <- gaussian_params

#scenario <- "stable"
#scenario <- "gradual"
scenario <- "abrupt"


if(scenario == "stable"){
  params <- rep(0.7, n_versions)
}else if(scenario == "gradual"){
  params <- seq(0.6, 0.9, length = n_versions)
}else if(scenario == "abrupt"){
  params <- c(rep(0.6, n_versions - 1), 0.9)
}



J <- length(params)
true_diff <- params[J] - params[1]



alpha_TEA_0 <- vector()
alpha_TEA_0.5 <- vector()
alpha_TEA_1 <- vector()
alpha_TEA_EB_0.5 <- vector()
alpha_TEA_EB_1 <- vector()
alpha_TEA_EB_epsilon <- vector()

#lambda <- 10
#lambda_vector <- (1:5)*10
lambda_vector <- (1:25)*2
for(lambda in lambda_vector){
  print(lambda)
  lambdas <- rep(lambda - 1, n_versions)
  
  CSD <- 1 #TEA without EB
  epsilon_const <- -1
  res <- posterior_sim(
    params = params,
    M = M,
    B = B,
    lambdas = lambdas,
    epsilon_const = epsilon_const,
    CSD = CSD,
    thres1 = thres1,
    model = model,
    model_params = model_params,
    coverage_interval = coverage_interval
  )
  rho <- res$rho
  alpha <- sum(rho > thres2)/B
  alpha_TEA_EB_epsilon <- c(alpha_TEA_EB_epsilon, alpha)
  
  
  
  epsilon_const <- 0
  res <- posterior_sim(
    params = params,
    M = M,
    B = B,
    lambdas = lambdas,
    epsilon_const = epsilon_const,
    CSD = CSD,
    thres1 = thres1,
    model = model,
    model_params = model_params,
    coverage_interval = coverage_interval
  )
  rho <- res$rho
  alpha <- sum(rho > thres2)/B
  alpha_TEA_0 <- c(alpha_TEA_0, alpha)
  
  epsilon_const <- 0.5
  res <- posterior_sim(
    params = params,
    M = M,
    B = B,
    lambdas = lambdas,
    epsilon_const = epsilon_const,
    CSD = CSD,
    thres1 = thres1,
    model = model,
    model_params = model_params,
    coverage_interval = coverage_interval
  )
  rho <- res$rho
  alpha <- sum(rho > thres2)/B
  alpha_TEA_0.5 <- c(alpha_TEA_0.5, alpha)
  
  epsilon_const <- 1
  res <- posterior_sim(
    params = params,
    M = M,
    B = B,
    lambdas = lambdas,
    epsilon_const = epsilon_const,
    CSD = CSD,
    thres1 = thres1,
    model = model,
    model_params = model_params,
    coverage_interval = coverage_interval
  )
  rho <- res$rho
  alpha <- sum(rho > thres2)/B
  alpha_TEA_1 <- c(alpha_TEA_1, alpha)
  
  
  CSD <- 0.2 #TEA with EB
  epsilon_const <- 0.5
  res <- posterior_sim(
    params = params,
    M = M,
    B = B,
    lambdas = lambdas,
    epsilon_const = epsilon_const,
    CSD = CSD,
    thres1 = thres1,
    model = model,
    model_params = model_params,
    coverage_interval = coverage_interval
  )
  rho <- res$rho
  alpha <- sum(rho > thres2)/B
  alpha_TEA_EB_0.5 <- c(alpha_TEA_EB_0.5, alpha)
  
  epsilon_const <- 1
  res <- posterior_sim(
    params = params,
    M = M,
    B = B,
    lambdas = lambdas,
    epsilon_const = epsilon_const,
    CSD = CSD,
    thres1 = thres1,
    model = model,
    model_params = model_params,
    coverage_interval = coverage_interval
  )
  rho <- res$rho
  alpha <- sum(rho > thres2)/B
  alpha_TEA_EB_1 <- c(alpha_TEA_EB_1, alpha)
  
  
  
  
}

t2 <- Sys.time()


cat("elapsed time <- ", round(difftime(t2, t1, units = "hours"), 2), "\n")
cat("elapsed time <- ", round(difftime(t2, t1, units = "mins"), 2), "\n")

cat("thres1 <- ", thres1, "\n")
cat("thres2 <- ", thres2, "\n")
cat("M <- ", M, "\n")
cat("B <- ", B, "\n")
cat("params <- c(", paste(params, collapse = ", "), ")\n")

cat("n_versions <- ", n_versions, "\n")
cat("model <- '", model, "'\n")
cat("scenario <- '", scenario, "'\n")
cat("lambda_vector <- c(", paste(lambda_vector, collapse = ", "), ")\n")
cat("alpha_TEA_0<- c(",paste(alpha_TEA_0, collapse = ", "), ")\n")
cat("alpha_TEA_0.5<- c(",paste(alpha_TEA_0.5, collapse = ", "), ")\n")
cat("alpha_TEA_1<- c(",paste(alpha_TEA_1, collapse = ", "), ")\n")
cat("alpha_TEA_EB_0.5<- c(",paste(alpha_TEA_EB_0.5, collapse = ", "), ")\n")
cat("alpha_TEA_EB_1<- c(",paste(alpha_TEA_EB_1, collapse = ", "), ")\n")
cat("alpha_TEA_EB_epsilon<- c(",paste(alpha_TEA_EB_epsilon, collapse = ", "), ")\n")

