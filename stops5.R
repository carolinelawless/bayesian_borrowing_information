remove(list = ls())
setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("functions5.R")

start_time <- Sys.time()
scenario <- "Gradual"


if(scenario == "Gradual"){
  params <- seq(0.6, 0.9, length = 10)
}else if(scenario == "Abrupt"){
  params <- c(rep(0.6, 9), 0.9)
}else if(scenario == "Drift"){
  params <- c(seq(0.6, 0.9, length = 5), seq(0.9, 0.6, length = 5))
}else if(scenario == "Stable"){
  params <- rep(0.7, 10)}



M <- 100
B <- 5000
a_theta <- 1
b_theta <- 1
a_eps <- 1
b_eps <- 1
CSD <- 0.1
thres <- 0.8
lambdas <- 1:50
#lambda <- 19



stops_TEA0 <- vector()
stops_TEA0.5 <- vector()
stops_TEA1 <- vector()
stops_EB <- vector()
stops_ATEA <- vector()
K <- length(params)

for(lambda in lambdas){
  print(lambda)
  res_TEA0 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0, CSD)
  res_TEA0.5 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 0.5, CSD)
  res_TEA1 <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, epsilon_type = "TEA", a_eps, b_eps, epsilon_const = 1, CSD)
  res_EB <- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "EB", a_eps, b_eps, epsilon_const, CSD)
  res_ATEA <- posterior_sim_binomial(params, M, B, lambdas, a_theta, b_theta, epsilon_type = "ATEA", a_eps, b_eps, epsilon_const, CSD)
  

  res <- res_TEA0
  k<- 0
  stop<- 0
  while(stop == 0 & k < K){
    k <- k + 1
    stat <- sum(res$thetas[[k]] - res$thetas[[1]] > CSD)/B
    if(stat > thres){
      stop <- 1
    }
  }
  stops_TEA0 <- c(stops_TEA0, k)

  res <- res_TEA0.5
  k<- 0
  stop<- 0
  while(stop == 0 & k < K){
    k <- k + 1
    stat <- sum(res$thetas[[k]] - res$thetas[[1]] > CSD)/B
    if(stat > thres){
      stop <- 1
    }
  }
  stops_TEA0.5 <- c(stops_TEA0.5, k)
  
  res <- res_TEA1
  k<- 0
  stop<- 0
  while(stop == 0 & k < K){
    k <- k + 1
    stat <- sum(res$thetas[[k]] - res$thetas[[1]] > CSD)/B
    if(stat > thres){
      stop <- 1
    }
  }
  stops_TEA1 <- c(stops_TEA1, k)
  
  res <- res_ATEA
  k<- 0
  stop<- 0
  while(stop == 0 & k < K){
    k <- k + 1
    stat <- sum(res$thetas[[k]] - res$thetas[[1]] > CSD)/B
    if(stat > thres){
      stop <- 1
    }
  }
  stops_ATEA <- c(stops_ATEA, k)
  
  res <- res_EB
  k<- 0
  stop<- 0
  while(stop == 0 & k < K){
    k <- k + 1
    stat <- sum(res$thetas[[k]] - res$thetas[[1]] > CSD)/B
    if(stat > thres){
      stop <- 1
    }
  }
  stops_EB <- c(stops_EB, k)
  

}
end_time <- Sys.time()

paste0("time elapsed =", end_time - start_time)
paste0("M =", M)
paste0("B = ", B)
paste0("a_theta =", a_theta)
paste0("b_theta =", b_theta)
paste0("a_eps =", a_eps)
paste0("b_eps =", b_eps)
paste0("CSD =", CSD)
paste0("thres =", thres)
cat("lambdas<- c(",paste(lambdas, collapse = ", "), ")", ")\n")


cat("stops_TEA0<- c(",paste(stops_TEA0, collapse = ", "), ")", ")\n")
cat("stops_TEA0.5<- c(",paste(stops_TEA0.5, collapse = ", "), ")", ")\n")
cat("stops_TEA1<- c(",paste(stops_TEA1, collapse = ", "), ")", ")\n")
cat("stops_ATEA<- c(",paste(stops_ATEA, collapse = ", "), ")", ")\n")
cat("stops_EB<- c(",paste(stops_EB, collapse = ", "), ")", ")\n")



# plot(lambdas, stops_TEA0, type = "l", xlab = expression(lambda), ylab = "Version")
# lines(lambdas, stops_TEA0.5, col = "purple")
# lines(lambdas, stops_TEA1, lty = 2)
# lines(lambdas, stops_EB, col = "blue")
# lines(lambdas, stops_ATEA, col = "red")
# legend(
#   "topright",
#   legend = c(
#     expression("TEA ("*epsilon*" = 0)"),
#     expression("TEA ("*epsilon*" = 0.5)"),
#     expression("TEA ("*epsilon*" = 1)"),
#     "Empirical Bayes",
#     "A-TEA"
#   ),
#   col = c("black", "purple", "black", "blue", "red"),
#   lty = c(1, 1, 2, 1, 1),
#   lwd = 2,
#   bty = "n"
# )

