remove(list = ls())


#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")

source("TEA_functions.R")


start_time <- Sys.time()

M <- 1e2 #maybe increase to 1e4
K <- 5

B <- 1e3
a_theta <- 1
b_theta <- 1

p_eps <- 0.5

B <- 1e3 #number of delta estimates
sampsize <- 5e2 #sample size for the posterior predictives
params_H0 <- rep(0.1, K)

powers_list <- list()



for(j in c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2)){
  print(j)
  n <- rep(j, K)
  
  diffs <- posterior_predictive_sim(params_H0, M, B, sampsize)
  abs_diffs <- abs(diffs)
  threshold <- quantile(abs_diffs, 0.95)
  threshold
  
  powers <- vector()
  for(delta in (0:89)/100){
    params_H1 <- c(rep(0.1, (K-1)), 0.1 + delta)
    diffs <- posterior_predictive_sim(params_H1, M, B, sampsize)
    abs_diffs <- abs(diffs)
    power <- length(which(abs_diffs >= threshold))/length(abs_diffs)
    powers <- c(powers, power)
  }
  powers
  powers_list[[length(powers_list)+1]] <- powers
}

end_time <- Sys.time()
elapsed <- round(end_time - start_time, 2)

print(paste0("elapsed time =", elapsed))

print(paste0("M =", M))
print(paste0("B =", B))
print(paste0("sampsize =", sampsize))
print(paste0("K =", K))

print(cat("y_10 <- c(", paste(powers_list[[1]], collapse = ", "), ")\n"))
print(cat("y_20 <- c(", paste(powers_list[[2]], collapse = ", "), ")\n"))
print(cat("y_50 <- c(", paste(powers_list[[3]], collapse = ", "), ")\n"))
print(cat("y_100 <- c(", paste(powers_list[[4]], collapse = ", "), ")\n"))
print(cat("y_200 <- c(", paste(powers_list[[5]], collapse = ", "), ")\n"))
print(cat("y_500 <- c(", paste(powers_list[[6]], collapse = ", "), ")\n"))









