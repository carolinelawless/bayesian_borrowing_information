remove(list = ls())


setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")


source("TEA_functions.R")
source("TEA_scenarios.R")


start_time <- Sys.time()

M <- 1e2
K <- 5


a_theta <- 1
b_theta <- 1

p_eps <- 0.5

B <- 1e3 #number of delta estimates
sampsize <- 1e2 #sample size for the posterior predictives

x <- (1:50)*4
x





params_list <- scenarios[1:8]
params_list <- scenarios[9:16]
params_list <- scenarios[17:20]
params_list <- scenarios[21:25]


n <- rep(50, K)
params_H0 <- rep(0.5, K)
diffs <- posterior_predictive_sim_diff(params_H0, M, B, sampsize)
abs_diffs <- abs(diffs)
threshold <- quantile(abs_diffs, 0.9)
threshold <- 0.18
threshold #0.18

powers_list <- list()

for(i in 1:length(params_list)){
  params_H1 <- params_list[[i]]
  K <- length(params_H1)
  
  powers <- vector()
  for(j in x){
    
    print(c(i, j))
    

    n <- rep(j, K)
    
    diffs <- posterior_predictive_sim_diff(params_H1, M, B, sampsize)
    abs_diffs <- abs(diffs)
    power <- length(which(abs_diffs >= threshold))/length(abs_diffs)
    powers <- c(powers, power)
  }
  powers_list[[length(powers_list)+1]] <- powers
}



end_time <- Sys.time()
elapsed <- round(end_time - start_time, 2)

print(paste0("elapsed time =", elapsed))

print(paste0("M =", M))
print(paste0("B =", B))
print(paste0("sampsize =", sampsize))
print(paste0("K =", K))
print(paste0("H0 =", params_H0))
print(paste0("H1 =", params_H1))






cat("y_sc1 <- c(", paste(powers_list[[1]], collapse = ", "), ")\n")
cat("y_sc2 <- c(", paste(powers_list[[2]], collapse = ", "), ")\n")
cat("y_sc3 <- c(", paste(powers_list[[3]], collapse = ", "), ")\n")
cat("y_sc4 <- c(", paste(powers_list[[4]], collapse = ", "), ")\n")
cat("y_sc5 <- c(", paste(powers_list[[5]], collapse = ", "), ")\n")
cat("y_sc6 <- c(", paste(powers_list[[6]], collapse = ", "), ")\n")
cat("y_sc7 <- c(", paste(powers_list[[7]], collapse = ", "), ")\n")
cat("y_sc8 <- c(", paste(powers_list[[8]], collapse = ", "), ")\n")

cat("y_sc9 <- c(", paste(powers_list[[1]], collapse = ", "), ")\n")
cat("y_sc10 <- c(", paste(powers_list[[2]], collapse = ", "), ")\n")
cat("y_sc11 <- c(", paste(powers_list[[3]], collapse = ", "), ")\n")
cat("y_sc12 <- c(", paste(powers_list[[4]], collapse = ", "), ")\n")
cat("y_sc13 <- c(", paste(powers_list[[5]], collapse = ", "), ")\n")
cat("y_sc14 <- c(", paste(powers_list[[6]], collapse = ", "), ")\n")
cat("y_sc15 <- c(", paste(powers_list[[7]], collapse = ", "), ")\n")
cat("y_sc16 <- c(", paste(powers_list[[8]], collapse = ", "), ")\n")

cat("y_sc17 <- c(", paste(powers_list[[1]], collapse = ", "), ")\n")
cat("y_sc18 <- c(", paste(powers_list[[2]], collapse = ", "), ")\n")
cat("y_sc19 <- c(", paste(powers_list[[3]], collapse = ", "), ")\n")
cat("y_sc20 <- c(", paste(powers_list[[4]], collapse = ", "), ")\n")

cat("y_sc21 <- c(", paste(powers_list[[1]], collapse = ", "), ")\n")
cat("y_sc22 <- c(", paste(powers_list[[2]], collapse = ", "), ")\n")
cat("y_sc23 <- c(", paste(powers_list[[3]], collapse = ", "), ")\n")
cat("y_sc24 <- c(", paste(powers_list[[4]], collapse = ", "), ")\n")
cat("y_sc25 <- c(", paste(powers_list[[5]], collapse = ", "), ")\n")
cat("y_sc26 <- c(", paste(powers_list[[6]], collapse = ", "), ")\n")
cat("y_sc27 <- c(", paste(powers_list[[7]], collapse = ", "), ")\n")
cat("y_sc28 <- c(", paste(powers_list[[8]], collapse = ", "), ")\n")
cat("y_sc29 <- c(", paste(powers_list[[9]], collapse = ", "), ")\n")






