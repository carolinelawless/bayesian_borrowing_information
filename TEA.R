remove(list = ls())


#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")

source("TEA_functions.R")


start_time <- Sys.time()

M <- 1e2
K <- 5


a_theta <- 1
b_theta <- 1

p_eps <- 0.5

B <- 1e3 #number of delta estimates
sampsize <- 1e2 #sample size for the posterior predictives
#params_H0 <- rep(0.1, K)

# powers_list <- list()
# for(j in c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2)){
#   print(j)
#   n <- rep(j, K)
# 
#   diffs <- posterior_predictive_sim_diff(params_H0, M, B, sampsize)
#   abs_diffs <- abs(diffs)
#   threshold <- quantile(abs_diffs, 0.95)
#   threshold
# 
#   powers <- vector()
#   for(delta in (0:89)/100){
#     params_H1 <- c(rep(0.1, (K-1)), 0.1 + delta)
#     diffs <- posterior_predictive_sim_diff(params_H1, M, B, sampsize)
#     abs_diffs <- abs(diffs)
#     power <- length(which(abs_diffs >= threshold))/length(abs_diffs)
#     powers <- c(powers, power)
#   }
#   powers
#   powers_list[[length(powers_list)+1]] <- powers
# }

x <- (1:50)*4
x

params_list <- list(
  c(rep(0.5, 4), 0.1),
  c(rep(0.5, 4), 0.2),
  c(rep(0.5, 4), 0.3),
  c(rep(0.5, 4), 0.4),
  c(rep(0.5, 4), 0.6),
  c(rep(0.5, 4), 0.7),
  c(rep(0.5, 4), 0.8),
  c(rep(0.5, 4), 0.9)
)

params_list <- list(
  rep(0.1, 5),
  rep(0.2, 5), 
  rep(0.3, 5),
  rep(0.4, 5),
  rep(0.5, 5),
  rep(0.6, 5),
  rep(0.7, 5),
  rep(0.8, 5),
  rep(0.9, 5)
)


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

#print(paste0("delta =", delta))
# cat("posterior_pred_delta_10 <- c(", paste(diffs_list[[1]], collapse = ", "), ")\n")
# cat("posterior_pred_delta_20 <- c(", paste(diffs_list[[2]], collapse = ", "), ")\n")
# cat("posterior_pred_delta_50 <- c(", paste(diffs_list[[3]], collapse = ", "), ")\n")
# cat("posterior_pred_delta_100 <- c(", paste(diffs_list[[4]], collapse = ", "), ")\n")
# cat("posterior_pred_delta_200 <- c(", paste(diffs_list[[5]], collapse = ", "), ")\n")
# cat("posterior_pred_delta_500 <- c(", paste(diffs_list[[6]], collapse = ", "), ")\n")




cat("y_0.1 <- c(", paste(powers_list[[1]], collapse = ", "), ")\n")
cat("y_0.2 <- c(", paste(powers_list[[2]], collapse = ", "), ")\n")
cat("y_0.3 <- c(", paste(powers_list[[3]], collapse = ", "), ")\n")
cat("y_0.4 <- c(", paste(powers_list[[4]], collapse = ", "), ")\n")
cat("y_0.5 <- c(", paste(powers_list[[5]], collapse = ", "), ")\n")
cat("y_0.6 <- c(", paste(powers_list[[6]], collapse = ", "), ")\n")
cat("y_0.7 <- c(", paste(powers_list[[7]], collapse = ", "), ")\n")
cat("y_0.8 <- c(", paste(powers_list[[8]], collapse = ", "), ")\n")
cat("y_0.9 <- c(", paste(powers_list[[9]], collapse = ", "), ")\n")


cat("y_sc1 <- c(", paste(powers_list[[1]], collapse = ", "), ")\n")
cat("y_sc2 <- c(", paste(powers_list[[2]], collapse = ", "), ")\n")
cat("y_sc3 <- c(", paste(powers_list[[3]], collapse = ", "), ")\n")
cat("y_sc4 <- c(", paste(powers_list[[4]], collapse = ", "), ")\n")
cat("y_sc5 <- c(", paste(powers_list[[5]], collapse = ", "), ")\n")
cat("y_sc6 <- c(", paste(powers_list[[6]], collapse = ", "), ")\n")
cat("y_sc7 <- c(", paste(powers_list[[7]], collapse = ", "), ")\n")
cat("y_sc8 <- c(", paste(powers_list[[8]], collapse = ", "), ")\n")

y_sc1 <- c( 0.65, 0.738, 0.803, 0.849, 0.867, 0.907, 0.936, 0.955, 0.953, 0.964, 0.973, 0.971, 0.986, 0.98, 0.984, 0.984, 0.984, 0.986, 0.991, 0.987, 0.993, 0.992, 0.994, 0.994, 0.992, 0.993, 0.993, 0.994, 0.996, 0.994, 0.999, 0.995, 0.997, 0.998, 1, 0.997, 0.996, 1, 0.997, 0.994, 0.997, 0.997, 0.998, 1, 0.994, 0.998, 0.997, 0.999, 0.997, 0.998 )
y_sc2 <- c( 0.559, 0.599, 0.592, 0.624, 0.683, 0.723, 0.712, 0.741, 0.768, 0.788, 0.784, 0.827, 0.842, 0.844, 0.868, 0.866, 0.885, 0.869, 0.874, 0.874, 0.882, 0.902, 0.874, 0.874, 0.905, 0.901, 0.905, 0.898, 0.916, 0.911, 0.918, 0.916, 0.915, 0.921, 0.932, 0.92, 0.906, 0.899, 0.922, 0.92, 0.92, 0.936, 0.933, 0.934, 0.929, 0.93, 0.93, 0.933, 0.947, 0.933 )
y_sc3 <- c( 0.47, 0.457, 0.443, 0.442, 0.449, 0.407, 0.442, 0.456, 0.457, 0.455, 0.477, 0.468, 0.512, 0.478, 0.506, 0.468, 0.475, 0.542, 0.513, 0.502, 0.504, 0.529, 0.533, 0.561, 0.533, 0.546, 0.554, 0.563, 0.569, 0.577, 0.585, 0.586, 0.559, 0.585, 0.576, 0.591, 0.575, 0.58, 0.616, 0.596, 0.582, 0.596, 0.621, 0.588, 0.601, 0.623, 0.611, 0.587, 0.618, 0.605 )
y_sc4 <- c( 0.449, 0.406, 0.353, 0.317, 0.306, 0.236, 0.234, 0.224, 0.22, 0.205, 0.201, 0.184, 0.191, 0.164, 0.183, 0.194, 0.167, 0.168, 0.167, 0.146, 0.155, 0.152, 0.163, 0.156, 0.149, 0.161, 0.159, 0.132, 0.129, 0.139, 0.17, 0.164, 0.137, 0.166, 0.143, 0.142, 0.173, 0.145, 0.16, 0.158, 0.153, 0.142, 0.147, 0.159, 0.151, 0.154, 0.129, 0.143, 0.147, 0.152 )
y_sc5 <- c( 0.486, 0.395, 0.34, 0.287, 0.277, 0.269, 0.247, 0.226, 0.226, 0.209, 0.199, 0.183, 0.191, 0.173, 0.164, 0.154, 0.175, 0.152, 0.171, 0.148, 0.164, 0.18, 0.175, 0.156, 0.171, 0.144, 0.178, 0.149, 0.158, 0.142, 0.174, 0.15, 0.151, 0.148, 0.144, 0.132, 0.136, 0.146, 0.118, 0.162, 0.146, 0.151, 0.144, 0.146, 0.138, 0.145, 0.149, 0.16, 0.14, 0.146 )
y_sc6 <- c( 0.498, 0.484, 0.436, 0.457, 0.421, 0.427, 0.428, 0.435, 0.456, 0.421, 0.452, 0.467, 0.451, 0.477, 0.451, 0.48, 0.505, 0.48, 0.482, 0.527, 0.525, 0.535, 0.525, 0.53, 0.546, 0.559, 0.492, 0.543, 0.529, 0.551, 0.572, 0.544, 0.562, 0.578, 0.542, 0.556, 0.552, 0.574, 0.563, 0.559, 0.558, 0.562, 0.554, 0.579, 0.592, 0.564, 0.576, 0.569, 0.573, 0.583 )
y_sc7 <- c( 0.553, 0.595, 0.593, 0.651, 0.644, 0.705, 0.734, 0.738, 0.755, 0.747, 0.774, 0.81, 0.823, 0.835, 0.844, 0.841, 0.869, 0.868, 0.876, 0.882, 0.881, 0.879, 0.882, 0.881, 0.9, 0.894, 0.887, 0.894, 0.893, 0.904, 0.906, 0.918, 0.903, 0.907, 0.927, 0.906, 0.916, 0.92, 0.927, 0.925, 0.914, 0.912, 0.91, 0.905, 0.913, 0.926, 0.93, 0.909, 0.921, 0.932 )
y_sc8 <- c( 0.638, 0.693, 0.777, 0.837, 0.883, 0.898, 0.928, 0.941, 0.963, 0.965, 0.976, 0.978, 0.976, 0.986, 0.981, 0.983, 0.993, 0.987, 0.983, 0.994, 0.993, 0.994, 0.992, 0.995, 0.992, 0.998, 0.999, 0.994, 0.998, 0.997, 0.996, 0.999, 0.996, 0.992, 0.995, 0.998, 0.999, 0.996, 0.999, 0.998, 1, 0.996, 0.999, 0.998, 0.996, 0.997, 0.998, 0.998, 0.997, 1 )
