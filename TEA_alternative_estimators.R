remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
source("TEA_functions.R")



M <- 1e2
B <- 1e3 #number of estimates

a_theta <- 1
b_theta <- 1
p_eps <- 0.5



params <- seq(0.5, 0.8, length = 5)
params <- rep(0.7, 5)

K <- length(params)

n <- rep(10, K)
B <- 100

theta_map <- matrix(nrow = B, ncol = K)
theta_bma <- matrix(nrow = B, ncol = K)
theta_mean <- matrix(nrow = B, ncol = K)

for(b in 1:B){
  V <- vector("list", K)
  for (k in 1:K) {
    V[[k]] <- rbinom(n[k], 1, params[k])
  }
  
  smc_out <- smc_sampler(V, M)
  theta_all <- smc_out$theta_all
  epsilon_all <- smc_out$epsilon_all
  
  estimates <- param_estimates(theta_all, M)
  theta_map[b,] <- estimates[[1]]
  theta_bma[b,] <- estimates[[2]]
  theta_mean[b,] <- sapply(theta_all, mean)
}


###
x <- 1:K

# ---- MAP summaries ----
map_mean <- colMeans(theta_map)
map_q10  <- apply(theta_map, 2, quantile, probs = 0.05)
map_q90  <- apply(theta_map, 2, quantile, probs = 0.95)

# ---- BMA summaries ----
bma_mean <- colMeans(theta_bma)
bma_q10  <- apply(theta_bma, 2, quantile, probs = 0.05)
bma_q90  <- apply(theta_bma, 2, quantile, probs = 0.95)

# ---- Posterior mean summaries ----
mean_mean <- colMeans(theta_mean)
mean_q10  <- apply(theta_mean, 2, quantile, probs = 0.05)
mean_q90  <- apply(theta_mean, 2, quantile, probs = 0.95)


###
par(mfrow = c(3,1), mar = c(4,4,2,1))

### MAP
plot(x, map_mean,
     ylim = c(0,1),
     pch = 19,
     ylab = "theta",
     xlab = "version",
     main = "MAP estimate")

arrows(x, map_q10, x, map_q90,
       angle = 90, code = 3, length = 0.05)

points(x, params, col = rgb(0,0,1,0.4), pch = 16, cex = 1.5)


### BMA
plot(x, bma_mean,
     ylim = c(0,1),
     pch = 19,
     ylab = "theta",
     xlab = "version",
     main = "BMA estimate")

arrows(x, bma_q10, x, bma_q90,
       angle = 90, code = 3, length = 0.05)

points(x, params, col = rgb(0,0,1,0.4), pch = 16, cex = 1.5)


### Posterior Mean
plot(x, mean_mean,
     ylim = c(0,1),
     pch = 19,
     ylab = "theta",
     xlab = "version",
     main = "Posterior Mean")

arrows(x, mean_q10, x, mean_q90,
       angle = 90, code = 3, length = 0.05)

points(x, params, col = rgb(0,0,1,0.4), pch = 16, cex = 1.5)

###

