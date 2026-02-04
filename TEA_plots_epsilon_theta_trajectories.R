##density plots of posteriors of theta and epsilon at each version

remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_functions.R")
source("TEA_scenarios.R")


M <- 1e2
B <- 1e3 #number of estimates

a_theta <- 1
b_theta <- 1
p_eps <- 0.5


#params <- scenarios[[8]] 

params <- rep(0.5, 2)
params <- rep(0.5, 5)
params <- scenarios[[25]]


K <- length(params)

n <- rep(10, K)
#n[2] <- 100



res <- posterior_sim(params, M, B)
thetas <- res[[1]]
epsilons <- res[[2]]


round(c(mean(thetas[[K]]), var(thetas[[K]])), 3) #00.500 0.012


par(mfrow = c(1,1))
x <- 1:K


#####

x <- 1:K
x_eps <- x - 0.1   # shift epsilon to the left

# Summaries
means_theta <- sapply(thetas, mean, na.rm = TRUE)
q10_theta   <- sapply(thetas, quantile, probs = 0.1, na.rm = TRUE)
q90_theta   <- sapply(thetas, quantile, probs = 0.9, na.rm = TRUE)

eps_means <- sapply(epsilons, mean, na.rm = TRUE)
eps_q10   <- sapply(epsilons, quantile, probs = 0.1, na.rm = TRUE)
eps_q90   <- sapply(epsilons, quantile, probs = 0.9, na.rm = TRUE)

# Base plot for theta
plot(x, means_theta,
     ylim = c(0, 1),
     pch = 16, col = "black",
     ylab = "Parameter value",
     xlab = "version",
     xaxt = "n")
axis(1, at = x, labels = x)

# Theta CI
arrows(x, q10_theta, x, q90_theta,
       angle = 90, code = 3, length = 0.05)

# True theta
points(x, params, pch = 16, col = rgb(0,0,1,0.4), cex = 1.8)

# Epsilon means (shifted)
points(x_eps[-1], eps_means[-1], pch = 4, col = "red", cex = 1.2)

# Epsilon CI (shifted)
arrows(x_eps[-1], eps_q10[-1], x_eps[-1], eps_q90[-1],
       angle = 90, code = 3, length = 0.05, col = "red")

# Legend
legend("topright",
       legend = c("theta posterior mean", "true theta", "epsilon posterior mean"),
       pch = c(16, 16, 4),
       col = c("black", rgb(0,0,1,0.4), "red"),
       bty = "n")


#####



# # Theta panel
# means_theta  <- sapply(thetas, mean, na.rm = TRUE)
# q10_theta    <- sapply(thetas, quantile, probs = 0.1, na.rm = TRUE)
# q90_theta    <- sapply(thetas, quantile, probs = 0.9, na.rm = TRUE)
# 
# 
# 
# plot(x, means_theta,
#      ylim = c(0, 1),
#      pch = 4,
#      cex = 1.4,
#      ylab = "theta",
#      xlab = "version",
#      xaxt = "n")   # suppress default x-axis
# 
# axis(1, at = x, labels = x)  # force integer ticks
# arrows(x0 = x, y0 = q10_theta, y1 = q90_theta, angle = 90, code = 3, length = 0.05, lty = 1, cex = 0.5)
# 
# points(x, params, pch = 16, col = rgb(0, 0, 1, 0.4), cex = 2)
# 
# 
# # Epsilon panel
# 
# 
# eps_means  <- sapply(epsilons, mean, na.rm = TRUE)
# eps_q10    <- sapply(epsilons, quantile, probs = 0.1, na.rm = TRUE)
# eps_q90    <- sapply(epsilons, quantile, probs = 0.9, na.rm = TRUE)
# 
# 
# plot(x, eps_means,
#      ylim = c(0, 1),
#      type = "n",  # empty plot
#      ylab = "epsilon",
#      xlab = "version",
#      xaxt = "n")
# 
# axis(1, at = x, labels = x)
# 
# points(x[-1], eps_means[-1], pch = 4, cex = 1.4)
# arrows(x[-1], eps_q10[-1], x[-1], eps_q90[-1],
#        angle = 90, code = 3, length = 0.05)
# 
# 
# 
# round(c(mean(thetas[[K]]), var(thetas[[K]])), 3) #0.845, 0.003
# 
