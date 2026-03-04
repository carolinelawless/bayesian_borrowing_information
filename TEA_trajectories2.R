##density plots of posteriors of theta and epsilon at each version

#####

#Binary setting

####

remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_functions.R")



M <- 1e2
B <- 1e3 #number of estimates

a_theta <- 1
b_theta <- 1
p_eps <- 0.5



params <- seq(0.5, 0.8, length = 5)


K <- length(params)

lambda <- 9




res <- posterior_sim(params, M, B)
thetas <- res[[1]]
epsilons <- res[[2]]





par(mfrow = c(1,1))


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

#Gaussian setting

####



remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_functions.R")



M <- 1e2
B <- 1e3 #number of estimates

a_theta <- 1
b_theta <- 1
p_eps <- 0.5



params <- seq(0.5, 0.8, length = 5)


K <- length(params)

lambda <- 9




res <- posterior_sim(params, M, B)
thetas <- res[[1]]
epsilons <- res[[2]]





par(mfrow = c(1,1))


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




