##density plots of posteriors of theta and epsilon at each version

remove(list = ls())

source("TEA_functions.R")



M <- 1e2
B <- 1e3 #number of estimates

a_theta <- 1
b_theta <- 1

p_eps <- 0.5

B <- 1e3 #number of estimates






params <- c(rep(0.5, 4), 0.1)
params <- 1:9/10
#params <- c(0.1, 0.5, 0.1, 0.5, 0.1, 0.5)
params <- c(0.1, 0.9)
params <- c(0.1, 0.1, 0.6, 0.1, 0.1)
K <- length(params)

n <- rep(20, K)
#n[2] <- 100



res <- posterior_sim(params, M, B)
thetas <- res[[1]]
epsilons <- res[[2]]


round(c(mean(thetas[[K]]), var(thetas[[K]])), 3) #0.845, 0.003


par(mfrow = c(2,1))
x <- 1:K

# Theta panel
means_theta  <- sapply(thetas, mean, na.rm = TRUE)
q10_theta    <- sapply(thetas, quantile, probs = 0.1, na.rm = TRUE)
q90_theta    <- sapply(thetas, quantile, probs = 0.9, na.rm = TRUE)



plot(x, means_theta,
     ylim = c(0, 1),
     pch = 4,
     cex = 1.4,
     ylab = "theta",
     xlab = "version",
     xaxt = "n")   # suppress default x-axis

axis(1, at = x, labels = x)  # force integer ticks
arrows(x0 = x, y0 = q10_theta, y1 = q90_theta, angle = 90, code = 3, length = 0.05, lty = 1, cex = 0.5)

points(x, params, pch = 16, col = rgb(0, 0, 1, 0.4), cex = 2)


# Epsilon panel


eps_means  <- sapply(epsilons, mean, na.rm = TRUE)
eps_q10    <- sapply(epsilons, quantile, probs = 0.1, na.rm = TRUE)
eps_q90    <- sapply(epsilons, quantile, probs = 0.9, na.rm = TRUE)


plot(x, eps_means,
     ylim = c(0, 1),
     type = "n",  # empty plot
     ylab = "epsilon",
     xlab = "version",
     xaxt = "n")

axis(1, at = x, labels = x)

points(x[-1], eps_means[-1], pch = 4, cex = 1.4)
arrows(x[-1], eps_q10[-1], x[-1], eps_q90[-1],
       angle = 90, code = 3, length = 0.05)



round(c(mean(thetas[[K]]), var(thetas[[K]])), 3) #0.845, 0.003

