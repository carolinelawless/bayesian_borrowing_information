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
params <- 1:10/10
params <- c(0.1, 0.5, 0.1, 0.5, 0.1, 0.5)


K <- length(params)

n <- rep(20, K)
#n[2] <- 100

res <- posterior_sim(params, M, B)
thetas <- res[[1]]
epsilons <- res[[2]]


par(mfrow = c(2,1))

# Theta panel
means  <- sapply(thetas, mean, na.rm = TRUE)
q10    <- sapply(thetas, quantile, probs = 0.1, na.rm = TRUE)
q90    <- sapply(thetas, quantile, probs = 0.9, na.rm = TRUE)
K <- length(thetas)
x <- 1:K

plot(x, means,
     ylim = c(0, 1),
     pch = 19,
     ylab = "theta",
     xlab = "version",
     xaxt = "n")   # suppress default x-axis

axis(1, at = x, labels = x)  # force integer ticks
arrows(x0 = x, y0 = q10, y1 = q90, angle = 90, code = 3, length = 0.05)
theta_true <- c(0.1, 0.1, 0.1, 0.1, 0.5)
points(x, params, pch = 16, col = rgb(0, 0, 1, 0.4), cex = 1.3)


# Epsilon panel


eps_means  <- sapply(epsilons, mean, na.rm = TRUE)
eps_q10    <- sapply(epsilons, quantile, probs = 0.1, na.rm = TRUE)
eps_q90    <- sapply(epsilons, quantile, probs = 0.9, na.rm = TRUE)


plot(x, eps_means,
     ylim = c(0, 1),
     pch = 19,
     ylab = "epsilon",
     xlab = "version",
     xaxt = "n")   # suppress default x-axis

axis(1, at = x, labels = x)  # force integer ticks
arrows(x0 = x, y0 = eps_q10, y1 = eps_q90, angle = 90, code = 3, length = 0.05)





