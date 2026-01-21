##density plots of posteriors of theta and epsilon at each version

remove(list = ls())


#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")

source("TEA_functions.R")




M <- 1e4 

a_theta <- 1
b_theta <- 1

p_eps <- 0.5

B <- 1e3 #number of delta estimates
sampsize <- 1e2 #sample size for the posterior predictives


params <- 1:5/10
K <- length(params)

n <- rep(20, K)


res <- posterior_sim(params, M, B)
thetas <- res[[1]]
epsilons <- res[[2]]


cols <- c("#08519C", "#3182BD", "#756BB1", "#9E9AC8", "#54278F")
cols <- c("#4E79A7", "#59A14F", "#9C755F", "#B07AA1", "#76B7B2")
cols <- c("#4E79A7", "#59A14F", "#9C755F", "#B07AA1", "#76B7B2")

# If K > length(cols), recycle safely or interpolate
if (K > length(cols)) {
  cols <- colorRampPalette(cols)(K)
}


dens_list <- lapply(thetas, density)
y_max <- max(sapply(dens_list, function(d) max(d$y)))
ylim  <- c(0, 1.1 * y_max)   # 10% headroom

# first density initializes the plot
plot(density(thetas[[1]]),
     col = cols[1],
     lwd = 2,
     xlim = c(0, 1),
     ylim = ylim,
     main = "Posterior of theta",
     xlab = "theta")

# remaining densities
for (k in 2:K) {
  lines(density(thetas[[k]]),
        col = cols[k],
        lwd = 2)
}


for (k in 1:K) {
  abline(v = params[[k]],
         col = cols[k],
         lty = 2,
         lwd = 2)
}


legend("topright",
       legend = paste0("v", 1:K),
       col    = cols,
       lwd    = 2,
       bty    = "n")




dens_list <- lapply(epsilons, density)
y_max <- max(sapply(dens_list, function(d) max(d$y)))
ylim  <- c(0, 1.1 * y_max)   # 10% headroom
# first density initializes the plot
plot(density(epsilons[[2]]),
     col = cols[2],
     lwd = 2,
     xlim = c(0, 1),
     ylim = ylim,
     main = "Posterior of theta",
     xlab = "theta")

# remaining densities
for (k in 3:K) {
  lines(density(epsilons[[k]]),
        col = cols[k],
        lwd = 2)
}




legend("topright",
       legend = paste0("epsilon", 2:K),
       col    = cols[2:K],
       lwd    = 2,
       bty    = "n")






