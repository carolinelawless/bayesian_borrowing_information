smc_sampler_gaussian <- function(V, M, mean_theta, sd_theta, sigma, p_eps) {
  # V: list of length K, continuous observations
  # M: number of particles
  # mean_theta, sd_theta: prior mean and sd for theta
  # sigma: known standard deviation of observations
  # p_eps: prior probability of borrowing
  
  K <- length(V)
  theta_all <- vector("list", K)
  epsilon_all <- vector("list", K)
  
  # Initial particles for first version
  theta <- rnorm(M, mean = mean_theta, sd = sd_theta)
  
  for (k in 1:K) {
    data1 <- V[[k]]
    
    if (k > 1) {
      epsilon <- rbinom(M, 1, p_eps)
      ind <- which(epsilon == 0)
      theta[ind] <- rnorm(length(ind), mean = mean_theta, sd = sd_theta)
    } else {
      epsilon <- rep(NA, M)  # no borrowing for first version
    }
    
    # Gaussian log-likelihoods
    logliks <- sapply(theta, function(theta_j) {
      sum(dnorm(data1, mean = theta_j, sd = sigma, log = TRUE))
    })
    
    # Normalize weights
    maxll <- max(logliks)
    w <- exp(logliks - maxll)
    w <- w / sum(w)
    
    # Resample
    samp <- sample(1:M, M, replace = TRUE, prob = w)
    theta <- theta[samp]
    epsilon <- epsilon[samp]
    
    theta_all[[k]] <- theta
    epsilon_all[[k]] <- epsilon
  }
  
  return(list(theta_all = theta_all, epsilon_all = epsilon_all))
}

posterior_sim_gaussian <- function(params, M, B, lambda, mean_theta, sd_theta, sigma, p_eps) {
  # params: true theta values for each version
  # M: number of particles in SMC
  # B: number of simulations
  # n: vector of sample sizes per version
  # mean_theta, sd_theta: prior for theta
  # sigma: SD of continuous observations
  # p_eps: prior probability of borrowing
  
  K <- length(params)
  
  # initialize lists to store simulated estimates
  thetas <- vector("list", K)
  epsilons <- vector("list", K)
  for (k in 1:K) {
    thetas[[k]] <- numeric(0)
    epsilons[[k]] <- numeric(0)
  }
  
  n <- vector(length = K)
  for (b in 1:B) {
    # simulate data for this iteration
    n <- rpois(K, lambda) + 1
    V <- lapply(1:K, function(k) rnorm(n[k], mean = params[k], sd = sigma))
    
    # run SMC for Gaussian model
    smc_out <- smc_sampler_gaussian(V, M, mean_theta, sd_theta, sigma, p_eps)
    theta_all <- smc_out$theta_all
    epsilon_all <- smc_out$epsilon_all
    

    
    for (k in 1:K) {
      # store estimates; you can switch between mean, MAP, or BMA
      thetas[[k]] <- c(thetas[[k]], mean(theta_all[[k]]))
      epsilons[[k]] <- c(epsilons[[k]], mean(epsilon_all[[k]], na.rm = TRUE))
    }
  }
  
  return(list(thetas = thetas, epsilons = epsilons))
}

plot_gaussian_trajectories <- function(thetas, epsilons, true_theta = NULL) {
  
  K <- length(thetas)
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
  
}
