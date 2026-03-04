smc_sampler_binomial <- function(V, M, a_theta, b_theta, p_eps) {
  # V: list of length K, binary observations
  # M: number of particles
  # a_theta, b_theta: Beta prior parameters
  # p_eps: prior probability of borrowing
  
  K <- length(V)
  theta_all <- vector("list", K)
  epsilon_all <- vector("list", K)
  
  # Initial particles for first version
  theta <- rbeta(M, a_theta, b_theta)
  
  for (k in 1:K) {
    data1 <- V[[k]]
    n1 <- length(data1)
    s1 <- sum(data1)
    
    if (k > 1) {
      epsilon <- rbinom(M, 1, p_eps)
      ind <- which(epsilon == 0)
      theta[ind] <- rbeta(length(ind), a_theta, b_theta)
    } else {
      epsilon <- rep(NA, M)  # no borrowing for first version
    }
    
    # Binomial log-likelihoods for each particle
    logliks <- sapply(theta, function(theta_j) {
      dbinom(s1, n1, theta_j, log = TRUE)
    })
    
    # Normalize weights
    maxll <- max(logliks)
    w <- exp(logliks - maxll)
    w <- w / sum(w)
    
    # Resample particles
    samp <- sample(1:M, M, replace = TRUE, prob = w)
    theta <- theta[samp]
    epsilon <- epsilon[samp]
    
    theta_all[[k]] <- theta
    epsilon_all[[k]] <- epsilon
  }
  
  return(list(theta_all = theta_all, epsilon_all = epsilon_all))
}


posterior_sim_binomial <- function(params, M, B, lambda, a_theta, b_theta, p_eps) {
  # params: true theta values for each version
  # M: number of particles
  # B: number of simulations
  # lambda: Poisson mean for sample size per version
  # a_theta, b_theta: Beta prior
  # p_eps: prior probability of borrowing
  
  K <- length(params)
  
  thetas <- vector("list", K)
  epsilons <- vector("list", K)
  for (k in 1:K) {
    thetas[[k]] <- numeric(0)
    epsilons[[k]] <- numeric(0)
  }
  
  for (b in 1:B) {
    n <- rpois(K, lambda) + 1
    V <- lapply(1:K, function(k) rbinom(n[k], 1, params[k]))
    
    smc_out <- smc_sampler_binomial(V, M, a_theta, b_theta, p_eps)
    theta_all <- smc_out$theta_all
    epsilon_all <- smc_out$epsilon_all
    
    for (k in 1:K) {
      thetas[[k]] <- c(thetas[[k]], mean(theta_all[[k]]))
      epsilons[[k]] <- c(epsilons[[k]], mean(epsilon_all[[k]], na.rm = TRUE))
    }
  }
  
  return(list(thetas = thetas, epsilons = epsilons))
}


plot_binomial_trajectories <- function(thetas, epsilons, true_theta = NULL) {
  K <- length(thetas)
  x <- 1:K
  x_eps <- x - 0.1  # shift epsilon slightly left
  
  # Summaries
  means_theta <- sapply(thetas, mean, na.rm = TRUE)
  q10_theta   <- sapply(thetas, quantile, probs = 0.1, na.rm = TRUE)
  q90_theta   <- sapply(thetas, quantile, probs = 0.9, na.rm = TRUE)
  
  eps_means <- sapply(epsilons, mean, na.rm = TRUE)
  eps_q10   <- sapply(epsilons, quantile, probs = 0.1, na.rm = TRUE)
  eps_q90   <- sapply(epsilons, quantile, probs = 0.9, na.rm = TRUE)
  
  # Theta panel
  plot(x, means_theta,
       ylim = c(0, 1),
       pch = 16, col = "black",
       ylab = "Parameter value",
       xlab = "Version",
       xaxt = "n")
  axis(1, at = x, labels = x)
  
  # Theta CI
  arrows(x, q10_theta, x, q90_theta,
         angle = 90, code = 3, length = 0.05)
  
  # True theta points
  if (!is.null(true_theta)) {
    points(x, true_theta, pch = 16, col = rgb(0,0,1,0.4), cex = 1.8)
  }
  
  # Epsilon panel (shifted)
  points(x_eps[-1], eps_means[-1], pch = 4, col = "red", cex = 1.2)
  arrows(x_eps[-1], eps_q10[-1], x_eps[-1], eps_q90[-1],
         angle = 90, code = 3, length = 0.05, col = "red")
  
  # Legend
  legend("topright",
         legend = c("theta posterior mean", "true theta", "epsilon posterior mean"),
         pch = c(16, 16, 4),
         col = c("black", rgb(0,0,1,0.4), "red"),
         bty = "n")
}
