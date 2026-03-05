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



smc_sampler_binomial <- function(V, M, a_theta, b_theta, p_eps) {
  # V: list of length K, continuous observations
  # M: number of particles
  # a_theta, b_theta: prior beta parameters for theta
  # p_eps: prior probability of borrowing
  
  K <- length(V)
  theta_all <- vector("list", K)
  epsilon_all <- vector("list", K)
  
  # Initial particles for first version
  theta <- rbeta(M, a_theta, b_theta)

  
  for (k in 1:K) {
    data1 <- V[[k]]
    
    if (k > 1) {
      epsilon <- rbinom(M, 1, p_eps)
      ind <- which(epsilon == 0)
      theta[ind] <- rbeta(length(ind), a_theta, b_theta)
    } else {
      epsilon <- rep(NA, M)  # no borrowing for first version
    }
    
    
    
  # Binomial log-likelihoods
    
    data1 <- V[[k]]
    n1 <- length(data1)
    s1 <- sum(data1)
    
    
    logliks <- dbinom(s1, n1, theta, log = TRUE)
    

    
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
  diffs <- numeric(B)
  
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
    
    theta1 <- theta_all[[1]]
    thetaK <- theta_all[[K]]
    diffs[b] <- abs(mean(thetaK) - mean(theta1))
    
    for (k in 1:K) {
      # store estimates; you can switch between mean, MAP, or BMA
      thetas[[k]] <- c(thetas[[k]], mean(theta_all[[k]]))
      epsilons[[k]] <- c(epsilons[[k]], mean(epsilon_all[[k]], na.rm = TRUE))
    }
  }
  
  return(list(thetas = thetas, epsilons = epsilons, diffs = diffs))
}


posterior_sim_binomial <- function(params, M, B, lambda, a_theta, b_theta, p_eps) {
  # params: true theta values for each version
  # M: number of particles in SMC
  # B: number of simulations
  # n: vector of sample sizes per version
  # a_theta, b_theta: Beta prior for theta
  # p_eps: prior probability of borrowing
  
  K <- length(params)
  
  # initialize lists to store simulated estimates
  thetas <- vector("list", K)
  epsilons <- vector("list", K)
  diffs <- numeric(B)
  
  for (k in 1:K) {
    thetas[[k]] <- numeric(0)
    epsilons[[k]] <- numeric(0)
  }
  
  n <- vector(length = K)
  for (b in 1:B) {
    # simulate data for this iteration
    n <- rpois(K, lambda) + 1
    
    V <- vector("list", K)
    for (k in 1:K) {
      V[[k]] <- rbinom(n[k], 1, params[k])
    }
    
    # run SMC for Gaussian model
    smc_out <- smc_sampler_binomial(V, M, a_theta, b_theta, p_eps)
    theta_all <- smc_out$theta_all
    epsilon_all <- smc_out$epsilon_all
    
    theta1 <- theta_all[[1]]
    thetaK <- theta_all[[K]]
    diffs[b] <- abs(mean(thetaK) - mean(theta1))
    
    for (k in 1:K) {
      # store estimates; you can switch between mean, MAP, or BMA
      thetas[[k]] <- c(thetas[[k]], mean(theta_all[[k]]))
      epsilons[[k]] <- c(epsilons[[k]], mean(epsilon_all[[k]], na.rm = TRUE))
    }
  }
  
  return(list(thetas = thetas, epsilons = epsilons, diffs = diffs))
}



plot_trajectories <- function(thetas, epsilons, params) {
  
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






compute_power_binomial <- function(params_list,
                                   lambdas,
                                   threshold,
                                   M, B,
                                   a_theta, b_theta,
                                   p_eps) {
  
  n_scen <- length(params_list)
  power <- vector("list", n_scen)
  
  for (s in seq_len(n_scen)) {
    
    params <- params_list[[s]]
    power_s <- numeric(length(lambdas))
    
    for (i in seq_along(lambdas)) {
      
      lambda <- lambdas[i]
      
      res <- posterior_sim_binomial(params, M, B,
                                    lambda, a_theta,
                                    b_theta, p_eps)
      
      power_s[i] <- mean(res$diffs > threshold)
    }
    
    power[[s]] <- power_s
  }
  
  return(power)
}


compute_power_gaussian <- function(params_list,
                                   lambdas,
                                   threshold,
                                   M, B,
                                   mean_theta, sd_theta, sigma,
                                   p_eps) {
  
  n_scen <- length(params_list)
  power <- vector("list", n_scen)
  
  for (s in seq_len(n_scen)) {
    
    params <- params_list[[s]]
    power_s <- numeric(length(lambdas))
    
    for (i in seq_along(lambdas)) {
      
      lambda <- lambdas[i]
      
      res <- posterior_sim_gaussian(params, M, B,
                                    lambda, mean_theta,
                                    sd_theta, sigma, p_eps)
      
      power_s[i] <- mean(res$diffs > threshold)
    }
    
    power[[s]] <- power_s
  }
  
  return(power)
}



plot_power_curves <- function(lambdas,
                              power,
                              labels,
                              hypothesis) {
  
  if (hypothesis == 0) {
    ylab <- "type-1 error"
    legend_pos <- "topright"
  } else {
    ylab <- "power"
    legend_pos <- "bottomright"
  }
  
  n_curves <- ncol(power)
  cols <- seq_len(n_curves)
  
  plot(lambdas, power[,1],
       type = "l",
       lwd = 2,
       col = cols[1],
       ylim = c(0,1),
       xlab = expression(lambda),
       ylab = ylab)
  
  if (n_curves > 1) {
    for (i in 2:n_curves) {
      lines(lambdas, power[,i],
            lwd = 2,
            col = cols[i])
    }
  }
  
  legend(legend_pos,
         legend = labels,
         col = cols,
         lwd = 2,
         bty = "n")
}





