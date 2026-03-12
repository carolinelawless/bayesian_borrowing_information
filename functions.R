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
       ylab = "θ / ε",
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
  legend("bottomright",
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



posterior_sim_naive_binomial <- function(params, M, B, lambda, a_theta, b_theta) {

  
  K <- length(params)
  
  thetas <- vector("list", K)
  diffs <- numeric(B)
  
  for (k in 1:K) {
    thetas[[k]] <- numeric(0)
  }
  
  for (b in 1:B) {
    
    # sample sizes
    n <- rpois(K, lambda) + 1
    
    # simulate data
    V <- vector("list", K)
    for (k in 1:K) {
      V[[k]] <- rbinom(n[k], 1, params[k])
    }
    

    
    # store means
    for(k in 1:K){
      y1 <- sum(V[[k]])
      n1 <- n[k]
      theta_post <- rbeta(M, a_theta + y1, b_theta + n1 - y1)
      thetas[[k]] <- c(thetas[[k]], mean(theta_post))
      if(k == 1){theta1_post <- theta_post}
      if(k == K){thetaK_post <- theta_post}
    }

    # difference
   diffs[b] <- abs(mean(thetaK_post) - mean(theta1_post))
  }
  
  return(list(
    thetas = thetas, 
    diffs = diffs
  ))
}


posterior_sim_naive_gaussian <- function(params, M, B, lambda, mean_theta, sd_theta, sigma) {

  
  K <- length(params)
  
  thetas <- vector("list", K)
  diffs <- numeric(B)
  
  for (k in 1:K) {
    thetas[[k]] <- numeric(0)
  }
  
  sigma0_sq <- sd_theta^2
  
  for (b in 1:B) {
    
    # sample sizes
    n <- rpois(K, lambda) + 1
    
    # simulate data
    V <- vector("list", K)
    for (k in 1:K) {
      V[[k]] <- rnorm(n[k], params[k], sigma)
    }
    
    # version 1
    x1 <- V[[1]]
    n1 <- length(x1)
    
    post_var1 <- 1 / (n1/sigma^2 + 1/sigma0_sq)
    post_mean1 <- post_var1 * (sum(x1)/sigma^2 + mean_theta/sigma0_sq)
    
    theta1_post <- rnorm(M, post_mean1, sqrt(post_var1))
    
    # version K
    xK <- V[[K]]
    nK <- length(xK)
    
    post_varK <- 1 / (nK/sigma^2 + 1/sigma0_sq)
    post_meanK <- post_varK * (sum(xK)/sigma^2 + mean_theta/sigma0_sq)
    
    thetaK_post <- rnorm(M, post_meanK, sqrt(post_varK))
    
    # store means
    thetas[[1]] <- c(thetas[[1]], mean(theta1_post))
    thetas[[K]] <- c(thetas[[K]], mean(thetaK_post))
    
    diffs[b] <- abs(mean(thetaK_post) - mean(theta1_post))
  }
  
  return(list(
    thetas = thetas,
    diffs = diffs
  ))
}


TEA_eval_binomial <- function(params, M, B, lambdas, a_theta, b_theta, p_eps){
  L <- length(lambdas)
  variance_TEA <- vector(length = L)
  bias_TEA <- vector(length = L)
  MSE_TEA <- vector(length = L)
  variance_naive <- vector(length = L)
  bias_naive <- vector(length = L)
  MSE_naive <- vector(length = L)
  
  for(i in 1:length(lambdas)){
    
    print(i)
    
    lambda <- lambdas[i]
    
    res_TEA <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, p_eps)
    variance_TEA[i] <- var(res_TEA$diffs)
    bias_TEA[i] <- abs(mean(res_TEA$diffs - abs(params[length(params)] - params[1])))
    MSE_TEA[i] <- variance_TEA[i] + bias_TEA[i]*bias_TEA[i]
    
    
    res_naive <- posterior_sim_naive_binomial(params, M, B, lambda, a_theta, b_theta)
    variance_naive[i] <- var(res_naive$diffs)
    bias_naive[i] <- abs(mean(res_naive$diffs - abs(params[length(params)] - params[1])))
    MSE_naive[i] <- variance_naive[i] + bias_naive[i]*bias_naive[i]
  }
  
  return(list(variance_TEA = variance_TEA, bias_TEA = bias_TEA, MSE_TEA = MSE_TEA, variance_naive = variance_naive, bias_naive = bias_naive, MSE_naive = MSE_naive))
}


TEA_eval_gaussian <- function(params, M, B, lambdas, mean_theta, sd_theta, sigma, p_eps){
  L <- length(lambdas)
  variance_TEA <- vector(length = L)
  bias_TEA <- vector(length = L)
  MSE_TEA <- vector(length = L)
  variance_naive <- vector(length = L)
  bias_naive <- vector(length = L)
  MSE_naive <- vector(length = L)
  
  for(i in 1:length(lambdas)){
    
    print(i)
    
    lambda <- lambdas[i]
    
    res_TEA <- posterior_sim_gaussian(params, M, B, lambda, mean_theta, sd_theta, sigma, p_eps)
    variance_TEA[i] <- var(res_TEA$diffs)
    bias_TEA[i] <- abs(mean(res_TEA$diffs - abs(params[length(params)] - params[1])))
    MSE_TEA[i] <- variance_TEA[i] + bias_TEA[i]*bias_TEA[i]
    
    
    res_naive <- posterior_sim_naive_gaussian(params, M, B, lambda, mean_theta, sd_theta, sigma)
    variance_naive[i] <- var(res_naive$diffs)
    bias_naive[i] <- abs(mean(res_naive$diffs - abs(params[length(params)] - params[1])))
    MSE_naive[i] <- variance_naive[i] + bias_naive[i]*bias_naive[i]
  }
  
  return(list(variance_TEA = variance_TEA, bias_TEA = bias_TEA, MSE_TEA = MSE_TEA, variance_naive = variance_naive, bias_naive = bias_naive, MSE_naive = MSE_naive))
}



plot_model_comparison <- function(lambdas,
                                  tea_values,
                                  naive_values,
                                  stat_name) {
  
  # basic checks
  if(length(lambdas) != length(tea_values) |
     length(lambdas) != length(naive_values)){
    stop("All input vectors must have the same length.")
  }
  
  # plot TEA
  plot(lambdas, tea_values,
       type = "l",
       lwd = 2,
       col = "1",
       xlab = expression(lambda),
       ylab = stat_name,
       ylim = range(c(tea_values, naive_values)))
  
  # add naive
  lines(lambdas, naive_values,
        lwd = 2,
        col = "2")
  
  # legend
  legend("topright",
         legend = c("TEA", "Naive"),
         col = c("blue", "red"),
         lwd = 2,
         bty = "n")
}




