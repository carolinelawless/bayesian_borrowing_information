smc_sampler <- function(V, M) {
  # V: list observations of length K
  # M: number of particles
  
  K <- length(V)
  theta_all <- list(K)   
  epsilon_all <- list(K)
  
  theta <- rbeta(M, a_theta, b_theta)

  for (k in 1:K) {
    
    epsilon <- rbinom(M, 1, p_eps)
    ind <- which(epsilon == 0)
    theta[ind] <- rbeta(length(ind), a_theta, b_theta)
    
    data1 <- V[[k]]
    n1 <- length(data1)
    s1 <- sum(data1)
    

    logliks <- dbinom(s1, n1, theta, log = TRUE)
    
    maxll <- max(logliks)
    w <- exp(logliks - maxll)
    w <- w / sum(w)
    
    samp <- sample(1:M, M, replace = TRUE, prob = w)
    theta <- theta[samp]
    epsilon <- epsilon[samp]
    
    
    theta_all[[k]] <- theta
    epsilon_all[[k]] <- epsilon
}
  
  return(list(theta_all = theta_all, epsilon_all = epsilon_all))
}


posterior_sim <- function(params, M, B){
  
  K <- length(params)
  
  # initialize list of theta vectors
  thetas <- vector("list", K)
  epsilons <- vector("list", K)
  
  for (k in 1:K) {
    thetas[[k]] <- numeric(0)
    epsilons[[k]] <- numeric(0)
  }
  
  for (b in 1:B) {
    
    V <- vector("list", K)
    for (k in 1:K) {
      V[[k]] <- rbinom(n[k], 1, params[k])
    }
    
    smc_out <- smc_sampler(V, M)
    theta_all <- smc_out$theta_all
    epsilon_all <- smc_out$epsilon_all
    
    for (k in 1:K) {
      thetas[[k]] <- c(thetas[[k]], mean(theta_all[[k]]))
      epsilons[[k]] <- c(epsilons[[k]], mean(epsilon_all[[k]]))
    }
  }
  
  return(list(thetas, epsilons))
}



posterior_predictive_sim_diff <- function(params, M, B, sampsize){
  diffs <- vector()
  for(b in 1:B){
    V <- list()
    for(k in 1:K){
      V[[k]] <- rbinom(n[k], 1, params[k])
    }
    smc_out <- smc_sampler(V, M)
    theta_all <- smc_out$theta_all
    theta1 <- theta_all[[1]]
    thetaK <- theta_all[[K]]

    samp1 <- sample(1:M, sampsize, replace = TRUE)
    sampk <- sample(1:M, sampsize, replace = TRUE)  
    theta1_samp <- theta1[samp1]
    thetaK_samp <- thetaK[sampk]
      
    y1_samp <- rbinom(sampsize, 1, theta1_samp)
    yK_samp <- rbinom(sampsize, 1, thetaK_samp)
      
    diffs[b] <- mean(yK_samp) - mean(y1_samp)  # e.g. current minus first
      
  }
  return(diffs)
}












