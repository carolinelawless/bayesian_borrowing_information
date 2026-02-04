smc_sampler <- function(V, M) {
  # V: list observations of length K
  # M: number of particles
  
  K <- length(V)
  theta_all <- list(K)   
  epsilon_all <- list(K)
  
  theta <- rbeta(M, a_theta, b_theta)

  for (k in 1:K) {
    
    if (k > 1) {
      epsilon <- rbinom(M, 1, p_eps)
      ind <- which(epsilon == 0)
      theta[ind] <- rbeta(length(ind), a_theta, b_theta)
    } else {
      epsilon <- rep(NA, M)
    }
    
    # epsilon <- rbinom(M, 1, p_eps)
    # ind <- which(epsilon == 0)
    # theta[ind] <- rbeta(length(ind), a_theta, b_theta)
    
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


naive_diff <- function(params, lambda, M, B){
  
  K <- length(params)
  diffs <- numeric(B)
  
  for(b in 1:B){
    n <- rpois(K, lambda) + 1
    
    s1 <- sum(rbinom(n[1], 1, params[1]))
    sK <- sum(rbinom(n[K], 1, params[K]))
    
    theta1 <- rbeta(M, a_theta + s1, b_theta + n[1] - s1)
    thetaK <- rbeta(M, a_theta + sK, b_theta + n[K] - sK)
    #diffs[b] <- abs(mean(thetaK) - mean(theta1))
    delta <- thetaK - theta1
    diffs[b] <- mean(delta) / sd(delta)
    
    
  }
  
  return(diffs)
}


tea_diff_mean <- function(params, lambda, M, B) {
  K <- length(params)
  diffs <- numeric(B)
  
  for (b in 1:B) {
    n <- rpois(K, lambda) + 1
    
    V <- lapply(1:K, function(k) rbinom(n[k], 1, params[k]))
    smc_out <- smc_sampler(V, M)
    theta1 <- smc_out$theta_all[[1]]
    thetaK <- smc_out$theta_all[[K]]
    
    #diffs[b] <- abs(mean(thetaK) - mean(theta1))   # posterior means of particles
    delta <- thetaK - theta1
    diffs[b] <- mean(delta) / sd(delta)
    
    
  }
  
  return(diffs)
}




naive_stopping_decision <- function(params, lambda, M, B){
  k <- 1
  new_trial <- 0
  while(new_trial == 0 & k < length(params)){
    k <- k+1
    params_subset <- params[1:k]
    diffs <- naive_diff(params_subset, lambda, M, B)
    if(sum(diffs > 0.2)/length(diffs) > 0.8 ){
      new_trial <- 1
    }
  }
  k
  return(k)
}


tea_stopping_decision <- function(params, lambda, M, B){
  k <- 1
  new_trial <- 0
  while(new_trial == 0 & k < length(params)){
    k <- k+1
    params_subset <- params[1:k]
    diffs <- tea_diff_mean(params_subset, lambda, M, B)
    if(sum(diffs > 0.2)/length(diffs) > 0.8 ){
      new_trial <- 1
    }
  }
  k
  return(k)
}

