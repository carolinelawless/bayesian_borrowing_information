smc_sampler <- function(V, M = 1e3) {
  # V: K x n matrix (rows = datasets/versions)
  # M: number of particles
  
  logsumexp <- function(x){
    m <- max(x)
    m + log(sum(exp(x - m)))
  }
  
  K <- nrow(V)
  mu_all <- vector("list", K)   # store particle populations AFTER resampling (posterior approx)
  w_all  <- vector("list", K)   # store pre-resampling normalized weights
  epsilons <- numeric(K - 1)
  
  # initial particles from prior Beta(1,1)
  mu <- rbeta(M, 1, 1)
  
  # ---- process dataset 1 ----
  data1 <- V[1, ]
  n1 <- length(data1)
  s1 <- sum(data1)
  
  loglik <- sapply(mu, function(m) sum(dbinom(data1, 1, m, log = TRUE)))
  maxll <- max(loglik)
  w <- exp(loglik - maxll)
  w <- w / sum(w)
  w_all[[1]] <- w                 # STORE pre-resampling weights
  
  ESS <- 1 / sum(w^2)
  if (ESS < M/2){
    mu <- mu[sample.int(M, M, replace = TRUE, prob = w)]
  }
  mu_all[[1]] <- mu

  
  # ---- loop k = 2..K ----
  for (k in 2:K) {
    
    data_k <- V[k, ]
    n_k <- length(data_k)
    s_k <- sum(data_k)
    
    # (1) predictive under posterior (particles mu represent posterior after k-1)
    loglik_post <- sapply(mu, function(m) sum(dbinom(data_k, 1, m, log = TRUE)))
    log_pred <- logsumexp(loglik_post) - log(M)
    
    # (2) baseline prior predictive
    log_p0 <- lbeta(1 + s_k, 1 + n_k - s_k) - lbeta(1, 1)
    
    # (3) epsilon
    epsilon <- 1 / (1 + exp(log_p0 - log_pred))
    epsilons[k - 1] <- epsilon
    
    # (4) forgetting U
    U <- rbinom(M, 1, epsilon)
    
    # (5) refresh U==0 from posterior of data_k
    idx_zero <- which(U == 0)
    if (length(idx_zero) > 0) {
      mu[idx_zero] <- rbeta(length(idx_zero), 1 + s_k, 1 + n_k - s_k)
    }
    
    # (6) importance weights (pre-resampling) for data_k
    loglik <- sapply(mu, function(m) sum(dbinom(data_k, 1, m, log = TRUE)))
    maxll <- max(loglik)
    w <- exp(loglik - maxll)
    w <- w / sum(w)
    w_all[[k]] <- w                 # STORE pre-resampling weights
    
    # (7) resample if ESS low
    ESS <- 1 / sum(w^2)
    if (ESS < M/2){
      mu <- mu[sample.int(M, M, replace = TRUE, prob = w)]
    }
    
    mu_all[[k]] <- mu
  }
  
  return(list(mu_all = mu_all, w_all = w_all, epsilons = epsilons, V = V))
}


simulate_S <- function(params, n, M = 1e3, n_sim = 1e2){
  S_vals <- numeric(n_sim)
  epsilons1 <- matrix(nrow = (length(params) - 1), ncol = n_sim)
  for(i in 1:n_sim){
    V <- sapply(params, function(p) rbinom(n, 1, p))
    V <- t(V)   # now V is K x n (row = dataset)
    result <- smc_sampler(V, M)
    eps <- result$epsilons
    S_vals[i] <- sum(eps)
    epsilons1[,i] <- eps
    print(i)
  }
  epsilons2 <- apply(epsilons1, 1, mean)
  return(list(S_vals, epsilons2))
}


posterior_predictive_diff <- function(smc_out, n, n_pred = 2e3) {
  # smc_out: output of smc_sampler(V, M)
  # n: dataset size (number of observations per version row in V)
  # n_pred: number of posterior predictive replicates
  
  mu_all <- smc_out$mu_all
  w_all  <- smc_out$w_all
  K <- length(mu_all)
  
  # posterior represented after dataset 1 is encoded in mu_all[[1]] (resampled),
  # pre-resample weights for that step are w_all[[1]]
  mu1 <- mu_all[[1]]
  w1  <- w_all[[1]]
  
  muK <- mu_all[[K]]
  wK  <- w_all[[K]]
  
  diffs <- numeric(n_pred)
  for (b in 1:n_pred) {
    # draw a particle index according to pre-resampling posterior weights
    i1 <- sample.int(length(mu1), size = 1, prob = w1)
    iK <- sample.int(length(muK), size = 1, prob = wK)
    
    theta1 <- mu1[i1]
    thetaK <- muK[iK]
    
    y1_rep <- rbinom(n, 1, theta1)
    yK_rep <- rbinom(n, 1, thetaK)
    
    diffs[b] <- mean(yK_rep) - mean(y1_rep)  # e.g. current minus first
  }
  
  return(diffs)
}

