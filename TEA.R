remove(list = ls())

smc_sampler <- function(params, M, n){
  
  # helper: stable log-sum-exp
  logsumexp <- function(x){
    m <- max(x)
    return(m + log(sum(exp(x - m))))
  }
  
  epsilons <- vector("numeric", length = length(params))
  K <- length(params)
  mu_all <- list()
  
  # simulate K datasets (each a vector of 100 Bernoulli trials) as in your original
  V <- sapply(params, function(p) rbinom(n, 1, p))
  V <- t(V)   # now V is K x 100 (row = dataset)
  
  # particle population
  k <- 1
  mu  <- rbeta(M, 1, 1)     # initialize particles from Beta(1,1) prior
  mu_all[[length(mu_all) + 1]] <- mu
  U <- rep(0, M)
  mu_mat <- matrix(mu, nrow = 1)   # store successive mu vectors (rows = steps)
  U_mat <- matrix(U, nrow = 1)
  
  while(k < K){
    k <- k + 1
    
    # data for transition (k-1) -> k
    data0 <- V[k-1, ]   # used only if needed for other diagnostics
    data1 <- V[k, ]     # new dataset (vector of 0/1)
    n1 <- length(data1)
    s1 <- sum(data1)
    
    # --- compute predictive probability of data1 under current particle approximation
    # compute log-likelihood of data1 given each particle mu_i
    loglik_particles <- sapply(mu, function(m) sum(dbinom(data1, 1, m, log = TRUE)))
    # log predictive (mixture average over particles): log( (1/M) * sum_i exp(loglik_i) )
    log_pred <- logsumexp(loglik_particles) - log(length(loglik_particles))
    
    # --- compute baseline (diffuse) predictive probability p0(data1)
    # use Beta(1,1) prior predictive (marginal likelihood)
    # log p0 = lbeta(1 + s1, 1 + n1 - s1) - lbeta(1,1)
    log_p0 <- lbeta(1 + s1, 1 + n1 - s1) - lbeta(1, 1)
    
    # epsilon = p_pred / (p_pred + p0)  computed stably from logs:
    # epsilon = 1 / (1 + exp(log_p0 - log_pred))
    epsilon <- 1 / (1 + exp(log_p0 - log_pred))
    epsilons[k] <- epsilon
    
    # draw U for each particle: inherit (1) with prob epsilon, else refresh (0)
    U <- rbinom(M, 1, epsilon)
    
    # --- innovation step for particles with U==0:
    # refresh them from posterior based on data1 alone (Beta(1 + s1, 1 + n1 - s1))
    idx_zero <- which(U == 0)
    if(length(idx_zero) > 0){
      mu_new <- rbeta(length(idx_zero), 1 + s1, 1 + n1 - s1)
      mu[idx_zero] <- mu_new
    }
    # if U==1, particles keep their mu values (inherit)
    
    # --- importance weighting by likelihood of the new data under each particle
    # (use log-likelihoods for stability)
    loglik_particles <- sapply(mu, function(m) sum(dbinom(data1, 1, m, log = TRUE)))
    # convert to non-negative weights stably
    maxll <- max(loglik_particles)
    w_unnorm <- exp(loglik_particles - maxll)
    w_sum <- sum(w_unnorm)
    if(w_sum == 0 || !is.finite(w_sum)){
      # fallback: uniform weights (avoid NaNs)
      w <- rep(1 / M, M)
    } else {
      w <- w_unnorm / w_sum
    }
    
    # Effective Sample Size
    ESS <- 1 / sum(w^2)
    
    # store population
    mu_mat <- rbind(mu_mat, mu)
    U_mat <- rbind(U_mat, U)
    
    # resample if ESS low
    if(ESS < M/2){
      samp <- sample(1:M, M, replace = TRUE, prob = w)
      mu <- mu[samp]
      U <- U[samp]
      # (mu_mat/U_mat are histories; we keep them as recorded before resampling step)
    }
    
    mu_all[[length(mu_all) + 1]] <- mu
  }
  
  return(list(mu_all = mu_all, V = V, epsilons = epsilons))
}

# ------------- run example --------------
M <- 1e3  # number of particles
n <- 1e3 # sample size

params <- c(0.3, 0.4, 0.5, 0.6, 0.7)
params <- rep(0.2, 5)
result <- smc_sampler(params, M, n)
mu_all <- result$mu_all
V <- result$V
epsilons <- result$epsilons
mu_means <- sapply(mu_all, mean)
mu_means
epsilons


# Example pseudocode for calibration
simulate_S <- function(params, M, n, n_sim=500){
  S_vals <- numeric(n_sim)
  for(i in 1:n_sim){
    result <- smc_sampler(params, M, n)
    eps <- result$epsilons[-1] # remove first zero
    S_vals[i] <- sum(eps)
    print(i)
  }
  return(S_vals)
}

# Simulate under H0 (same params for all datasets)
S_H0 <- simulate_S(rep(0.5, 5), M, n)
threshold <- quantile(S_H0, 0.05)  # 5% quantile for 5% false positive rate

# Simulate under H1 (last param different)
params_H1 <- c(0.5, 0.5, 0.5, 0.5, 0.7)  # last dataset different
S_H1 <- simulate_S(params_H1, M, n)

# You can check how often S_H1 < threshold to estimate power
length(which(S_H1 < threshold))/length(S_H1)





