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



posterior_predictive_diff <- function(smc_out, B) {
  # smc_out: output of smc_sampler(V, M, lambda)
  # n: dataset size (number of observations per version row in V)
  # npred: number of posterior predictive replicates
  
  theta_all <- smc_out$theta_all

  K <- length(theta_all)
  M <- length(theta_all[[1]])

  theta1 <- theta_all[[1]]

  
  thetaK <- theta_all[[K]]

  
  diffs <- numeric(B)
  for (b in 1:B) {
    # draw a particle index according to pre-resampling posterior weights
    samp1 <- sample(1:M, 1, replace = TRUE)
    sampK <- sample(1:M, 1, replace = TRUE)
    
    theta1_samp <- theta1[samp1]
    thetaK_samp <- thetaK[sampK]
    
    y1_rep <- rbinom(n[K], 1, theta1_samp)
    yK_rep <- rbinom(n[K], 1, thetaK_samp)
    
    diffs[b] <- mean(yK_rep) - mean(y1_rep)  # e.g. current minus first
  }
  
  return(diffs)
}



# simulate_eps_samples <- function(params, n, M, B){  
#   J <- length(params)  
#   epsilons_vector <- vector()
#   
#   for(b in 1:B){
#     print(b)
#     V <- list()
#     for(j in 1:J){
#       V[[j]] <- rbinom(n[j], 1, params[j])
#     }
#     
#     smc_res <- smc_sampler(V, M)
#     epsilon_all <- smc_res$epsilon_all
#     w_all <- smc_res$w_all
#     
#     epsilons <- epsilon_all[[J]]
#     weights <- w_all[[J]]
#     
#     samp <- sample(epsilons, size = M, replace = TRUE, prob = weights)
#     epsilons_vector <- c(epsilons_vector, mean(samp))
#   }
# 
#   return(epsilons_vector)
# 
# }







