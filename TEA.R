remove(list = ls())


M <- 1e4  #number of particles
alpha <- beta <- 1 # uninformative prior for the beta

params <- c(0.3, 0.4, 0.5, 0.6, 0.7)
true_response <- params[length(params)]

n_versions <- length(params)


#K <- n_versions
epsilon <- 1
epsilons <- rep(epsilon, n_versions)

smc_sampler <- function(params, epsilons){
  
  
  mu_all <- list()
  v1 <- rbinom(1e1, 1, params[1])
  v2 <- rbinom(1e1, 1, params[2])
  v3 <- rbinom(1e1, 1, params[3])
  v4 <- rbinom(1e1, 1, params[4])
  v5 <- rbinom(1e1, 1, params[5])
  V <- rbind(v1, v2, v3, v4, v5)
  
  K <- nrow(V)
  k = 0
  mu  <- rbeta(M, 1, 1)
  mu_all[[length(mu_all) + 1]] <- mu
  U <- rep(0, M)
  mu_mat <- mu
  U_mat <- U
  
  while(k < K){
    k <- k+1
    U <- rbinom(M, 1, epsilons[k])
    
    mu2 <- rbeta(length(which(U == 0)), 1, 1)
    mu[which(U == 0)] <- mu2
    data <- V[k,]
    w <- sapply(mu, function(m) prod(dbinom(data, 1, m)))
    s <- sum(w)
    w <- w/s
    ESS <- 1/sum(w^2)
    mu_mat <- rbind(mu_mat, mu)
    U_mat <- rbind(U_mat, U)
    if(ESS < M/2){
      samp <- sample(1:M, M, replace = TRUE, prob = w)
      mu_mat <- mu_mat[,samp]
      U_mat <- U_mat[, samp]
      mu <- mu[samp]
      U <- U[samp]
      
    }
    mu_all[[length(mu_all) + 1]] <- mu
  }
  #return(mu_all)
  return(list(mu_all = mu_all, V = V))
}


res <- function(params, epsilons, credible_interval, reps, true_response){
  count <- 0
  bias <- vector()
  mse <- vector()
  for(rep in 1:reps){
    mu_all <- smc_sampler(params, epsilons)[[1]]
    posterior_smc <- mu_all[[length(mu_all)]]
    if(quantile(posterior_smc, (1 - credible_interval)/2) < true_response & quantile(posterior_smc, 1 - (1 - credible_interval)/2) > true_response){
      count <- count + 1
    }
    bias <- c(bias, mean(posterior_smc) - true_response)
    mse1 <- (mean(posterior_smc) - true_response)^2
    mse <- c(mse, mse1)
  }
  coverage_rate <- count/reps
  return(list(bias, coverage_rate, mse))
}

credible_interval <- 0.95
reps <- 5e2

result <- res(params, epsilons, credible_interval, reps, true_response)
bias <- mean(result[[1]]) #-0.196 (epsilon = 1) #-0.037 (epsilon = 0)
cr <- result[[2]] #0.52 (epsilon = 1) #0.936 (epsilon = 0)
mse <- mean(result[[3]]) #0.048 (epsilon = 1) #0.017 (epsilon = 0)

#sanity check

# Flatten V to get total counts
V_vec <- as.vector(V)
n_total <- length(V_vec)
y_total <- sum(V_vec)

# Exact Beta posterior
alpha_post <- 1 + y_total
beta_post  <- 1 + n_total - y_total
posterior_exact <- rbeta(1e5, alpha_post, beta_post)

# Compare SMC posterior (last mu in mu_all)
posterior_smc <- mu_all[[length(mu_all)]]

# Check means
cat("Exact posterior mean: ", mean(posterior_exact), "\n")
cat("SMC posterior mean:   ", mean(posterior_smc), "\n")

# Compare densities visually
hist(posterior_smc, breaks=50, probability=TRUE, col=rgb(1,0,0,0.5),
     main="SMC vs Exact Posterior", xlab="mu")
lines(density(posterior_exact), col="blue", lwd=2)
legend("topright", legend=c("SMC","Exact Beta"), fill=c(rgb(1,0,0,0.5),"blue"))









