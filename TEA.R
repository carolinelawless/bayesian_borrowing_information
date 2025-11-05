remove(list = ls())


smc_sampler <- function(params, epsilons, M){
  
  K <- length(params)
  mu_all <- list()
  
  V <- sapply(params, function(p) rbinom(100, 1, p))
  V <- t(V)
  

  k = 1
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


res <- function(params, epsilons, credible_interval, reps, true_response, M){
  count <- 0
  bias <- vector()
  mse <- vector()
  for(rep in 1:reps){
    mu_all <- smc_sampler(params, epsilons, M)[[1]]
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


M <- 1e3  #number of particles
params <- c(0.3, 0.4, 0.5, 0.6, 0.7)
true_response <- params[length(params)]

epsilon <- 0
epsilons <- rep(epsilon, length((params)))

credible_interval <- 0.95
reps <- 5e2

result <- res(params, epsilons, credible_interval, reps, true_response, M)
bias <- mean(result[[1]]) #-0.196 (epsilon = 1) #-0.037 (epsilon = 0)
cr <- result[[2]] #0.52 (epsilon = 1) #0.936 (epsilon = 0)
mse <- mean(result[[3]]) #0.048 (epsilon = 1) #0.017 (epsilon = 0)



bias
cr
mse


#simulations


# --- Simulation setup ---
M <- 1e3
params <- c(0.3, 0.4, 0.5, 0.6, 0.7)
true_response <- params[length(params)]
credible_interval <- 0.95
reps <- 3e2  # you can increase later

# --- Range of epsilon values to test ---
eps_grid <- seq(0, 1, by = 0.25)

# --- Run simulation for each epsilon value ---
library(dplyr)

results <- lapply(eps_grid, function(eps) {
  epsilons <- rep(eps, length(params))
  res_out <- res(params, epsilons, credible_interval, reps, true_response, M)
  data.frame(
    epsilon = eps,
    bias = mean(res_out[[1]]),
    coverage = res_out[[2]],
    mse = mean(res_out[[3]])
  )
})

results_df <- bind_rows(results)
print(results_df)





#sanity check



smc_sample <- smc_sampler(params, epsilons, M)
mu_all <- smc_sample[[1]]
V <- smc_sample[[2]]


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









