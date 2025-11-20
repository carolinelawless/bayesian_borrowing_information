#########################################################################
# Hierarchical commensurate prior for multiple historical binomial datasets
#########################################################################

remove(list = ls())
library(rjags)
library(coda)
library(tictoc)

params <- (3:7)/10
n <- 1000


model_string <- "
model {

  for (j in 1:J) {
    for (i in 1:n) {
      y[j,i] ~ dbern(params[j])
    }
  }

  for (j in 1:J) {
    logit(params[j]) <- theta[j]
  }


  theta[1] ~ dnorm(0, 0.01)

  for (j in 2:J) {
    theta[j] ~ dnorm(theta[j-1], tau[j])
  }


  for (j in 2:J) {
    sigma[j] ~ dunif(0, 5)
    tau[j] <- pow(sigma[j], -2)
  }

}
"

jags_sampler <- function(params, n, model_string){
  
  J <- length(params)
  
  # --- Simulate data ---
  y <- matrix(NA, nrow = J, ncol = n)
  for(j in 1:J){
    y[j, ] <- rbinom(n, 1, params[j])
  }

  
  # --- JAGS data list ---
  data_jags <- list(
    y = y,
    n = n,
    J = J
  )
  


# --- Build model ---
model <- jags.model(
  textConnection(model_string),
  data = data_jags,
  n.chains = 3,
  n.adapt = 1000
)

# Burn-in
update(model, 1000)

# Sample
samples <- coda.samples(
  model,
  variable.names = c(
    "params",
    "theta",
    "tau"
  ),
  n.iter = 10000
)

samples_mat <- as.matrix(samples)
return(samples_mat)
}









simulate_S <- function(params, n, n_sim = 10){
  S_vals <- numeric(n_sim)
  J <- length(params)
  for(i in 1:n_sim){
    print(i)
    samples_mat <- jags_sampler(params, n, model_string)
    S <- 0 
    for(j in 2:J){
      S <- S + mean(samples_mat[, paste0("tau[",j,"]")])
    }
    S_vals[i] <- S
  }

  return(S_vals)
}

t1 <- proc.time()
# Simulate under H0 (same params for all datasets)
S_H0 <- simulate_S(rep(0.5, 5), n)
threshold <- quantile(S_H0, 0.05)  # 5% quantile for 5% false positive rate

# Simulate under H1 (last param different)
params_H1 <- c(0.5, 0.5, 0.5, 0.5, 0.7)  # last dataset different
S_H1 <- simulate_S(params_H1, n)

# You can check how often S_H1 < threshold to estimate power
length(which(S_H1 < threshold))/length(S_H1)
t2 <- proc.time()
t2 - t1


#Bayesian predictive 

predictive_similarity <- function(samples_mat, delta = 0.05, n_pred = 1000) {
  # samples_mat: the MCMC sample matrix returned by jags_sampler()
  # delta: similarity tolerance (e.g., 0.05)
  # n_pred: number of predictive draws to approximate P(similarity)
  
  # Extract posterior samples for theta[J]
  theta_cols <- grep("^theta\\[", colnames(samples_mat))
  theta_samples <- samples_mat[, theta_cols, drop = FALSE]
  
  J <- ncol(theta_samples)  # number of historical datasets
  
  theta_1 <- theta_samples[, 1]   # latent parameter for first historical version
  
  # Extract tau[J+1] if present, else create a prior for predictive variance
  tau_cols <- grep("^tau\\[", colnames(samples_mat))
  tau_samples <- samples_mat[, tau_cols, drop = FALSE]
  
  # If your model does NOT include tau[J+1] (e.g., only tau[2:J]),
  # then use the last tau as commensurate precision for predictive step.
  tau_for_pred <- tau_samples[, ncol(tau_samples)]
  
  # Number of posterior draws
  S <- length(theta_J)
  
  # For posterior predictive, draw theta_{J+1} from:
  # theta_{J+1} ~ Normal(theta_J, tau_for_pred^{-1})
  theta_pred <- rnorm(S, mean = theta_J, sd = 1 / sqrt(tau_for_pred))
  
  # Convert theta to probabilities
  p_1     <- plogis(theta_1)
  p_Jplus <- plogis(theta_pred)
  
  # Compute predictive probability of similarity
  similarity_prob <- mean(abs(p_Jplus - p_1) < delta)
  
  return(similarity_prob)
}



