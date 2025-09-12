#########################################################################
# Location commensurate power prior (LCPP)
# CL 10/09/2025
#########################################################################
#test
#test3

remove(list = ls())

library(rjags)
library(coda)

# ------------------------------
# Data setup
# ------------------------------


n0 <- 10
n <- 10
mu0 <- 0
sigma0 <- 1
mu <- 0
sigma <- 1



g <- function(tau){
  return(max(log(tau), 1))
}






# ------------------------------
# JAGS model
# ------------------------------
model_string <- "
model {
  # --- Prior on commensurability parameter tau ---
  logtau ~ dt(0, pow(30,-2), 1)   # Cauchy(0,30) prior on logtau
  tau <- exp(logtau)

  # --- Define g(tau) = max(log(tau), 1) ---
  logtau_clip <- logtau * step(logtau - 1) + 1 * (1 - step(logtau - 1))
  gtau <- logtau_clip   # equivalent to max(logtau, 1)

  # --- Power prior weight ---
  alpha ~ dbeta(gtau, 1)

  # --- Prior on sampling variance ---
  tau_sigma ~ dgamma(1.0E-3, 1.0E-3)   # vague prior on precision
  sigma <- 1 / sqrt(tau_sigma)

  # --- Commensurate prior for mu ---
  precision <- 1 / (1/tau + Y0_sd*Y0_sd/alpha*n0)
  mu ~ dnorm(Y0_mean, precision)

  # --- Likelihood for current data ---
  for (i in 1:n) {
    Y[i] ~ dnorm(mu, tau_sigma)
  }
}
"
results <- vector()

nsim <- 1e2


for(i in 1:nsim){
  
  print(paste0("iteration ", i))
  
  # Historical sample 
  Y0 <- rnorm(n0, mu0, sigma0)
  Y0_mean <- mean(Y0)
  Y0_sd <- sqrt(var(Y0))
  
  # Current sample 
  Y  <- rnorm(n, mu, sigma)
  
  # JAGS data list
  data_jags <- list(
    Y0_mean = Y0_mean,
    Y0_sd   = Y0_sd,
    Y       = Y,
    n       = length(Y),
    n0      = n0
  )
  model <- jags.model(textConnection(model_string),
                      data = data_jags,
                      n.chains = 3,
                      n.adapt = 1000)
  
  update(model, 1000)  # burn-in
  samples <- coda.samples(model,
                          variable.names = c("mu", "alpha", "tau", "sigma"),
                          n.iter = 5000)
  samples_mat <- as.matrix(samples)
  if(quantile(samples_mat[, "mu"], 0.1) <=0){
    res <- "accept"
  }else{res <- "reject"}
  results <- c(results, res)
}
results

sum(results == "reject")/length(results)



