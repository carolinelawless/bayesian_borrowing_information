#########################################################################
# Hierarchical commensurate prior for multiple historical binomial datasets
#########################################################################

remove(list = ls())
library(rjags)
library(coda)
library(parallel)




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

#jags_sampler <- function(params, n, model_string){
jags_sampler <- function(y, n, model_string){
  
  J <- nrow(y)
  
  # # --- Simulate data ---
  # y <- matrix(NA, nrow = J, ncol = n)
  # for(j in 1:J){
  #   y[j, ] <- rbinom(n, 1, params[j])
  # }

  
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



# --- Simulate data ---
params <- rep(0.5, 5)
J <- length(params)
n <- 1e2
y <- matrix(NA, nrow = J, ncol = n)

for(j in 1:J){
  y[j, ] <- rbinom(n, 1, params[j])
}






simulate_S <- function(params, n, n_sim = 1e3){
  S_vals <- numeric(n_sim)
  J <- length(params)
  for(i in 1:n_sim){
    print(i)
    y <- matrix(NA, nrow = J, ncol = n)
    for(j in 1:J){
      y[j, ] <- rbinom(n, 1, params[j])
    }
    samples_mat <- jags_sampler(y, n, model_string)
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
threshold #20442 21775 17194 1991 (1e3 simulations) 19311 (1e3 sim)
t2 <- proc.time()
t2 - t1


# Simulate under H1 (last param different)
params_H1 <- c(0.5, 0.5, 0.5, 0.5, 0.7)  # last dataset different
S_H1 <- simulate_S(params_H1, n)

# You can check how often S_H1 < threshold to estimate power
length(which(S_H1 < threshold))/length(S_H1) #0.124 0.132
t2 <- proc.time()
t2 - t1


# Simulate under H1 
params_H1 <- c(0.5, 0.56, 0.6, 0.65, 0.7)  
S_H1 <- simulate_S(params_H1, n)

# You can check how often S_H1 < threshold to estimate power
length(which(S_H1 < threshold))/length(S_H1) #0.11
t2 <- proc.time()
t2 - t1

# Simulate under H1 (last param different)
params_H1 <- c(0.5, 0.5, 0.5, 0.5, 0.9)  # last dataset different
S_H1 <- simulate_S(params_H1, n)

# You can check how often S_H1 < threshold to estimate power
length(which(S_H1 < threshold))/length(S_H1) #0.147 0.141
t2 <- proc.time()
t2 - t1


# Simulate under H1 
params_H1 <- c(0.5, 0.6, 0.7, 0.8, 0.9)  
S_H1 <- simulate_S(params_H1, n)

# You can check how often S_H1 < threshold to estimate power
length(which(S_H1 < threshold))/length(S_H1) #0.371
t2 <- proc.time()
t2 - t1


#Bayesian predictive 


posterior_predictive_diff <- function(samples_mat, n, n_pred = 2e3) {
  p_cols <- grep("^params\\[", colnames(samples_mat))
  p_samples <- samples_mat[, p_cols, drop = FALSE]
  p1 <- p_samples[,1]
  pK <- p_samples[,ncol(p_samples)]
  diffs <- numeric(n_pred)
  for(b in 1:n_pred){
    
    p1_samp <- sample(p1, 1)
    pK_samp <- sample(pK, 1)
    
    y1_rep <- rbinom(n, 1, p1_samp)
    yK_rep <- rbinom(n, 1, pK_samp)
    diffs[b] <- mean(yK_rep) - mean(y1_rep)
  }
  return(diffs)
}


#to do:
# speed up (parallel? optimise?)
# write-up results. simulation set-up, scenarios, results, discussion





