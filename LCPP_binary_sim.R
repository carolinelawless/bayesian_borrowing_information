#########################################################################
# Function: Run hierarchical commensurate (Markov) prior model
#########################################################################

library(rjags)
library(coda)

run_commensurate_model <- function(p_hist, p_current, n = 10,
                                   credible_interval = 0.95,
                                   reps = 100, n_iter = 5000) {
  
  # Storage
  bias_vec <- c()
  mse_vec <- c()
  coverage_count <- 0
  
  # Loop over repeated simulations
  for (rep in 1:reps) {
    # --- Simulate data ---
    y_hist <- matrix(NA, nrow = length(p_hist), ncol = n)
    for (j in 1:length(p_hist)) {
      y_hist[j, ] <- rbinom(n, 1, p_hist[j])
    }
    y_current <- rbinom(n, 1, p_current)
    
    # --- Data for JAGS ---
    data_jags <- list(
      y_hist = y_hist,
      y_current = y_current,
      n = n,
      J = length(p_hist)
    )
    
    # --- JAGS model ---
    model_string <- "
    model {
      # Likelihood
      for (j in 1:J) {
        for (i in 1:n) {
          y_hist[j, i] ~ dbern(p_hist[j])
        }
      }
      for (i in 1:n) {
        y_current[i] ~ dbern(p_current)
      }

      # Logistic link
      for (j in 1:J) {
        logit(p_hist[j]) <- theta_hist[j]
      }
      logit(p_current) <- theta_current

      # Markov structure over theta
      theta_hist[1] ~ dnorm(0, 0.01)      # diffuse prior for first version
      for (j in 2:J) {
        theta_hist[j] ~ dnorm(theta_hist[j-1], tau_theta)
      }

      # Current version depends on last historical theta
      theta_current ~ dnorm(theta_hist[J], tau_current)

      # Priors for precisions (hyperparameters)
      tau_theta <- pow(sigma_theta, -2)
      tau_current <- pow(sigma_current, -2)
      sigma_theta ~ dunif(0, 5)
      sigma_current ~ dunif(0, 5)
    }
    "
    
    # --- Run model ---
    model <- jags.model(textConnection(model_string), data = data_jags,
                        n.chains = 3, n.adapt = 500)
    update(model, 500)
    
    samples <- coda.samples(model,
                            variable.names = c("p_current"),
                            n.iter = n_iter)
    
    samples_mat <- as.matrix(samples)
    posterior_current <- samples_mat[, "p_current"]
    
    # --- Metrics ---
    post_mean <- mean(posterior_current)
    post_ci <- quantile(posterior_current, c((1 - credible_interval) / 2,
                                             1 - (1 - credible_interval) / 2))
    
    bias_vec <- c(bias_vec, post_mean - p_current)
    mse_vec <- c(mse_vec, (post_mean - p_current)^2)
    
    if (p_current >= post_ci[1] && p_current <= post_ci[2]) {
      coverage_count <- coverage_count + 1
    }
  }
  
  # --- Summarize ---
  coverage_rate <- coverage_count / reps
  bias <- mean(bias_vec)
  mse <- mean(mse_vec)
  
  return(list(bias = bias,
              coverage = coverage_rate,
              mse = mse))
}




set.seed(123)

p_hist <- c(0.3, 0.4, 0.5, 0.6)
p_current <- 0.7

result <- run_commensurate_model(
  p_hist = p_hist,
  p_current = p_current,
  n = 10,
  credible_interval = 0.95,
  reps = 50,       # increase to 500 later
  n_iter = 2000
)

print(result)


