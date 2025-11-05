#########################################################################
# Hierarchical commensurate prior for multiple historical binomial datasets
#########################################################################

remove(list = ls())
library(rjags)
library(coda)

# --- Simulate data ---
n <- 10

p_hist <- c(0.3, 0.4, 0.5, 0.6)   # 4 historical datasets
p_current <- 0.7

# Simulate historical binary data
y_hist <- matrix(NA, nrow = length(p_hist), ncol = n)
for(j in 1:length(p_hist)){
  y_hist[j, ] <- rbinom(n, 1, p_hist[j])
}

apply(y_hist, 1, sum)

# Simulate current binary data
y_current <- rbinom(n, 1, p_current)

# --- Data for JAGS ---
data_jags <- list(
  y_hist = y_hist,
  y_current = y_current,
  n = n,
  J = length(p_hist)
)

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
                    n.chains = 3, n.adapt = 1000)
update(model, 1000)
samples <- coda.samples(model, variable.names = c("p_current", "p_hist",
                                                  "theta_current", "theta_hist"),
                        n.iter = 10000)



samples_mat <- as.matrix(samples)


mean(samples_mat[, "p_hist[1]"])
mean(samples_mat[, "p_hist[2]"])
mean(samples_mat[, "p_hist[3]"])
mean(samples_mat[, "p_hist[4]"])
mean(samples_mat[, "p_current"])

