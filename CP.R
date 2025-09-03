#########################################################################
# Hierarchical commensurate prior with tau prior and a0 ~ Beta(tau,1)
#########################################################################

remove(list = ls())

library(rjags)
library(coda)

# --- Historical data ---
Y0 <- 5      # historical mean
n0 <- 20       # historical sample size
sigma0 <- 2    # known sd for historical data

# --- Current data ---
Y <- 15      # current mean
n <- 20        # current sample size
sigma <- 2     # known sd for current data

# --- Compute base precisions ---
prec_curr <- n / sigma^2
prec_hist_base <- n0 / sigma0^2   # precision of historical mean (before a0 scaling)

# --- Data for JAGS ---
data_jags <- list(
  Y = Y,
  Y0 = Y0,
  prec_curr = prec_curr,
  prec_hist_base = prec_hist_base
)

data_jags <- list(
  Y0 = Y0,
  prec_hist_base = prec_hist_base
)

# --- JAGS model ---
model_string <- "
model {
  # Prior on tau (precision between theta_curr and theta_hist)
  tau ~ dgamma(0.1, 0.1)

  # Power prior weight a0 with Beta(tau,1) prior
  a0 ~ dbeta(tau, 1)

  # Historical theta (precision scaled by a0)
  # theta_hist ~ dnorm(Y0, prec_hist_base * a0)
  theta_hist ~ dnorm(Y0, prec_hist_base * a0 + 1.0E-6)


  # Current theta follows a normal centered at historical theta with precision tau
  theta_curr ~ dnorm(theta_hist, tau)

  # Likelihood for current data
  #Y ~ dnorm(theta_curr, prec_curr)
}
"

# --- Run model ---
model <- jags.model(textConnection(model_string), data = data_jags, n.chains = 3, n.adapt = 1000)
update(model, 1000)
samples <- coda.samples(model, variable.names = c("theta_hist", "theta_curr", "tau", "a0"), n.iter = 5000)



# Extract samples as a matrix
samples_mat <- as.matrix(samples)

head(samples_mat)
mean(samples_mat[,"theta_hist"])
mean(samples_mat[,"theta_curr"])
mean(samples_mat[, "tau"])
mean(samples_mat[, "a0"])
mean(rbeta(1e3, mean(samples_mat[, "tau"]), 1))

mean(samples_mat[, "tau"]/(1 + samples_mat[, "tau"]))


