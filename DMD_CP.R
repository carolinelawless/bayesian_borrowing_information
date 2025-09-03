#########################################################################
# Hierarchical commensurate prior 
#########################################################################

remove(list = ls())

library(rjags)
library(coda)

N <- 20
# --- Historical data ---
H <- 3      # number of versions
Y0 <- c(2,3,4)      # historical mean
n0 <- rep(N, H)       # historical sample size for each version
sigma0 <- rep(2, H)    # known sd for historical data (maybe put a prior on this)

# --- Current data ---
Y <- 5    # current mean
n <- N       # current sample size
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


# --- JAGS model ---
model_string <- "
model {
  # Prior on tau (precision between theta_curr and theta_hist)
  #tau ~ dgamma(0.1, 0.1)

  logtau ~ dnorm(log(5), 10)   
  tau <- exp(logtau)


  # Power prior weight a0 with Beta(tau,1) prior
  a0 ~ dbeta(tau, 1)

  # Historical theta (precision scaled by a0)
  theta_hist ~ dnorm(Y0, prec_hist_base * a0)
  # theta_hist ~ dnorm(Y0, prec_hist_base * a0 + 1.0E-6)


  # Current theta follows a normal centered at historical theta with precision tau
  theta_curr ~ dnorm(theta_hist, tau)

  # Likelihood for current data
  Y ~ dnorm(theta_curr, prec_curr)
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



