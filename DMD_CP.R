#########################################################################
# Hierarchical commensurate prior with 3 historical datasets + current
# Each historical dataset likelihood is raised to a power a0[j]
#########################################################################

remove(list = ls())

library(rjags)
library(coda)

# ------------------------------
# Data setup
# ------------------------------
N <- 100
sigma <- 2

# Historical sample means
Y1 <- 2
Y2 <- 2
Y3 <- 10

# Current sample mean
Y  <- 12

# Precisions (all datasets same size and sigma for simplicity)
prec_hist_base <- N / sigma^2
prec_curr <- N / sigma^2

# JAGS data list
data_jags <- list(
  Y1 = Y1,
  Y2 = Y2,
  Y3 = Y3,
  Y  = Y,
  prec_hist_base = prec_hist_base,
  prec_curr = prec_curr
)

# ------------------------------
# JAGS model
# ------------------------------
model_string <- "
model {
  # Priors for tau and a0
  for (i in 1:3) {
    logtau[i] ~ dnorm(log(5), 1.0)    # weakly informative prior on commensurability
    tau[i] <- exp(logtau[i])
    a0[i] ~ dbeta(tau[i], 1)          # power prior weight for each dataset
  }

  # Historical 1
  theta[1] ~ dnorm(0, 1.0E-6)                     # vague prior
  Y1 ~ dnorm(theta[1], prec_hist_base * a0[1])    # downweighted likelihood

  # Historical 2
  theta[2] ~ dnorm(theta[1], tau[1])              # commensurate link
  Y2 ~ dnorm(theta[2], prec_hist_base * a0[2])    # downweighted likelihood

  # Historical 3
  theta[3] ~ dnorm(theta[2], tau[2])              # commensurate link
  Y3 ~ dnorm(theta[3], prec_hist_base * a0[3])    # downweighted likelihood

  # Current
  theta_curr ~ dnorm(theta[3], tau[3])            # commensurate link
  Y ~ dnorm(theta_curr, prec_curr)                # full likelihood
}
"

# ------------------------------
# Run the model
# ------------------------------
model <- jags.model(textConnection(model_string),
                    data = data_jags,
                    n.chains = 3,
                    n.adapt = 1000)

update(model, 1000)  # burn-in

samples <- coda.samples(model,
                        variable.names = c("theta", "theta_curr", "tau", "a0"),
                        n.iter = 5000)

# ------------------------------
# Posterior summaries
# ------------------------------
samples_mat <- as.matrix(samples)

head(samples_mat)

cat("Posterior means:\n")
cat("theta[1]:", mean(samples_mat[, 'theta[1]']), "\n")
cat("theta[2]:", mean(samples_mat[, 'theta[2]']), "\n")
cat("theta[3]:", mean(samples_mat[, 'theta[3]']), "\n")
cat("theta_curr:", mean(samples_mat[, 'theta_curr']), "\n\n")

cat("tau means:\n")
print(colMeans(samples_mat[, grep('tau', colnames(samples_mat))]))

cat("a0 means:\n")
print(colMeans(samples_mat[, grep('a0', colnames(samples_mat))]))

