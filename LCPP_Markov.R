#########################################################################
# Location commensurate power prior (LCPP)
# CL 10/09/2025
#########################################################################

remove(list = ls())

library(rjags)
library(coda)


model_string <- "
model {
  
  precision0 ~ dgamma(1.0E-3, 1.0E-3)
  theta0 ~ dnorm(0, 1.0E-6)  # vague prior
  
  for (j in 1:3) {
    logtau[j] ~ dt(0, pow(30,-2), 1)   # Cauchy(0,30)
    tau[j] <- exp(logtau[j]) + 1.0E-6
    
    # g(tau) = max(logtau, 1)
    logtau_clip[j] <- logtau[j] * step(logtau[j] - 1) + 1 * (1 - step(logtau[j] - 1))
    gtau[j] <- logtau_clip[j]
    
    alpha[j] ~ dbeta(gtau[j], 1)
    
    precision[j] ~ dgamma(1.0E-3, 1.0E-3)
  }
  
  # Sequential prior on theta
  theta[1] ~ dnorm(theta0, tau[1])
  for (j in 2:3) {
    theta[j] ~ dnorm(theta[j-1], tau[j])
  }
  
  # Likelihoods
  for(i in 1:n0) {
    Y0[i] ~ dnorm(theta0, precision0)
  }
  for(j in 1:3) {
    for(i in 1:n[j]) {
      Y[j,i] ~ dnorm(theta[j], precision[j])
    }
  }
}
"





sd_theta <- 1
sigma <- 1
theta0 <- 0
n0 <- 1e1
n <- rep(1e1, 3)

# Data
theta1 <- rnorm(1, theta0, sd_theta)
theta2 <- rnorm(1, theta1, sd_theta)
theta3 <- rnorm(1, theta2, sd_theta)
Y0 <- rnorm(n0, theta0, sigma)
Y1 <- rnorm(n[1], theta1, sigma)
Y2 <- rnorm(n[2], theta2, sigma)
Y3 <- rnorm(n[3], theta3, sigma)
Y <- rbind(Y1, Y2, Y3)




# JAGS data list
data_jags <- list(
  Y0 = Y0,
  Y = Y,
  n0 = n0,
  n = n
)

model <- jags.model(textConnection(model_string),
                    data = data_jags,
                    n.chains = 3,
                    n.adapt = 1000)

update(model, 1000)  # burn-in
samples <- coda.samples(model,
                        variable.names = c("theta3", "alpha", "tau", "sigma"),
                        n.iter = 5000)
samples_mat <- as.matrix(samples)
