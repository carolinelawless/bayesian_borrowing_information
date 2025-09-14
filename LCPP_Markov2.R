#########################################################################
# Hierarchical commensurate prior, multiple historical datasets
#########################################################################

remove(list = ls())

library(rjags)
library(coda)

n <- 20
sigmah1 <- 1
sigmah2 <- 1
sigmah3 <- 1
sigmac <- 1

# historical data
hist1 = rnorm(n, 1, sigmah1)
hist2 = rnorm(n, 1.5, sigmah2)
hist3 = rnorm(n, 2, sigmah3)
histmean1 = mean(hist1)
histvar1 = var(hist1)
histmean2 = mean(hist2)
histvar2 = var(hist2)
histmean3 = mean(hist3)
histvar3 = var(hist3)

# current data
y = rnorm(n, 2.5, sigmac)


# --- Data for JAGS ---
data_jags <- list(
  y = y,
  n = n,
  n_hist = rep(n, 3),
  mean_hist = c(histmean1, histmean2, histmean3),
  var_hist = c(histvar1, histvar2, histvar3)
)


# --- JAGS model ---
model_string <- "
model {
  for(j in 1:3){
      logtau[j] ~ dnorm(log(5), 10)   
      tau[j] <- exp(logtau[j])
      alpha[j] ~ dbeta(tau[j], 1)
      w[j] <- 1/(var_hist[j]/(alpha[j]*n_hist[j]) + 1/tau[j])
      num[j] <- w[j]*mean_hist[j]
    }
  mu_hist <- sum(num[]) / sum(w[])
  v_hist <- 1 / sum(w[])

  # Prior for mu centered at historical weighted mean
  mu ~ dnorm(mu_hist, 1/v_hist)

  # Likelihood for current data
  for(i in 1:n){
      y[i] ~ dnorm(mu, 1/sigmac^2)
  }

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
