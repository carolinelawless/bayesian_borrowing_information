#########################################################################
# Hierarchical commensurate prior, multiple historical datasets, Markov AR(1) model on the commensurability parameters
#########################################################################

remove(list = ls())


library(rjags)
library(coda)

n <- 200
sigmah1 <- 1
sigmah2 <- 1
sigmah3 <- 1
sigmac <- 1

# historical data
hist1 = rnorm(n, 1, sigmah1)
hist2 = rnorm(n, 1.5, sigmah2)
hist3 = rnorm(n, 2, sigmah3)
meanhist1 = mean(hist1)
varhist1 = var(hist1)
meanhist2 = mean(hist2)
varhist2 = var(hist2)
meanhist3 = mean(hist3)
varhist3 = var(hist3)

# current data
y = rnorm(n, 2.5, sigmac)


# --- Data for JAGS ---
data_jags <- list(
  y = y,
  n = n,
  n_hist = rep(n, 3),
  mean_hist = c(meanhist1, meanhist2, meanhist3),
  var_hist = c(varhist1, varhist2, varhist3)
)


# --- JAGS model ---
model_string <- "
model {
  a ~ dnorm(0, 0.001) # vague
  b ~ dunif(0, 1) # ensure stationarity
  c ~ dunif(0, 1) # positive sd

  logtau[1] ~ dnorm(log(5), 0.1)     # prior for first version
  for (j in 2:3){
    logtau[j] ~ dnorm(a + b*logtau[j-1], 1/c^2)
  }

  for(j in 1:3){
      tau[j] <- exp(logtau[j])
      alpha[j] ~ dbeta(tau[j], 1)
      w[j] <- 1/(var_hist[j]/(alpha[j]*n_hist[j]) + 1/tau[j])
      num[j] <- w[j]*mean_hist[j]
    }
  mu_hist <- sum(num[]) / sum(w[])
  v_hist <- 1 / sum(w[])

  # Prior for mu centered at historical weighted mean
  mu ~ dnorm(mu_hist, 1/v_hist)

  prec ~ dgamma(0.001, 0.001) 
  # Likelihood for current data
  for(i in 1:n){
      y[i] ~ dnorm(mu, prec)
  }

}
"

# --- Run model ---
model <- jags.model(textConnection(model_string), data = data_jags, n.chains = 3, n.adapt = 1000)
update(model, 1000)


samples <- coda.samples(model, variable.names = c("tau", "alpha", "prec", "a", "b", "c", "mu"), n.iter = 5000)



# Extract samples as a matrix
samples_mat <- as.matrix(samples)

head(samples_mat)

c(meanhist1, meanhist2, meanhist3, mean(y))
mean(samples_mat[, "tau[1]"])
mean(samples_mat[, "tau[2]"])
mean(samples_mat[, "tau[3]"])
mean(samples_mat[, "a"])
mean(samples_mat[, "b"])
mean(samples_mat[, "c"])

mean(samples_mat[, "mu"]) #2.572
var(samples_mat[, "mu"]) #0.006

quantile(samples_mat[, "mu"], 0.995) - quantile(samples_mat[, "mu"], 0.005) #0.41


c(mean(samples_mat[,"alpha[1]"]), mean(samples_mat[,"alpha[2]"]), mean(samples_mat[,"alpha[3]"])) # 0.21, 0.34, 0.47

