remove(list = ls())
mu_all <- list()
M <- 100 #number of particles

v1 <- round(rnorm(5, 1, 1))
v2 <- round(rnorm(5, 2, 1))
v3 <- round(rnorm(5, 3, 1))
v4 <- round(rnorm(5, 4, 1))
v5 <- round(rnorm(5, 5, 1))
V <- rbind(v1, v2, v3, v4, v5)
K <- nrow(V)
epsilon <- rep(0.5, K)

k = 0
mu  <- round(runif(M, 0, 5))
mu_all[[length(mu_all) + 1]] <- mu
U <- rep(0, M)
mu_mat <- mu
U_mat <- U

while(k < K){
  k <- k+1
  U <- rbinom(M, 1, epsilon[k])
  
  mu2 <- round(runif(length(which(U == 0)), 0, 5))
  mu[which(U == 0)] <- mu2
  data <- V[k,]
  w <- vector(length = length(mu))
  for(i in 1:length(mu)){
    likelihoods <- dnorm(data, mu[i], 1)
    w[i] <- prod(likelihoods)
  }
  s <- sum(w)
  w <- w/s
  ESS <- 1/sum(w^2)
  mu_mat <- rbind(mu_mat, mu)
  U_mat <- rbind(U_mat, U)
  if(ESS < M/2){
    samp <- sample(1:M, M, replace = TRUE, prob = w)
    mu_mat <- mu_mat[,samp]
    U_mat <- U_mat[, samp]
    mu <- mu[samp]
    U <- U[samp]

  }
  mu_all[[length(mu_all) + 1]] <- mu
}

U_mat
mu_mat
mu_all


