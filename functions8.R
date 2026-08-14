smc_sampler_binomial <- function(
    V = V,
    M = M,
    a_theta = a_theta,
    b_theta = b_theta,
    epsilon_const = epsilon_const,
    CSD = CSD
){
  
  
  K <- length(V)
  theta_all <- vector("list", K)
  epsilon_all <- vector("list", K)
  
  # first version uses baseline prior
  theta <- rbeta(M, a_theta, b_theta)
  
  for (k in 1:K) {
    data1 <- V[[k]]
    
    n1 <- length(data1)
    s1 <- sum(data1)
    
    
    logliks0 <- dbinom(
      s1,
      n1,
      theta,
      log = TRUE
    )
    
    logliks1 <- dbinom(
      s1,
      n1,
      pmax(theta - CSD, 0.01),
      log = TRUE
    )
    
    logliks2 <- dbinom(
      s1,
      n1,
      pmin(theta + CSD, 0.99),
      log = TRUE
    )
    
    logliks_borrow <- logliks0
    logliks_alternative <- pmax(logliks1, logliks2)
    epsilon <- ifelse(logliks_borrow > logliks_alternative, epsilon_const, 0)
    
    
    
    z <- rbinom(M, 1, epsilon)
    ind <- which(z == 0)
    
    if (length(ind) > 0) {
      theta[ind] <- rbeta(
        length(ind),
        a_theta,
        b_theta
      )
      
    }
    
    logliks <- dbinom(
      s1,
      n1,
      theta,
      log = TRUE
    )
    
    maxll <- max(logliks)
    w <- exp(logliks - maxll)
    w <- w / sum(w)
    samp <- sample(
      1:M,
      M,
      replace = TRUE,
      prob = w
    )
    
    theta <- theta[samp]
    if (k > 1) {
      epsilon <- epsilon[samp]
    }
    
    theta_all[[k]] <- theta
    epsilon_all[[k]] <- epsilon
    
  } 
  
  return(
    list(
      theta_all = theta_all,
      epsilon_all = epsilon_all
    )
  )
  
}



posterior_sim_binomial <- function(
    params = params,
    M = M,
    B = B,
    lambdas = lambdas,
    a_theta = a_theta,
    b_theta = b_theta,
    epsilon_const = epsilon_const,
    CSD = CSD
) {
  
  K <- length(params)
  thetas <- vector("list", K)
  epsilons <- vector("list", K)
  diffs <- numeric(B)
  
  for (k in 1:K) {
    thetas[[k]] <- numeric(0)
    epsilons[[k]] <- numeric(0)
  }
  
  
  for (b in 1:B) {
    
    n <- rpois(K, lambdas) + 1
    V <- vector("list", K)
    
    for (k in 1:K) {
      V[[k]] <- rbinom(
        n[k],
        1,
        params[k]
      )
    }
    
    smc_out <- smc_sampler_binomial(
      V = V,
      M = M,
      a_theta = a_theta,
      b_theta = b_theta,
      epsilon_const = epsilon_const,
      CSD = CSD
    )
    
    theta_all <- smc_out$theta_all
    epsilon_all <- smc_out$epsilon_all
    
    # -------------------------
    # Difference metric
    # -------------------------
    
    theta1 <- theta_all[[1]]
    thetaK <- theta_all[[K]]
    
    diffs[b] <- 
      mean(thetaK) -
      mean(theta1)
    
    
    # -------------------------
    # Store summaries
    # -------------------------
    
    for (k in 1:K) {
      
      thetas[[k]] <- c(
        thetas[[k]],
        mean(theta_all[[k]])
      )
      
      if (k > 1) {
        
        epsilons[[k]] <- c(
          epsilons[[k]],
          mean(epsilon_all[[k]])
        )
        
      }
    }
  }
  
  return(
    list(
      thetas = thetas,
      epsilons = epsilons,
      diffs = diffs
    )
  )
}



plot_trajectories <- function(thetas, epsilons, params, epsilon_type) {
  
  K <- length(thetas)
  x <- 1:K
  x_eps <- x - 0.1   # shift epsilon to the left
  
  # Summaries
  means_theta <- sapply(thetas, mean, na.rm = TRUE)
  q10_theta   <- sapply(thetas, quantile, probs = 0.1, na.rm = TRUE)
  q90_theta   <- sapply(thetas, quantile, probs = 0.9, na.rm = TRUE)
  
  eps_means <- sapply(epsilons, mean, na.rm = TRUE)
  eps_q10   <- sapply(epsilons, quantile, probs = 0.1, na.rm = TRUE)
  eps_q90   <- sapply(epsilons, quantile, probs = 0.9, na.rm = TRUE)
  
  
  
  # Base plot for theta
  plot(x, means_theta,
       ylim = c(0, 1),
       pch = 16, 
       col = "black",
       #cex.lab = 3,
       ylab = "θ / ε",
       xlab = "Version",
       main = epsilon_type,
       xaxt = "n")
  axis(1, at = x, labels = x)
  
  # Theta CI
  arrows(x, q10_theta, x, q90_theta,
         angle = 90, code = 3, length = 0.05)
  
  # True theta
  points(x, params, pch = 16, col = rgb(0,0,1,0.4), cex = 1.8)
  
  # Epsilon means (shifted)
  points(x_eps[-1], eps_means[-1], pch = 4, col = "red", cex = 1.2)
  
  # Epsilon CI (shifted)
  arrows(x_eps[-1], eps_q10[-1], x_eps[-1], eps_q90[-1],
         angle = 90, code = 3, length = 0.05, col = "red")
  
  # # Legend
  # legend("bottomright",
  #        legend = c("theta posterior mean", "true theta", "epsilon posterior mean"),
  #        pch = c(16, 16, 4),
  #        col = c("black", rgb(0,0,1,0.4), "red"),
  #        cex = 1.5,
  #        bty = "n")
  # 
}




