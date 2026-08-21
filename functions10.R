smc_sampler_gaussian <- function(
    V = V,
    M = M,
    mean_theta = mean_theta,
    sd_theta = sd_theta,
    sd_obs = sd_obs,
    epsilon_const = epsilon_const,
    CSD = CSD
){
  
  
  K <- length(V)
  
  theta_all <- vector("list", K)
  epsilon_all <- vector("list", K)
  
  
  loglik_normal <- function(theta, data, sd_obs) {
    
    n <- length(data)
    sum_y <- sum(data)
    sum_y2 <- sum(data^2)
    
    loglik <- (
      -n / 2 * log(2 * pi * sd_obs^2)
      -
        1 / (2 * sd_obs^2) *
        (sum_y2 - 2 * theta * sum_y + n * theta^2)
    )
    
    return(loglik)
  }
  

  
  theta <- rnorm(
    M,
    mean = mean_theta,
    sd = sd_theta
  )

  
  for (k in 1:K) {
    
    data1 <- V[[k]]
    
    logliks0 <- loglik_normal(
      theta = theta,
      data = data1,
      sd_obs = sd_obs
    )
    
    logliks1 <- loglik_normal(
      theta = theta - CSD,
      data = data1,
      sd_obs = sd_obs
    )
    
    
    logliks2 <- loglik_normal(
      theta = theta + CSD,
      data = data1,
      sd_obs = sd_obs
    )
    
    logliks_borrow <- logliks0
    
    logliks_alternative <- pmax(
      logliks1,
      logliks2
    )
    
    epsilon <- ifelse(
      logliks_borrow > logliks_alternative,
      epsilon_const,
      0
    )
    
    
    z <- rbinom(
      M,
      size = 1,
      prob = epsilon
    )
    
    ind <- which(z == 0)
    
    if (length(ind) > 0) {
      
      theta[ind] <- rnorm(
        length(ind),
        mean = mean_theta,
        sd = sd_theta
      )
      
    }
    
    logliks <- loglik_normal(
      theta = theta,
      data = data1,
      sd_obs = sd_obs
    )
    
    maxll <- max(logliks)
    
    w <- exp(
      logliks - maxll
    )
    
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


posterior_sim <- function(
    params = params,
    M = M,
    B = B,
    lambdas = lambdas,
    epsilon_const = epsilon_const,
    CSD = CSD, 
    thres1 = thres1, 
    model = model,
    model_params = model_params
) {
  
  K <- length(params)
  
  thetas <- vector("list", K)
  epsilons <- vector("list", K)
  rho <- numeric(B)
  
  for (k in 1:K) {
    thetas[[k]] <- numeric(0)
    epsilons[[k]] <- numeric(0)
  }
  
  
  for (b in 1:B) {
    
    n <- rpois(K, lambdas) + 1
    V <- vector("list", K)
    
    
    # --------------------------------------------------
    # Generate data
    # --------------------------------------------------
    
    for (k in 1:K) {
      
      if (model == "binomial") {
        
        V[[k]] <- rbinom(
          n[k],
          size = 1,
          prob = params[k]
        )
        
      } else if (model == "gaussian") {
        
        V[[k]] <- rnorm(
          n[k],
          mean = params[k],
          sd = model_params$sd_obs
        )
        
      } else {
        
        stop(
          "model must be either 'binomial' or 'gaussian'"
        )
        
      }
    }
    
    
    # --------------------------------------------------
    # Run appropriate SMC sampler
    # --------------------------------------------------
    
    if (model == "binomial") {
      
      smc_out <- smc_sampler_binomial(
        V = V,
        M = M,
        a_theta = model_params$a_theta,
        b_theta = model_params$b_theta,
        epsilon_const = epsilon_const,
        CSD = CSD
      )
      
    } else if (model == "gaussian") {
      
      smc_out <- smc_sampler_gaussian(
        V = V,
        M = M,
        mean_theta = model_params$mean_theta,
        sd_theta = model_params$sd_theta,
        sd_obs = model_params$sd_obs,
        epsilon_const = epsilon_const,
        CSD = CSD
      )
      
    }
    
    
    # --------------------------------------------------
    # Extract SMC output
    # --------------------------------------------------
    
    theta_all <- smc_out$theta_all
    epsilon_all <- smc_out$epsilon_all
    
    
    # --------------------------------------------------
    # Posterior probability of difference exceeding
    # threshold
    # --------------------------------------------------
    
    rho[b] <- sum(
      theta_all[[K]] - theta_all[[1]] > thres1
    ) / M
    
    
    # --------------------------------------------------
    # Store posterior means
    # --------------------------------------------------
    
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
  
  
  # --------------------------------------------------
  # Return
  # --------------------------------------------------
  
  return(
    list(
      thetas = thetas,
      epsilons = epsilons,
      rho = rho
    )
  )
  
}


plot_estimates <- function(thetas, epsilons, params, lambda) {
  
  K <- length(thetas)
  x <- 1:K
  
  means_theta <- sapply(thetas, mean, na.rm = TRUE)
  q10_theta   <- sapply(thetas, quantile, probs = 0.1, na.rm = TRUE)
  q90_theta   <- sapply(thetas, quantile, probs = 0.9, na.rm = TRUE)
  
  
  prop <- function(eps){
    proportion <- length(which(eps !=0))/length(eps)
    return(proportion)
  }
  proportions <- sapply(epsilons, prop)
  proportions <- as.vector(proportions)
  
  # Base plot for theta
  plot(x, means_theta,
       ylim = c(0, 1),
       pch = 16,
       col = "black",
       #cex.lab = 3,
       ylab = "θ",
       xlab = "Version",
       #main = paste0("λ = ", lambda),
       xaxt = "n")
  axis(1, at = x, labels = x)
  
  # Theta CI
  arrows(x, q10_theta, x, q90_theta,
         angle = 90, code = 3, length = 0.05)
  
  # True theta
  points(x, params, pch = 16, col = rgb(0,0,1,0.4), cex = 2)
  
  points(x[-1], proportions[-1], pch = 4, col = "red")
  
  legend("bottomleft",
         legend = c("θ posterior mean", "true theta", "η posterior mean"),
         pch = c(16, 16, 4),
         col = c("black", rgb(0,0,1,0.4), "red"),
         cex = 1.5,
         bty = "n"
  )
  
}


plot_rho <- function(n_versions, scenario, lambda_vector,
                     proba_TEA_0,
                     proba_TEA_0.5,
                     proba_TEA_1,
                     proba_TEA_EB_0.5,
                     proba_TEA_EB_1) {
  
  
  if(scenario == "stable"){
    plot(lambda_vector, proba_TEA_0, type = "l",
         col = "black",
         xlab = "λ",
         ylab = "α",
         main = paste0(n_versions," versions, ", scenario, " scenario"),
         ylim = c(0, 1)
    )
  }else{
    plot(lambda_vector, proba_TEA_0, type = "l",
         col = "black",
         xlab = "λ",
         ylab = "Power",
         main = paste0(n_versions," versions, ", scenario, " scenario"),
         ylim = c(0, 1)
    )
  }
  
  
  
  lines(lambda_vector, proba_TEA_0.5, col = "blue")
  lines(lambda_vector, proba_TEA_1, col = "blue", lty = 2)
  lines(lambda_vector, proba_TEA_EB_0.5, col = "red")
  lines(lambda_vector, proba_TEA_EB_1, col = "red", lty = 2)
  
  if(scenario == "stable"){
    legend("topright",
           legend = c("ε = 0 (no borrowing)", "ε = 0.5 without EB", "ε = 1 without EB (full borrowing)", "ε = 0.5 with EB", "ε = 1 with EB"),
           col = c("black", "blue", "blue", "red", "red"),
           lty = c(1, 1, 2, 1, 2))
  }else{
    legend("bottomright",
           legend = c("ε = 0 (no borrowing)", "ε = 0.5 without EB", "ε = 1 without EB (full borrowing)", "ε = 0.5 with EB", "ε = 1 with EB"),
           col = c("black", "blue", "blue", "red", "red"),
           lty = c(1, 1, 2, 1, 2))
  }
  
  
}
