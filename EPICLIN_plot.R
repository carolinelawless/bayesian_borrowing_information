library(ggplot2)
library(dplyr)


source("functions.R")


params_gradual <- seq(0.6, 0.9, length = 20)   # true means per version
params_abrupt <- c(rep(0.6, 19), 0.9)
params_drift <- c(seq(0.6, 0.9, length = 10), seq(0.9, 0.6, length = 10))
params_stable <- rep(0.7, 20)


M <- 100
B <- 1000
lambda <- 9




a_theta <- 1
b_theta <- 1
mean_theta <- 0.5                          # prior mean
sd_theta   <- 0.5                          # prior sd
sigma      <- 0.5                          # measurement noise
p_eps <- 0.5

summarize_posterior <- function(thetas, epsilons, params, scenario_name) {
  
  K <- length(thetas)
  x <- 1:K
  
  # theta summaries
  theta_df <- data.frame(
    version = x,
    mean = sapply(thetas, mean, na.rm = TRUE),
    q10  = sapply(thetas, quantile, probs = 0.1, na.rm = TRUE),
    q90  = sapply(thetas, quantile, probs = 0.9, na.rm = TRUE),
    truth = params,
    type = "theta"
  )
  
  # epsilon summaries
  eps_df <- data.frame(
    version = x[-1] - 0.1,
    mean = sapply(epsilons[-1], mean, na.rm = TRUE),
    q10  = sapply(epsilons[-1], quantile, probs = 0.1, na.rm = TRUE),
    q90  = sapply(epsilons[-1], quantile, probs = 0.9, na.rm = TRUE),
    truth = NA,
    type = "epsilon"
  )
  
  bind_rows(theta_df, eps_df) %>%
    mutate(scenario = scenario_name)
}


# gradual
res_gradual <- posterior_sim_binomial(
  params_gradual, M, B, lambda,
  a_theta, b_theta, p_eps
)

df_gradual <- summarize_posterior(
  res_gradual$thetas,
  res_gradual$epsilons,
  params_gradual,
  "Gradual"
)

# abrupt
res_abrupt <- posterior_sim_binomial(
  params_abrupt, M, B, lambda,
  a_theta, b_theta, p_eps
)

df_abrupt <- summarize_posterior(
  res_abrupt$thetas,
  res_abrupt$epsilons,
  params_abrupt,
  "Abrupt"
)

# drift
res_drift <- posterior_sim_binomial(
  params_drift, M, B, lambda,
  a_theta, b_theta, p_eps
)

df_drift <- summarize_posterior(
  res_drift$thetas,
  res_drift$epsilons,
  params_drift,
  "Drift"
)

# stable
res_stable <- posterior_sim_binomial(
  params_stable, M, B, lambda,
  a_theta, b_theta, p_eps
)

df_stable <- summarize_posterior(
  res_stable$thetas,
  res_stable$epsilons,
  params_stable,
  "Stable"
)


plot_df <- bind_rows(
  df_gradual,
  df_abrupt,
  df_drift,
  df_stable
)


plot_df$scenario <- factor(
  plot_df$scenario,
  levels = c("Gradual", "Abrupt", "Drift", "Stable")
)


ggplot(plot_df,
       aes(x = version,
           y = mean,
           color = type)) +
  
  geom_point(size = 2) +
  
  geom_errorbar(aes(ymin = q10,
                    ymax = q90),
                width = 0.05) +
  
  # true theta only for theta rows
  geom_point(
    data = subset(plot_df, type == "theta"),
    aes(y = truth),
    inherit.aes = FALSE,
    color = "blue",
    alpha = 0.5,
    size = 2.5,
    x = subset(plot_df, type == "theta")$version
  ) +
  
  #facet_wrap(~scenario, nrow = 2) +
  facet_wrap(~scenario, nrow = 1) +
  
  scale_color_manual(values = c(
    theta = "black",
    epsilon = "red"
  )) +
  
  ylim(0, 1) +
  
  labs(
    x = "Version",
    y = expression(theta ~ "/" ~ epsilon),
    color = ""
  ) +
  
  theme_bw(base_size = 16) +
  
  theme(
    strip.text = element_text(face = "bold", size = 16),
    #legend.position = "bottom",
    legend.position = "none",
    panel.grid.minor = element_blank()
  )


lambdas <- 1:50 * 4
col_blue <- rgb(45, 118, 255, maxColorValue = 255)
col_purple <- rgb(145, 70, 255, maxColorValue = 255)
col1 <- rgb(31, 119, 180, maxColorValue = 255)
col2 <- rgb(174, 199, 232, maxColorValue = 255)

plot_power_curves <- function(lambdas,
                              power,
                              labels,
                              hypothesis) {
  
  if (hypothesis == 0) {
    ylab <- "erreur de type I"
    legend_pos <- "topright"
  } else {
    ylab <- "puissance"
    legend_pos <- "bottomright"
  }
  
  n_curves <- ncol(power)
  cols <- seq_len(n_curves)
  cols[1] <-  col1
  cols[2] <- col2
  
  
  par(mar = c(5, 6, 4, 2) + 0.1)
  
  
  
  plot(lambdas, power[,1],
       type = "l",
       lwd = 5,
       col = cols[1],
       cex.lab = 3,
       cex.axis = 1.5,
       ylim = c(0,1),
       xlab = expression(lambda),
       ylab = ylab)
  
  if (n_curves > 1) {
    for (i in 2:n_curves) {
      lines(lambdas, power[,i],
            lwd = 5,
            col = cols[i])
    }
  }
  
  legend(legend_pos,
         legend = labels,
         col = cols,
         lwd = 5,
         cex = 3,
         bty = "n")
}

p_binomial_gradual <- c( 0.598, 0.673, 0.685, 0.733, 0.767, 0.766, 0.792, 0.804, 0.823, 0.848, 0.847, 0.881, 0.861, 0.884, 0.892, 0.883, 0.904, 0.921, 0.904, 0.918, 0.927, 0.924, 0.94, 0.942, 0.95, 0.948, 0.943, 0.961, 0.963, 0.96, 0.963, 0.945, 0.969, 0.967, 0.964, 0.967, 0.975, 0.979, 0.971, 0.974, 0.983, 0.985, 0.977, 0.984, 0.979, 0.983, 0.99, 0.983, 0.99, 0.986 )
p_binomial_jump <- c( 0.476, 0.529, 0.552, 0.607, 0.653, 0.674, 0.706, 0.757, 0.78, 0.816, 0.824, 0.853, 0.849, 0.867, 0.874, 0.895, 0.894, 0.91, 0.929, 0.931, 0.932, 0.931, 0.93, 0.936, 0.932, 0.945, 0.936, 0.954, 0.958, 0.965, 0.966, 0.964, 0.972, 0.968, 0.976, 0.984, 0.967, 0.965, 0.977, 0.98, 0.982, 0.98, 0.982, 0.986, 0.99, 0.984, 0.989, 0.983, 0.99, 0.99 )
p_binomial_drift <- c( 0.357, 0.256, 0.196, 0.16, 0.116, 0.09, 0.08, 0.058, 0.073, 0.043, 0.045, 0.033, 0.015, 0.017, 0.013, 0.01, 0.015, 0.008, 0.009, 0.004, 0.003, 0.005, 0.004, 0.004, 0.001, 0.006, 0.002, 0, 0.001, 0, 0, 0, 0, 0.001, 0, 0, 0, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
p_binomial_stable <- c( 0.284, 0.231, 0.158, 0.115, 0.092, 0.068, 0.06, 0.037, 0.029, 0.017, 0.012, 0.011, 0.017, 0.01, 0.005, 0.006, 0.005, 0.002, 0.002, 0.001, 0.002, 0.001, 0, 0.001, 0.001, 0, 0, 0.001, 0, 0, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )

plot_power_curves(lambdas, cbind(p_binomial_gradual, p_binomial_jump), c("gradual", "abrupt"), 1)
plot_power_curves(lambdas, cbind(p_binomial_drift, p_binomial_stable), c("drift", "stable"), 0)






