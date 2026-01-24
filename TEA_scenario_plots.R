remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_scenarios.R")
par(mfrow = c(1,1))

################
x <- 1:5
cols <- c("#4E79A7", "#59A14F", "#9C755F", "#B07AA1", "#76B7B2")

################

params_list <- scenarios[1:8]
n_scenarios <- length(params_list)

if (n_scenarios > length(cols)) {
  cols <- colorRampPalette(cols)(n_scenarios)
}

plot(x, scenarios[[1]], 
     type = "p", pch = 16, cex = 1.3,
     col = cols[1],
     ylim = c(0, 1),
     xlab = "version",
     ylab = "true θ")

for(i in 2:length(params_list)){
  points(x, params_list[[i]], col = cols[i], pch = 16, cex = 1.3)
}

legend("topleft",
       legend = c("θ = (0.1, 0.1, 0.1, 0.1, 0.2)", "θ = (0.1, 0.1, 0.1, 0.1, 0.3)", "θ = (0.1, 0.1, 0.1, 0.1, 0.4)", "θ = (0.1, 0.1, 0.1, 0.1, 0.5)", "θ = (0.1, 0.1, 0.1, 0.1, 0.6)", "θ = (0.1, 0.1, 0.1, 0.1, 0.7)", "θ = (0.1, 0.1, 0.1, 0.1, 0.8)", "θ = (0.1, 0.1, 0.1, 0.1, 0.9)"),
       col = cols, pch = 16, cex = 0.8)

################

params_list <- scenarios[9:16]
n_scenarios <- length(params_list)

if (n_scenarios > length(cols)) {
  cols <- colorRampPalette(cols)(n_scenarios)
}

plot(x, params_list[[1]], 
     type = "p", pch = 16, cex = 1.3,
     col = cols[1],
     ylim = c(0, 1),
     xlab = "version",
     ylab = "true θ")

for(i in 2:length(params_list)){
  points(x, params_list[[i]], col = cols[i], pch = 16, cex = 1.3)
}

legend("topleft",
       legend = c("θ = (0.1, 0.125, 0.15, 0.175, 0.2)", "θ = (0.1, 0.15, 0.2, 0.25, 0.3)", "θ = (0.1, 0.175, 0.25, 0.325, 0.4)", "θ = (0.1, 0.2, 0.3, 0.4, 0.5)", "θ = (0.1, 0.225, 0.35, 0.475, 0.6)", "θ = (0.1, 0.25, 0.4, 0.55, 0.7)", "θ = (0.1, 0.275, 0.45, 0.625, 0.8)", "θ = (0.1, 0.3, 0.5, 0.7, 0.9)"),
       col = cols, pch = 16, cex = 0.8)

################

params_list <- scenarios[17:20]
n_scenarios <- length(params_list)

if (n_scenarios > length(cols)) {
  cols <- colorRampPalette(cols)(n_scenarios)
}

plot(x, params_list[[1]], 
     type = "p", pch = 16, cex = 1.3,
     col = cols[1],
     ylim = c(0, 1),
     xlab = "version",
     ylab = "true θ")

for(i in 2:length(params_list)){
  points(x, params_list[[i]], col = cols[i], pch = 16, cex = 1.3)
}

legend("topleft",
       legend = c("θ = (0.1, 0.1, 0.9, 0.1, 0.1)", "θ = (0.1, 0.5, 0.9, 0.5, 0.1)", "θ = (0.1, 0.1, 0.1, 0.9, 0.1)", "θ = (0.1, 0.4, 0.7, 0.9, 0.1)"),
       col = cols, pch = 16, cex = 0.8)

################

params_list <- scenarios[21:25]

n_scenarios <- length(params_list)

if (n_scenarios > length(cols)) {
  cols <- colorRampPalette(cols)(n_scenarios)
}

plot(x, params_list[[1]], 
     type = "p", pch = 16, cex = 1.3,
     col = cols[1],
     ylim = c(0, 1),
     xlab = "version",
     ylab = "true θ")

for(i in 2:length(params_list)){
  points(x, params_list[[i]], col = cols[i], pch = 16, cex = 1.3)
}

legend("topleft",
       legend = c("θ = (0.1, 0.1, 0.1, 0.1, 0.1)", "θ = (0.2, 0.2, 0.2, 0.2, 0.2)", "θ = (0.3, 0.3, 0.3, 0.3, 0.3)", "θ = (0.4, 0.4, 0.4, 0.4, 0.4)", "θ = (0.5, 0.5, 0.5, 0.5, 0.5)"),
       col = cols, pch = 16, cex = 0.8)



